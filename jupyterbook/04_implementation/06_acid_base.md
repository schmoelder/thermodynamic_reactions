---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.0
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

(implementation-acid-base)=
# pH and Acid-Base Speciation

Proton-transfer reactions are always fast relative to chromatographic transport
(@acid-base), so they are handled as equilibrium constraints throughout this library.
This chapter implements the buffer machinery from Part 3: the `pKa()` factory,
polyprotic speciation, and activity corrections via the Davies equation.

```{code-cell} ipython3
import numpy as np

from reactions.api import (
    C_REF,
    Species, Component,
    ActivityCoefficientDavies,
    EquilibriumConstant, EquilibriumConstantVantHoff,
    ThermodynamicReaction, ReactionModel,
    IonicStrengthBackground,
    pKa,
)
from reactions.solver import solve_equilibrium

# Reusable components shared across examples
proton    = Component("proton",    [Species("H+",   charge=+1)])
hydroxide = Component("hydroxide", [Species("OH-",  charge=-1)])
water     = Component("water",     [Species("H2O",  charge=0, is_solvent=True)])
```

## The pKa factory

`pKa()` converts a p$K_a$ value measured at 25 °C into an equilibrium constant,
so the caller never has to compute $K_a = 10^{-\text{p}K_a}$ by hand.
Without `dH` the constant is temperature-independent:

```{code-cell} ipython3
k_acetic   = pKa(4.756)
k_ammonium = pKa(9.245)

print(f"Acetic acid  Ka(25 °C) = {k_acetic.K(298.15):.3e}  "
      f"pKa = {-np.log10(k_acetic.K(298.15)):.3f}")
print(f"Ammonium     Ka(25 °C) = {k_ammonium.K(298.15):.3e}  "
      f"pKa = {-np.log10(k_ammonium.K(298.15)):.3f}")
```

Passing `dH` (J/mol) adds van't Hoff temperature dependence; the factory
back-calculates $\Delta S^\circ$ so that $K(298.15\,\text{K}) = 10^{-\text{p}K_a}$
exactly.
The temperature-dependent form and its effect on $K(T)$ are demonstrated in
@implementation-equilibrium.

## Polyprotic buffers: phosphate

Phosphoric acid dissociates in three steps spanning the full biochemical pH range:

| Step | Reaction | p$K_a$ (25 °C) |
|------|----------|----------------|
| 1 | H₃PO₄ ⇌ H₂PO₄⁻ + H⁺ | 2.148 |
| 2 | H₂PO₄⁻ ⇌ HPO₄²⁻ + H⁺ | 7.198 |
| 3 | HPO₄²⁻ ⇌ PO₄³⁻ + H⁺ | 12.350 |

Each step is one `ThermodynamicReaction`.
All four phosphate forms belong to a single `Component` because total phosphate is
conserved.

```{code-cell} ipython3
pKa1, pKa2, pKa3 = 2.148, 7.198, 12.350

phosphate = Component("phosphate", [
    Species("H3PO4",  charge=0),
    Species("H2PO4-", charge=-1),
    Species("HPO4-2", charge=-2),
    Species("PO4-3",  charge=-3),
])

model_phos = ReactionModel(
    components=[phosphate, proton, hydroxide, water],
    reactions=[
        ThermodynamicReaction("H3PO4 <-> H2PO4- + H+",  mode="equil",
                              equilibrium_constant=pKa(pKa1)),
        ThermodynamicReaction("H2PO4- <-> HPO4-2 + H+", mode="equil",
                              equilibrium_constant=pKa(pKa2)),
        ThermodynamicReaction("HPO4-2 <-> PO4-3 + H+",  mode="equil",
                              equilibrium_constant=pKa(pKa3)),
        ThermodynamicReaction("H2O <-> H+ + OH-",        mode="equil",
                              equilibrium_constant=EquilibriumConstant(1e-14)),
    ],
    T=298.15,
)
print(model_phos)
```

```{code-cell} ipython3
def phosphate_fractions(pH, Ka1, Ka2, Ka3):
    h = 10.0 ** (-pH)
    D = h**3 + Ka1*h**2 + Ka1*Ka2*h + Ka1*Ka2*Ka3
    return h**3/D, Ka1*h**2/D, Ka1*Ka2*h/D, Ka1*Ka2*Ka3/D

Ka1 = 10**(-pKa1); Ka2 = 10**(-pKa2); Ka3 = 10**(-pKa3)
c_tot_phos = 100.0   # mol/m³

print(f"{'pH':>6}  {'H3PO4 ana':>12}  {'H3PO4 num':>12}  {'error':>10}")
for pH in [4.0, 7.2, 10.0]:
    a0, a1, a2, a3 = phosphate_fractions(pH, Ka1, Ka2, Ka3)
    H  = 10.0**(-pH) * C_REF
    OH = 1e-14 * C_REF**2 / H
    c0 = {
        "H3PO4":  max(a0*c_tot_phos, 1e-10),
        "H2PO4-": max(a1*c_tot_phos, 1e-10),
        "HPO4-2": max(a2*c_tot_phos, 1e-10),
        "PO4-3":  max(a3*c_tot_phos, 1e-10),
        "H+": H, "OH-": OH,
    }
    c_eq = solve_equilibrium(model_phos, c0, T=298.15)
    ana  = a0 * c_tot_phos
    print(f"{pH:>6.1f}  {ana:>12.4f}  {c_eq['H3PO4']:>12.4f}  {abs(c_eq['H3PO4']-ana):>10.2e}")
```

## Ionic strength corrections

In real buffers, activity coefficients shift the apparent p$K_a$ values at nonzero
ionic strength (see @nonidealities, @acid-base).
`ActivityCoefficientDavies` implements the Davies equation, valid up to
$I \approx 0.5\ \text{mol/L}$, and is passed per reaction:

```{code-cell} ipython3
model_phos_davies = ReactionModel(
    components=[phosphate, proton, hydroxide, water],
    reactions=[
        ThermodynamicReaction("H3PO4 <-> H2PO4- + H+",  mode="equil",
                              equilibrium_constant=pKa(pKa1),
                              activity_coefficient=ActivityCoefficientDavies()),
        ThermodynamicReaction("H2PO4- <-> HPO4-2 + H+", mode="equil",
                              equilibrium_constant=pKa(pKa2),
                              activity_coefficient=ActivityCoefficientDavies()),
        ThermodynamicReaction("HPO4-2 <-> PO4-3 + H+",  mode="equil",
                              equilibrium_constant=pKa(pKa3),
                              activity_coefficient=ActivityCoefficientDavies()),
        ThermodynamicReaction("H2O <-> H+ + OH-",        mode="equil",
                              equilibrium_constant=EquilibriumConstant(1e-14),
                              activity_coefficient=ActivityCoefficientDavies()),
    ],
    ionic_strength=IonicStrengthBackground(I_bg=150.0),
    T=298.15,
)
```

At physiological ionic strength ($I = 150\ \text{mol/m}^3$), phosphate is particularly
sensitive because its species carry charges 0, $-1$, $-2$, and $-3$:

```{code-cell} ipython3
pH  = 7.2
a0, a1, a2, a3 = phosphate_fractions(pH, Ka1, Ka2, Ka3)
H   = 10.0**(-pH) * C_REF
OH  = 1e-14 * C_REF**2 / H
c0  = {
    "H3PO4":  max(a0*c_tot_phos, 1e-10),
    "H2PO4-": max(a1*c_tot_phos, 1e-10),
    "HPO4-2": max(a2*c_tot_phos, 1e-10),
    "PO4-3":  max(a3*c_tot_phos, 1e-10),
    "H+": H, "OH-": OH,
}
eq_ideal  = solve_equilibrium(model_phos,        c0, T=298.15)
eq_davies = solve_equilibrium(model_phos_davies, c0, T=298.15)

print(f"At pH 7.2, I = 150 mol/m³ vs ideal (mol/m³):")
print(f"{'Species':>10}  {'ideal':>10}  {'Davies':>10}  {'shift':>10}")
for sp in ["H3PO4", "H2PO4-", "HPO4-2", "PO4-3"]:
    print(f"{sp:>10}  {eq_ideal[sp]:>10.3f}  {eq_davies[sp]:>10.3f}  "
          f"{eq_davies[sp]-eq_ideal[sp]:>+10.3f}")
```

---

The next chapter combines multiple equilibria to compute pH curves and buffer
capacity, including ionic strength coupling (@implementation-buffer).
