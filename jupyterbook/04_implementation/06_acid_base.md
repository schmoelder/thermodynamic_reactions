---
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

(implementation-acid-base)=
# pH and Acid-Base Speciation

Proton-transfer reactions are always fast relative to chromatographic transport (@acid-base) and are treated throughout as algebraic constraints.
The `pKa` factory, water autoionisation, and Davies activity corrections each introduced in earlier chapters now combine: solving a realistic buffer system requires all three simultaneously.

```{code-cell} ipython3
import numpy as np

from reactions.api import (
    Species, Component,
    ActivityCoefficientDavies,
    EquilibriumConstant, EquilibriumConstantVantHoff,
    ThermodynamicReaction, ReactionModel,
    IonicStrengthBackground,
    pKa,
)
from reactions.solver import solve_equilibrium

C_REF = 1000.0  # mol/m³

proton    = Component("proton",    [Species("H+",   charge=+1)])
hydroxide = Component("hydroxide", [Species("OH-",  charge=-1)])
water     = Component("water",     [Species("H2O",  charge=0)])
```


## The pKa factory

`pKa()` converts a p$K_a$ value measured at 25 °C into an equilibrium constant.
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


## Monoprotic buffer: acetic acid

A monoprotic weak acid establishes the pattern used throughout.
Acetic acid dissociates as $\ce{HAc <=> Ac- + H+}$; water autoionisation is included as a second constraint so that pH is exact near neutrality:

```{code-cell} ipython3
acetic = Component("acetic", [Species("HAc", charge=0), Species("Ac-", charge=-1)])

model_acetic = ReactionModel(
    components=[acetic, proton, hydroxide, water],
    reactions=[
        ThermodynamicReaction(
            "HAc <-> Ac- + H+",
            mode="equil",
            equilibrium_constant=k_acetic,
        ),
        ThermodynamicReaction(
            "H2O <-> H+ + OH-",
            mode="equil",
            equilibrium_constant=EquilibriumConstant(1e-14),
        ),
    ],
    T=298.15,
)

c_tot = 100.0   # mol/m³ = 100 mM total acetate
c0    = {"HAc": c_tot, "Ac-": 1e-6, "H+": 1e-4, "OH-": 1e-10 * C_REF}
c_eq  = solve_equilibrium(model_acetic, c0, T=298.15, prescribed={"H2O": C_REF})

pH    = -np.log10(c_eq["H+"] / C_REF)
pKa_v = -np.log10(k_acetic.K(298.15))
pH_HH = pKa_v + np.log10(c_eq["Ac-"] / c_eq["HAc"])

print(f"pH (numerical)          = {pH:.4f}")
print(f"pH (Henderson-Hasselbalch) = {pH_HH:.4f}")
```

Water appears in the model so that the autoionisation constraint $\ce{H2O <=> H+ + OH-}$ is stoichiometrically balanced.
The `prescribed={"H2O": C_REF}` argument fixes $c_{\ce{H2O}} = c^\circ = 1000\,\mathrm{mol/m}^3$ throughout the solve: no balance equation is written for it, and each constraint it appears in receives activity $a_{\ce{H2O}} = c_{\ce{H2O}}/c^\circ = 1$.
The value $1000\,\mathrm{mol/m}^3$ is not the physical concentration of water ($\approx 55\,\mathrm{mol/L}$); it is the standard-state reference that gives unit solvent activity — the same approximation already built into the tabulated $K_w = 10^{-14}$.
In dilute aqueous solutions the mole fraction of water is very close to one, so $a_{\ce{H2O}} \approx 1$ regardless; setting $c_{\ce{H2O}} = c^\circ$ makes this exact within the model.
With $a_{\ce{H2O}} = 1$ the equilibrium condition reduces to $a_{\ce{H+}}\,a_{\ce{OH-}} = K_w$, where $K_w = 10^{-14}$ is the dimensionless water-ion product at 25 °C: the conventional $\mathrm{p}K_w = 14$ expressed as an equilibrium constant with activities referenced to $c^\circ = 1\,\mathrm{mol/L}$.
This is a prototype convention: in CADET-Core the transport solver owns the full concentration vector and passes it unchanged to the reaction model at each step, so the reaction model never prescribes any species directly.

`prescribed` is not limited to solvents.
Fixing a non-solvent species, for example holding $c_{\ce{H+}}$ at a target value, turns `solve_equilibrium` into a speciation calculator at a given proton activity.
Sweeping a sequence of prescribed proton activities maps out equilibrium composition as a function of pH, equivalent to a pH-stat measurement; @implementation-buffer uses this pattern to compute buffer capacity numerically.

The numerical result matches Henderson-Hasselbalch (@acid-base) because the ideal-dilute approximation holds at 100 mM.
Activity corrections break this agreement at finite ionic strength, as shown below.


## Polyprotic buffers: phosphate

Phosphoric acid dissociates in three steps spanning the full biochemical pH range:

| Step | Reaction | p$K_a$ (25 °C) |
|------|----------|----------------|
| 1 | $\ce{H3PO4 <=> H2PO4- + H+}$ | 2.148 |
| 2 | $\ce{H2PO4- <=> HPO4^2- + H+}$ | 7.198 |
| 3 | $\ce{HPO4^2- <=> PO4^3- + H+}$ | 12.350 |

Each step is one `ThermodynamicReaction`.
All four phosphate forms belong to a single `Component` because total phosphate is conserved.
Species names include conventional charge notation for readability; the `charge` argument is authoritative for all activity calculations:

```{code-cell} ipython3
pKa1, pKa2, pKa3 = 2.148, 7.198, 12.350

phosphate = Component("phosphate", [
    Species("H3PO4",  charge=0),
    Species("H2PO4-", charge=-1),
    Species("HPO4-2", charge=-2),
    Species("PO4-3",  charge=-3),
])

model_phosphate_ideal = ReactionModel(
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
print(model_phosphate_ideal)
```

The analytic speciation fractions serve as both initial guesses for the Newton solver and a verification reference:

```{code-cell} ipython3
def phosphate_fractions(pH, Ka1, Ka2, Ka3):
    h = 10.0 ** (-pH)
    D = h**3 + Ka1*h**2 + Ka1*Ka2*h + Ka1*Ka2*Ka3
    return h**3/D, Ka1*h**2/D, Ka1*Ka2*h/D, Ka1*Ka2*Ka3/D

Ka1 = 10**(-pKa1); Ka2 = 10**(-pKa2); Ka3 = 10**(-pKa3)
c_tot_phosphate = 100.0   # mol/m³

print(f"{'pH':>6}  {'H3PO4 ana':>12}  {'H3PO4 num':>12}  {'error':>10}")
for pH in [4.0, 7.2, 10.0]:
    a0, a1, a2, a3 = phosphate_fractions(pH, Ka1, Ka2, Ka3)
    H  = 10.0**(-pH) * C_REF
    OH = 1e-14 * C_REF**2 / H
    c0 = {
        "H3PO4":  max(a0*c_tot_phosphate, 1e-10),
        "H2PO4-": max(a1*c_tot_phosphate, 1e-10),
        "HPO4-2": max(a2*c_tot_phosphate, 1e-10),
        "PO4-3":  max(a3*c_tot_phosphate, 1e-10),
        "H+": H, "OH-": OH,
    }
    c_eq = solve_equilibrium(model_phosphate_ideal, c0, T=298.15, prescribed={"H2O": C_REF})
    ana  = a0 * c_tot_phosphate
    print(f"{pH:>6.1f}  {ana:>12.4f}  {c_eq['H3PO4']:>12.4f}  {abs(c_eq['H3PO4']-ana):>10.2e}")
```


## Ionic strength corrections

Activity coefficients shift the apparent p$K_a$ values at nonzero ionic strength (@nonidealities, @acid-base).
`ActivityCoefficientDavies` is passed per reaction; ionic strength is recomputed from the current concentration state and passed into $\gamma_i$ before each residual evaluation, so activity corrections enter the nonlinear solve consistently (@implementation-activity):

```{code-cell} ipython3
model_phosphate_davies = ReactionModel(
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

At physiological ionic strength ($I = 150\ \mathrm{mol/m}^3$), phosphate is particularly sensitive because its species carry charges 0, $-1$, $-2$, and $-3$:

```{code-cell} ipython3
pH  = 7.2
a0, a1, a2, a3 = phosphate_fractions(pH, Ka1, Ka2, Ka3)
H   = 10.0**(-pH) * C_REF
OH  = 1e-14 * C_REF**2 / H
c0  = {
    "H3PO4":  max(a0*c_tot_phosphate, 1e-10),
    "H2PO4-": max(a1*c_tot_phosphate, 1e-10),
    "HPO4-2": max(a2*c_tot_phosphate, 1e-10),
    "PO4-3":  max(a3*c_tot_phosphate, 1e-10),
    "H+": H, "OH-": OH,
}
eq_ideal  = solve_equilibrium(model_phosphate_ideal,  c0, T=298.15, prescribed={"H2O": C_REF})
eq_davies = solve_equilibrium(model_phosphate_davies, c0, T=298.15, prescribed={"H2O": C_REF})

print(f"At pH 7.2, I = 150 mol/m³ vs ideal (mol/m³):")
print(f"{'Species':>10}  {'ideal':>10}  {'Davies':>10}  {'shift':>10}")
for sp in ["H3PO4", "H2PO4-", "HPO4-2", "PO4-3"]:
    print(f"{sp:>10}  {eq_ideal[sp]:>10.3f}  {eq_davies[sp]:>10.3f}  "
          f"{eq_davies[sp]-eq_ideal[sp]:>+10.3f}")
```

---

The next chapter combines multiple equilibria to compute pH curves and buffer
capacity, including ionic strength coupling (@implementation-buffer).
