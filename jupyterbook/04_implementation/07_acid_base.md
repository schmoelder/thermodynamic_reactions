---
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

(implementation-acid-base)=
# pH and Acid-Base Speciation

Proton-transfer reactions are always fast relative to chromatographic transport (@acid-base) and are treated throughout as algebraic constraints.
The `pKa()` factory (@implementation-practical), water autoionisation, and Davies activity corrections now combine into a complete buffer model: solving a realistic system requires all three simultaneously.

```{code-cell} ipython3
import numpy as np

from reactions.api import (
    Species, Component,
    ActivityCoefficientDavies,
    ThermodynamicReaction, ReactionModel,
    IonicStrengthBackground,
    pKa,
    water, H_plus, OH_minus,
)
from reactions.solver import solve_equilibrium
```


## Monoprotic buffer: acetic acid

A monoprotic weak acid establishes the pattern used throughout.
Acetic acid dissociates as $\ce{HAc <=> Ac- + H+}$; water autoionisation is included as a second constraint so that pH is exact near neutrality:

```{code-cell} ipython3
acetic = Component("acetic", [Species("HAc", charge=0), Species("Ac-", charge=-1)])
k_acetic = pKa(4.756)

model_acetic = ReactionModel(
    components=[acetic, H_plus, OH_minus, water],
    reactions=[
        ThermodynamicReaction(
            "HAc <-> Ac- + H+",
            mode="equil",
            equilibrium_constant=k_acetic,
        ),
        ThermodynamicReaction(
            "H2O <-> H+ + OH-",
            mode="equil",
            equilibrium_constant=pKa(14.00),
        ),
    ],
    T=298.15,
)

c_tot = 100.0   # mol/m³ = 100 mM total acetate
c0    = {"HAc": c_tot, "Ac-": 1e-6, "H+": 1e-4, "OH-": 1e-7}
c_eq  = solve_equilibrium(model_acetic, c0, T=298.15, prescribed={"H2O": water.c_ref})

pH    = -np.log10(c_eq["H+"] / H_plus.c_ref)
pKa_v = -np.log10(k_acetic.K(298.15))
pH_HH = pKa_v + np.log10(c_eq["Ac-"] / c_eq["HAc"])

print(f"pH (numerical)          = {pH:.4f}")
print(f"pH (Henderson-Hasselbalch) = {pH_HH:.4f}")
```

Water appears in the model so that the autoionisation constraint $\ce{H2O <=> H+ + OH-}$ is stoichiometrically balanced.
The pre-built `water` component carries `molar_mass=0.018015` and `density=1000.0`, so `water.c_ref` $= \rho/M \approx 55{,}509\ \mathrm{mol/m^3}$ is the pure-component molar concentration of liquid water.
Passing `prescribed={"H2O": water.c_ref}` fixes $c_{\ce{H2O}}$ at that value throughout the solve, giving $a_{\ce{H2O}} = c/c^\circ = 1$; no balance equation is written for it.
With $a_{\ce{H2O}} = 1$ the equilibrium condition reduces to $a_{\ce{H+}}\,a_{\ce{OH-}} = K_w = 10^{-14}$, consistent with the dilute-solution convention under which $K_w$ is tabulated.
In production the `Solution` class (@implementation-practical) handles solvent prescriptions automatically; the explicit call here makes the mechanism transparent.

`prescribed` is the mechanism for parametric sweeps: fixing $c_{\ce{H+}}$ at a target value turns `solve_equilibrium` into a speciation calculator at a given proton activity.
Sweeping a sequence of fixed proton activities maps out equilibrium composition as a function of pH, equivalent to a pH-stat measurement; @implementation-buffer uses this pattern to compute buffer capacity numerically.

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
    components=[phosphate, H_plus, OH_minus, water],
    reactions=[
        ThermodynamicReaction("H3PO4 <-> H2PO4- + H+",  mode="equil",
                              equilibrium_constant=pKa(pKa1)),
        ThermodynamicReaction("H2PO4- <-> HPO4-2 + H+", mode="equil",
                              equilibrium_constant=pKa(pKa2)),
        ThermodynamicReaction("HPO4-2 <-> PO4-3 + H+",  mode="equil",
                              equilibrium_constant=pKa(pKa3)),
        ThermodynamicReaction("H2O <-> H+ + OH-",        mode="equil",
                              equilibrium_constant=pKa(14.00)),
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
    H  = 10.0**(-pH) * H_plus.c_ref
    OH = 1e-14 * H_plus.c_ref**2 / H
    c0 = {
        "H3PO4":  max(a0*c_tot_phosphate, 1e-10),
        "H2PO4-": max(a1*c_tot_phosphate, 1e-10),
        "HPO4-2": max(a2*c_tot_phosphate, 1e-10),
        "PO4-3":  max(a3*c_tot_phosphate, 1e-10),
        "H+": H, "OH-": OH,
    }
    c_eq = solve_equilibrium(model_phosphate_ideal, c0, T=298.15, prescribed={"H2O": water.c_ref})
    ana  = a0 * c_tot_phosphate
    print(f"{pH:>6.1f}  {ana:>12.4f}  {c_eq['H3PO4']:>12.4f}  {abs(c_eq['H3PO4']-ana):>10.2e}")
```


## Ionic strength corrections

Activity coefficients shift the apparent p$K_a$ values at nonzero ionic strength (@nonidealities, @acid-base).
`ActivityCoefficientDavies` is passed per reaction; ionic strength is recomputed from the current concentration state and passed into $\gamma_i$ before each residual evaluation, so activity corrections enter the nonlinear solve consistently (@implementation-activity):

```{code-cell} ipython3
model_phosphate_davies = ReactionModel(
    components=[phosphate, H_plus, OH_minus, water],
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
                              equilibrium_constant=pKa(14.00),
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
H   = 10.0**(-pH) * H_plus.c_ref
OH  = 1e-14 * H_plus.c_ref**2 / H
c0  = {
    "H3PO4":  max(a0*c_tot_phosphate, 1e-10),
    "H2PO4-": max(a1*c_tot_phosphate, 1e-10),
    "HPO4-2": max(a2*c_tot_phosphate, 1e-10),
    "PO4-3":  max(a3*c_tot_phosphate, 1e-10),
    "H+": H, "OH-": OH,
}
eq_ideal  = solve_equilibrium(model_phosphate_ideal,  c0, T=298.15, prescribed={"H2O": water.c_ref})
eq_davies = solve_equilibrium(model_phosphate_davies, c0, T=298.15, prescribed={"H2O": water.c_ref})

print(f"At pH 7.2, I = 150 mol/m³ vs ideal (mol/m³):")
print(f"{'Species':>10}  {'ideal':>10}  {'Davies':>10}  {'shift':>10}")
for sp in ["H3PO4", "H2PO4-", "HPO4-2", "PO4-3"]:
    print(f"{sp:>10}  {eq_ideal[sp]:>10.3f}  {eq_davies[sp]:>10.3f}  "
          f"{eq_davies[sp]-eq_ideal[sp]:>+10.3f}")
```

## Apparent pKw and the shift of neutral pH

Activity corrections shift more than the apparent pKa values of buffer acids — they shift neutral pH itself.
At the neutral point, the proton and hydroxide activities are equal: $a_{\ce{H+}} = a_{\ce{OH-}}$.
Substituting into the $K_w$ constraint with $a_{\ce{H2O}} = 1$:

$$
K_w = a_{\ce{H+}}\,a_{\ce{OH-}} = a_{\ce{H+}}^2 = 10^{-14}
\quad\Rightarrow\quad
a_{\ce{H+}}^\text{neut} = 10^{-7}.
$$

Writing $a_{\ce{H+}} = \gamma_\pm\,c_{\ce{H+}}/c^\circ$, the neutral concentration is $c_{\ce{H+}}^\text{neut} = c^\circ \cdot 10^{-7} / \gamma_\pm$, so:

$$
\mathrm{pH}_\text{neut}
= -\log_{10}\!\frac{c_{\ce{H+}}^\text{neut}}{c^\circ}
= 7 + \log_{10}\gamma_\pm.
$$

Since $\gamma_\pm < 1$ at nonzero ionic strength, $\mathrm{pH}_\text{neut} < 7$.
The Davies equation gives the shift directly (@fig-pKw-shift):

```{code-cell} ipython3
:tags: [remove-cell]

A_davies = 0.509   # Davies constant at 25 °C, I in mol/L

def log10_gamma_davies(I_mol_m3, z=1):
    I = I_mol_m3 / 1000.0
    sI = np.sqrt(I)
    return -A_davies * z**2 * (sI / (1 + sI) - 0.3 * I)

def log10_gamma_dh(I_mol_m3, z=1):
    I = I_mol_m3 / 1000.0
    sI = np.sqrt(I)
    return -A_davies * z**2 * sI / (1 + sI)

I_range = np.linspace(0, 500, 400)
pH_dav = 7.0 + np.array([log10_gamma_davies(I) for I in I_range])
pH_dh  = 7.0 + np.array([log10_gamma_dh(I)     for I in I_range])
```

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-pKw-shift

import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(6, 3.5))
ax.plot(I_range, pH_dav, color="C0", label="Davies")
ax.plot(I_range, pH_dh,  color="C1", label="Debye–Hückel", ls="--")
ax.axhline(7.0, color="gray", lw=0.8, ls=":", label="ideal ($\\gamma = 1$)")
ax.axvline(150,  color="gray", lw=0.6, ls=":")
pH_150 = 7.0 + log10_gamma_davies(150)
ax.annotate(
    f"pH = {pH_150:.2f}  at  $I = 150$",
    xy=(150, pH_150),
    xytext=(210, pH_150 + 0.025),
    fontsize=8,
    arrowprops=dict(arrowstyle="->", lw=0.8),
)
ax.set_xlabel(r"ionic strength $I$ [mol/m³]")
ax.set_ylabel(r"$\mathrm{pH}_\mathrm{neut}$")
ax.legend(fontsize=8)
fig.tight_layout()
```

```{figure} #cell-pKw-shift
:name: fig-pKw-shift

Neutral pH as a function of ionic strength from the Davies (solid) and Debye--Hückel (dashed) equations.
The dotted lines mark physiological ionic strength ($I = 150\ \mathrm{mol/m^3}$) and ideal neutral pH (7.00).
Davies predicts a shift of 0.12 pH units at physiological conditions; Debye--Hückel overestimates the shift at high $I$ because it lacks the $-0.3I$ empirical correction.
```

At physiological ionic strength ($I = 150\ \mathrm{mol/m^3}$), Davies predicts neutral pH $\approx 6.88$ — a 0.12 unit shift.
A solution prepared to pH 7.00 in high-salt media is slightly acidic relative to the true equal-activity midpoint.

This effect requires two ingredients simultaneously.
The $K_w$ constraint gives neutral pH its algebraic content: without $\ce{OH-}$ as a tracked species the condition $a_{\ce{H+}} = a_{\ce{OH-}}$ has nothing to enforce it.
The Davies correction makes $\gamma_\pm \neq 1$: without it $\log_{10}\gamma_\pm = 0$ and neutral pH is always 7.0 regardless of ionic strength.
Neither alone produces the shift; both are necessary.

In CADET's original implementation, neither $\ce{OH-}$ nor ionic strength corrections were included — a consistent approximation at a coarser level.
Adding Davies corrections to buffer equilibria without also tracking $\ce{OH-}$ and enforcing $K_w$ would have been the problematic intermediate state: ionic strength effects applied selectively, with no mechanism to propagate the shift to the proton--hydroxide balance.

---

The next chapter combines multiple equilibria to compute pH curves and buffer
capacity, including ionic strength coupling (@implementation-buffer).
