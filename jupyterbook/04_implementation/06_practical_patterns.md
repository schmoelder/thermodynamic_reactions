---
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

(implementation-practical)=
# Practical Patterns

This chapter collects self-contained idioms and pre-built components that reduce boilerplate in common modelling tasks.
Each section covers one tool or one pattern with a minimal example; the chapter is a reference rather than a narrative.
The chapter breaks the conceptual rhythm of the part intentionally; readers may navigate directly to the section they need without reading sequentially.


## Verifying conservation with `check_conservation`

`ReactionModel.check_conservation()` verifies that the sum of species within each multi-species component is a conserved moiety of the reaction network — that is, the total amount of each component is unchanged by the reactions.
The check operates on the left null space of the stoichiometric matrix (@implementation-source-term) and works for kinetic, equilibrium, and mixed models alike.

```{code-cell} ipython3
import numpy as np

from reactions.species import Component, Species
from reactions.model import ReactionModel
from reactions.reaction import ThermodynamicReaction
from reactions.equilibrium import pKa

acetic_acid_comp = Component(
    "acetic_acid",
    [Species("HAc", charge=0), Species("Ac-", charge=-1)],
)
H_plus_comp = Component("H+", [Species("H+", charge=1)])

model_check = ReactionModel(
    components=[acetic_acid_comp, H_plus_comp],
    reactions=[
        ThermodynamicReaction(
            "HAc <-> Ac- + H+",
            mode="equil",
            equilibrium_constant=pKa(4.756),
        ),
    ],
)

for report in model_check.check_conservation():
    print(
        f"{report.component.name}: conserved={report.conserved}, residual={report.residual:.2e}"
    )
```

`acetic_acid` is conserved because the reaction moves protons between species within the component without creating or destroying the acetate moiety.
`H+` is a single-species component and is not checked.
Run `check_conservation()` after constructing any model with polyprotic or multi-form components to detect missing reactions or stoichiometric errors before running a solve.


## Mixed equil/kinetic models

A `ReactionModel` may contain reactions in both `mode="equil"` and `mode="kinetic"` simultaneously.
Equilibrium reactions are enforced as algebraic constraints at each time step; kinetic reactions drive the dynamics.
This is the natural pattern when some steps are fast (proton transfer, fast isomerisation) and others are slow (enzymatic conversion, slow conformational change).

```{code-cell} ipython3
from reactions.equilibrium import EquilibriumConstant
from reactions.rate import RateConstantFixed
from reactions.solver import simulate

model_mixed = ReactionModel(
    components=[Component("A"), Component("B"), Component("C")],
    reactions=[
        ThermodynamicReaction(
            "A <-> B",
            mode="kinetic",
            equilibrium_constant=EquilibriumConstant(K_eq=2.0),
            rate_constant=RateConstantFixed(kf_value=0.05),
        ),
        ThermodynamicReaction(
            "B <-> C",
            mode="equil",
            equilibrium_constant=EquilibriumConstant(K_eq=3.0),
        ),
    ],
)

result_mixed = simulate(model_mixed, c0={"A": 1.0}, t_span=(0, 100.0))
print(
    f"c_B / c_C  = {result_mixed['B'][-1] / result_mixed['C'][-1]:.4f}  (expect {1 / 3.0:.4f})"
)
```

The `B⇌C` constraint is satisfied at every time step: $c_B/c_C = 1/K_{BC} = 1/3$ is maintained throughout.
`A⇌B` drives the approach to the overall equilibrium on the kinetic timescale.


## Initial guess strategy for `solve_equilibrium`

`solve_equilibrium` uses Newton's method in log-space.
Convergence is guaranteed only if the initial guess is close enough to the solution; a guess that places a concentration near zero or several orders of magnitude from equilibrium may fail or converge slowly.

For acid-base systems, the Henderson-Hasselbalch equation (@acid-base) provides a reliable starting guess.
At a target pH, the conjugate-acid/conjugate-base ratio follows:

$$
\frac{[\text{A}^-]}{[\text{HA}]} = 10^{\,\mathrm{pH} - \mathrm{p}K_a}.
$$

```{code-cell} ipython3
# Good initial guess for an acetic acid system at pH ≈ 4.5
pKa_val = 4.756
c_total = 100.0  # mol/m³  total acetate
C_REF = 1000.0  # mol/m³
pH_guess = 4.5

ratio = 10.0 ** (pH_guess - pKa_val)
c_Ac_neg = c_total * ratio / (1.0 + ratio)
c_HAc = c_total - c_Ac_neg
c_H = 10.0 ** (-pH_guess) * C_REF
c_OH = 1e-14 * C_REF**2 / c_H

c0_good = {"HAc": c_HAc, "Ac-": c_Ac_neg, "H+": c_H, "OH-": c_OH}
print(f"HAc={c_HAc:.1f}, Ac-={c_Ac_neg:.1f}, H+={c_H:.3e}, OH-={c_OH:.3e}  mol/m³")
```

For simple A⇌B with known K, use $c_B^0 = K\,c_\text{tot}/(1+K)$.
When the system has multiple equilibria, solve them sequentially to bootstrap a consistent starting point.

```{admonition} Note
:class: note

`solve_equilibrium` operates in log-space ($\ln c$ as the solver variable), so it requires all initial concentrations to be strictly positive.
Any species with a physically expected near-zero concentration should be set to a small positive value (e.g. $10^{-12}\,c^\circ$) rather than zero.
```


## Pre-defined components: `common.py`

`reactions.common` provides `Component` instances for frequently used buffer species and the solvent.
These are standalone objects with no reaction dependency; import them directly without constructing new `Component` instances.
The acid components carry only `name` and `charge` on each species: no `molar_mass`, `density`, or `heat_capacity`.
This is intentional: for dilute aqueous equilibria, solute molar masses are not needed (volume is accounted for by the solvent), and thermal properties require application-specific values.
If a simulation needs physical properties on an acid species, construct a `Species` explicitly rather than patching the pre-defined component.

```{code-cell} ipython3
from reactions.common import (
    water,
    H_plus,
    OH_minus,
    acetic_acid,
    phosphate,
    citric_acid,
    tris,
    hepes,
    mops,
)

print(f"water c_ref  = {water.species[0].c_ref:.0f} mol/m³")
print(f"H+  charge   = {H_plus.species[0].charge:+d}")
print(f"OH- charge   = {OH_minus.species[0].charge:+d}")
print(f"phosphate species: {[s.name for s in phosphate.species]}")
```

The `water` component carries `molar_mass = 0.018015 kg/mol` and `density = 1000.0 kg/m³`, giving $c_\text{ref} = \rho/M \approx 55\,509\ \text{mol/m}^3$ for pure water.

```{admonition} Energy balance
:class: tip

`water` does not include `heat_capacity`.
To use the coupled energy balance (@implementation-energy-balance), attach it explicitly:

    from reactions.species import Species
    water_with_Cp = Component(
        "water",
        [Species("H2O", charge=0, molar_mass=0.018015, density=1000.0, heat_capacity=75.3)],
    )
```


## The `pKa()` factory

`pKa(value)` converts a thermodynamic pK$_a$ to an `EquilibriumConstant` such that $K(298.15\ \text{K}) = 10^{-\text{p}K_a}$ exactly.
An optional `dH` argument adds van't Hoff temperature dependence (@equilibrium-temperature) with the standard-state enthalpy $\Delta H^\circ$ [J/mol]:

```{code-cell} ipython3
from reactions.equilibrium import pKa

K_acetic = pKa(4.756)  # temperature-independent
K_water = pKa(14.00, dH=55800.0)  # van't Hoff: pKw shifts with T
K_tris = pKa(8.072, dH=-47450.0)  # ΔH < 0 → pKa drops with T

T_amb = 298.15
T_warm = 310.0

print(f"K_acetic(298.15 K) = {K_acetic.K(T_amb):.4e}  (expect {10**-4.756:.4e})")
print(f"pKw(298.15 K)      = {-np.log10(K_water.K(T_amb)):.4f}  (expect 14.0000)")
print(
    f"pKw(310.0 K)       = {-np.log10(K_water.K(T_warm)):.4f}  (shifted by van't Hoff)"
)
print(f"pKa_Tris(298.15 K) = {-np.log10(K_tris.K(T_amb)):.4f}")
print(f"pKa_Tris(310.0 K)  = {-np.log10(K_tris.K(T_warm)):.4f}  (−0.03/°C expected)")
```

The back-calculated $\Delta S^\circ$ ensures $K(T_\text{ref}) = 10^{-\text{p}K_a}$ exactly, regardless of the sign or magnitude of `dH`.


## Reaction factories

`reactions.common` provides one factory function per buffer system, each returning a list of `ThermodynamicReaction` instances.
The `*factory()` unpacking idiom passes the list directly into `reactions=[...]`:

```{code-cell} ipython3
from reactions.common import (
    autoionisation,
    acetic_acid_equilibria,
    phosphate_equilibria,
    citric_acid_equilibria,
    tris_equilibria,
    hepes_equilibria,
    mops_equilibria,
)
from reactions.activity import ActivityCoefficientDavies

davies = ActivityCoefficientDavies()

reactions = [
    *autoionisation(activity_coefficient=davies),  # pKw = 14.00, ΔH° = +55.8 kJ/mol
    *acetic_acid_equilibria(activity_coefficient=davies),  # pKa = 4.756
    *phosphate_equilibria(activity_coefficient=davies),  # pKa = 2.148 / 7.198 / 12.35
]
print(f"Reactions assembled: {len(reactions)}")
for rxn in reactions:
    print(f"  {rxn}")
```

All pK$_a$ values are at 298.15 K.
Factories accept `activity_coefficient=None` (ideal) or any `ActivityCoefficientBase` instance.
The temperature-sensitive systems are:

| Factory                    | pKa      | ΔH° [kJ/mol] | dpKa/dT [pH/°C]  |
|----------------------------|----------|--------------|------------------|
| `autoionisation()`         | 14.00    | +55.8        | −0.033           |
| `tris_equilibria()`        | 8.072    | −47.45       | +0.028           |
| All others                 | —        | —            | ≈ 0              |


## The `Solution` class

`Solution` converts a buffer recipe (solvent + solute concentrations) to the `c0` and `prescribed` dicts required by `solve_equilibrium` and `simulate`.
Solvent species are set to their pure-component concentration scaled by volume fraction; their entries are automatically included in both `c0` and `prescribed`.

```{code-cell} ipython3
from reactions.formulation import Solution

sol = Solution(
    water,
    solutes={"HAc": 100.0, "Ac-": 0.0, "H+": 1e-4, "OH-": 1e-10},
)

print("c0       :", sol.c0)
print("prescribed:", sol.prescribed)
```

For mixed solvents, pass a dict of `{Component: volume_fraction}` pairs (fractions must sum to 1).
The `prescribed` dict from `Solution` is passed directly to `solve_equilibrium(..., prescribed=sol.prescribed)`.


## pH-stat pattern

Prescribing $[\text{H}^+]$ at a fixed target fixes the pH throughout an equilibrium solve or simulation.
This is the pH-stat pattern: combine `prescribed={"H+": c_H}` with a sweep over a target pH array to compute speciation or rates as a function of pH.

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-phstat-speciation

import matplotlib.pyplot as plt
from reactions.common import (
    acetic_acid,
    H_plus,
    OH_minus,
    water,
    acetic_acid_equilibria,
    autoionisation,
)
from reactions.ionic import IonicStrengthIdeal
from reactions.plots import setup_figure
from reactions.solver import solve_equilibrium

model_stat = ReactionModel(
    components=[acetic_acid, H_plus, OH_minus, water],
    reactions=[
        *acetic_acid_equilibria(),
        *autoionisation(),
    ],
    ionic_strength=IonicStrengthIdeal(),
)

pH_targets = np.linspace(3.0, 9.0, 40)
f_Ac_neg = []
H_CREF = H_plus.species[0].c_ref
WATER_CREF = water.species[0].c_ref

for pH in pH_targets:
    c_H = 10.0 ** (-pH) * H_CREF
    c_OH = 1e-14 * H_CREF**2 / c_H
    c0 = {"HAc": 50.0, "Ac-": 50.0, "H+": c_H, "OH-": c_OH}
    c_eq = solve_equilibrium(
        model_stat,
        c0,
        prescribed={"H+": c_H, "H2O": WATER_CREF},
    )
    total = c_eq["HAc"] + c_eq["Ac-"]
    f_Ac_neg.append(c_eq["Ac-"] / total)

pKa_val = 4.756
fig, ax = setup_figure()
ax.plot(pH_targets, f_Ac_neg, color="C0", lw=2.0, label=r"$f_{\mathrm{Ac}^-}$")
ax.plot(pH_targets, [1 - f for f in f_Ac_neg], color="C1", lw=2.0, label=r"$f_{\mathrm{HAc}}$")
ax.axvline(pKa_val, color="gray", ls=":", lw=1.0)
ax.text(pKa_val + 0.1, 0.52, rf"$\mathrm{{p}}K_a = {pKa_val}$", fontsize=9, color="gray")
ax.set_xlabel("pH")
ax.set_ylabel("Speciation fraction")
ax.set_xlim(3, 9)
ax.set_ylim(-0.02, 1.05)
ax.legend(fontsize=10)
fig.tight_layout()
```

```{figure} #cell-phstat-speciation
:name: fig-phstat-speciation

Speciation of acetic acid computed by the pH-stat pattern: $[\text{H}^+]$ is prescribed at each target pH and the equilibrium solve determines the remaining concentrations.
The result is identical to the Bjerrum diagram (@fig-bjerrum-acetic), confirming that prescribing pH and solving for speciation is equivalent to the analytical Henderson-Hasselbalch expression.
```

The `prescribed` dict overrides the equilibrium solve for the listed species: their concentrations are held fixed at the specified value while all others are solved freely.
The pattern applies equally to `simulate(..., prescribed={"H+": c_H})` for time-dependent simulations at a controlled pH.

---

The following chapter applies these patterns to a complete acid-base speciation problem, where the `pKa()` factory, reaction factories, `Solution`, and `prescribed` are used together to model buffer equilibria from scratch (@implementation-acid-base).
