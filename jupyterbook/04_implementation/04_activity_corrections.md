---
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

(implementation-activity)=
# Activity and Non-ideal Corrections

The previous chapter enforced $k_r = k_f / K(T)$ but still used ideal activities $a_i = c_i / c^\circ$ with $\gamma_i = 1$.
For charged species at finite ionic strength, activity coefficients deviate significantly from unity, shifting the apparent equilibrium.
The Debye-Hückel and Davies models for $\gamma_i(I, z_i)$ are established in @nonidealities; this chapter wires them into the reaction model.

The measurable consequence is a shift in the apparent pKa (@speciation-buffers):

$$
\text{pKa}^{\text{app}} = \text{pKa} + \sum_i \nu_i \log_{10} \gamma_i.
$$

For a monovalent acid HA $\rightleftharpoons$ A$^-$ + H$^+$, the ions A$^-$ and H$^+$ both have $\gamma_i < 1$ at finite $I$ while the neutral species has $\gamma_\text{HA} = 1$, so the apparent pKa is always lower than the thermodynamic value.
The size of this shift is the figure of merit for the activity correction.

## Ionic strength computation

The ionic strength $I = \frac{1}{2}\sum_i c_i z_i^2$ must be evaluated from the current composition before $\gamma_i$ can be computed.
The library provides three strategies depending on how much of the ionic environment is modelled explicitly:

```{code-cell} ipython3
import numpy as np
from reactions.api import (
    C_REF,
    Species, Component,
    IonicStrengthIdeal, IonicStrengthBackground, IonicStrengthFixed,
    ActivityCoefficientDebyeHuckel,
    ActivityCoefficientDavies,
    ActivityCoefficientCustom,
    EquilibriumConstant,
    ThermodynamicReaction, ReactionModel,
    pKa,
)
from reactions.solver import solve_equilibrium

proton    = Component("proton",    [Species("H+",  charge=+1)])
hydroxide = Component("hydroxide", [Species("OH-", charge=-1)])
water     = Component("water",     [Species("H2O", charge=0, is_solvent=True)])
chloride  = Component("chloride",  [Species("Cl-", charge=-1)])
sodium    = Component("sodium",    [Species("Na+", charge=+1)])
```

**`IonicStrengthIdeal`** computes $I$ from the dynamic species in the model.
For a 1:1 electrolyte such as NaCl this gives $I = c_\text{NaCl}$:

```{code-cell} ipython3
model_nacl = ReactionModel(
    components=[proton, hydroxide, water, chloride, sodium],
    reactions=[
        ThermodynamicReaction(
            "H2O <-> H+ + OH-",
            mode="equil",
            equilibrium_constant=EquilibriumConstant(1e-14),
        ),
    ],
    ionic_strength=IonicStrengthIdeal(),
    T=298.15,
)

c = np.array([1e-4, 1e-4, 150.0, 150.0])   # H+, OH-, Cl-, Na+ [mol/m³]
state = model_nacl.make_state(c)
print(f"I = {state.I / 1000:.4f} mol/L   (expected 0.1500 for 150 mM NaCl)")
```

**`IonicStrengthBackground`** adds a fixed background $I_\text{bg}$ to the species-computed value, representing an unmodelled salt such as the NaCl gradient in ion-exchange chromatography.

**`IonicStrengthFixed`** holds $I$ at a constant value regardless of species concentrations, useful for sensitivity studies where the ionic environment is externally controlled.

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-ionic-strength

import matplotlib.pyplot as plt

c_nacl = np.linspace(0, 300, 300)
I_bg   = 50.0   # mol/m³

def eval_I_ideal(c_salt):
    c_vec = np.array([1e-4, 1e-4, c_salt, c_salt])
    return model_nacl.make_state(c_vec).I

I_ideal      = np.array([eval_I_ideal(c) for c in c_nacl])
I_background = I_ideal + I_bg
I_fixed_val  = 150.0 * np.ones_like(c_nacl)

fig, ax = plt.subplots()
ax.plot(c_nacl / 10, I_ideal      / 10, label="IonicStrengthIdeal")
ax.plot(c_nacl / 10, I_background / 10, label=f"IonicStrengthBackground ($I_{{bg}}$ = {I_bg/10:.0f} mM)")
ax.plot(c_nacl / 10, I_fixed_val  / 10, label="IonicStrengthFixed (150 mM)", ls="--")
ax.set_xlabel(r"$c_\mathrm{NaCl}$ [mM]")
ax.set_ylabel(r"$I$ [mM]")
ax.legend()
fig.tight_layout()
```

```{figure} #cell-ionic-strength
:name: fig-ionic-strength

Ionic strength vs NaCl concentration for the three models.
`IonicStrengthIdeal` tracks the solution composition exactly.
`IonicStrengthBackground` shifts the curve upward by a fixed $I_\text{bg}$, representing an unmodelled background salt.
`IonicStrengthFixed` is constant by construction.
```

## Activity coefficients and the self-consistent loop

The Debye-Hückel and Davies models for $\gamma_i(I, z_i)$ are documented in @nonidealities (@fig-activity).
Both are passed as the `activity_coefficient` argument of `ThermodynamicReaction`;
the solver evaluates them at the current $I$ during every Newton step.

The pKa shift is not a one-time correction: equilibrium concentrations determine $I$, which determines $\gamma$, which shifts the apparent equilibrium, which changes the concentrations.
`solve_equilibrium` closes this loop when $Q = K$ is satisfied at concentrations consistent with their own ionic strength.

The suppression is measurable: for a monovalent weak acid the apparent pKa falls $\approx 0.25$--$0.45$ units at 150 mM, typical of IEX loading conditions (@fig-activity-pka):

```{code-cell} ipython3
acid = Component("acid", [Species("HA", charge=0), Species("A-", charge=-1)])

pKa_thermo = 7.0
model = ReactionModel(
    components=[acid, proton, hydroxide, water],
    reactions=[
        ThermodynamicReaction(
            "HA <-> A- + H+",
            mode="equil",
            equilibrium_constant=pKa(pKa_thermo),
            activity_coefficient=ActivityCoefficientDavies(),
        ),
        ThermodynamicReaction(
            "H2O <-> H+ + OH-",
            mode="equil",
            equilibrium_constant=EquilibriumConstant(1e-14),
            activity_coefficient=ActivityCoefficientDavies(),
        ),
    ],
    ionic_strength=IonicStrengthBackground(I_bg=150.0),
    T=298.15,
)

c0    = {"HA": 0.5, "A-": 0.5, "H+": 10**(-pKa_thermo) * C_REF, "OH-": 10**(-7) * C_REF}
c_eq  = solve_equilibrium(model, c0, T=298.15)
pH_eq = -np.log10(c_eq["H+"] / C_REF)
A_frac   = c_eq["A-"] / (c_eq["HA"] + c_eq["A-"])
pKa_app  = pH_eq - np.log10(A_frac / (1 - A_frac))
print(f"apparent pKa at 150 mM: {pKa_app:.3f}  (thermodynamic: {pKa_thermo})")
```

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-activity-pka

I_bg_range = np.array([0, 10, 25, 50, 100, 150, 200, 300, 400, 500])
c0_loop    = {"HA": 1e-3, "A-": 1e-3, "H+": 10**(-pKa_thermo) * C_REF, "OH-": 10**(-7) * C_REF}

def make_model(activity_coeff, I_bg):
    return ReactionModel(
        components=[acid, proton, hydroxide, water],
        reactions=[
            ThermodynamicReaction(
                "HA <-> A- + H+",
                mode="equil",
                equilibrium_constant=pKa(pKa_thermo),
                activity_coefficient=activity_coeff,
            ),
            ThermodynamicReaction(
                "H2O <-> H+ + OH-",
                mode="equil",
                equilibrium_constant=EquilibriumConstant(1e-14),
                activity_coefficient=activity_coeff,
            ),
        ],
        ionic_strength=IonicStrengthBackground(I_bg=I_bg),
        T=298.15,
    )

pKa_dh  = []
pKa_dav = []

for I_bg_mM in I_bg_range:
    I_bg = float(I_bg_mM)
    for results, coeff in [(pKa_dh, ActivityCoefficientDebyeHuckel()),
                           (pKa_dav, ActivityCoefficientDavies())]:
        model  = make_model(coeff, I_bg)
        c_eq   = solve_equilibrium(model, c0_loop, T=298.15)
        pH_eq  = -np.log10(c_eq["H+"] / C_REF)
        A_frac = c_eq["A-"] / (c_eq["HA"] + c_eq["A-"])
        results.append(pH_eq - np.log10(A_frac / (1 - A_frac)))

fig, ax = plt.subplots()
ax.plot(I_bg_range, pKa_dh,  label="Debye-Hückel", marker="o", ms=4)
ax.plot(I_bg_range, pKa_dav, label="Davies",        marker="s", ms=4)
ax.axhline(pKa_thermo, color="gray", lw=0.8, ls="--",
           label=f"thermodynamic pKa = {pKa_thermo}")
ax.axvline(100, color="gray", lw=0.6, ls=":")
ax.set_xlabel(r"$I$ [mM]")
ax.set_ylabel(r"apparent $\mathrm{p}K_a$")
ax.legend()
fig.tight_layout()
```

```{figure} #cell-activity-pka
:name: fig-activity-pka

Apparent pKa of a monovalent weak acid (thermodynamic pKa = 7, dilute solution) as a function of background ionic strength.
At zero added salt the curve starts at the thermodynamic value; as ionic strength increases both A$^-$ and H$^+$ are stabilised by the ionic atmosphere, lowering the apparent pKa.
At 150 mM (typical IEX loading) the suppression reaches $\approx 0.25$--$0.35$ units; ignoring this shift introduces the same error into any pH calculation that uses the thermodynamic pKa directly.
```

For activity models outside the Debye-Hückel family, `ActivityCoefficientCustom` accepts any callable `(state, charges) -> np.ndarray`, slotting into any `ThermodynamicReaction` without changing the rest of the framework.

---

With activities and ionic strength in place, the next chapter applies Arrhenius kinetics to show how $k_f(T)$ and $k_r(T) = k_f(T)/K(T)$ enforce thermodynamic consistency across temperature (@implementation-kinetics).
