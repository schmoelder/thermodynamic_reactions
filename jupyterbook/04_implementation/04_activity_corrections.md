---
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

(implementation-activity)=
# Activity and Non-ideal Corrections

The previous chapters used ideal activities $a_i = c_i / c^\circ$ with $\gamma_i = 1$.
Replacing concentrations with activities $a_i = \gamma_i c_i / c^\circ$ in every occurrence of $Q$ and $K$ is the full non-ideal correction; all non-ideality enters through $\gamma_i$.
The sources of non-ideality vary: long-range electrostatic interactions between ions (dominant in aqueous electrolytes at finite ionic strength), excluded-volume effects at high concentration, and specific ion interactions such as salting effects.
The choice of activity model determines which sources are captured.
Activity coefficients for real solutions and the Debye-Hückel / Davies models are established in @nonidealities; this chapter wires them into the reaction model.

## Activity models

`ThermodynamicReaction` accepts an `activity_coefficient` argument.
`ActivityCoefficientCustom` accepts any callable `(state, charges) -> np.ndarray`, allowing any activity model without modifying the framework:

```{code-cell} ipython3
import numpy as np
from reactions.api import (
    C_REF,
    Species, Component,
    IonicStrengthIdeal, IonicStrengthBackground, IonicStrengthFixed,
    ActivityCoefficientDebyeHuckel,
    ActivityCoefficientDavies,
    ActivityCoefficientCustom,
    PhysicalState,
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

For aqueous electrolyte solutions the two built-in models cover the relevant range: `ActivityCoefficientDebyeHuckel` is accurate to $I \approx 100\ \text{mol/m}^3$; `ActivityCoefficientDavies` extends validity to $\approx 500\ \text{mol/m}^3$ (@fig-activity).


## Ionic strength

For the Debye-Hückel and Davies models, $\gamma_i$ depends on the ionic strength $I = \frac{1}{2}\sum_i c_i z_i^2 / c^\circ$, which must be evaluated from the current composition before $\gamma_i$ can be computed.
Three strategies cover different modelling choices:

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

**`IonicStrengthBackground`** adds a fixed $I_\text{bg}$ to the species-computed value, representing an unmodelled salt such as the NaCl gradient in ion-exchange chromatography.

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


## The self-consistent loop

Activity corrections introduce a coupling: the equilibrium composition determines the activity coefficients, which shift the apparent equilibrium, which changes the composition.
`solve_equilibrium` closes this loop: Newton iteration continues until $Q = K$ is satisfied at concentrations consistent with their own $\gamma_i$.


## pH and activity

The thermodynamic definition of pH uses the activity of $\ce{H+}$, not its molar concentration:

$$
\text{pH} = -\log_{10} a_{\ce{H+}} = -\log_{10}\!\left(\gamma_{\ce{H+}}\,\frac{c_{\ce{H+}}}{c^\circ}\right).
$$

At finite ionic strength $\gamma_{\ce{H+}} < 1$, so pH exceeds $-\log_{10}(c_{\ce{H+}}/c^\circ)$.
This is not a correction applied on top of pH; it is the definition that pH measurements (glass electrode, calibrated against standard buffers) actually track.

```{code-cell} ipython3
dav  = ActivityCoefficientDavies()
c_H  = 1e-4     # mol/m³  (= 10⁻⁷ mol/L, nominal pH 7)
I    = 150.0    # mol/m³  (= 150 mM)

state   = PhysicalState(c=np.array([c_H]), T=298.15, I=I)
gamma_H = dav.activity(state, np.array([1.0]))[0]

print(f"γ(H+) at 150 mM    = {gamma_H:.4f}")
print(f"pH from c(H+)      = {-np.log10(c_H / C_REF):.4f}")
print(f"pH from a(H+)      = {-np.log10(gamma_H * c_H / C_REF):.4f}")
```


## Apparent pKa

The same activity correction shifts the observable pKa.
For $\ce{HA <=> A- + H+}$, substituting activities into $Q = K$ and taking logarithms gives (@speciation-buffers):

$$
\text{pKa}^\text{app} = \text{pKa} + \sum_i \nu_i \log_{10} \gamma_i.
$$

The neutral species $\ce{HA}$ has $\gamma_\text{HA} \approx 1$, while $\ce{A-}$ and $\ce{H+}$ both have $\gamma_i < 1$ at finite $I$, so $\text{pKa}^\text{app} < \text{pKa}$ at any non-zero ionic strength.
At 150 mM, typical of IEX loading conditions, the suppression reaches $\approx 0.25$--$0.35$ units (@fig-activity-pka):

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

c0       = {"HA": 0.5, "A-": 0.5, "H+": 10**(-pKa_thermo) * C_REF, "OH-": 10**(-7) * C_REF}
c_eq     = solve_equilibrium(model, c0, T=298.15)
pH_eq    = -np.log10(c_eq["H+"] / C_REF)
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

Apparent pKa of a monovalent weak acid (thermodynamic pKa = 7) as a function of background ionic strength.
At zero added salt the curve starts at the thermodynamic value; as $I$ increases, $\ce{A-}$ and $\ce{H+}$ are stabilised by the ionic atmosphere, lowering the apparent pKa.
At 150 mM the suppression reaches $\approx 0.25$--$0.35$ units; ignoring this introduces the same error into any pH calculation that uses the thermodynamic pKa directly.
```

---

With activities and ionic strength in place, the next chapter applies Arrhenius kinetics to show how $k_f(T)$ and $k_r(T) = k_f(T)/K(T)$ enforce thermodynamic consistency across temperature (@implementation-kinetics).
