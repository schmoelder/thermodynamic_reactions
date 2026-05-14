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
The sources of non-ideality vary: long-range electrostatic interactions between ions (dominant in aqueous electrolytes at finite ionic strength), excluded-volume effects at high concentration, and specific ion interactions.
The choice of activity model determines which sources are captured.
Activity coefficients and the Debye-Hückel/Davies models are derived in @nonidealities; this chapter wires them into the reaction model.


## General activity corrections

`ThermodynamicReaction` accepts an `activity_coefficient` argument.
The contract is a callable `(state, charges) -> np.ndarray` that returns $\gamma_i$ for each species given the current physical state and charge numbers.
`ActivityCoefficientCustom` wraps any such callable, making any activity model available without modifying the framework.

A constant activity coefficient already illustrates the effect.
For $\ce{A <=> B}$ with ideal $K = 4$, setting $\gamma_\text{B} = 0.7$ while $\gamma_\text{A} = 1$ reduces $a_\text{B}$ relative to $c_\text{B}$, so the solver must increase $c_\text{B}$ to satisfy $Q = K$.
The apparent equilibrium constant in concentration space becomes $K_c = K\,\gamma_\text{A}/\gamma_\text{B} = 4/0.7 \approx 5.71$:

```{code-cell} ipython3
import numpy as np
from reactions.api import (
    Species,
    Component,
    IonicStrengthIdeal,
    IonicStrengthBackground,
    IonicStrengthFixed,
    ActivityCoefficientDebyeHuckel,
    ActivityCoefficientDavies,
    ActivityCoefficientCustom,
    PhysicalState,
    EquilibriumConstant,
    RateConstantFixed,
    ThermodynamicReaction,
    ReactionModel,
    pKa,
)
from reactions.solver import solve_equilibrium, simulate

C_REF = 1000.0  # mol/m³

K_thermo = 4.0
gamma_B = 0.7
c_tot = 1000.0  # mol/m³

comp_ab = Component("ab", [Species("A"), Species("B")])

model_ideal = ReactionModel(
    components=[comp_ab],
    reactions=[
        ThermodynamicReaction(
            "A <-> B",
            mode="equil",
            equilibrium_constant=EquilibriumConstant(K_thermo),
        ),
    ],
)

model_custom = ReactionModel(
    components=[comp_ab],
    reactions=[
        ThermodynamicReaction(
            "A <-> B",
            mode="equil",
            equilibrium_constant=EquilibriumConstant(K_thermo),
            activity_coefficient=ActivityCoefficientCustom(
                lambda state, charges: np.array([1.0, gamma_B])
            ),
        ),
    ],
)

c_ideal = solve_equilibrium(model_ideal, c0={"A": c_tot / 2, "B": c_tot / 2})
c_custom = solve_equilibrium(model_custom, c0={"A": c_tot / 2, "B": c_tot / 2})
```

```{code-cell} ipython3
:tags: [remove-cell]

K_c_expected = K_thermo / gamma_B  # apparent K in concentration space
print(f"Ideal:   c_B/c_A = {c_ideal['c'].sel(species='B').item() / c_ideal['c'].sel(species='A').item():.4f}  (K = {K_thermo})")
print(
    f"Custom:  c_B/c_A = {c_custom['c'].sel(species='B').item() / c_custom['c'].sel(species='A').item():.4f}  (K_c = K/γ_B = {K_c_expected:.4f})"
)
```

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-custom-gamma

import matplotlib.pyplot as plt
from reactions.plots import setup_figure, COLORS

gamma_B_range = np.linspace(0.3, 1.5, 300)
K_c_range = K_thermo / gamma_B_range
fB_range = K_c_range / (1 + K_c_range)
fB_ideal = K_thermo / (1 + K_thermo)

fig, ax = setup_figure()
ax.plot(
    gamma_B_range,
    fB_range,
    color="C0",
    label=r"non-ideal ($\gamma_\mathrm{A}=1$, $K=4$)",
)
ax.axhline(
    fB_ideal,
    color="gray",
    lw=0.8,
    ls="--",
    label=f"ideal ($\\gamma_i = 1$,  $c_B/c_{{tot}} = {fB_ideal:.2f}$)",
)
ax.axvline(1.0, color="gray", lw=0.6, ls=":")
ax.plot(
    gamma_B,
    K_thermo / gamma_B / (1 + K_thermo / gamma_B),
    "o",
    color="C0",
    ms=7,
    zorder=5,
    label=rf"example $\gamma_\mathrm{{B}} = {gamma_B}$",
)
ax.set_xlabel(r"$\gamma_\mathrm{B}$")
ax.set_ylabel(r"$c_\mathrm{B}^\mathrm{eq}\,/\,c_\mathrm{tot}$")
ax.legend()
fig.tight_layout()
```

```{figure} #cell-custom-gamma
:name: fig-custom-gamma

Equilibrium fraction of B as a function of $\gamma_\text{B}$ for $\ce{A <=> B}$ with $K = 4$ and $\gamma_\text{A} = 1$.
The dashed line is the ideal result ($\gamma_i = 1$); $\gamma_\text{B} < 1$ stabilises B and shifts the equilibrium toward the product, while $\gamma_\text{B} > 1$ destabilises it.
The filled circle marks the example value $\gamma_\text{B} = 0.7$.
```

The ratio $c_\text{B}/c_\text{A}$ shifts from $K = 4$ to $K/\gamma_\text{B} \approx 5.71$, confirming that $\gamma_i \ne 1$ is absorbed into the equilibrium composition rather than the equilibrium constant.
The same shift appears in kinetic simulations: both systems relax toward equilibrium at the same rate, but converge to different long-time compositions (@fig-custom-gamma-kinetics):

```{code-cell} ipython3
model_ideal_kin = ReactionModel(
    components=[comp_ab],
    reactions=[
        ThermodynamicReaction(
            "A <-> B",
            mode="kinetic",
            equilibrium_constant=EquilibriumConstant(K_thermo),
            rate_constant=RateConstantFixed(2000.0),
        ),
    ],
)

model_custom_kin = ReactionModel(
    components=[comp_ab],
    reactions=[
        ThermodynamicReaction(
            "A <-> B",
            mode="kinetic",
            equilibrium_constant=EquilibriumConstant(K_thermo),
            rate_constant=RateConstantFixed(2000.0),
            activity_coefficient=ActivityCoefficientCustom(
                lambda state, charges: np.array([1.0, gamma_B])
            ),
        ),
    ],
)

result_ideal = simulate(model_ideal_kin, c0={"A": c_tot, "B": 0.0}, t_span=(0, 2.0))
result_custom = simulate(model_custom_kin, c0={"A": c_tot, "B": 0.0}, t_span=(0, 2.0))
```

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-custom-gamma-kinetics

fB_ideal_eq = K_thermo / (1 + K_thermo)
fB_custom_eq = (K_thermo / gamma_B) / (1 + K_thermo / gamma_B)

fig, ax = setup_figure()
ax.plot(result_ideal.coords["time"], result_ideal["c"].sel(species="B"), color="C0", label="ideal ($\\gamma_i = 1$)")
ax.plot(
    result_custom.coords["time"],
    result_custom["c"].sel(species="B"),
    color="C1",
    label=rf"non-ideal ($\gamma_\mathrm{{B}} = {gamma_B}$)",
)
ax.axhline(fB_ideal_eq * c_tot, color="C0", lw=0.8, ls="--")
ax.axhline(fB_custom_eq * c_tot, color="C1", lw=0.8, ls="--")
ax.set_xlabel("time [s]")
ax.set_ylabel(r"$c_\mathrm{B}$ [mol/m³]")
ax.legend()
fig.tight_layout()
```

```{figure} #cell-custom-gamma-kinetics
:name: fig-custom-gamma-kinetics

Concentration of B over time for ideal and non-ideal ($\gamma_\text{B} = 0.7$, $\gamma_\text{A} = 1$) cases with $K = 4$.
Both trajectories start from pure A and follow the same initial slope; the non-ideal system converges to a higher equilibrium concentration (dashed lines) because $\gamma_\text{B} < 1$ stabilises B relative to A.
```

Any callable that maps `(state, charges)` to an array of activity coefficients slots in identically; the solver loop is unchanged.


## Activity models for aqueous electrolytes

For aqueous electrolyte solutions the two built-in models cover the relevant range.
`ActivityCoefficientDebyeHuckel` implements the Debye-Hückel limiting law, accurate to $I \approx 100\ \mathrm{mol/m}^3$; `ActivityCoefficientDavies` adds an empirical linear correction that extends validity to $\approx 500\ \mathrm{mol/m}^3$ (@nonidealities, @fig-activity):

```{code-cell} ipython3
proton = Component("proton", [Species("H+", charge=+1)])
hydroxide = Component("hydroxide", [Species("OH-", charge=-1)])
water = Component(
    "water", [Species("H2O", charge=0, molar_mass=0.018015, density=1000.0)]
)
chloride = Component("chloride", [Species("Cl-", charge=-1)])
sodium = Component("sodium", [Species("Na+", charge=+1)])
```

Species names are labels only; the `charge` argument is authoritative for all electrostatic calculations.
Both scale as $z_i^2$: a divalent ion at the same ionic strength receives a correction four times larger than a monovalent one.


## Ionic strength

For the Debye-Hückel and Davies models, $\gamma_i$ depends on the ionic strength $I = \frac{1}{2}\sum_i z_i^2 c_i$ [mol/m³], which must be evaluated from the current composition before $\gamma_i$ can be computed.
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

c = np.array([1e-4, 1e-4, 1000.0, 150.0, 150.0])  # H+, OH-, H2O, Cl-, Na+ [mol/m³]
state = model_nacl.make_state(c)
print(f"I = {state.I / 1000:.4f} mol/L   (expected 0.1500 for 150 mM NaCl)")
```

**`IonicStrengthBackground`** adds a fixed $I_\text{bg}$ to the species-computed value, representing an unmodelled salt such as the NaCl gradient in ion-exchange chromatography.

**`IonicStrengthFixed`** holds $I$ at a constant value regardless of species concentrations, useful for sensitivity studies where the ionic environment is externally controlled.

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-ionic-strength

import matplotlib.pyplot as plt
from reactions.plots import setup_figure, COLORS

c_nacl = np.linspace(0, 300, 300)
I_bg = 50.0  # mol/m³


def eval_I_ideal(c_salt):
    c_vec = np.array([1e-4, 1e-4, 1000.0, c_salt, c_salt])
    return model_nacl.make_state(c_vec).I


I_ideal = np.array([eval_I_ideal(c) for c in c_nacl])
I_background = I_ideal + I_bg
I_fixed_val = 150.0 * np.ones_like(c_nacl)

fig, ax = setup_figure()
ax.plot(c_nacl / 10, I_ideal / 10, label="IonicStrengthIdeal")
ax.plot(
    c_nacl / 10,
    I_background / 10,
    label=f"IonicStrengthBackground ($I_{{bg}}$ = {I_bg / 10:.0f} mM)",
)
ax.plot(c_nacl / 10, I_fixed_val / 10, label="IonicStrengthFixed (150 mM)", ls="--")
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
dav = ActivityCoefficientDavies()
c_H = 1e-4  # mol/m³  (= 10⁻⁷ mol/L, nominal pH 7)
I = 150.0  # mol/m³  (= 150 mM)

state = PhysicalState(c=np.array([c_H]), T=298.15, I=I)
gamma_H = dav.activity(state, np.array([1.0]))[0]

print(f"γ(H+) at 150 mM    = {gamma_H:.4f}")
print(f"pH from c(H+)      = {-np.log10(c_H / C_REF):.4f}")
print(f"pH from a(H+)      = {-np.log10(gamma_H * c_H / C_REF):.4f}")
```


## Apparent pKa

The same activity correction shifts the observable pKa.
From the general result (@speciation-buffers), $\text{pKa}^\text{app} = \text{pKa} + \sum_i \nu_i \log_{10} \gamma_i$, and for a neutral acid ($z = 0$) at physiological ionic strength the shift is approximately $-0.24$.
Since $\gamma_{\ce{HA}} \approx 1$ while $\gamma_{\ce{A-}}, \gamma_{\ce{H+}} < 1$ at finite $I$, the sum is negative and $\text{pKa}^\text{app} < \text{pKa}$ at any non-zero ionic strength.
At 150 mM, typical of IEX loading conditions, the suppression reaches $\approx 0.25-0.35$ units (@fig-activity-pka):

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

c0 = {
    "HA": 0.5,
    "A-": 0.5,
    "H+": 10 ** (-pKa_thermo) * C_REF,
    "OH-": 10 ** (-7) * C_REF,
}
c_eq = solve_equilibrium(model, c0, T=298.15, prescribed={"H2O": water.c_ref})
pH_eq = -np.log10(c_eq["c"].sel(species="H+").item() / C_REF)
A_frac = c_eq["c"].sel(species="A-").item() / (c_eq["c"].sel(species="HA").item() + c_eq["c"].sel(species="A-").item())
pKa_app = pH_eq - np.log10(A_frac / (1 - A_frac))
print(f"apparent pKa at 150 mM: {pKa_app:.3f}  (thermodynamic: {pKa_thermo})")
```

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-activity-pka

I_bg_range = np.array([0, 10, 25, 50, 100, 150, 200, 300, 400, 500])
c0_loop = {
    "HA": 1e-3,
    "A-": 1e-3,
    "H+": 10 ** (-pKa_thermo) * C_REF,
    "OH-": 10 ** (-7) * C_REF,
}


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


pKa_dh = []
pKa_dav = []

for I_bg_mM in I_bg_range:
    I_bg = float(I_bg_mM)
    for results, coeff in [
        (pKa_dh, ActivityCoefficientDebyeHuckel()),
        (pKa_dav, ActivityCoefficientDavies()),
    ]:
        model_i = make_model(coeff, I_bg)
        c_eq = solve_equilibrium(
            model_i, c0_loop, T=298.15, prescribed={"H2O": water.c_ref}
        )
        pH_eq = -np.log10(c_eq["c"].sel(species="H+").item() / C_REF)
        A_frac = c_eq["c"].sel(species="A-").item() / (c_eq["c"].sel(species="HA").item() + c_eq["c"].sel(species="A-").item())
        results.append(pH_eq - np.log10(A_frac / (1 - A_frac)))

fig, ax = setup_figure()
ax.plot(I_bg_range, pKa_dh, label="Debye-Hückel", marker="o", ms=4)
ax.plot(I_bg_range, pKa_dav, label="Davies", marker="s", ms=4)
ax.axhline(
    pKa_thermo, color="gray", lw=0.8, ls="--", label=f"thermodynamic pKa = {pKa_thermo}"
)
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


## Temperature dependence of $A$

The Debye-Hückel parameter $A$ is not a constant: it depends on the dielectric constant $\varepsilon_r$ of the solvent and temperature through (@nonidealities)

$$
A \propto (\varepsilon_r\, T)^{-3/2}.
$$

For water, $\varepsilon_r$ decreases from 88 at 0 °C to 56 at 100 °C, so $A$ grows above ambient temperature: ions interact more strongly because thermal motion weakens the dielectric screening.
Passing `epsilon_r=None` (the default) uses the stored 25 °C value $A = 0.509$ at all temperatures and emits a `UserWarning` when $|T - 298.15| > 5\,\mathrm{K}$.

The library provides `_water_epsilon_r`, a {cite:t}`malmberg1956` correlation valid from 0 to 100 °C, as a ready-made callable:

```{code-cell} ipython3
import warnings
from reactions.api import _water_epsilon_r

charges = np.array([1.0, -1.0])
I = 150.0  # mol/m³

dav_fixed = ActivityCoefficientDavies()
dav_Tdep = ActivityCoefficientDavies(epsilon_r=_water_epsilon_r)

for T in [298.15, 310.0, 330.0]:
    state = PhysicalState(c=np.array([1e-4, 1e-4]), T=T, I=I)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        g_fixed = dav_fixed.activity(state, charges)[0]
    g_Tdep = dav_Tdep.activity(state, charges)[0]
    print(f"T = {T:.2f} K:  γ (fixed A) = {g_fixed:.4f},  γ (T-dep A) = {g_Tdep:.4f}")
```

At 298 K the two agree by construction; above 298 K the temperature-dependent path gives lower $\gamma$ (stronger non-ideality) as $A$ grows.
A scalar `epsilon_r` value overrides the 25 °C constant with a fixed dielectric environment, for example `epsilon_r=33.0` for methanol at 25 °C.

---

With activities and ionic strength in place, the next chapter collects the pre-built components, reaction factories, and convenience patterns that the acid-base and buffer chapters use directly (@implementation-practical).
