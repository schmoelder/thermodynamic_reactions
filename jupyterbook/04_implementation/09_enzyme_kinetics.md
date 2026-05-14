---
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

(implementation-enzyme)=
# Enzyme Kinetics

Mass action kinetics are unbounded in substrate concentration in the model form: the rate is linear or polynomial in concentration with no ceiling.
For enzyme-catalyzed reactions the rate saturates because the enzyme active site becomes limiting; under quasi-steady-state or rapid equilibrium assumptions, this produces Michaelis-Menten kinetics: the same finite-site saturation that produces the Langmuir isotherm at equilibrium (@adsorption, @saturation).

This chapter is an orthogonal extension: it depends on the kinetic framework introduced in @implementation-kinetics, but is independent of the acid-base and activity chemistry in @implementation-activity through @implementation-buffer.
`EnzymaticReaction` slots into the same $\mathbf{S}\boldsymbol{\varphi}$ framework as `ThermodynamicReaction`; only the rate closure $\varphi$ changes.
Three reaction classes cover the library: `ThermodynamicReaction` enforces thermodynamic consistency ($k^r = k^f/K$ at every evaluation), `MassActionReaction` applies a linear rate law with a free reverse constant, and `EnzymaticReaction` accepts a custom rate closure.
`EnzymaticReaction` does not enforce $k^r = k^f/K$: the rate closure is supplied by the user and may represent an empirical fit that is not thermodynamically reversible in the strict sense (@implementation-kinetics).
This is appropriate for enzyme catalysis, where the Michaelis-Menten and Hill forms describe net forward flux under physiological conditions and the reverse reaction is typically negligible or absorbed into the effective parameters.
The Michaelis-Menten form is not an arbitrary empirical fit: it emerges from a finite-site binding equilibrium (the enzyme-substrate complex) when that complex is not tracked as an explicit species — the same occupancy constraint that produces the Langmuir isotherm at equilibrium (@adsorption).
The coarse-graining is in the elimination of the intermediate, not in the physics.


## Saturation rate laws

`MichaelisMenten` and `HillRate` implement the two standard saturation forms:

$$
v_\text{MM} = \frac{V_\text{max}\,[S]}{K_m + [S]},
\qquad
v_\text{Hill} = \frac{V_\text{max}\,[S]^n}{K_m^n + [S]^n}.
$$

$V_\text{max} = k_\text{cat} E_\text{tot}$ is set by the catalytic rate constant and total enzyme concentration; $K_m$ is the substrate concentration at half-maximal rate.
At low substrate concentration both reduce to a linear (mass-action) rate $v \approx (V_\text{max}/K_m)\,[S]$; at saturation ($[S] \gg K_m$) the rate approaches $V_\text{max}$ (@fig-saturation).
The Hill coefficient $n > 1$ introduces effective cooperativity: the response is sigmoidal and steeper than the Michaelis-Menten curve.

```{code-cell} ipython3
import numpy as np
from reactions.api import (
    Component,
    EnzymaticReaction,
    MichaelisMenten,
    HillRate,
    ReactionModel,
)
from reactions.solver import simulate

Vmax = 1.0  # mol/(m³·s)
Km = 200.0  # mol/m³

model_mm = ReactionModel(
    components=[Component("S"), Component("P")],
    reactions=[
        EnzymaticReaction(
            "S -> P",
            rate=MichaelisMenten(Vmax=Vmax, Km=Km, substrate="S"),
        ),
    ],
)

model_hill = ReactionModel(
    components=[Component("S"), Component("P")],
    reactions=[
        EnzymaticReaction(
            "S -> P",
            rate=HillRate(Vmax=Vmax, Km=Km, n=3, substrate="S"),
        ),
    ],
)
```

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-saturation

import matplotlib.pyplot as plt

S_vals = np.linspace(0, 1000, 400)
v_mm = Vmax * S_vals / (Km + S_vals)
v_hill = Vmax * S_vals**3 / (Km**3 + S_vals**3)
v_lin = (Vmax / Km) * S_vals

fig, ax = plt.subplots()
ax.plot(S_vals, v_mm, label="Michaelis-Menten")
ax.plot(S_vals, v_hill, label=r"Hill ($n = 3$)")
ax.plot(
    S_vals, v_lin, ls=":", color="gray", label=r"linear ($k^f = V_\mathrm{max}/K_m$)"
)
ax.axhline(Vmax, color="gray", lw=0.8, ls="--", label=r"$V_\mathrm{max}$")
ax.axvline(Km, color="gray", lw=0.8, ls="-.", label=r"$K_m$")
ax.set_xlabel(r"$[S]\ [\mathrm{mol/m^3}]$")
ax.set_ylabel(r"$v\ [\mathrm{mol/(m^3 \cdot s)}]$")
ax.set_ylim(0, 1.5)
ax.legend()
fig.tight_layout()
```

```{figure} #cell-saturation
:name: fig-saturation

Michaelis-Menten and Hill ($n=3$) rate laws compared to the linear (mass-action) limit, for $V_\text{max} = 1\ \mathrm{mol/(m^3 \cdot s)}$ and $K_m = 200\ \mathrm{mol/m^3}$.
Both saturate at $V_\text{max}$ for $[S] \gg K_m$; the Hill curve is sigmoidal due to cooperativity.
The dotted line is the linear approximation valid at low substrate concentration.
```


## Kinetic simulation

The rate law $\varphi$ feeds into the same source term $\mathbf{S}\boldsymbol{\varphi}$ as all other reaction types.
For a single reaction, $v \equiv \varphi_j$ in CADET notation; the stoichiometric matrix distributes $\varphi$ across species.
For $\ce{S -> P}$ with $\nu_S = -1$, $\nu_P = +1$ this gives $\dot{c}_S = -\varphi$ and $\dot{c}_P = +\varphi$.

```{code-cell} ipython3
S0 = 600.0  # mol/m³  (initial substrate)

result_mm = simulate(model_mm, c0={"S": S0}, t_span=(0, 2000.0))
result_hill = simulate(model_hill, c0={"S": S0}, t_span=(0, 2000.0))
```

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-enzyme-sim

fig, axes = plt.subplots(1, 2, figsize=(9, 3.8), sharey=True)

for ax, result, label in [
    (axes[0], result_mm, "Michaelis-Menten"),
    (axes[1], result_hill, r"Hill ($n = 3$)"),
]:
    ax.plot(result.coords["time"], result["c"].sel(species="S"), label="S (substrate)", color="C0")
    ax.plot(result.coords["time"], result["c"].sel(species="P"), label="P (product)", color="C1")
    ax.set_xlabel("time [s]")
    ax.set_title(label)

axes[0].set_ylabel(r"concentration [mol/m³]")
axes[0].legend()
fig.tight_layout()
```

```{figure} #cell-enzyme-sim
:name: fig-enzyme-sim

Substrate depletion and product formation for Michaelis-Menten (left) and Hill $n=3$ (right)
kinetics with $V_\text{max} = 1\ \mathrm{mol/(m^3 \cdot s)}$, $K_m = 200\ \mathrm{mol/m}^3$,
$[S]_0 = 600\ \mathrm{mol/m}^3$.
The Hill model shows a sigmoidal initial phase: when $[S] \gg K_m$ effective cooperativity maintains a high rate near $V_\text{max}$; depletion accelerates once $[S]$ drops
toward $K_m$, then slows again in the linear regime below $K_m$.
```

At high substrate ($[S] \gg K_m$) the rate is near $V_\text{max}$ and the trajectory is nearly linear in time.
Below $K_m$ the rate drops into the linear regime and the approach to full conversion slows.
The Hill model reaches full conversion more abruptly because the cooperative response maintains a high rate until $[S]$ falls close to $K_m$, then drops steeply.

## Combined model: pH-dependent enzyme

The three-state protonation model (@enzyme-ph-activity) gives $f_\text{active}$ as the middle Bjerrum fraction of the enzyme active site.
The coupling between proton-transfer equilibria and saturating kinetics is natural in the `ReactionModel` framework: both reaction types occupy the same model, driven by the same $\mathbf{S}\boldsymbol{\varphi}$ source term, without any additional machinery.

```{code-cell} ipython3
from reactions.api import CustomRate
from reactions.common import H_plus, OH_minus, water, autoionization

C_REF = 1000.0  # mol/m³  # todo: let's check all uses of C_REF again. do we still *need* it?
pKa1, pKa2 = 5.5, 8.5
Ka1_d = 10.0 ** (-pKa1)
Ka2_d = 10.0 ** (-pKa2)
Vmax_max = 1.0  # mol/(m³·s)  peak rate
Km_enz = 200.0  # mol/m³


def ph_enzyme_rate(state, species_index):
    a_H = state.c[species_index["H+"]] / C_REF  # todo: do we need activity coefficient? should we have an "activity" pseudo state?
    S = state.c[species_index["S"]]
    denom = a_H**2 + Ka1_d * a_H + Ka1_d * Ka2_d
    return Vmax_max * (Ka1_d * a_H / denom) * S / (Km_enz + S)


model_ph = ReactionModel(
    components=[Component("S"), Component("P"), H_plus, OH_minus, water],
    reactions=[
        EnzymaticReaction("S -> P", rate=CustomRate(fn=ph_enzyme_rate)),
        *autoionization(),
    ],
)
```

The profile (@fig-ph-bell) shows that $f_\text{active}$ enters as a multiplicative prefactor on $V_\text{max}$.
Simulations at pH 4 (well below $\text{p}K_{a1}$, $f_\text{active} \approx 0.03$), pH 7 (optimal, $f_\text{active} \approx 1$), and pH 9 (above $\text{p}K_{a2}$, $f_\text{active} \approx 0.24$) sample three distinct regions of the bell.

Simulations at three prescribed pH values use the pH-stat pattern (@implementation-practical):

```{code-cell} ipython3
S0 = 600.0
t_span = (0, 3000.0)

results_ph = {}
for pH in [4.0, 7.0, 9.0]:
    c_H = 10.0 ** (-pH) * C_REF
    c_OH = 1e-14 * C_REF**2 / c_H
    results_ph[pH] = simulate(
        model_ph,
        c0={"S": S0, "H+": c_H, "OH-": c_OH},
        t_span=t_span,
        prescribed={"H+": c_H, "H2O": C_REF},
    )
```

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-ph-enzyme-sim

fig, ax = plt.subplots(figsize=(5, 3.5))
for pH, result in results_ph.items():
    ax.plot(result.coords["time"], result["c"].sel(species="S"), label=f"pH {pH:.0f}")
ax.set_xlabel("time [s]")
ax.set_ylabel(r"$c_S\ [\mathrm{mol/m^3}]$")
ax.legend()
fig.tight_layout()
```

```{figure} #cell-ph-enzyme-sim
:name: fig-ph-enzyme-sim

Substrate depletion at pH 4.0, 7.0, and 9.0 ($V_\text{max} = 1\ \mathrm{mol/(m^3 \cdot s)}$, $K_m = 200\ \mathrm{mol/m^3}$, $[S]_0 = 600\ \mathrm{mol/m^3}$).
At pH 7.0 (optimal) conversion is rapid; at pH 9.0 the active fraction is $\approx 0.24$ and the rate is substantially reduced; at pH 4.0 the enzyme is nearly inactive ($f_\text{active} \approx 0.03$) and the substrate barely depletes on this timescale.
```

At the optimal pH, $f_\text{active} \approx 1$ and the enzyme operates near $V_\text{max}$; the difference between pH 9.0 and pH 4.0 reflects the asymmetry of the bell: pH 4.0 is well below $\text{p}K_{a1} = 5.5$, suppressing activity far more than the equivalent distance above $\text{p}K_{a2} = 8.5$.
A `CustomRate` makes the rate closure available to the solver with no changes to the surrounding framework: the proton-transfer equilibrium (autoionization) and the enzymatic kinetics are co-resident in the same `ReactionModel` and are handled uniformly.

---

The reaction framework covers mass-action kinetics (@implementation-kinetics), thermodynamic equilibria (@implementation-equilibrium), and saturation kinetics — all as composable elements of a single `ReactionModel`.
The next chapter defines the integration contract that CADET-Core expects: the residual and Jacobian interface that links this prototype library to the production solver (@implementation-interface).
