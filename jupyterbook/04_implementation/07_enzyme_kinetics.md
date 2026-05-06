---
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

(implementation-enzyme)=
# Enzyme Kinetics

Mass action kinetics are unbounded in substrate concentration in the model form: the rate is linear or polynomial in concentration with no ceiling.
For enzyme-catalysed reactions the rate saturates because the enzyme active site becomes limiting; under quasi-steady-state or rapid equilibrium assumptions, this produces Michaelis-Menten kinetics: the same finite-site saturation that produces the Langmuir isotherm at equilibrium (@adsorption, @saturation).

This chapter is an orthogonal extension: it depends on the kinetic framework introduced in @implementation-kinetics, but is independent of the acid-base and activity chemistry in @implementation-activity through @implementation-buffer.
`EnzymaticReaction` slots into the same $\mathbf{S}\boldsymbol{\varphi}$ framework as `ThermodynamicReaction`; only the rate closure $\varphi$ changes.
This is a third closure class: `ThermodynamicReaction` closes via an equilibrium constraint, `MassActionReaction` via a linear rate law, and `EnzymaticReaction` via a saturating nonlinear closure.


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

Vmax = 1.0    # mol/(m³·s)
Km   = 200.0  # mol/m³

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
v_mm   = Vmax * S_vals / (Km + S_vals)
v_hill = Vmax * S_vals**3 / (Km**3 + S_vals**3)
v_lin  = (Vmax / Km) * S_vals

fig, ax = plt.subplots()
ax.plot(S_vals, v_mm,   label="Michaelis-Menten")
ax.plot(S_vals, v_hill, label=r"Hill ($n = 3$)")
ax.plot(S_vals, v_lin,  ls=":", color="gray", label=r"linear ($k_f = V_\mathrm{max}/K_m$)")
ax.axhline(Vmax, color="gray", lw=0.8, ls="--", label=r"$V_\mathrm{max}$")
ax.axvline(Km,   color="gray", lw=0.8, ls="-.", label=r"$K_m$")
ax.set_xlabel(r"$[S]\ [\mathrm{mol/m^3}]$")
ax.set_ylabel(r"$v\ [\mathrm{mol/(m^3 \cdot s)}]$")
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
S0 = 600.0   # mol/m³  (initial substrate)

result_mm   = simulate(model_mm,   c0={"S": S0}, t_span=(0, 2000.0))
result_hill = simulate(model_hill, c0={"S": S0}, t_span=(0, 2000.0))
```

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-enzyme-sim

fig, axes = plt.subplots(1, 2, figsize=(9, 3.8), sharey=True)

for ax, result, label in [
    (axes[0], result_mm,   "Michaelis-Menten"),
    (axes[1], result_hill, r"Hill ($n = 3$)"),
]:
    ax.plot(result.t, result["S"], label="S (substrate)", color="C0")
    ax.plot(result.t, result["P"], label="P (product)",   color="C1")
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

---

The kinetic framework is now complete: rate constants can be temperature-independent (`RateConstantFixed`), Arrhenius (`RateConstantArrhenius`), or replaced by saturation kinetics (`MichaelisMenten`, `HillRate`).
The next chapter defines the integration contract between this library and CADET-Core: the residual and Jacobian that the solver calls at every step (@implementation-interface).
