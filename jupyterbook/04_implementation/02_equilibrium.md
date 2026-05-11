---
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

(implementation-equilibrium)=
# Equilibrium and Thermodynamic Consistency

Part 3 established the equilibrium condition $Q = K$ and its temperature dependence via the van't Hoff equation (@equilibrium, @equilibrium-temperature).
`ThermodynamicReaction` is the physically constrained formulation: given $K(T)$, it enforces $\ln Q(\mathbf{a}, T) = \ln K(T)$ as a residual constraint, and derives $k^r(T) = k^f(T)/K(T)$ whenever a rate constant is also supplied.
`MassActionReaction`, by contrast, treats $k^r$ as a free parameter with no thermodynamic grounding (@implementation-source-term).
This chapter works through the equilibrium mode, leaving kinetic mode for the next chapter.


## $\ce{A <=> B}$: equilibrium with a fixed $K$

The simplest case is a single reversible reaction.
For $\ce{A <=> B}$ with ideal activities the equilibrium condition $Q = K$ gives (@equilibrium)

$$
c_\text{B}^\text{eq} = \frac{K}{1+K}\,c_\text{tot}, \qquad
c_\text{A}^\text{eq} = \frac{1}{1+K}\,c_\text{tot}.
$$

As $K$ increases, A converts more completely to B; for $K \gg 1$ the reaction goes
to completion (@fig-eq-ab).
`ThermodynamicReaction(mode="equil")` enforces this constraint algebraically; `solve_equilibrium` finds the composition by Newton iteration on the residual $\ln Q - \ln K$ in log-activity space.

```{code-cell} ipython3
import numpy as np
from reactions.api import (
    Component,
    EquilibriumConstant,
    EquilibriumConstantVantHoff,
    EquilibriumConstantVantHoffCp,
    ThermodynamicReaction,
    ReactionModel,
)
from reactions.solver import simulate, solve_equilibrium

K_val = 4.0
c_tot = 1000.0  # mol/m³

model_eq = ReactionModel(
    components=[Component("A"), Component("B")],
    reactions=[
        ThermodynamicReaction(
            "A <-> B",
            mode="equil",
            equilibrium_constant=EquilibriumConstant(K_val),
        ),
    ],
)

c_eq = solve_equilibrium(
    model_eq,
    c0={"A": c_tot / 2, "B": c_tot / 2},
)
```

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-eq-ab

import matplotlib.pyplot as plt

print(f"c_A = {c_eq['A']:.4f} mol/m³")
print(f"c_B = {c_eq['B']:.4f} mol/m³")
print(f"c_B / c_A = {c_eq['B'] / c_eq['A']:.6f}   (K = {K_val})")

K_vals = np.linspace(0, 10, 300)
fA = 1 / (1 + K_vals)
fB = K_vals / (1 + K_vals)

fig, ax = plt.subplots()
ax.plot(K_vals, fA, label=r"$c_\mathrm{A}^\mathrm{eq}/c_\mathrm{tot}$", color="C0")
ax.plot(K_vals, fB, label=r"$c_\mathrm{B}^\mathrm{eq}/c_\mathrm{tot}$", color="C1")
ax.axvline(K_val, color="gray", lw=0.8, ls=":", label=f"K = {K_val:.0f} (example)")
ax.set_xlabel(r"$K$")
ax.set_ylabel("equilibrium fraction")
ax.legend()
fig.tight_layout()
```

```{figure} #cell-eq-ab
:name: fig-eq-ab

Equilibrium fractions of A and B as a function of $K$ for $\ce{A <=> B}$.
At $K = 1$ both species are equally present; at the example value $K = 4$ (dotted line)
80 % of material ends up as B.
```


## Temperature dependence of $K$

When temperature varies, a fixed ratio $k^f/k^r$ calibrated at one temperature diverges from the true $K(T)$; `ThermodynamicReaction` tracks it exactly at every temperature (@fig-eq-drift).

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-eq-drift

dH = -20e3  # J/mol  (exothermic)
dS = -50.0  # J/(mol·K)
K_vH_drift = EquilibriumConstantVantHoff(dH=dH, dS=dS)

T_ref = 298.15
K_ref = K_vH_drift.K(T_ref)

T_sweep = np.linspace(270, 370, 200)
K_true = np.array([K_vH_drift.K(T) for T in T_sweep])
K_fixed = np.full_like(T_sweep, K_ref)

fig, ax = plt.subplots()
ax.plot(T_sweep - 273.15, K_true, color="C0", label=r"$K(T)$ [van't Hoff]")
ax.plot(
    T_sweep - 273.15,
    K_fixed,
    color="C1",
    ls="--",
    label=f"$k^f/k^r$ [fixed at {T_ref - 273.15:.0f} °C]",
)
ax.fill_between(T_sweep - 273.15, K_true, K_fixed, alpha=0.15, color="C3")
ax.axvline(
    T_ref - 273.15, color="gray", lw=0.8, ls=":", label="calibration temperature"
)
ax.set_xlabel(r"$T$ [°C]")
ax.set_ylabel(r"$K$")
ax.legend()
fig.tight_layout()
```

```{figure} #cell-eq-drift
:name: fig-eq-drift

Equilibrium constant $K(T)$ (van't Hoff, $\Delta H^\circ = -20\ \mathrm{kJ/mol}$) alongside the fixed ratio $k^f/k^r$ calibrated at 25 °C.
The two curves coincide only at the calibration temperature; the shaded region shows how the divergence grows away from it.
`ThermodynamicReaction` with `EquilibriumConstantVantHoff` keeps the ratio on the $K(T)$ curve at every temperature.
```

The van't Hoff equation (@equilibrium-temperature) gives $K(T)$ from standard thermodynamic data,

$$
\ln K(T) = -\frac{\Delta H^\circ}{RT} + \frac{\Delta S^\circ}{R}.
$$

For an exothermic reaction ($\Delta H^\circ < 0$) heating decreases $K$, shifting the equilibrium toward the reactant (@fig-eq-vanthoff).

```{code-cell} ipython3
dH = -20e3  # J/mol  (exothermic)
dS = -50.0  # J/(mol·K)

K_vH = EquilibriumConstantVantHoff(dH=dH, dS=dS)
```

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-eq-vanthoff

print(f"{'T (°C)':>8}  {'K':>8}  {'c_B/c_tot':>10}")
for T in [280.0, 298.15, 320.0, 350.0]:
    K = K_vH.K(T)
    print(f"{T - 273.15:>8.1f}  {K:>8.4f}  {K / (1 + K):>10.4f}")

T_arr = np.linspace(270, 370, 300)
K_arr = np.array([K_vH.K(T) for T in T_arr])
fB_arr = K_arr / (1 + K_arr)

fig, axes = plt.subplots(1, 2, figsize=(8, 3.5))

axes[0].plot(T_arr - 273.15, K_arr)
axes[0].set_xlabel("T (°C)")
axes[0].set_ylabel(r"$K$")
axes[0].set_title("Equilibrium constant")

axes[1].plot(T_arr - 273.15, fB_arr)
axes[1].set_xlabel("T (°C)")
axes[1].set_ylabel(r"$c_\mathrm{B}^\mathrm{eq}\,/\,c_\mathrm{tot}$")
axes[1].set_title("Equilibrium composition")

fig.tight_layout()
```

```{figure} #cell-eq-vanthoff
:name: fig-eq-vanthoff

Temperature dependence of $K$ (left) and equilibrium fraction of B (right) for $\Delta H^\circ = -20\ \mathrm{kJ/mol}$, $\Delta S^\circ = -50\ \mathrm{J\,mol^{-1}\,K^{-1}}$.
Heating shifts the equilibrium towards A, consistent with Le Chatelier's principle.
```

Passing `EquilibriumConstantVantHoff` to `ThermodynamicReaction` makes the equilibrium composition temperature-aware: `solve_equilibrium` calls $K(T)$ at the specified temperature.

```{code-cell} ipython3
model_vH = ReactionModel(
    components=[Component("A"), Component("B")],
    reactions=[
        ThermodynamicReaction(
            "A <-> B",
            mode="equil",
            equilibrium_constant=K_vH,
        ),
    ],
)

c_320 = solve_equilibrium(
    model_vH,
    c0={"A": c_tot / 2, "B": c_tot / 2},
    T=320.0,
)
```

```{code-cell} ipython3
:tags: [remove-cell]

K_320 = K_vH.K(320.0)
print(f"K(320 K)  = {K_320:.6f}")
print(
    f"c_B / c_A = {c_320['B'] / c_320['A']:.6f}   (error: {abs(c_320['B'] / c_320['A'] - K_320):.2e})"
)
```


## Heat capacity correction

For reactions where $\Delta H^\circ$ itself varies with temperature, `EquilibriumConstantVantHoffCp` applies the Kirchhoff relation (@equilibrium-temperature):

$$
\Delta H^\circ(T) = \Delta H^\circ_\text{ref} + \Delta C_p\,(T - T_\text{ref}),
\qquad
\Delta S^\circ(T) = \Delta S^\circ_\text{ref} + \Delta C_p\ln\frac{T}{T_\text{ref}}.
$$

A large positive $\Delta C_p$ means $\Delta H^\circ(T)$ itself increases with temperature.
If $\Delta C_p$ is large enough, an initially exothermic reaction becomes endothermic at high $T$, reversing the direction of the $K(T)$ curve, a behavior the two-parameter van't Hoff model cannot capture:

```{code-cell} ipython3
from reactions.api import EquilibriumConstantCustom, EquilibriumConstantTabulated

K_cp = EquilibriumConstantVantHoffCp(dH=dH, dS=dS, dCp=+200.0, T_ref=298.15)
```

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-eq-cp

T_wide = np.linspace(250, 450, 400)
K_vH_arr = np.array([K_vH.K(T) for T in T_wide])
K_cp_arr = np.array([K_cp.K(T) for T in T_wide])

fig, ax = plt.subplots()
ax.plot(1000 / T_wide, np.log(K_vH_arr), label="VantHoff (straight line)")
ax.plot(
    1000 / T_wide,
    np.log(K_cp_arr),
    "--",
    label=r"VantHoffCp ($\Delta C_p = +200\ \mathrm{J\,mol^{-1}\,K^{-1}}$)",
)
ax.axvline(1000 / 298.15, color="gray", lw=0.8, ls=":", label=r"$T_\mathrm{ref}$")
ax.set_xlabel(r"$1000\,/\,T\ [\mathrm{K}^{-1}]$")
ax.set_ylabel(r"$\ln K$")
ax.legend()
fig.tight_layout()
```

```{figure} #cell-eq-cp
:name: fig-eq-cp

Van't Hoff plot ($\ln K$ vs $1/T$) comparing `EquilibriumConstantVantHoff` and `EquilibriumConstantVantHoffCp` ($\Delta H^\circ = -20\ \mathrm{kJ/mol}$, $\Delta C_p = +200\ \mathrm{J\,mol^{-1}\,K^{-1}}$).
`VantHoff` is a straight line by construction; `VantHoffCp` curves because $\Delta H^\circ(T)$ changes with temperature, eventually reversing sign and turning the slope around.
```

The Kirchhoff correction matters when $\Delta C_p$ is large or the temperature range is wide; for most aqueous reactions at near-ambient conditions `EquilibriumConstantVantHoff` is adequate (@fig-eq-cp).
For empirically fitted $K(T)$, `EquilibriumConstantCustom` accepts any callable `(T: float) -> float`, and `EquilibriumConstantTabulated` interpolates linearly through measured $(T, K)$ pairs:

```{code-cell} ipython3
# Fitted exponential from experiment
K_exp = EquilibriumConstantCustom(lambda T: 2.5 * np.exp(-1500 / T))

# Measured data points
K_tab = EquilibriumConstantTabulated(
    T_data=[280.0, 298.15, 320.0, 350.0],
    K_data=[5.2, 4.0, 2.8, 1.7],
)
```

Both slot into `ThermodynamicReaction` identically to the van't Hoff forms.

---

Equilibrium mode determines *where* the system ends up; it says nothing about *how fast* it gets there.
The next chapter introduces kinetic mode, where a rate constant sets the relaxation timescale while $k^r(T) = k^f(T)/K(T)$ remains enforced throughout (@implementation-kinetics).
