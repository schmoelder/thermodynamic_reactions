---
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

(implementation-equilibrium)=
# Thermodynamic Consistency and K(T)

The previous chapter showed that treating $k_f$ and $k_r$ as independent, temperature-independent parameters breaks the link between kinetics and thermodynamics: their ratio no longer tracks the equilibrium constant $K(T)$ as temperature changes.
`ThermodynamicReaction` restores this consistency by taking $K(T)$ as the primary input and deriving $k_r(T) = k_f(T)/K(T)$ at every temperature, rather than treating $k_r$ as a free parameter.
In the fast-reaction limit, the kinetic formulation transitions naturally to equilibrium mode, where the ODE for the dependent species is replaced by the algebraic constraint $\ln Q_j - \ln K_j(T) = 0$, enforced at each grid point.

## A⇌B: equilibrium with a fixed K

The simplest case is a single reversible reaction.
For $\text{A} \rightleftharpoons \text{B}$ with ideal activities the equilibrium condition $Q = K$ gives (@equilibrium)

$$
c_\text{B}^\text{eq} = \frac{K}{1+K}\,c_\text{tot}, \qquad
c_\text{A}^\text{eq} = \frac{1}{1+K}\,c_\text{tot}.
$$

As $K$ increases, A converts more completely to B; for $K \gg 1$ the reaction goes
to completion (@fig-eq-ab).
`ThermodynamicReaction(mode="equil")` enforces this constraint algebraically; `solve_equilibrium` finds the composition by Newton iteration in log space.

```{code-cell} ipython3
import numpy as np
from reactions.api import (
    Component,
    EquilibriumConstant,
    EquilibriumConstantVantHoff,
    EquilibriumConstantVantHoffCp,
    MassActionReaction,
    ThermodynamicReaction,
    ReactionModel,
)
from reactions.solver import simulate, solve_equilibrium

K_val = 4.0
c_tot = 1000.0    # mol/m³

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

Equilibrium fractions of A and B as a function of $K$ for A⇌B.
At $K = 1$ both species are equally present; at the example value $K = 4$ (dotted line)
80 % of material ends up as B.
```

## Temperature dependence of $K$

`ThermodynamicReaction` with `EquilibriumConstantVantHoff` keeps $k_r(T) = k_f(T)/K(T)$ exact at every temperature.
`MassActionReaction` with a fixed ratio is correct only at the calibration temperature; elsewhere the ratio drifts from $K(T)$, with the error growing with distance from that point (@fig-eq-drift).

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-eq-drift

dH = -20e3    # J/mol  (exothermic)
dS = -50.0    # J/(mol·K)
K_vH_drift = EquilibriumConstantVantHoff(dH=dH, dS=dS)

T_ref   = 298.15
K_ref   = K_vH_drift.K(T_ref)

T_sweep = np.linspace(270, 370, 200)
K_true  = np.array([K_vH_drift.K(T) for T in T_sweep])
K_fixed = np.full_like(T_sweep, K_ref)

fig, ax = plt.subplots()
ax.plot(T_sweep - 273.15, K_true,  color="C0", label=r"$K(T)$ [van't Hoff]")
ax.plot(T_sweep - 273.15, K_fixed, color="C1", ls="--", label=f"$k_f/k_r$ [fixed at {T_ref - 273.15:.0f} °C]")
ax.fill_between(T_sweep - 273.15, K_true, K_fixed, alpha=0.15, color="C3")
ax.axvline(T_ref - 273.15, color="gray", lw=0.8, ls=":", label="calibration temperature")
ax.set_xlabel(r"$T$ [°C]")
ax.set_ylabel(r"$K$")
ax.legend()
fig.tight_layout()
```

```{figure} #cell-eq-drift
:name: fig-eq-drift

Equilibrium constant $K(T)$ (van't Hoff, $\Delta H^\circ = -20\ \text{kJ/mol}$) vs the fixed ratio $k_f/k_r$ calibrated at 25 °C.
The shaded region is the error incurred by `MassActionReaction` at temperatures other than the calibration point.
`ThermodynamicReaction` with `EquilibriumConstantVantHoff` tracks $K(T)$ exactly by recomputing $k_r = k_f / K(T)$ at every evaluation.
```

The van't Hoff equation (@equilibrium-temperature) relates $K$ to the standard enthalpy and entropy of reaction,

$$
\ln K(T) = -\frac{\Delta H^\circ}{RT} + \frac{\Delta S^\circ}{R}.
$$

`EquilibriumConstantVantHoff` evaluates this at any temperature.
For an exothermic reaction ($\Delta H^\circ < 0$) heating decreases $K$, shifting the equilibrium back towards the reactant (@fig-eq-vanthoff).

```{code-cell} ipython3
dH = -20e3    # J/mol  (exothermic)
dS = -50.0    # J/(mol·K)

K_vH = EquilibriumConstantVantHoff(dH=dH, dS=dS)
```

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-eq-vanthoff

print(f"{'T (°C)':>8}  {'K':>8}  {'c_B/c_tot':>10}")
for T in [280.0, 298.15, 320.0, 350.0]:
    K = K_vH.K(T)
    print(f"{T - 273.15:>8.1f}  {K:>8.4f}  {K / (1 + K):>10.4f}")

T_arr  = np.linspace(270, 370, 300)
K_arr  = np.array([K_vH.K(T) for T in T_arr])
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

Temperature dependence of $K$ (left) and equilibrium fraction of B (right) for $\Delta H^\circ = -20\ \text{kJ/mol}$, $\Delta S^\circ = -50\ \text{J mol}^{-1}\text{K}^{-1}$.
Heating shifts the equilibrium towards A, consistent with Le Chatelier's principle.
```

Passing `EquilibriumConstantVantHoff` to `ThermodynamicReaction` makes the equilibrium composition temperature-aware: `solve_equilibrium` calls `K(T)` at the specified temperature.

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
print(f"c_B / c_A = {c_320['B'] / c_320['A']:.6f}   (error: {abs(c_320['B'] / c_320['A'] - K_320):.2e})")
```

## Heat capacity correction

For reactions where $\Delta H^\circ$ itself varies with temperature, `EquilibriumConstantVantHoffCp` applies the Kirchhoff relation (@equilibrium-temperature):

$$
\Delta H^\circ(T) = \Delta H^\circ_\text{ref} + \Delta C_p\,(T - T_\text{ref}),
\qquad
\Delta S^\circ(T) = \Delta S^\circ_\text{ref} + \Delta C_p\ln\frac{T}{T_\text{ref}}.
$$

A large positive $\Delta C_p$ means $\Delta H^\circ(T)$ itself increases with temperature.
If $\Delta C_p$ is large enough, an initially exothermic reaction becomes endothermic at high $T$, reversing the direction of the $K(T)$ curve, a behaviour the two-parameter van't Hoff model cannot capture:

```{code-cell} ipython3
K_cp = EquilibriumConstantVantHoffCp(dH=dH, dS=dS, dCp=+200.0, T_ref=298.15)
```

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-eq-cp

T_wide   = np.linspace(250, 450, 400)
K_vH_arr = np.array([K_vH.K(T) for T in T_wide])
K_cp_arr = np.array([K_cp.K(T) for T in T_wide])

fig, ax = plt.subplots()
ax.plot(1000 / T_wide, np.log(K_vH_arr), label="VantHoff (straight line)")
ax.plot(1000 / T_wide, np.log(K_cp_arr), "--", label=r"VantHoffCp ($\Delta C_p = +200\ \mathrm{J\,mol^{-1}\,K^{-1}}$)")
ax.axvline(1000 / 298.15, color="gray", lw=0.8, ls=":", label=r"$T_\mathrm{ref}$")
ax.set_xlabel(r"$1000\,/\,T\ [\mathrm{K}^{-1}]$")
ax.set_ylabel(r"$\ln K$")
ax.legend()
fig.tight_layout()
```

```{figure} #cell-eq-cp
:name: fig-eq-cp

Van't Hoff plot ($\ln K$ vs $1/T$) comparing `EquilibriumConstantVantHoff` and `EquilibriumConstantVantHoffCp` ($\Delta H^\circ = -20\ \text{kJ/mol}$, $\Delta C_p = +200\ \text{J mol}^{-1}\text{K}^{-1}$).
`VantHoff` is a straight line by construction; `VantHoffCp` curves because $\Delta H^\circ(T)$ changes with temperature, eventually reversing sign and turning the slope around.
```

The Kirchhoff correction matters when $\Delta C_p$ is large or the temperature range is wide; for most aqueous reactions at near-ambient conditions `EquilibriumConstantVantHoff` is adequate (@fig-eq-cp).

---

With $K(T)$ enforced, ideal activities ($\gamma_i = 1$) are still assumed.
The next chapter replaces concentrations with activities $a_i = \gamma_i c_i / c^\circ$, completing the non-ideal picture before temperature-dependent kinetics are introduced (@implementation-activity).
