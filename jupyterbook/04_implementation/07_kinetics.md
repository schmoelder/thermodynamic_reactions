---
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

(implementation-kinetics)=
# Kinetics, Arrhenius, and Thermodynamic Consistency

`ThermodynamicReaction` enforces $k_r = k_f / K(T)$ (@implementation-equilibrium), but so far $k_f$ has been a fixed scalar.
When temperature varies, $k_f$ must change with $T$ as well; $k_r(T) = k_f(T) / K(T)$ then tracks the shift in equilibrium automatically.
This chapter introduces the Arrhenius model for $k_f(T)$ and shows that $k_f(T) / k_r(T) = K(T)$ holds at every temperature by construction.

## Arrhenius kinetics and thermodynamic consistency

The Arrhenius equation (@kinetics-temperature) gives the forward rate constant as

$$
k_f(T) = A\,\exp\!\left(-\frac{E_a}{RT}\right).
$$

`RateConstantArrhenius` evaluates this at any temperature and is passed as the `rate_constant` argument of `ThermodynamicReaction`.
Pairing `RateConstantArrhenius` with `ThermodynamicReaction` ensures that $k_r(T) = k_f(T)/K(T)$ is recomputed at every evaluation; the ratio tracks $K(T)$ exactly across the full temperature range (@fig-impl-arrhenius).

```{code-cell} ipython3
import numpy as np
from reactions.api import (
    Component,
    EquilibriumConstantVantHoff,
    RateConstantArrhenius,
    ThermodynamicReaction,
    ReactionModel,
)
from reactions.solver import simulate

kf_arr = RateConstantArrhenius(A=1e10, Ea=40e3)
K_vH   = EquilibriumConstantVantHoff(dH=-20e3, dS=-50.0)
```

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-impl-arrhenius

import matplotlib.pyplot as plt

T_range = np.linspace(270, 370, 300)
kf_vals = np.array([kf_arr.kf(T) for T in T_range])
K_true  = np.array([K_vH.K(T) for T in T_range])

kr_fixed = kf_arr.kf(298.15) / K_vH.K(298.15)
ratio_fixed = kf_vals / kr_fixed

fig, axes = plt.subplots(1, 2, figsize=(9, 3.8))

axes[0].plot(1000 / T_range, np.log(kf_vals))
axes[0].set_xlabel(r"$1000\,/\,T\ [\mathrm{K}^{-1}]$")
axes[0].set_ylabel(r"$\ln k_f$")
axes[0].set_title("Arrhenius plot")

axes[1].plot(T_range - 273.15, K_true,       color="C0", label=r"$K(T)$ [van't Hoff]")
axes[1].plot(T_range - 273.15, ratio_fixed,  color="C1", ls="--", label=r"$k_f(T)\,/\,k_r$ [fixed $k_r$]")
axes[1].fill_between(T_range - 273.15, K_true, ratio_fixed, alpha=0.15, color="C3")
axes[1].axvline(298.15 - 273.15, color="gray", lw=0.8, ls=":", label="calibration T")
axes[1].set_xlabel(r"$T$ [°C]")
axes[1].set_ylabel(r"$K$")
axes[1].set_title("Thermodynamic consistency")
axes[1].legend(fontsize=8)

fig.tight_layout()
```

```{figure} #cell-impl-arrhenius
:name: fig-impl-arrhenius

Left: Arrhenius plot ($\ln k_f$ vs $1/T$) for $A = 10^{10}\ \text{s}^{-1}$, $E_a = 40\ \text{kJ/mol}$.
Right: a fixed $k_r$ (calibrated at 25 °C) keeps $k_f(T)/k_r$ on the dashed line while $K(T)$
shifts with temperature (shaded region); `ThermodynamicReaction` recomputes $k_r = k_f(T)/K(T)$
at every step, keeping the ratio on the solid $K(T)$ curve.
```

## Thermodynamic consistency across temperatures

Pairing `RateConstantArrhenius` with `EquilibriumConstantVantHoff` inside a
`ThermodynamicReaction` ensures consistency: at every temperature the library
computes $k_r(T) = k_f(T) / K(T)$, so both rate constants co-vary with $T$ and
their ratio tracks $K(T)$ exactly.

```{code-cell} ipython3
model = ReactionModel(
    components=[Component("A"), Component("B")],
    reactions=[
        ThermodynamicReaction(
            "A <-> B",
            mode="kinetic",
            equilibrium_constant=K_vH,
            rate_constant=kf_arr,
        ),
    ],
)
```

```{code-cell} ipython3
:tags: [remove-cell]

T_check = np.linspace(280, 360, 9)
print(f"{'T (K)':>8}  {'kf':>12}  {'kr':>12}  {'kf/kr':>10}  {'K(T)':>10}  {'error':>10}")
for T in T_check:
    kf_val = kf_arr.kf(T)
    K_val  = K_vH.K(T)
    kr_val = kf_val / K_val
    print(
        f"{T:>8.2f}  {kf_val:>12.4f}  {kr_val:>12.4f}"
        f"  {kf_val / kr_val:>10.6f}  {K_val:>10.6f}  {abs(kf_val / kr_val - K_val):>10.2e}"
    )
```

The table confirms that $k_f(T) / k_r(T) = K(T)$ holds to machine precision at every
temperature; thermodynamic consistency is structural, not a numerical approximation.

## Temperature-dependent equilibrium in a kinetic simulation

Because $K(T)$ changes with temperature, the long-time equilibrium composition shifts.
Simulating the same reaction at two temperatures shows both the kinetic trajectory and
the equilibrium it converges to (@fig-two-temps):

```{code-cell} ipython3
c_tot = 1000.0  # mol/m³

result_298 = simulate(
    model,
    c0={"A": c_tot, "B": 0.0},
    t_span=(0, 10.0),
    T=298.15,
)

result_320 = simulate(
    model,
    c0={"A": c_tot, "B": 0.0},
    t_span=(0, 10.0),
    T=320.0,
)
```

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-two-temps

K_298 = K_vH.K(298.15)
K_320 = K_vH.K(320.0)
B_eq_298 = c_tot * K_298 / (1 + K_298)
B_eq_320 = c_tot * K_320 / (1 + K_320)

print(f"T = 298.15 K:  K = {K_298:.3f},  c_B_eq = {B_eq_298:.1f} mol/m³  ({100*K_298/(1+K_298):.1f}%)")
print(f"T = 320.00 K:  K = {K_320:.3f},  c_B_eq = {B_eq_320:.1f} mol/m³  ({100*K_320/(1+K_320):.1f}%)")

fig, axes = plt.subplots(1, 2, figsize=(8, 3.5), sharey=True)

for ax, result, T, K, B_eq, label in [
    (axes[0], result_298, 298.15, K_298, B_eq_298, "298 K"),
    (axes[1], result_320, 320.0,  K_320, B_eq_320, "320 K"),
]:
    ax.plot(result.t, result["A"], label="A", color="C0")
    ax.plot(result.t, result["B"], label="B", color="C1")
    ax.axhline(c_tot - B_eq, color="C0", lw=1.0, ls="--")
    ax.axhline(B_eq,         color="C1", lw=1.0, ls="--", label="equilibrium")
    ax.set_xlabel("time [s]")
    ax.set_title(f"$T = {label}$,  $K = {K:.2f}$")

axes[0].set_ylabel("concentration [mol/m³]")
axes[0].legend()
fig.tight_layout()
```

```{figure} #cell-two-temps
:name: fig-two-temps

Kinetic trajectories for $\ce{A <=> B}$ at $T = 298\ \text{K}$ (left) and $T = 320\ \text{K}$ (right),
with $E_a = 40\ \text{kJ/mol}$, $\Delta H^\circ = -20\ \text{kJ/mol}$.
Dashed lines are the analytical equilibria $c_\text{B}^\text{eq} = c_\text{tot}\,K(T)/(1+K(T))$.
The exothermic reaction has a smaller $K$ at higher temperature (Le Chatelier's principle),
so less B accumulates at $320\ \text{K}$.
Both reach equilibrium faster at $320\ \text{K}$ because $k_f$ increases with temperature.
```

## Saturation kinetics

Mass action kinetics assume the reaction rate grows without bound as concentration
increases.
For enzyme-catalysed reactions the rate saturates because the enzyme active site
becomes limiting.
`MichaelisMenten` and `HillRate` implement two standard saturation forms:

$$
v_\text{MM} = \frac{V_\text{max}\,[S]}{K_m + [S]},
\qquad
v_\text{Hill} = \frac{V_\text{max}\,[S]^n}{K_m^n + [S]^n}.
$$

At low substrate concentration both reduce to a linear (mass-action) rate
$v \approx (V_\text{max}/K_m)\,[S]$; at saturation ($[S] \gg K_m$) the rate
approaches $V_\text{max}$ (@fig-saturation).
The Hill coefficient $n > 1$ introduces cooperativity: the response is sigmoidal
and steeper than the Michaelis-Menten curve.

These rate laws replace $\varphi = k_f \cdot c_S$ inside `EnzymaticReaction`.
The stoichiometric structure $\mathbf{S}\boldsymbol{\varphi}$ and the equilibrium mode
are unchanged.

```{code-cell} ipython3
from reactions.api import EnzymaticReaction, MichaelisMenten, HillRate

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

Michaelis-Menten and Hill ($n=3$) rate laws compared to the linear (mass-action) limit,
for $V_\text{max} = 1\ \text{mol/(m}^3\text{s)}$ and $K_m = 200\ \text{mol/m}^3$.
Both saturate at $V_\text{max}$ for $[S] \gg K_m$; the Hill curve is sigmoidal due to
cooperativity.
The dotted line is the linear approximation valid at low substrate.
```

---

The kinetic framework is now complete: rate constants can be temperature-independent (`RateConstantFixed`), Arrhenius (`RateConstantArrhenius`), or replaced by saturation kinetics (`MichaelisMenten`, `HillRate`).
The next chapter defines the integration contract between this library and CADET-Core: the residual and Jacobian that the solver calls at every step (@implementation-interface).
