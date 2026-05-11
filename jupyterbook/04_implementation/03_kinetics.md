---
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

(implementation-kinetics)=
# Kinetic Mode and Arrhenius Rate Constants

Equilibrium mode enforces $Q = K(T)$ algebraically and determines the composition without any rate information.
Kinetic mode adds a rate constant $k^f$ that sets the relaxation timescale: the system evolves as an ODE and converges to the same long-time equilibrium, but at a finite speed controlled by $k^f$.
The constraint $k^r(T) = k^f(T)/K(T)$ is enforced in both modes.


## Kinetic mode with a fixed rate constant

`ThermodynamicReaction(mode="kinetic")` requires a `rate_constant` argument.
`RateConstantFixed` supplies a constant $k^f$ independent of temperature; the simplest choice when the rate is known at one operating condition.

The rate constant in `ThermodynamicReaction` has units $\text{mol}/(\text{m}^3\cdot\mathrm{s})$ because the reaction flux is evaluated in terms of dimensionless activities $a_i = c_i / c^\circ$ rather than concentrations directly.
For a first-order reaction at ideal conditions this is related to the mass-action rate constant by $k^{f,\text{thermo}} = k^{f,\text{MA}} \cdot c^\circ$ (@implementation-source-term, @mass-action-law).
The factor $c^\circ = 1000\ \mathrm{mol/m^3}$ reflects the standard-state choice: `MassActionReaction` implicitly uses $c^\circ = 1\ \mathrm{mol/m^3}$, so the same physical rate requires a 1000-fold larger $k^f$ in `ThermodynamicReaction`.
`MassActionReaction` is the $\gamma_i = 1$, free-$k^r$ limit of `ThermodynamicReaction` (@implementation-source-term): equivalent at fixed temperature when $k^r$ is set to $k^f/K$ by hand, but without the mechanism to enforce this as temperature changes.

```{code-cell} ipython3
import numpy as np
from reactions.api import (
    Component,
    EquilibriumConstant,
    EquilibriumConstantVantHoff,
    RateConstantFixed,
    RateConstantArrhenius,
    ThermodynamicReaction,
    ReactionModel,
)
from reactions.solver import simulate, solve_equilibrium

K_val = 4.0
c_tot = 1000.0    # mol/m³
kf_thermo = 2000.0  # mol/(m³·s); τ ≈ C_REF / (kf * (1 + 1/K)) ≈ 0.4 s

comp_a = Component("A")
comp_b = Component("B")

model_kin = ReactionModel(
    components=[comp_a, comp_b],
    reactions=[
        ThermodynamicReaction(
            "A <-> B",
            mode="kinetic",
            equilibrium_constant=EquilibriumConstant(K_val),
            rate_constant=RateConstantFixed(kf_value=kf_thermo),
        ),
    ],
)

model_eq = ReactionModel(
    components=[comp_a, comp_b],
    reactions=[
        ThermodynamicReaction(
            "A <-> B",
            mode="equil",
            equilibrium_constant=EquilibriumConstant(K_val),
        ),
    ],
)

result_kin = simulate(model_kin, c0={"A": c_tot}, t_span=(0, 3.0))
c_eq = solve_equilibrium(model_eq, c0={"A": c_tot / 2, "B": c_tot / 2})
```

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-kin-fixed

import matplotlib.pyplot as plt

B_eq = c_tot * K_val / (1 + K_val)
A_eq = c_tot - B_eq

print(f"solve_equilibrium: c_A = {c_eq['A']:.4f}, c_B = {c_eq['B']:.4f}")
print(f"simulate long-time: c_A = {result_kin['A'][-1]:.4f}, c_B = {result_kin['B'][-1]:.4f}")

fig, ax = plt.subplots()
ax.plot(result_kin.t, result_kin["A"], label="A", color="C0")
ax.plot(result_kin.t, result_kin["B"], label="B", color="C1")
ax.axhline(A_eq, color="C0", lw=1.0, ls="--")
ax.axhline(B_eq, color="C1", lw=1.0, ls="--", label="equilibrium")
ax.set_xlabel("time [s]")
ax.set_ylabel(r"concentration [mol/m³]")
ax.legend()
fig.tight_layout()
```

```{figure} #cell-kin-fixed
:name: fig-kin-fixed

Kinetic simulation of $\ce{A <=> B}$ with $K = 4$ and $k^f = 2000\ \mathrm{mol/(m^3 \cdot s)}$.
Dashed lines are the equilibrium concentrations from `solve_equilibrium`; the trajectory converges to the same values.
The rate constant sets the timescale; $K$ alone sets the endpoint.
```

The equilibrium composition is identical whether computed via `mode="equil"` or reached by time integration in `mode="kinetic"`, because both enforce the same fixed point $\mathbf{f}_\text{react}(\mathbf{c}) = 0$.
The distinction is purely about whether the constraint is imposed algebraically (fast-reaction limit) or resolved dynamically.


## Arrhenius rate constant

When temperature varies, $k^f$ must change with $T$ as well.
The Arrhenius equation (@kinetics-temperature) gives

$$
k^f(T) = A\,\exp\!\left(-\frac{E_a}{RT}\right).
$$

`RateConstantArrhenius` evaluates this at any temperature and is passed as the `rate_constant` argument of `ThermodynamicReaction`.
With $K(T)$ from `EquilibriumConstantVantHoff`, the derived $k^r(T) = k^f(T)/K(T)$ is recomputed at every evaluation, so the ratio tracks $K(T)$ exactly across the full temperature range (@fig-impl-arrhenius).

```{code-cell} ipython3
kf_arr = RateConstantArrhenius(A=1e10, Ea=40e3)
K_vH   = EquilibriumConstantVantHoff(dH=-20e3, dS=-50.0)
```

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-impl-arrhenius

T_range = np.linspace(270, 370, 300)
kf_vals = np.array([kf_arr.kf(T) for T in T_range])
K_true  = np.array([K_vH.K(T) for T in T_range])

kr_fixed = kf_arr.kf(298.15) / K_vH.K(298.15)
ratio_fixed = kf_vals / kr_fixed

fig, axes = plt.subplots(1, 2, figsize=(9, 3.8))

axes[0].plot(1000 / T_range, np.log(kf_vals))
axes[0].set_xlabel(r"$1000\,/\,T\ [\mathrm{K}^{-1}]$")
axes[0].set_ylabel(r"$\ln k^f$")
axes[0].set_title("Arrhenius plot")

axes[1].plot(T_range - 273.15, K_true,       color="C0", label=r"$K(T)$ [van't Hoff]")
axes[1].plot(T_range - 273.15, ratio_fixed,  color="C1", ls="--", label=r"$k^f(T)\,/\,k^r$ [fixed $k^r$]")
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

Left: Arrhenius plot ($\ln k^f$ vs $1/T$) for $A = 10^{10}\ \mathrm{s}^{-1}$, $E_a = 40\ \mathrm{kJ/mol}$.
Right: a fixed $k^r$ (calibrated at 25 °C) keeps $k^f(T)/k^r$ on the dashed line while $K(T)$
shifts with temperature (shaded region); `ThermodynamicReaction` recomputes $k^r = k^f(T)/K(T)$
at every step, keeping the ratio on the solid $K(T)$ curve.
```


## Consistency across temperatures

Pairing `RateConstantArrhenius` with `EquilibriumConstantVantHoff` inside a `ThermodynamicReaction` ensures that at every temperature the library computes $k^r(T) = k^f(T)/K(T)$, so both rate constants co-vary with $T$ and their ratio tracks $K(T)$ exactly.

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

The table confirms that $k^f(T)/k^r(T) = K(T)$ holds to machine precision at every temperature, because $k^r(T)$ is defined as $k^f(T)/K(T)$: thermodynamic consistency is structural, not a numerical approximation.


## Temperature-dependent equilibrium in a kinetic simulation

Because $K(T)$ changes with temperature, the long-time equilibrium composition shifts.
Simulating the same reaction at two temperatures shows both the kinetic trajectory and the equilibrium it converges to (@fig-two-temps):

```{code-cell} ipython3
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

axes[0].set_ylabel(r"concentration [mol/m³]")
axes[0].legend()
fig.tight_layout()
```

```{figure} #cell-two-temps
:name: fig-two-temps

Kinetic trajectories for $\ce{A <=> B}$ at $T = 298\ \text{K}$ (left) and $T = 320\ \text{K}$ (right),
with $E_a = 40\ \mathrm{kJ/mol}$, $\Delta H^\circ = -20\ \mathrm{kJ/mol}$.
Dashed lines are the analytical equilibria $c_\text{B}^\text{eq} = c_\text{tot}\,K(T)/(1+K(T))$.
The exothermic reaction has a smaller $K$ at higher temperature (Le Chatelier's principle),
so less B accumulates at $320\ \text{K}$.
Both reach equilibrium faster at $320\ \text{K}$ because $k^f$ increases with temperature.
```

```{admonition} Limitation: elementary kinetics assumed
:class: warning

`ThermodynamicReaction(mode="kinetic")` with `RateConstantFixed` or `RateConstantArrhenius` assumes the overall rate law is elementary in activities: reaction order equals stoichiometric coefficient.
For complex mechanisms, the apparent rate constant can depend on concentration (Michaelis-Menten is the canonical example), and the apparent macroscopic rate constants may not satisfy $k^f/k^r = K$ (@mass-action-law).
In those cases, use `EnzymaticReaction` with a custom rate closure (@implementation-enzyme).
Thermodynamic consistency at the overall reaction level cannot be enforced automatically; it must be verified against the microscopic mechanism.
```

## Prescribed temperature programme

The two-temperature comparison runs each trajectory at a fixed temperature.
Passing a callable `T(t)` to `simulate` lets the temperature vary continuously within a single integration: the solver evaluates $k^f(T(t))$, $K(T(t))$, and $k^r(T(t)) = k^f(T(t))/K(T(t))$ at every step, so the instantaneous equilibrium composition shifts along with the temperature.

```{code-cell} ipython3
T_start, T_end, t_end = 298.15, 320.0, 10.0
T_ramp = lambda t: T_start + (T_end - T_start) * t / t_end  # noqa: E731

result_ramp = simulate(
    model,
    c0={"A": c_tot, "B": 0.0},
    t_span=(0, t_end),
    T=T_ramp,
)
```

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-T-ramp

T_inst = result_ramp.T_profile
K_inst = np.array([K_vH.K(T) for T in T_inst])
B_inst = c_tot * K_inst / (1 + K_inst)

fig, axes = plt.subplots(1, 2, figsize=(9, 3.5))

axes[0].plot(result_ramp.t, result_ramp["A"], color="C0", label="A")
axes[0].plot(result_ramp.t, result_ramp["B"], color="C1", label="B")
axes[0].plot(result_ramp.t, B_inst,           color="C1", ls="--", lw=1.0,
             label=r"$c_B^\mathrm{eq}(T(t))$")
axes[0].plot(result_ramp.t, c_tot - B_inst,   color="C0", ls="--", lw=1.0)
axes[0].set_xlabel("time [s]")
axes[0].set_ylabel(r"concentration [mol/m³]")
axes[0].legend(fontsize=8)

axes[1].plot(result_ramp.t, T_inst - 273.15, color="C3")
axes[1].set_xlabel("time [s]")
axes[1].set_ylabel(r"$T$ [°C]")

fig.tight_layout()
```

```{figure} #cell-T-ramp
:name: fig-T-ramp

Kinetic simulation under a linear temperature ramp from $298\ \text{K}$ to $320\ \text{K}$ over $10\ \text{s}$.
Left: concentration trajectories; dashed lines track the instantaneous equilibrium
$c_\text{B}^\text{eq}(T(t)) = c_\text{tot}\,K(T(t))/(1+K(T(t)))$.
Right: the prescribed temperature programme stored in `result_ramp.T_profile`.
The exothermic reaction loses product B as temperature rises, consistent with Le Chatelier's principle.
```

The trajectory tracks the moving equilibrium closely because $k^f$ is large relative to the ramp rate.
A slower rate constant or a faster ramp would introduce a visible lag between the solid and dashed curves, showing that temperature-induced equilibrium shifts are rate-limited by kinetics, not instantaneous.


## Structure of a reaction model

Three chapters have now introduced three API classes — `MassActionReaction`, `ThermodynamicReaction(mode="equil")`, and `ThermodynamicReaction(mode="kinetic")` — that all reduce to the same source term $\mathbf{S}\boldsymbol{\varphi}$ but differ in what constrains $\varphi$.
A reaction model is a composition of five independent layers:

| Layer | What it specifies |
| ----- | ----------------- |
| **Stoichiometry** | $\nu_{ij}$, kinetic exponents $e^f_{ij}$, $e^r_{ij}$: which species participate and how |
| **Equilibrium relation** | $K(T)$: the admissible long-run composition |
| **Kinetic closure** | $k^f(T)$ and rate form: the relaxation timescale |
| **Activity model** | $\gamma_i(c, I, T)$: replaces concentrations with activities |
| **State dependencies** | how $T$ and $I$ enter the layers above |

The three current classes occupy different positions along the first three layers:

| Class | Equilibrium relation | Kinetic closure |
| ----- | -------------------- | --------------- |
| `MassActionReaction` | none — $k^r$ is a free parameter | mass-action polynomial |
| `ThermodynamicReaction(mode="kinetic")` | $K(T)$ enforces $k^r = k^f/K$ automatically | Arrhenius, polynomial, tabulated, or custom |
| `ThermodynamicReaction(mode="equil")` | $K(T)$ as an algebraic constraint | none — fast-reaction limit |
| `EnzymaticReaction` *(introduced in @implementation-enzyme)* | none — reverse rate absent | custom saturating closure (MM, Hill, user-supplied) |

The key distinction is whether a class carries an explicit equilibrium relation.
`MassActionReaction` and `EnzymaticReaction` are siblings at the stoichiometry-plus-closure level: neither enforces thermodynamic consistency between forward and reverse rates.
`MassActionReaction` is not a deficient form — it is permissive, empirical, and appropriate when thermodynamic parameters are unavailable or when the two rate constants are treated as independently fitted quantities.
`ThermodynamicReaction` adds the equilibrium relation as an explicit layer; in kinetic mode that relation derives $k^r$ automatically so consistency is maintained across all temperatures, and in equilibrium mode the kinetic closure disappears entirely (the fast-reaction limit).

The activity model and state dependency layers are orthogonal: they extend any of the four classes without changing stoichiometry, equilibrium relation, or kinetic closure.
The remaining chapters of this part each extend exactly one layer.

---

The kinetic and equilibrium frameworks are now complete for homogeneous reactions with ideal activities.
An optional next chapter extends the solver to treat temperature as a dynamic state, coupling the energy balance to the reaction kinetics (@implementation-energy-balance).
The following chapter extends the activity term $a_i = \gamma_i c_i / c^\circ$ to non-ideal solutions, where the activity coefficient $\gamma_i \neq 1$ shifts the apparent equilibrium composition (@implementation-activity).
