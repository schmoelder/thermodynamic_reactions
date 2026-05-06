---
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

(implementation-rate-laws)=
# Rate Laws and the Mass Action Principle

The previous chapter introduced $\mathbf{f}_\text{react} = \mathbf{S}\boldsymbol{\varphi}$ as the general structure of the reaction source term.
This chapter works out what $\boldsymbol{\varphi}$ looks like for the simplest class of rate laws: mass action kinetics.

Mass action kinetics express each reaction flux as a product of concentrations raised to their kinetic orders.
For an elementary step the order equals the stoichiometric coefficient; for empirical rate laws it is a fitted parameter.
`MassActionReaction` implements this directly on concentrations, with $k_f$ and $k_r$ as independent inputs.

## Irreversible first-order: A→B

The simplest case is a single species converting irreversibly,

$$
\text{A} \xrightarrow{k_f} \text{B}, \qquad \varphi = k_f c_\text{A}.
$$

The source terms are $\dot{c}_\text{A} = -\varphi$ and $\dot{c}_\text{B} = +\varphi$, giving the analytical solution (@kinetics)

$$
c_\text{A}(t) = c_{\text{A},0}\,e^{-k_f t}, \qquad c_\text{B}(t) = c_{\text{A},0} - c_\text{A}(t).
$$

```{code-cell} ipython3
import numpy as np
from reactions.api import Component, MassActionReaction, ReactionModel
from reactions.solver import simulate

a = Component("A")
b = Component("B")

model = ReactionModel(
    components=[a, b],
    reactions=[MassActionReaction("A -> B", kf=1.0)],
)

result = simulate(
    model,
    c0={"A": 1000.0},
    t_span=(0, 6.0),
)
```

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-rate-irrev

import matplotlib.pyplot as plt

c_A_exact = 1000.0 * np.exp(-1.0 * result.t)
c_B_exact = 1000.0 - c_A_exact
print(f"max error A: {np.max(np.abs(result['A'] - c_A_exact)):.2e} mol/m³")
print(f"max error B: {np.max(np.abs(result['B'] - c_B_exact)):.2e} mol/m³")
print(f"conservation (A+B): {np.max(np.abs(result['A'] + result['B'] - 1000.0)):.2e} mol/m³")

fig, ax = plt.subplots()
ax.plot(result.t, result['A'], label="A (simulated)")
ax.plot(result.t, result['B'], label="B (simulated)")
ax.plot(result.t, c_A_exact, "k--", lw=1, label="A (analytical)")
ax.plot(result.t, c_B_exact, "k--", lw=1, label="B (analytical)")
ax.set_xlabel("time [s]")
ax.set_ylabel("concentration [mol/m³]")
ax.legend()
fig.tight_layout()
```

```{figure} #cell-rate-irrev
:name: fig-rate-irrev

Irreversible first-order reaction A→B with $k_f = 1\ \text{s}^{-1}$ and $c_{\text{A},0} = 1000\ \text{mol/m}^3$.
Dashed lines are the analytical solutions $c_\text{A}(t) = c_{\text{A},0}\,e^{-k_f t}$ and $c_\text{B}(t) = c_{\text{A},0}(1 - e^{-k_f t})$; simulated trajectories are indistinguishable at this scale.
```

## Consecutive reactions: A→B→C

When a product is itself a reactant, two reactions share the same `ReactionModel`.
B is an intermediate: it is produced by the first reaction and consumed by the second,

$$
\text{A} \xrightarrow{k_1} \text{B} \xrightarrow{k_2} \text{C},
\qquad \varphi_1 = k_1 c_\text{A}, \quad \varphi_2 = k_2 c_\text{B}.
$$

The source terms are $\dot{c}_\text{A} = -\varphi_1$, $\dot{c}_\text{B} = \varphi_1 - \varphi_2$, $\dot{c}_\text{C} = \varphi_2$.
For $k_1 \neq k_2$ the analytical solutions are

$$
c_\text{A}(t) = c_{\text{A},0}\,e^{-k_1 t}, \qquad
c_\text{B}(t) = c_{\text{A},0}\,\frac{k_1}{k_2 - k_1}\!\left(e^{-k_1 t} - e^{-k_2 t}\right).
$$

B rises from zero, reaches a maximum at $t^* = \ln(k_2/k_1)/(k_2 - k_1)$, then decays as C accumulates.

```{code-cell} ipython3
k1, k2 = 1.0, 3.0

model_consec = ReactionModel(
    components=[Component("A"), Component("B"), Component("C")],
    reactions=[
        MassActionReaction("A -> B", kf=k1),
        MassActionReaction("B -> C", kf=k2),
    ],
)

result_consec = simulate(
    model_consec,
    c0={"A": 1000.0},
    t_span=(0, 5.0),
)
```

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-rate-consec

c_A_exact_c = 1000.0 * np.exp(-k1 * result_consec.t)
c_B_exact_c = 1000.0 * k1 / (k2 - k1) * (np.exp(-k1 * result_consec.t) - np.exp(-k2 * result_consec.t))
c_C_exact_c = 1000.0 - c_A_exact_c - c_B_exact_c
print(f"max error A: {np.max(np.abs(result_consec['A'] - c_A_exact_c)):.2e} mol/m³")
print(f"max error B: {np.max(np.abs(result_consec['B'] - c_B_exact_c)):.2e} mol/m³")
print(f"max error C: {np.max(np.abs(result_consec['C'] - c_C_exact_c)):.2e} mol/m³")
print(f"conservation (A+B+C): {np.max(np.abs(result_consec['A'] + result_consec['B'] + result_consec['C'] - 1000.0)):.2e} mol/m³")

fig, ax = plt.subplots()
ax.plot(result_consec.t, result_consec['A'], label="A (simulated)")
ax.plot(result_consec.t, result_consec['B'], label="B (simulated)")
ax.plot(result_consec.t, result_consec['C'], label="C (simulated)")
ax.plot(result_consec.t, c_A_exact_c, "k--", lw=1, label="A (analytical)")
ax.plot(result_consec.t, c_B_exact_c, "k--", lw=1, label="B (analytical)")
ax.plot(result_consec.t, c_C_exact_c, "k--", lw=1, label="C (analytical)")
ax.set_xlabel("time [s]")
ax.set_ylabel("concentration [mol/m³]")
ax.legend()
fig.tight_layout()
```

```{figure} #cell-rate-consec
:name: fig-rate-consec

Consecutive first-order reactions A→B→C with $k_1 = 1\ \text{s}^{-1}$, $k_2 = 3\ \text{s}^{-1}$.
The intermediate B rises to a maximum at $t^* = \ln(k_2/k_1)/(k_2 - k_1) \approx 0.55\ \text{s}$ then decays as C accumulates.
Dashed lines are analytical solutions; $c_\text{C}$ follows by conservation $c_\text{C} = c_{\text{A},0} - c_\text{A} - c_\text{B}$.
```

## Second-order: A+B→C

When two species must collide, the flux is proportional to both concentrations,

$$
\text{A} + \text{B} \xrightarrow{k_f} \text{C}, \qquad \varphi = k_f c_\text{A} c_\text{B}.
$$

For equal initial concentrations $c_{\text{A},0} = c_{\text{B},0} = c_0$ the analytical solution is

$$
c_\text{A}(t) = c_\text{B}(t) = \frac{c_0}{1 + k_f c_0\, t}.
$$

```{code-cell} ipython3
kf = 1e-3   # m³ mol⁻¹ s⁻¹
c0 = 500.0  # mol/m³

model_2nd = ReactionModel(
    components=[Component("A"), Component("B"), Component("C")],
    reactions=[MassActionReaction("A + B -> C", kf=kf)],
)

result_2nd = simulate(
    model_2nd,
    c0={"A": c0, "B": c0},
    t_span=(0, 10.0 / (kf * c0)),
)
```

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-rate-2nd

c_A_exact = c0 / (1.0 + kf * c0 * result_2nd.t)
c_C_exact = c0 - c_A_exact
print(f"max error A: {np.max(np.abs(result_2nd['A'] - c_A_exact)):.2e} mol/m³")
print(f"max error B: {np.max(np.abs(result_2nd['B'] - c_A_exact)):.2e} mol/m³")
print(f"max error C: {np.max(np.abs(result_2nd['C'] - c_C_exact)):.2e} mol/m³")
print(f"conservation (A+C): {np.max(np.abs(result_2nd['A'] + result_2nd['C'] - c0)):.2e} mol/m³")

fig, ax = plt.subplots()
ax.plot(result_2nd.t, result_2nd['A'], label="A (simulated)")
ax.plot(result_2nd.t, result_2nd['B'], label="B (simulated)")
ax.plot(result_2nd.t, result_2nd['C'], label="C (simulated)")
ax.plot(result_2nd.t, c_A_exact, "k--", lw=1, label="A = B (analytical)")
ax.plot(result_2nd.t, c_C_exact, "k--", lw=1, label="C (analytical)")
ax.set_xlabel("time [s]")
ax.set_ylabel("concentration [mol/m³]")
ax.legend()
fig.tight_layout()
```

```{figure} #cell-rate-2nd
:name: fig-rate-2nd

Second-order irreversible reaction A+B→C with $k_f = 10^{-3}\ \text{m}^3\text{mol}^{-1}\text{s}^{-1}$ and equal initial concentrations $c_0 = 500\ \text{mol/m}^3$.
Because $c_{\text{A},0} = c_{\text{B},0}$, both reactants follow the same analytical solution $c_\text{A}(t) = c_\text{B}(t) = c_0 / (1 + k_f c_0 t)$; product C follows by conservation $c_\text{C} = c_0 - c_\text{A}$.
```

## Reversible: A⇌B

Adding $k_r > 0$ makes the reaction reversible,

$$
\text{A} \underset{k_r}{\stackrel{k_f}{\rightleftharpoons}} \text{B},
\qquad \varphi = k_f c_\text{A} - k_r c_\text{B}.
$$

At long times $\varphi \to 0$ and the system reaches equilibrium at
$c_{\text{B},\text{eq}} / c_{\text{A},\text{eq}} = k_f / k_r$.
The approach is exponential with relaxation time $\tau = 1/(k_f + k_r)$:

$$
c_\text{A}(t) = c_{\text{A},\text{eq}} + (c_{\text{A},0} - c_{\text{A},\text{eq}})\,e^{-(k_f+k_r)t},
\qquad c_\text{B}(t) = c_{\text{A},0} - c_\text{A}(t).
$$

```{code-cell} ipython3
kf, kr = 2.0, 0.5

model_rev = ReactionModel(
    components=[Component("A"), Component("B")],
    reactions=[MassActionReaction("A <-> B", kf=kf, kr=kr)],
)

result_rev = simulate(
    model_rev,
    c0={"A": 1000.0, "B": 0.0},
    t_span=(0, 5.0),
)
```

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-rate-rev

K    = kf / kr
A_eq = 1000.0 / (1 + K)
B_eq = 1000.0 - A_eq
c_A_exact = A_eq + (1000.0 - A_eq) * np.exp(-(kf + kr) * result_rev.t)
c_B_exact = 1000.0 - c_A_exact
print(f"K = kf/kr = {K:.1f}")
print(f"simulated  A_eq = {result_rev['A'][-1]:.2f} mol/m³")
print(f"analytical A_eq = {A_eq:.2f} mol/m³")
print(f"max error A: {np.max(np.abs(result_rev['A'] - c_A_exact)):.2e} mol/m³")
print(f"max error B: {np.max(np.abs(result_rev['B'] - c_B_exact)):.2e} mol/m³")

fig, ax = plt.subplots()
ax.plot(result_rev.t, result_rev['A'], label="A (simulated)")
ax.plot(result_rev.t, result_rev['B'], label="B (simulated)")
ax.plot(result_rev.t, c_A_exact, "k--", lw=1, label="A (analytical)")
ax.plot(result_rev.t, c_B_exact, "k--", lw=1, label="B (analytical)")
ax.axhline(A_eq, color="C0", lw=0.8, ls=":")
ax.axhline(B_eq, color="C1", lw=0.8, ls=":", label="equilibrium")
ax.set_xlabel("time [s]")
ax.set_ylabel("concentration [mol/m³]")
ax.legend()
fig.tight_layout()
```

```{figure} #cell-rate-rev
:name: fig-rate-rev

Reversible first-order reaction A⇌B with $k_f = 2\ \text{s}^{-1}$, $k_r = 0.5\ \text{s}^{-1}$ ($K = k_f/k_r = 4$).
Dashed lines are the analytical solutions; dotted lines mark the equilibrium concentrations $c_{\text{A},\text{eq}} = c_{\text{A},0}/(1+K)$ and $c_{\text{B},\text{eq}} = K\,c_{\text{A},0}/(1+K)$.
```

---

At fixed temperature `MassActionReaction` is exact: $k_f/k_r$ sets the equilibrium composition directly.
When temperature changes, $K(T)$ shifts with it; the next chapter introduces `ThermodynamicReaction`, which derives $k_r(T) = k_f(T)/K(T)$ at every evaluation so the long-time limit always tracks the correct equilibrium (@implementation-equilibrium).
