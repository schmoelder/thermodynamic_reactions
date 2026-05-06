---
kernelspec:
  name: python3
  display_name: Python 3 (ipykernel)
---

(reaction-gibbs-energy)=
# Reaction Gibbs Energy and Driving Force

The extent of reaction $\xi$ (@reaction-coordinates) parameterises the composition.
Since all $n_i$ are determined by $\xi$, the Gibbs energy becomes a function $G(\xi)$ at fixed $T$ and $P$, and the question of which direction a reaction proceeds becomes a calculus problem.


## Reaction Gibbs energy

The Gibbs energy at composition $\xi$ is:
$$
G(\xi) = \sum_i \mu_i(\xi)\, n_i(\xi).
$$

The thermodynamic driving force is the slope of this curve:
$$
\Delta_r G \equiv \frac{dG}{d\xi} = \sum_i \nu_i \mu_i.
$$

The sign of $\Delta_r G$ determines the direction of spontaneous change.
$\Delta_r G < 0$ drives the reaction forward; $\Delta_r G > 0$ drives it in reverse.
Equilibrium is the stationary point where $\Delta_r G = 0$.


## Geometric interpretation

For $\mathrm{A} \rightleftharpoons \mathrm{B}$, taking $V = 1/c^\circ$ with $n_{\mathrm{A},0} = 1\,\mathrm{mol}$ so that $a_i = n_i/(V c^\circ) = n_i/\mathrm{mol}$, inserting
$a_\mathrm{A} = 1 - \xi$ and $a_\mathrm{B} = \xi$ into $G = \sum_i \mu_i n_i$ gives:

$$
\frac{G(\xi) - \mu_\mathrm{A}^\circ}{RT}
= \frac{\Delta_r G^\circ}{RT}\,\xi
+ (1-\xi)\ln(1-\xi) + \xi\ln\xi
$$

This separates into a linear term controlled by $\Delta_r G^\circ$ and a nonlinear mixing contribution from the entropy of mixing.
The derivative of $G(\xi)$ reproduces the reaction Gibbs energy, as shown in @fig-g-xi.

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-g-xi

import numpy as np
import matplotlib.pyplot as plt

xi = np.linspace(0.001, 0.999, 600)
delta = -1.0

def g_norm(xi, delta):
    return delta * xi + (1 - xi) * np.log(1 - xi) + xi * np.log(xi)

def drg_norm(xi, delta):
    return delta + np.log(xi / (1 - xi))

def xi_eq(delta):
    K = np.exp(-delta)
    return K / (1 + K)

def add_tangent(ax, xi_pt, delta, half_width=0.10, color="mediumpurple"):
    g_pt = g_norm(xi_pt, delta)
    slope = drg_norm(xi_pt, delta)
    x_tan = np.array([xi_pt - half_width, xi_pt + half_width])
    y_tan = g_pt + slope * (x_tan - xi_pt)
    ax.plot(x_tan, y_tan, color=color, linewidth=1.8, zorder=4)
    ax.plot(xi_pt, g_pt, "o", color="white", markersize=7, zorder=5,
            markeredgecolor="#1c4f8a", markeredgewidth=1.5)

g = g_norm(xi, delta)
xi_e = xi_eq(delta)
g_e = g_norm(xi_e, delta)
xi_left = 0.22
xi_right = 0.95

fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(xi, g, color="#1c4f8a", linewidth=2.5, zorder=3)
ax.plot([0.0, 1.04], [0.0, 0.0], "--", color="gray", linewidth=1.0)
ax.plot([0.0, 1.04], [delta, delta], "--", color="gray", linewidth=1.0)

add_tangent(ax, xi_left, delta)
add_tangent(ax, xi_right, delta)

ax.plot([xi_e - 0.10, xi_e + 0.10], [g_e, g_e],
        color="mediumpurple", linewidth=1.8, zorder=4)
ax.plot(xi_e, g_e, "o", color="white", markersize=10, zorder=6,
        markeredgecolor="#1c4f8a", markeredgewidth=2)
ax.axvline(xi_e, color="gray", linestyle=":", linewidth=1.0)

ax.set_xlabel(r"Extent of reaction, $\xi$")
ax.set_ylabel("Gibbs energy, $G$")
ax.set_xticks([])

fig.tight_layout()
```

```{figure} #cell-g-xi
:name: fig-g-xi

Gibbs energy $G(\xi)$ for $\mathrm{A} \rightleftharpoons \mathrm{B}$ with $\Delta_r G^\circ < 0$.
Tangent lines at two off-equilibrium compositions have nonzero slope ($\Delta_r G \ne 0$); the tangent at the minimum is horizontal ($\Delta_r G = 0$).
```

```{admonition} Geometric interpretation
:class: note
$G(\xi)$ defines an energy landscape along the reaction coordinate: the equilibrium state corresponds to the minimum, where $\Delta_r G = 0$.
Away from the minimum, the slope $\mathrm{d}G/\mathrm{d}\xi = \Delta_r G$ determines the direction of spontaneous change: negative values drive the reaction forward, positive values drive it backward.
The standard Gibbs energy change $\Delta_r G^\circ$ sets the tilt of the landscape and therefore the position of the minimum through the equilibrium constant $K$.
The mixing entropy contribution introduces curvature, ensuring that for reversible reactions with finite $K$, the equilibrium state lies between the pure reactants and pure products.
```


## Reaction quotient

Substituting $\mu_i = \mu_i^\circ + RT \ln a_i$ into $\Delta_r G = \sum_i \nu_i \mu_i$:
$$
\Delta_r G
= \sum_i \nu_i \mu_i^\circ + RT \sum_i \nu_i \ln a_i
= \Delta_r G^\circ + RT \ln Q
$$

where the **reaction quotient** is:
$$
Q = \prod_i a_i^{\nu_i}
$$

$Q$ measures the instantaneous composition.
When $Q < K$ the system is below equilibrium and $\Delta_r G < 0$; when $Q > K$ the reaction has overshot and $\Delta_r G > 0$.

## Le Chatelier's principle

The $Q$/$K$ framework gives a quantitative account of Le Chatelier's principle (@thermodynamic-potentials).
Any perturbation to an intensive variable drives the system toward a new equilibrium that partially opposes the change.
The three conjugate pairs of the fundamental relation each produce one class of perturbation.

**Composition ($\mu_i$, $n_i$):** a perturbation in species amount shifts $Q$ but not $K$.
Adding a reactant reduces $Q$ below $K$, making $\Delta_r G < 0$ and driving the reaction forward; adding a product raises $Q$ above $K$ and reverses it.
Removing a species has the opposite effect.

**Temperature ($T$, $S$):** a temperature change shifts $K$ itself, not $Q$.
For an exothermic reaction ($\Delta_r H^\circ < 0$), raising $T$ decreases $K$, shifting equilibrium toward reactants; the system absorbs heat to partially oppose the imposed temperature rise.
The quantitative form is the van't Hoff equation (@equilibrium-temperature).

**Pressure ($P$, $V$):** a pressure increase shifts equilibrium toward the side with fewer moles of gas, reducing total volume and opposing the compression.
Reactions in solution are nearly pressure-independent because partial molar volumes $\bar{V}_i$ are small.

---

The equilibrium condition $\Delta_r G = 0$ defines the equilibrium constant $K$ and fixes the composition at which the reaction stops.
The next chapter derives $K$ and connects it to tabulated thermodynamic data.
