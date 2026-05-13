---
kernelspec:
  name: python3
  display_name: Python 3 (ipykernel)
---

(adsorption)=
# Adsorption Equilibrium

Adsorption is the partitioning of solute between a mobile phase and a stationary phase.
In chromatography, this partitioning determines retention: solute adsorbs onto the stationary phase, slows down relative to the mobile phase, then desorbs and continues through the column.
The process is a quasi-reaction,

$$\ce{A_{mobile} <=> A_{adsorbed}},$$

with equilibrium constant

$$K_\text{ads} = \exp\!\left(-\frac{\Delta_\text{ads} G^\circ}{RT}\right),$$

where $\Delta_\text{ads} G^\circ$ is the standard Gibbs energy of adsorption.
This is the same relation established in @equilibrium; adsorption is one more instance of $\Delta_r G^\circ = -RT \ln K$.


## The finite-site constraint

What distinguishes adsorption from a homogeneous reaction is a structural constraint: the stationary phase has a fixed total density of adsorption sites $q_\text{max}$.
At any loading, each site is either occupied by a solute molecule or free.
Defining the occupancy $\theta = q / q_\text{max}$, the constraint reads

$$\theta + \theta_\text{free} = 1,$$

where $\theta_\text{free} = 1 - \theta$ is the fraction of free sites.
This bound limits loading: as $\theta \to 1$, no free sites remain regardless of solution concentration.
This finite-site constraint has no counterpart in stoichiometric equilibria, where concentrations are free to range without an imposed ceiling.


## The Langmuir isotherm

{cite:t}`langmuir1918` derived this isotherm by treating adsorption as a surface reaction between mobile solute and a free site $\ce{S}$:

$$\ce{A_{mobile} + S_{free} <=> AS_{occupied}}$$

The equilibrium condition is $\mu_A + \mu_S = \mu_{AS}$.
For the mobile solute in ideal dilute solution, $\mu_A = \mu^\circ_A + RT \ln(c/c^\circ)$ (from @chemical-potential).
For the surface species, the chemical potential depends on the fraction of sites in each state: $\mu_S = \mu^\circ_S + RT \ln(1-\theta)$ and $\mu_{AS} = \mu^\circ_{AS} + RT \ln \theta$.

Combining these under the equilibrium condition gives the mass-action form of $Q = K$:

$$K_\text{ads} = \frac{\theta}{(1 - \theta)\,(c/c^\circ)}.$$

Solving for $\theta$ and substituting $q = q_\text{max} \theta$ yields the **Langmuir isotherm**:

$$q = \frac{q_\text{max}\, K_\text{ads}\, c}{1 + K_\text{ads}\, c},$$

where $c^\circ = 1\,\mathrm{mol/L}$ has been absorbed into $K_\text{ads}$.

At low loading ($K_\text{ads} c \ll 1$) the denominator approaches unity and the isotherm is linear: $q \approx q_\text{max} K_\text{ads} c \equiv Hc$, with $H = q_\text{max} K_\text{ads}$ the Henry coefficient.
At high loading ($K_\text{ads} c \gg 1$) the isotherm saturates at $q_\text{max}$, because all sites are occupied.
The saturation is an equilibrium consequence of the finite-site constraint, not a kinetic phenomenon; $q_\text{max}$ is a thermodynamic capacity.
@fig-langmuir shows the isotherm for two values of $K_\text{ads}$, with Henry tangents marking the linear limit.

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-langmuir

import numpy as np
import matplotlib.pyplot as plt

q_max = 1.0
c = np.linspace(0, 5.0, 500)
K_vals = [0.5, 2.0]
colors = ["#1c4f8a", "#b05c2e"]
labels = [
    r"$K_\mathrm{ads} = 0.5\ \mathrm{L\,mol^{-1}}$",
    r"$K_\mathrm{ads} = 2.0\ \mathrm{L\,mol^{-1}}$",
]

fig, ax = plt.subplots(figsize=(6, 4))

for K, color, lbl in zip(K_vals, colors, labels):
    q = q_max * K * c / (1.0 + K * c)
    ax.plot(c, q, color=color, linewidth=2.5, label=lbl, zorder=3)
    c_h = np.array([0.0, 0.75 / K])
    ax.plot(c_h, q_max * K * c_h, "--", color=color, linewidth=1.4, alpha=0.75, zorder=2)

ax.axhline(q_max, color="gray", linewidth=1.0, linestyle=":", zorder=1)
ax.text(0.12, q_max + 0.04, r"$q_\mathrm{max}$", fontsize=10, color="gray", va="bottom")

ax.set_xlabel(r"Concentration, $c$ / mol$\,$L$^{-1}$")
ax.set_ylabel(r"Loading, $q$")
ax.set_xlim(0, 5)
ax.set_ylim(0, 1.18)
ax.set_yticks([0, 0.5, 1.0])
ax.set_yticklabels(["$0$", r"$q_\mathrm{max}/2$", r"$q_\mathrm{max}$"])
ax.legend(fontsize=9, framealpha=0.9, loc="lower right")

fig.tight_layout()
```

```{figure} #cell-langmuir
:name: fig-langmuir

The Langmuir isotherm at two values of $K_\mathrm{ads}$.
At low loading each curve follows its Henry tangent $q = Hc$ (dashed), where $H = q_\mathrm{max} K_\mathrm{ads}$.
At high loading both saturate at $q_\mathrm{max}$ (dotted), regardless of concentration.
A larger $K_\mathrm{ads}$ shifts the saturation onset to lower concentrations.
```

```{admonition} Intuition
:class: tip

The Langmuir isotherm behaves like seats on a train.
In the Henry region the train is nearly empty: every additional passenger finds a seat, and the number seated is proportional to the number on the platform ($q \propto c$).
As loading increases, fewer seats remain and each additional molecule is less likely to adsorb.
At saturation the train is full: $q = q_\text{max}$ regardless of how many more molecules are in solution.
A large $K_\text{ads}$ means the train fills at lower platform concentration.
```

Adsorption sits between phase equilibrium (@phase-equilibrium) and reaction equilibrium (@equilibrium): like phase equilibrium it describes partitioning between two regions; like reaction equilibrium it is governed by a $K$ derived from $\Delta G^\circ$.
The occupancy constraint is what sets it apart from both.

The same finite-site logic extends to multiple competing species sharing one site pool (@multicomponent), and reappears in kinetics as Michaelis-Menten saturation when a finite pool of enzyme active sites is the catalytic resource (@saturation).

---

The next chapter generalises the single-component isotherm to competing species, ion exchange, and Donnan equilibrium, all instances of the same $\mu$-equality framework applied to a shared finite-site pool (@multicomponent).
