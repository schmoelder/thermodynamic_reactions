---
kernelspec:
  name: python3
  display_name: Python 3 (ipykernel)
---

(maxwell-boltzmann)=
# The Maxwell-Boltzmann Distribution

The previous chapter established that the equilibrium distribution is the one that maximizes $\mathcal{H} = \ln \Omega$.
There are two complementary routes to finding it: a symmetry argument due to Maxwell, and a direct entropy maximisation.
Both arrive at the same answer.

## Maxwell's argument: isotropy and independence

Maxwell arrived at the equilibrium distribution in 1860 by symmetry alone, without any reference to entropy.

Consider the three velocity components $v_x$, $v_y$, $v_z$ of a single particle.
In an isotropic gas with no preferred direction, the joint distribution of $(v_x, v_y, v_z)$ can depend only on the speed $v = \sqrt{v_x^2 + v_y^2 + v_z^2}$, not on direction.
At the same time, there is no physical coupling between different directions: the components should be **statistically independent**, so the joint distribution factorises:

$$
f(v_x, v_y, v_z) = g(v_x)\, g(v_y)\, g(v_z)
$$

for some single-variable distribution $g$.
Isotropy and independence together determine $g$ uniquely.
The only function satisfying $g(v_x)\, g(v_y)\, g(v_z) = \phi(v_x^2 + v_y^2 + v_z^2)$ is the **Gaussian**:

$$
g(v_x) \propto \exp\!\left(-\frac{v_x^2}{2\sigma^2}\right)
$$

So the equilibrium velocity in each direction is normally distributed with mean zero and variance $\sigma^2$.

**From velocity components to speed.**
With three independent Gaussian components, converting to the speed $v$ accounts for the surface area of the sphere of radius $v$ in velocity space ($4\pi v^2$):

$$
p(v) \propto v^2\, \exp\!\left(-\frac{v^2}{2\sigma^2}\right)
$$

This is the **Maxwell speed distribution** (see @fig-mb-derivation, center panel).

**From speed to energy.**
Setting $\varepsilon = \frac{1}{2}mv^2$, so $v = \sqrt{2\varepsilon/m}$, the Jacobian $dv/d\varepsilon = 1/(mv)$ gives:

$$
p(\varepsilon) \propto v^2\, e^{-v^2/2\sigma^2} \cdot \frac{1}{mv}
= \frac{v}{m}\, e^{-v^2/2\sigma^2}
\propto \varepsilon^{1/2}\, e^{-\varepsilon / m\sigma^2}
$$

The $\sqrt{\varepsilon}$ prefactor is purely geometric: it counts how much velocity-space volume is available at energy $\varepsilon$.
The result is a **Gamma distribution** with shape $\frac{3}{2}$ and scale $m\sigma^2$.

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-mb-derivation

import numpy as np
import matplotlib.pyplot as plt
from reactions.plots import setup_figure, COLORS
from scipy import stats

sigma = 1.0  # kT/m = 1 in arbitrary units

fig, axes = setup_figure(1, 3, figsize=(12, 3.5))

v_x = np.linspace(-4, 4, 400)
axes[0].plot(v_x, stats.norm(scale=sigma).pdf(v_x), "C0")
axes[0].set_xlabel(r"$v_x$  [a.u.]")
axes[0].set_ylabel("probability density")
axes[0].set_title(r"component: $\mathcal{N}(0,\,\sigma^2)$")

v = np.linspace(0, 5, 400)
axes[1].plot(v, stats.maxwell(scale=sigma).pdf(v), "C1")
axes[1].set_xlabel(r"speed $v$  [a.u.]")
axes[1].set_title(r"speed: $p(v) \propto v^2\,e^{-v^2/2\sigma^2}$")

eps = np.linspace(0, 8, 400)
axes[2].plot(eps, stats.gamma(a=3 / 2, scale=sigma**2).pdf(eps), "C2")
axes[2].set_xlabel(r"energy $\varepsilon$  [a.u.]")
axes[2].set_title(
    r"energy: $p(\varepsilon) \propto \varepsilon^{1/2}\,e^{-\varepsilon/m\sigma^2}$"
)

for ax in axes[1:]:
    ax.set_ylabel("")

plt.suptitle(
    r"From Gaussian velocity components to the Maxwell–Boltzmann distribution", y=1.02
)
fig.tight_layout()
```

```{figure} #cell-mb-derivation
:name: fig-mb-derivation

The three-step transformation from Gaussian velocity components to the Maxwell–Boltzmann energy distribution.
Left: each component $v_x$ is normally distributed.
Centre: the speed $v = |\mathbf{v}|$ picks up a $v^2$ factor from the spherical shell volume, giving an asymmetric distribution peaked away from zero.
Right: the energy $\varepsilon = \frac{1}{2}mv^2$ acquires the $\sqrt{\varepsilon}$ Jacobian factor, yielding a Gamma$(3/2,\, m\sigma^2)$ distribution.
```

## Maximising entropy with constraints

The Lagrange multiplier approach confirms and extends Maxwell's result: it identifies $\sigma^2$ and shows that the Gaussian is the *unique* maximum-entropy distribution.

The previous chapter showed that $\mathcal{H} = -N\sum_i p_i \ln p_i$ where $p_i = n_i/N$.
Since $N$ is fixed, maximising $\mathcal{H}$ reduces to maximising $-\sum_i p_i \ln p_i$ subject to two constraints: normalisation $\sum_i p_i = 1$ and fixed mean energy $\sum_i p_i \varepsilon_i = \langle\varepsilon\rangle$.
One Lagrange multiplier per constraint gives:

$$
\mathcal{L} = -\sum_i p_i \ln p_i
  - \alpha\!\left(\sum_i p_i - 1\right)
  - \beta\!\left(\sum_i p_i \varepsilon_i - \langle\varepsilon\rangle\right)
$$

Differentiating with respect to $p_j$ and setting to zero:

$$
\frac{\partial \mathcal{L}}{\partial p_j} = 0
\quad\Longrightarrow\quad
-\ln p_j - 1 - \alpha - \beta \varepsilon_j = 0
\quad\Longrightarrow\quad
p_j \propto e^{-\beta \varepsilon_j}
$$

The equilibrium probability falls off **exponentially** with energy.
The multiplier $\alpha$ is fixed by $\sum p_i = 1$; $\beta$ remains to be identified.

Using the result from the previous chapter, $\partial S/\partial U = 1/T$ with $S = k_B \mathcal{H}$, and evaluating this derivative at the entropy maximum gives $\beta = 1/k_BT$, so $m\sigma^2 = k_BT$, i.e. $\sigma^2 = k_BT/m$.
The energy distribution is:

$$
p(\varepsilon) \propto \varepsilon^{1/2}\, e^{-\varepsilon / k_B T}
$$

a Gamma distribution with shape $\frac{3}{2}$ and scale $k_BT$, consistent with what both approaches give.

## The Boltzmann factor

The result $n_j \propto e^{-\varepsilon_j / k_B T}$ is the **Boltzmann factor** (@fig-mb).
It is one of the central results of statistical mechanics: in any system in thermal equilibrium at temperature $T$, states with higher energy are exponentially less probable, with the energy scale set by $k_BT$.

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-mb

import numpy as np
import matplotlib.pyplot as plt
from reactions.plots import setup_figure, COLORS
from scipy import stats

x_max = 8.0
eps = np.linspace(0, x_max, 500)
kT_vals = [0.5, 1.0, 2.0, 3.0]
colors = ["C0", "C1", "C2", "C3"]

fig, ax = setup_figure()

for kT, color in zip(kT_vals, colors):
    pdf = stats.gamma(a=3 / 2, scale=kT).pdf(eps)
    ax.plot(
        eps,
        pdf,
        color=color,
        linewidth=2.5,
        label=f"$k_BT = {kT}$",
    )

ax.set_xlabel(r"energy $\varepsilon$ [a.u.]")
ax.set_ylabel("probability density")
ax.set_xlim(0, x_max)
ax.set_ylim(bottom=0)
ax.legend()
fig.tight_layout()
```

```{figure} #cell-mb
:name: fig-mb

Maxwell-Boltzmann energy distribution at four temperatures.
Higher $k_BT$ shifts the peak to larger energies and flattens the distribution, while the decay rate $1/k_BT$ decreases, reflecting the Lagrange multiplier on the energy constraint.
```

## Consequences of the distribution

**Mean energy.**
The mean kinetic energy per particle is:

$$
\langle \varepsilon \rangle = \int_0^\infty \varepsilon\, p(\varepsilon)\, d\varepsilon = \frac{3}{2} k_B T
$$

This is the equipartition theorem for a monatomic ideal gas: each of the three translational degrees of freedom contributes $\frac{1}{2} k_B T$.
It also confirms $\sigma^2 = k_BT/m$, consistent with Maxwell's argument.

**Characteristic speeds.**
The speed distribution has three natural landmarks that appear frequently in textbooks and in later derivations:

$$
v_p = \sqrt{\frac{2k_BT}{m}}, \qquad
\langle v \rangle = \sqrt{\frac{8k_BT}{\pi m}}, \qquad
v_\text{rms} = \sqrt{\frac{3k_BT}{m}}
$$

$v_p$ is the mode (peak of $p(v)$); $\langle v \rangle$ is the mean; $v_\text{rms}$ is the root-mean-square speed, whose square links directly to the mean energy: $\frac{1}{2}mv_\text{rms}^2 = \langle \varepsilon \rangle = \frac{3}{2}k_BT$.
All three scale as $\sqrt{T}$, so doubling the temperature increases every characteristic speed by $\sqrt{2}$.
The corresponding energy distribution has mode $\varepsilon_\text{mode} = \frac{1}{2}k_BT$ and mean $\langle \varepsilon \rangle = \frac{3}{2}k_BT$ (@fig-mb-landmarks).

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-mb-landmarks

import numpy as np
import matplotlib.pyplot as plt
from reactions.plots import setup_figure, COLORS
from scipy import stats

kT = 1.0
m = 1.0
sigma = np.sqrt(kT / m)

v_p = np.sqrt(2 * kT / m)
v_mean = np.sqrt(8 * kT / (np.pi * m))
v_rms = np.sqrt(3 * kT / m)

E_a = 3.0 * kT
eps_mode = 0.5 * kT
eps_mean = 1.5 * kT

fig, axes = setup_figure(1, 2)

v = np.linspace(0, 5, 400)
axes[0].plot(v, stats.maxwell(scale=sigma).pdf(v), "C0", lw=2)
for xv, col, lbl in [
    (v_p, "C1", r"$v_p = \sqrt{2k_BT/m}$"),
    (v_mean, "C2", r"$\langle v\rangle = \sqrt{8k_BT/\pi m}$"),
    (v_rms, "C3", r"$v_\mathrm{rms} = \sqrt{3k_BT/m}$"),
]:
    axes[0].axvline(xv, color=col, lw=1.5, label=lbl)
axes[0].set_xlabel(r"speed $v$  [a.u.]")
axes[0].set_ylabel("probability density")
axes[0].set_title("Speed distribution")
axes[0].legend(fontsize=9)

eps = np.linspace(0, 8, 400)
pdf_e = stats.gamma(a=3 / 2, scale=kT).pdf(eps)
axes[1].plot(eps, pdf_e, "C0", lw=2)
axes[1].axvline(
    eps_mode, color="C1", lw=1.5, label=r"$\varepsilon_\mathrm{mode} = k_BT/2$"
)
axes[1].axvline(
    eps_mean, color="C3", lw=1.5, label=r"$\langle\varepsilon\rangle = 3k_BT/2$"
)
axes[1].axvline(E_a, color="C4", lw=1.5, ls="--", label=r"$E_a$ (activation energy)")
axes[1].fill_between(
    eps, pdf_e, where=(eps >= E_a), alpha=0.25, color="C4", label="_nolegend_"
)
axes[1].set_xlabel(r"energy $\varepsilon$  [a.u.]")
axes[1].set_title("Energy distribution")
axes[1].legend(fontsize=9)

fig.tight_layout()
```

```{figure} #cell-mb-landmarks
:name: fig-mb-landmarks

Characteristic values of the Maxwell-Boltzmann distribution.
Left: the three standard speed markers $v_p < \langle v\rangle < v_\text{rms}$, all scaling as $\sqrt{T}$.
Right: the corresponding energy landmarks, plus an activation energy threshold $E_a$ (dashed); the shaded tail is the fraction of particles that can overcome it.
```

**Why this shape and not another.**
The Maxwell-Boltzmann distribution is the *unique* distribution that maximizes $\mathcal{H}$ subject to fixed $N$ and $U$.
Any other distribution has lower $\mathcal{H}$ and is driven toward this shape by collisions.
The simulation in @collisions is a direct demonstration: starting from a spike ($\mathcal{H} = 0$), the system explores microstates and converges to the maximum-entropy distribution.

**The high-energy tail.**
Even though most particles sit near $\langle \varepsilon \rangle$, the exponential decay is slow enough that there is always a population at energies far above the mean.
For a chemical reaction to occur, colliding particles typically need at least a minimum energy $E_a$, the **activation energy**, to overcome the barrier separating reactants from products.
The fraction of particles with $\varepsilon > E_a$ is proportional to $e^{-E_a / k_BT}$: small when $E_a \gg k_BT$, but never exactly zero, and growing rapidly with temperature.
The Arrhenius equation $k = A\, e^{-E_a/k_BT}$ captures this dependence and is derived in @kinetics.

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-mb-ea

import numpy as np
import matplotlib.pyplot as plt
from reactions.plots import setup_figure, COLORS
from scipy import stats
from scipy.special import gammaincc

kT = 1.0
x_max = 10.0
eps = np.linspace(0, x_max, 500)
pdf = stats.gamma(a=3 / 2, scale=kT).pdf(eps)

Ea_vals = [2.0, 4.0, 6.0]
colors_ea = ["C1", "firebrick", "C3"]

fig, ax = setup_figure()
ax.plot(eps, pdf, color="steelblue", linewidth=2.5, zorder=3)

for Ea, color in zip(Ea_vals, colors_ea):
    mask = eps >= Ea
    eps_tail = np.concatenate([[Ea], eps[mask]])
    pdf_tail = np.concatenate([[stats.gamma(a=3 / 2, scale=kT).pdf(Ea)], pdf[mask]])
    frac = gammaincc(3 / 2, Ea / kT)
    ax.fill_between(
        eps_tail,
        pdf_tail,
        alpha=0.35,
        color=color,
        label=rf"$E_a = {Ea:.0f}\,k_BT$  ({frac * 100:.1f}%)",
    )
    ax.axvline(Ea, color=color, linestyle="--", linewidth=1.2, alpha=0.7)

ax.set_xlabel(r"energy $\varepsilon$ [a.u.]  ($k_BT = 1$)")
ax.set_ylabel("probability density")
ax.set_xlim(0, x_max)
ax.set_ylim(bottom=0)
ax.legend(loc="upper right")
fig.tight_layout()
```

```{figure} #cell-mb-ea
:name: fig-mb-ea

Fraction of particles that can overcome an activation energy barrier $E_a$ (shaded tails, $k_BT = 1$).
At $E_a = 2\,k_BT$ roughly 26% of particles are reactive; at $E_a = 6\,k_BT$ the fraction drops to 0.4%.
The tail fraction $\propto e^{-E_a/k_BT}$: each additional $k_BT$ of barrier reduces the reactive population by a factor of $e$.
```
