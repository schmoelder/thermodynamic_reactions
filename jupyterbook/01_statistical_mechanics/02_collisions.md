---
kernelspec:
  name: python3
  display_name: Python 3 (ipykernel)
---

(collisions)=
# Random Collisions of Particles

Elastic collisions are the mechanism by which a gas explores its available microstates.
This chapter traces what happens to the energy distribution as collisions accumulate, starting from the most ordered possible initial condition.

## Starting point: a perfectly ordered gas

To make the question concrete, start from the most ordered initial state imaginable: $N$ particles all moving at exactly the same speed $v_0$, half to the left and half to the right.
Every particle carries identical kinetic energy $\varepsilon_0 = \frac{1}{2}mv_0^2$.

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-initial-state

import numpy as np
import matplotlib.pyplot as plt
from reactions.plots import setup_figure, COLORS

rng = np.random.default_rng(42)

N = 1000
v0 = 1.0
v = np.full(N, v0, dtype=float)
v[: N // 2] *= -1

fig, axes = setup_figure(1, 2)

axes[0].hist(v, bins=30, color="C0", edgecolor="white")
axes[0].set_xlabel("velocity $v$")
axes[0].set_ylabel("count")
axes[0].set_title("Initial state: all particles at $|v| = v_0$")

KE = 0.5 * v**2
axes[1].hist(KE, bins=30, color="C1", edgecolor="white")
axes[1].set_xlabel(r"kinetic energy $\varepsilon = \frac{1}{2}mv^2$  [a.u.]")
axes[1].set_ylabel("count")
axes[1].set_title("Initial energy distribution: a single spike")

fig.tight_layout()
```

```{figure} #cell-initial-state
:name: fig-initial-state

Initial state of $N = 1000$ particles all moving at the same speed $v_0$.
Left: the velocity distribution is two spikes at $\pm v_0$.
Right: the energy distribution is a single spike at $\varepsilon_0 = \frac{1}{2}mv_0^2$.
This is not a stable configuration: collisions will immediately begin spreading energy away from this spike.
```

As shown in @fig-initial-state, this is not a stable situation.
In any real gas, particles interact: they collide, exchange energy, and continuously move from one microstate to another.

## Elastic collisions

Each elastic collision conserves total momentum and kinetic energy: $\varepsilon_1' + \varepsilon_2' = \varepsilon_1 + \varepsilon_2$.
The individual energies are free to change, so each collision genuinely redistributes energy between the two particles.
The model randomises the direction of the relative velocity uniformly over all angles while keeping its magnitude fixed.
This preserves both conservation laws exactly: the centre-of-mass velocity $\mathbf{v}_\text{cm} = (\mathbf{v}_1 + \mathbf{v}_2)/2$ is unchanged (momentum conserved, component by component), and the relative speed $|\mathbf{v}_1 - \mathbf{v}_2|$ is unchanged (kinetic energy conserved).
Repeated many times across all particle pairs, this gradually spreads the initial spike across a broad range of energies.

Collisions immediately begin spreading energy away from the initial spike.
After many collisions the system reaches a **steady state** in which the distribution no longer changes on average, even though individual particles continue to exchange energy.
@fig-collisions shows the transition to a steady state.

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-collisions

import numpy as np
import matplotlib.pyplot as plt
from reactions.plots import setup_figure, COLORS
from scipy import stats


def simulate_collisions_with_snapshots(
    N: int,
    v0: float,
    step_counts: list[int],
    rng: np.random.Generator,
) -> dict[int, np.ndarray]:
    """Run an elastic-collision simulation and return energy snapshots.

    Particles start with 3-D velocities of magnitude v0 in random directions.
    Each collision randomises the relative-velocity direction while conserving
    both total momentum (centre-of-mass velocity unchanged) and kinetic energy
    (relative speed unchanged).

    Returns a dict mapping each value in step_counts to an array of per-particle
    kinetic energies (½|v|²) captured after exactly that many collisions.
    """
    directions = rng.normal(size=(N, 3))
    directions /= np.linalg.norm(directions, axis=1, keepdims=True)
    v = directions * v0
    snapshots: dict[int, np.ndarray] = {}
    done = 0
    for target in sorted(step_counts):
        for _ in range(target - done):
            i, j = rng.choice(N, size=2, replace=False)
            v_cm = (v[i] + v[j]) / 2
            v_rel_mag = np.linalg.norm(v[i] - v[j])
            new_dir = rng.normal(size=3)
            new_dir /= np.linalg.norm(new_dir)
            v[i] = v_cm + 0.5 * v_rel_mag * new_dir
            v[j] = v_cm - 0.5 * v_rel_mag * new_dir
        done = target
        snapshots[target] = 0.5 * np.sum(v**2, axis=1).copy()
    return snapshots


rng = np.random.default_rng(42)
N = 2000
v0 = 1.0
step_counts = [0, 600, 4_000, 50_000]

snapshots = simulate_collisions_with_snapshots(N, v0, step_counts, rng)

kT = v0**2 / 3
eps_max = 3.0
eps_range = np.linspace(0, eps_max, 300)
mb_pdf = stats.gamma(a=3 / 2, scale=kT).pdf(eps_range)

y_max = mb_pdf.max() * 7
labels = ["0 collisions", "600 collisions", "4 000 collisions", "50 000 collisions"]

fig, axes = setup_figure(1, 4, sharey=True)

for ax, n, label in zip(axes, step_counts, labels):
    ax.hist(
        snapshots[n],
        bins=31,
        range=(0, eps_max),
        density=True,
        color="steelblue",
        alpha=0.85,
        edgecolor="white",
    )
    ax.fill_between(eps_range, mb_pdf, alpha=0.4, color="firebrick")
    ax.plot(
        eps_range, mb_pdf, color="firebrick", linewidth=2, label="Maxwell-Boltzmann"
    )
    ax.set_xlabel(r"kinetic energy $\varepsilon$ [a.u.]")
    ax.set_ylim(0, y_max)
    ax.set_xlim(0, eps_max)
    ax.set_title(label)

axes[0].set_ylabel("probability density")
axes[-1].legend(loc="upper right")
fig.tight_layout()
```

```{figure} #cell-collisions
:name: fig-collisions

Evolution of the energy distribution under elastic collisions ($N = 2000$ particles, initial speed $v_0$).
The distribution begins as a spike at $\varepsilon = \frac{1}{2}v_0^2$ and converges to the Maxwell-Boltzmann shape (red) as collisions accumulate.
The initial spike clips the y-axis; the scale is set to show the equilibrium distribution clearly.
```

The simulation reveals the equilibrium distribution empirically: regardless of where the system starts, collisions always drive it to the same smooth, asymmetric shape.
Most particles carry modest energy; very energetic ones are rare, but the tail never reaches zero.

The starting configuration does not matter: any initial distribution converges to the same outcome.
Why is this particular shape the unique attractor?
That is what the next chapters answer: the next chapter introduces a way to count microscopic arrangements, and @maxwell-boltzmann uses it to show that the MB distribution is the only one that maximises entropy under fixed energy and particle number.

