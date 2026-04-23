---
kernelspec:
  name: python3
  display_name: Python 3 (ipykernel)
---

(particles)=
# Particles in a Box

To build intuition for thermodynamics from the ground up, start with the simplest possible system:
$N$ identical particles confined to a box, each with some velocity $\mathbf{v}_i$.

```{figure} figures/translational_motion.gif
:name: fig-translational-motion
:width: 400px
:align: center

Helium atoms in a box at 1,950 atm, slowed down two trillion fold from room temperature.
The size of the atoms relative to their spacing is shown to scale.
The temperature of the gas is proportional to the average kinetic energy of the particles.
Source: [Wikimedia Commons](https://en.wikipedia.org/wiki/File:Translational_motion.gif), public domain.
```

Without considering interactions between particles, the only energy each particle carries is kinetic:

$$
\varepsilon_i = \frac{1}{2} m v_i^2
$$

The total internal energy of the system is then just the sum over all particles:

$$
U = \sum_{i=1}^{N} \varepsilon_i = \sum_{i=1}^{N} \frac{1}{2} m v_i^2
$$

This is the **internal energy**, the total microscopic kinetic energy stored in the motion of the particles.
It is a fixed quantity for an isolated system: the box has rigid walls, no heat flows in or out, and no work is done.
At any instant, the system is described by the full set of velocities $\{v_1, v_2, \ldots, v_N\}$, a point in a very high-dimensional space.
This full description is called a **microstate**.
For a system of $N = 10^{23}$ particles, tracking the microstate is hopeless in practice.
The question then is: given that $U$ is fixed, what is the *distribution* of velocities?
That question is what the next chapter begins to answer.

