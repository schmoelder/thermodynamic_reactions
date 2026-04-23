---
kernelspec:
  name: python3
  display_name: Python 3 (ipykernel)
---

(pressure)=
# Pressure and the Ideal Gas Law

## Pressure as momentum transfer

Pressure is force per unit area, and force is rate of momentum transfer.
So the pressure on a wall is the rate at which gas particles deposit momentum on it.

Consider a cubic box with side length $\ell$ and volume $V = \ell^3$.
Focus on one wall perpendicular to the $x$-axis.
A particle with $x$-velocity $v_x > 0$ hits the wall and bounces back elastically: its $x$-momentum changes from $mv_x$ to $-mv_x$, depositing $2mv_x$ of momentum per collision.

The particle travels distance $2\ell$ between successive hits on the same wall, so its collision frequency is $v_x / 2\ell$.
The force it exerts on the wall is therefore:

$$
f_i = 2mv_{x,i} \cdot \frac{v_{x,i}}{2\ell} = \frac{mv_{x,i}^2}{\ell}
$$

Summing over all $N$ particles and dividing by the wall area $\ell^2$:

$$
P = \frac{1}{\ell^2} \sum_{i=1}^N \frac{mv_{x,i}^2}{\ell} = \frac{Nm}{V}\langle v_x^2 \rangle
$$

## Connecting to temperature

By symmetry, motion is equally distributed among all three directions: $\langle v_x^2 \rangle = \frac{1}{3}\langle v^2 \rangle$.
From equipartition (@maxwell-boltzmann):

$$
\frac{1}{2}m\langle v^2 \rangle = \frac{3}{2}k_B T
\quad\Longrightarrow\quad
m\langle v_x^2 \rangle = k_B T
$$

Substituting:

$$
P = \frac{N k_B T}{V}
$$

Pressure arises purely from the kinetic energy of the particles.
No intermolecular forces, no potential energy: just particles bouncing off walls.
Computing $Nm\langle v_x^2\rangle/V$ directly from the equilibrium velocities of @collisions confirms $P = Nk_BT/V$ numerically.

## The ideal gas law

In chemistry it is more convenient to count in moles.
With $n = N/N_A$ moles and defining the **gas constant**:

$$
R = N_A k_B = 6.022 \times 10^{23} \times 1.380 \times 10^{-23} \;\text{J/K} = 8.314 \;\text{J/(mol K)}
$$

the pressure equation becomes the ideal gas law:

$$
PV = nRT
$$

$R$ is not a new constant; it is $k_B$ expressed per mole rather than per particle.
The law holds whenever particles do not interact: no intermolecular forces, no excluded volume.

The next chapter traces how the same result was discovered experimentally over two centuries, before anyone knew about molecules, and examines where the law begins to fail.

