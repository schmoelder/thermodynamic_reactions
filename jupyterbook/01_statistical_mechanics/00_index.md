---
kernelspec:
  name: python3
  display_name: Python 3 (ipykernel)
---

(statistical-mechanics)=
# Part 1: From Statistical Mechanics to Thermodynamics

The central objects are particles and energy.
Starting from a box of identical particles with fixed total energy, the question is: how is that energy distributed among the particles at equilibrium?
The answer requires a way to count how many microscopic arrangements are consistent with a given macroscopic distribution.
That counting function is entropy.

From entropy, temperature emerges as the equilibrium condition between two systems.
Maximising entropy subject to fixed energy and particle number yields the Boltzmann factor: the probability that a particle occupies a state of energy $\varepsilon$ is proportional to $e^{-\varepsilon / k_B T}$.
Equipartition follows directly, and pressure, derived from molecular impacts on a wall, gives the ideal gas law.
This part is optional; Part 2 begins from the ideal gas law as an empirical entry point.

```{admonition} Key results
:class: tip

**Boltzmann factor.**
In thermal equilibrium at temperature $T$, the probability of a particle being in a state with energy $\varepsilon$ is

$$p(\varepsilon) \propto e^{-\varepsilon / k_B T}.$$

Higher energy states are exponentially less likely; the energy scale is set by $k_B T$, where $k_B = 1.380\times10^{-23}$ J/K converts the molecular energy scale to the macroscopic temperature scale, just as $N_A$ converts particle counts to moles.

**Equipartition.**
Each translational degree of freedom carries average kinetic energy $\frac{1}{2}k_BT$.
For a particle moving in 3D: $\langle \varepsilon \rangle = \frac{3}{2}k_BT$.

**Ideal gas law.**
For $N$ non-interacting particles in volume $V$:

$$PV = Nk_BT = nRT,$$

where $n$ is the amount in moles and $R = N_A k_B = 8.314$ J/(mol K) is the gas constant.
```

The ideal gas law is the first macroscopic expression of the chemical potential for non-interacting particles.
In Part 2 it re-emerges as the defining property of the ideal reference state for $\mu_i$, and the logarithmic concentration dependence $\mu_i = \mu_i^\circ + RT \ln (c_i / c^\circ)$ follows directly from it.
