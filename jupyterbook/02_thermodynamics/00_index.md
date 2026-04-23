---
kernelspec:
  name: python3
  display_name: Python 3 (ipykernel)
---

(thermodynamics)=
# Part 2: Thermodynamic Fundamentals

The ideal gas law $PV = nRT$ is the empirical entry point.
From it, the laws of thermodynamics establish the constraints on how energy flows and how equilibrium is reached.
Combining the first and second laws yields the fundamental relation $dU = TdS - PdV + \sum_i \mu_i\, dn_i$, from which the Gibbs energy $G$ emerges as the natural potential at constant temperature and pressure: at equilibrium, $G$ is minimised.
This part can be read independently of Part 1; @ideal-gas introduces $PV = nRT$ from historical experiments, and the statistical foundation of entropy is not required.

The central quantity is the **chemical potential** $\mu_i = \partial G / \partial n_i$: the Gibbs energy change per mole of species $i$ added at fixed $T$, $P$, and composition.
For an ideal solution, $\mu_i$ depends logarithmically on concentration.
Real solutions require the use of activity $a_i = \gamma_i c_i / c^\circ $, where $c^\circ$ is the standard concentration.
The activity coefficient $\gamma_i$ accounts for non-ideal interactions, approaching 1 as $c_i \to 0$.
In @reactions, $\mu_i$ is used directly to define reaction driving forces and equilibrium constants, ensuring thermodynamic consistency in kinetic models.


```{admonition} Key results
:class: note

**Second law of thermodynamics.**
The entropy of an isolated system never decreases:

$$dS \geq 0$$

At constant $T$ and $P$, maximising the entropy of the universe is equivalent to minimising the Gibbs energy of the system.
Equilibrium is the state where $G$ is minimal and $dG = 0$ for all possible changes.

**Gibbs energy and spontaneity.**
The Gibbs energy $G = H - TS$ determines the direction of change at constant $T$ and $P$:

$$\Delta G = \Delta H - T\Delta S < 0 \quad \text{(spontaneous)}$$

**Chemical potential.**
The chemical potential is the partial molar Gibbs energy:

$$\mu_i = \left.\frac{\partial G}{\partial n_i}\right|_{T,P,n_{j \neq i}}$$

Equilibrium requires $\mu_i$ to be uniform across all phases and equal for all species connected by a reaction.

**Concentration dependence (ideal solution).**
For an ideal solution, $\mu_i$ depends logarithmically on concentration:

$$\mu_i = \mu_i^\circ + RT \ln \frac{c_i}{c^\circ}$$

**Activity.**
Real solutions replace concentration with activity $a_i = \gamma_i c_i / c^\circ$, where $\gamma_i \to 1$ as $c_i \to 0$:

$$\mu_i = \mu_i^\circ(T) + RT \ln a_i$$
```
