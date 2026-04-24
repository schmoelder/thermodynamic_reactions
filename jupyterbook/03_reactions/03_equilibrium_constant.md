---
kernelspec:
  name: python3
  display_name: Python 3 (ipykernel)
---

(equilibrium)=
# Chemical Equilibrium and the Equilibrium Constant

The reaction Gibbs energy $\Delta_r G = \Delta_r G^\circ + RT \ln Q$ (@reaction-gibbs-energy) vanishes at the composition where the reaction stops.
That condition defines the equilibrium constant and connects the observable equilibrium position to standard thermodynamic data.


## Equilibrium condition

At equilibrium $\Delta_r G = 0$, so:
$$
\Delta_r G^\circ = -RT \ln K
$$

where the **equilibrium constant** is the value of $Q$ at equilibrium:
$$
K = Q\big|_{\mathrm{eq}} = \prod_i a_i^{\nu_i}\bigg|_{\mathrm{eq}}
$$

For $\mathrm{A} + 2\mathrm{B} \rightleftharpoons \mathrm{C}$:
$$
K = \frac{a_\mathrm{C}}{a_\mathrm{A} a_\mathrm{B}^2}
$$

In dilute ideal solution $a_i = c_i/c^\circ$, recovering the familiar concentration form.

```{admonition} Sign and direction
:class: note

$\Delta_r G^\circ < 0$ gives $K > 1$: products are favoured at equilibrium.
$\Delta_r G^\circ > 0$ gives $K < 1$: reactants are favoured.
$\Delta_r G^\circ$ fixes the equilibrium position; $\Delta_r G$ determines the direction of evolution at any given composition.
```


## Standard thermodynamic data

$\Delta_r G^\circ$ is computed from tabulated formation data:
$$
\Delta_r G^\circ = \sum_i \nu_i \Delta_f G^\circ_i
= \sum_i \nu_i \Delta_f H^\circ_i - T \sum_i \nu_i S^\circ_i
$$

The standard enthalpy of formation $\Delta_f H^\circ_i$ and the standard molar entropy $S^\circ_i$ are tabulated at $T = 298.15\,\mathrm{K}$ and $P^\circ = 1\,\mathrm{bar}$.
The entropy values are absolute because the third law provides a reference: $S \to 0$ as $T \to 0$.

This gives $K$ at 298.15 K.
The temperature dependence of $K$ is developed in @equilibrium-temperature.

---

Reactions are not new physics: $K$ is the composition at which $\mu$-equality is satisfied.
The equilibrium constant encodes where the system must end up; it says nothing about how fast it gets there.
The next chapter applies $K$ immediately to a practical case: adsorption of molecules onto a surface.
