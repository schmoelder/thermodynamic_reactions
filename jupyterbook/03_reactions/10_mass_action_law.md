---
kernelspec:
  name: python3
  display_name: Python 3 (ipykernel)
---

(mass-action-law)=
# The Mass Action Law

The rate laws in @kinetics describe individual reaction mechanisms.
The mass action law (MAL) provides a general rate expression for reversible reactions that is guaranteed to reproduce the correct thermodynamic equilibrium.

## Rate expression

For a reversible reaction $\sum_i \nu_i \text{A}_i = 0$, the MAL assigns exponents equal to the stoichiometric coefficients:

$$
r = k_f \prod_{\text{reactants}} c_i^{|\nu_i|} - k_r \prod_{\text{products}} c_j^{\nu_j}
$$

At equilibrium $r = 0$, so:

$$
\frac{k_f}{k_r} = \frac{\prod_{\text{products}} c_j^{\nu_j}}{\prod_{\text{reactants}} c_i^{|\nu_i|}}\bigg|_\text{eq} = K
$$

The ratio of rate constants equals the equilibrium constant.
This **thermodynamic consistency condition** is a constraint that any kinetic model must satisfy to be consistent with @equilibrium.
Specifying $K$ (from thermodynamics) and either $k_f$ or $k_r$ fully determines both.

The MAL is exact for elementary reactions, where stoichiometric coefficients reflect the actual molecularity of the collision.
For complex mechanisms it is an approximation, but one that is thermodynamically consistent by construction: regardless of the path, the steady state reached by MAL coincides with the thermodynamic equilibrium.

```{important}
**MAL vs. empirical power laws.**
Users sometimes specify reaction orders that differ from stoichiometry, for instance when fitting an observed half-order experimentally.
The resulting power law $r = k_f \prod c_i^{n_i} - k_r \prod c_j^{m_j}$ with $n_i \neq |\nu_i|$ is no longer the mass action law.
At equilibrium $r = 0$ gives $k_f/k_r = \prod c_j^{m_j} / \prod c_i^{n_i}$, which depends on composition and is not equal to $K$ in general.
Thermodynamic consistency is lost: the model does not guarantee the correct steady state.
Use empirical orders when they fit data; use MAL when thermodynamic consistency matters.
```

```{admonition} Equilibrium is independent of mechanism
:class: note

$K$ is fixed by $\Delta_r G^\circ = -RT \ln K$, which depends only on the standard chemical potentials of the species, not on how the reaction proceeds.
A catalyst, a solvent change, or a complex multi-step mechanism alters the rate at which equilibrium is reached, but not the equilibrium composition itself.
The consistency condition $k_f/k_r = K$ is exact for elementary steps.
For a complex mechanism, the apparent overall rate constants can be concentration-dependent (Michaelis-Menten is the canonical example), and their ratio need not equal $K$ globally.
The relation holds only for each microscopic elementary step individually.
```

@fig-mal illustrates the consistency condition for $\ce{A <=> B}$ with $K = 3$: the forward and reverse rates cross at exactly the thermodynamic equilibrium extent.

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-mal

import numpy as np
import matplotlib.pyplot as plt

K   = 3.0
k_f = 3.0
k_r = k_f / K

xi    = np.linspace(0, 1, 400)
c_A   = 1 - xi
c_B   = xi
r_f   = k_f * c_A
r_r   = k_r * c_B
r_net = r_f - r_r
xi_eq = K / (1 + K)

fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(xi, r_f,   color="C0",     linewidth=2.0, label=r"$r_f = k_f[\mathrm{A}]$")
ax.plot(xi, r_r,   color="C1",     linewidth=2.0, label=r"$r_r = k_r[\mathrm{B}]$")
ax.plot(xi, r_net, color="#1c4f8a",linewidth=2.0, linestyle="--",
        label=r"$r_\mathrm{net}$")
ax.axvline(xi_eq, color="gray", linestyle=":", linewidth=1.2)
ax.axhline(0, color="gray", linewidth=0.8)
ax.text(xi_eq + 0.02, k_f * 0.9, r"$\xi_\mathrm{eq}$, $r_\mathrm{net}=0$",
        fontsize=9, color="gray", va="top")
ax.set_xlabel(r"Extent of reaction $\xi$ / $\xi_\mathrm{max}$", fontsize=11)
ax.set_ylabel("Rate [a.u.]", fontsize=11)
ax.legend(fontsize=9)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
fig.tight_layout()
```

```{figure} #cell-mal
:name: fig-mal

Forward rate $r_f$, reverse rate $r_r$, and net rate for $\ce{A <=> B}$ with $K = k_f/k_r = 3$.
The net rate is zero at $\xi_\text{eq} = K/(1+K) = 0.75$, the equilibrium extent predicted by @equilibrium.
Setting $k_f/k_r = K$ is what makes the kinetic model land on the correct thermodynamic equilibrium.
```

---

The next chapter asks how the rate constant $k$ itself depends on temperature, connecting the Arrhenius equation back to the Maxwell-Boltzmann distribution and the barrier-enthalpy consequence of the consistency condition $k_f/k_r = K(T)$.
