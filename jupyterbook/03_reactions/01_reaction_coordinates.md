---
kernelspec:
  name: python3
  display_name: Python 3 (ipykernel)
---

(reaction-coordinates)=
# Reaction Coordinates and Stoichiometry

@thermodynamic-potentials established that the Gibbs energy $G$ is minimized at equilibrium under constant $T$ and $P$.
@chemical-potential showed that $G = \sum_i \mu_i n_i$ and that $\mu_i = \mu_i^\circ + RT \ln a_i$ encodes how $G$ responds to changes in composition.
Applying these results to a chemical reaction requires first describing how the composition itself changes.


## Stoichiometric coefficients

When a reaction occurs, the amounts of all species change simultaneously and in fixed ratios.
The key observation is that these changes are not independent: all composition changes are locked to a single degree of freedom.
The molar ratios are encoded by the **stoichiometric coefficients** $\nu_i$: positive for products, negative for reactants.

Take the reaction $\ce{A + 2B <=> C}$ as a concrete example.
If the system produces a small amount $d\xi$ of $\mathrm{C}$, then by stoichiometry:

$$dn_\mathrm{C} = +d\xi, \qquad dn_\mathrm{B} = -2\,d\xi, \qquad dn_\mathrm{A} = -d\xi$$

The single variable $\xi$ drives all three changes.

For a general reaction written as
$$
\sum_i \nu_i \mathrm{A}_i = 0,
$$
this becomes:
$$
dn_i = \nu_i\, d\xi \qquad \text{for all } i
$$

where $\nu_i > 0$ for products and $\nu_i < 0$ for reactants.


## Extent of reaction

Integrating from initial amounts $n_{i,0}$ gives:
$$
n_i(\xi) = n_{i,0} + \nu_i \xi
$$

The variable $\xi$ is the **extent of reaction** (SI unit: mol).
An increase of $1\,\mathrm{mol}$ corresponds to one stoichiometric event of the reaction as written.
The physically allowed range of $\xi$ is determined by $n_i(\xi) \ge 0$ for all species.

The composition of the system at any moment is therefore a single point on a one-dimensional manifold in composition space, parameterised entirely by $\xi$.
All the thermodynamic machinery developed for multicomponent systems now reduces to a function of one variable.

---

The next chapter uses $\xi$ to express the Gibbs energy as a function of composition and identifies the thermodynamic driving force for the reaction.
