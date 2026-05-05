---
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

(implementation-source-term)=
# The Reaction Source Term

CADET solves the transport equation

$$
\frac{\partial \mathbf{c}}{\partial t}
= -\nabla \cdot \mathbf{J}(\mathbf{c})
+ \mathbf{f}_\text{react}(\mathbf{c}, T)
$$

at every spatial discretisation point along the column.
The transport term $-\nabla \cdot \mathbf{J}$ represents convection and dispersion;
the source term $\mathbf{f}_\text{react}(\mathbf{c}, T)$ encodes all chemical reactions and is provided by this library.
Parts 1--3 established the thermodynamic and kinetic foundations that determine $\mathbf{f}_\text{react}$; Part 4 translates this theory into a computational form compatible with CADET-Core.

## Two modes for $\mathbf{f}_\text{react}$

The source term takes two complementary forms depending on the timescale of reactions relative to transport.

**Kinetic mode** applies when reactions evolve on timescales comparable to transport.
Each reaction $j$ contributes a reaction flux $\varphi_j(\mathbf{c}, T)$, and the net source term for species $i$ is

$$
f_{\text{react},i} = \sum_j \nu_{ij}\, \varphi_j(\mathbf{c}, T).
$$

The residual contribution is $F_i = \dot{c}_i - f_{\text{react},i}$, integrated as an ordinary differential equation in the transport PDE.

**Equilibrium mode** applies when reactions are fast relative to transport and can be treated as instantaneous (@kinetics).
One differential equation is replaced by an algebraic constraint enforcing thermodynamic equilibrium,

$$
\ln Q(\mathbf{c}, T) - \ln K(T) = 0,
$$

turning the transport PDE into a partial differential-algebraic equation (PDAE).
Constraints are enforced independently at each spatial discretisation point.

## Canonical reaction formalism

Both modes fit within a single unified structure:

$$
\frac{\partial \mathbf{c}}{\partial t}
= -\nabla \cdot \mathbf{J}(\mathbf{c})
+ \mathbf{S}\,\boldsymbol{\varphi}(\mathbf{c}, T),
$$

where $\mathbf{S}$ is the stoichiometric matrix and $\boldsymbol{\varphi}$ collects the reaction fluxes.
Each reaction $j$ contributes

$$
\varphi_j = k_j(T)\, \prod_i a_i^{\nu_{ij}},
$$

with activities $a_i = \gamma_i(\mathbf{c}, I)\, c_i / c^\circ$.
In CADET the kinetic exponents can differ from the stoichiometric coefficients $\nu_{ij}$; the notation here uses $\nu_{ij}$ for both, with the understanding that they coincide for elementary steps and may be fitted independently for empirical rate laws.
Thermodynamic equilibrium enters as

$$
\sum_i \nu_{ij}\, \mu_i = 0
\quad \Leftrightarrow \quad
Q_j(\mathbf{a}) = K_j(T).
$$

Kinetic and equilibrium formulations are different closures of the same structure: kinetic systems evolve explicitly in time; equilibrium systems impose algebraic constraints; thermodynamic consistency links both through $k_f(T)/k_r(T) = K(T)$.

## Building blocks

Three objects cover all cases.

**`Component`** is the primary abstraction for model construction.
It represents a named chemical entity and groups the species that entity can exist as, such as different protonation or charge states.
Keyword arguments are forwarded to the underlying `Species`, so charge and solvent status can be set directly:

```{code-cell} ipython3
from reactions.api import Component, Species

a      = Component("A")                   # shorthand: one species, charge 0
proton = Component("H+",  charge=+1)
water  = Component("H2O", is_solvent=True)
```

An explicit species list is needed when a component spans multiple charge states:

```{code-cell} ipython3
phosphate = Component("phosphate", [
    Species("H3PO4",  charge=0),
    Species("H2PO4-", charge=-1),
    Species("HPO4-2", charge=-2),
    Species("PO4-3",  charge=-3),
])
```

```{admonition} Components are a library concept
:class: note

CADET-Core operates on a flat list of species concentrations and has no notion of components.
The `Component` abstraction exists to simplify model construction and enable automation, for example in pKa-based reaction generation.
```

**`ReactionModel`** assembles components and reactions into the system passed to CADET-Core.
It exposes `residual(c, c_dot, T)` and `jacobian(c, c_dot, T)`, which implement the source term contract at each solver step (@implementation-interface).

---

The next chapter builds the rate laws that populate $\mathbf{f}_\text{react}$, starting from irreversible reactions and extending to fully reversible multi-species systems.
