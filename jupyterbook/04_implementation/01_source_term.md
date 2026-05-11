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

at every spatial discretization point along the column.
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

The residual contribution $F_i = \dot{c}_i - f_{\text{react},i}$ is included as an ODE residual in the transport solver.

**Equilibrium mode** applies when reactions are fast relative to transport and can be treated as instantaneous (@kinetics).
One differential equation is replaced by an algebraic constraint enforcing thermodynamic equilibrium,

$$
\ln Q(\mathbf{a}, T) - \ln K(T) = 0,
$$

turning the transport PDE into a partial differential-algebraic equation (PDAE).
Constraints are enforced independently at each spatial discretization point.
This is a timescale-separation limit (fast reactions), not a limit in the value of $K$.

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
\varphi_j = k_j^f(T)\prod_i a_i^{e_{ij}^f}
           - k_j^r(T)\prod_i a_i^{e_{ij}^r},
$$

with activities $a_i = \gamma_i(\mathbf{c}, I)\, c_i / c^\circ$ and kinetic exponents $e_{ij}^f, e_{ij}^r \geq 0$.
In CADET, $e_{ij}^f$, $e_{ij}^r$ and the stoichiometric coefficients $\nu_{ij}$ are independent parameters.
For an elementary reaction $e_{ij}^f = |\nu_{ij}|$ for reactants (zero for products) and $e_{ij}^r = |\nu_{ij}|$ for products (zero for reactants), but they may differ for empirical rate laws.
Thermodynamic equilibrium enters as

$$
\sum_i \nu_{ij}\, \mu_i = 0
\quad \Leftrightarrow \quad
Q_j(\mathbf{a}) = K_j(T).
$$

Kinetic and equilibrium modes are two closures of the same stoichiometric structure, selected by timescale separation: kinetic mode resolves the relaxation trajectory; equilibrium mode imposes $Q(\mathbf{a}, T) = K(T)$ directly.
Thermodynamic consistency links both through $k^f(T)/k^r(T) = K(T)$.

## Building blocks

Three objects cover all cases.

**`Component`** is the primary user-facing object: a named chemical entity.
Most components map to a single chemical form and need no further configuration:

```{code-cell} ipython3
from reactions.api import Component

a      = Component("A")
proton = Component("H+", charge=+1)
```

When a component spans multiple charge states, such as the protonation forms of the same molecule, it takes an explicit list of **`Species`**, each representing one chemical form with its own charge:

```{code-cell} ipython3
from reactions.api import Component, Species

phosphate = Component("phosphate", [
    Species("H3PO4",  charge=0),
    Species("H2PO4-", charge=-1),
    Species("HPO4-2", charge=-2),
    Species("PO4-3",  charge=-3),
])
```

Each `Species` carries a `c_ref` attribute (default 1000 mol/m³) that sets the standard-state concentration $c^\circ_i$ in the activity $a_i = \gamma_i c_i / c^\circ_i$.
The default corresponds to 1 mol/L, the conventional aqueous standard state; a different value is appropriate for gas-phase species or non-aqueous systems where the natural reference concentration differs.
Three further attributes are optional and only required when the energy balance is active (@implementation-energy-balance): `molar_mass` (kg/mol), `density` (kg/m³), and `heat_capacity` (J/(mol·K)).

```{admonition} Components are a library concept
:class: note

CADET-Core operates on a flat list of species concentrations and has no notion of components.
The `Component` abstraction exists to simplify model construction and enable automation, for example in pKa-based reaction generation.
```

**`ReactionModel`** assembles components and reactions into the system passed to CADET-Core.
It exposes `residual(c, c_dot, T)` and `jacobian(c, c_dot, T)`, which implement the source term contract at each solver step (@implementation-interface).
The `reactions.solver` module wraps these into a standalone prototype (`simulate()` for kinetic trajectories and `solve_equilibrium()` for steady states) used throughout Part 4 for development and verification outside CADET-Core.

**Reaction classes** populate a `ReactionModel` and implement one of two constraint structures developed in Part 3.
Stoichiometric reactions (a single extent $\xi$ relates all species changes; equilibrium means $Q = K$) are covered by two closures: `ThermodynamicReaction`, which enforces thermodynamic consistency $k^r(T) = k^f(T)/K(T)$, and `MassActionReaction`, which treats $k^r$ as a free parameter.
Finite-site reactions (a conserved pool of sites partitioned between occupied and free states, producing saturation kinetics) are covered by `EnzymaticReaction`.
Phase equilibria (equality of chemical potentials across coexisting phases) are developed in Part 2 but are not implemented as reaction closures here.

## `MassActionReaction` as starting point

`MassActionReaction` is the existing reaction interface in CADET.
For a reversible step it takes $k^f$ and $k^r$ as independent inputs:

```python
MassActionReaction("A <-> B", kf=2.0, kr=0.5)
```

At fixed temperature `MassActionReaction` is exact: $k^f/k^r$ directly sets the equilibrium composition, here $c_\text{B}/c_\text{A} = k^f/k^r = 4$.
This ratio is a free parameter: nothing in the interface requires it to equal the equilibrium constant $K$ derived from $\Delta_r G^\circ = -RT \ln K$ (@equilibrium).
When temperature changes, $K(T)$ shifts; `MassActionReaction` has no mechanism to follow it.
Structurally, `MassActionReaction` is the unit-activity, free-$k^r$ limit of the stoichiometric closure: equivalent to `ThermodynamicReaction` with $\gamma_i = 1$ and $k^r$ treated as an independent parameter.
`ThermodynamicReaction` restores thermodynamic consistency by enforcing $k^r(T) = k^f(T)/K(T)$ at every evaluation, so the long-time limit always tracks the correct equilibrium.

---

The next chapter introduces `ThermodynamicReaction` in equilibrium mode, connecting the implementation directly to the equilibrium thermodynamics of Part 3 (@implementation-equilibrium).
