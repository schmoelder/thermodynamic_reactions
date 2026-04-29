---
kernelspec:
  name: python3
  display_name: Python 3 (ipykernel)
---

(implementation)=
# Part 4: Implementation in CADET

Parts 1--3 built the thermodynamic and kinetic theory.
This part translates it into working code using the `reactions` library.
CADET solves a transport equation with a reaction source term $\mathbf{r}(\mathbf{c}, T)$ at every spatial point along the column; the flux term handles convection and dispersion, and this library provides $\mathbf{r}$.

The part develops one reaction model, extended one feature at a time:

$$
\varphi = k_f \prod_i c_i^{e_i} - k_r \prod_j c_j^{e_j}
\quad\xrightarrow{k_r = k_f/K(T)}\quad
\varphi(a_i, T)
\quad\xrightarrow{\text{fast reaction}}\quad
\ln Q = \ln K(T).
$$

Each step addresses one specific limitation of the previous.
`MassActionReaction` and `ThermodynamicReaction` are not two different tools; they are the same model with different features enabled.
Under ideal conditions ($\gamma = 1$, $c^\circ = 1\ \text{M}$) and at fixed temperature they are numerically equivalent; the only difference is whether $k_r = k_f/K(T)$ is enforced.

```{admonition} Install
:class: note

Install the library once from the project root:

    pip install -e .
```

```{admonition} Key design decisions
:class: tip

**Kinetic source term.**
The reaction contribution to $\partial c_i / \partial t$ is $r_i = \sum_j \nu_{ij}\, \varphi_j(\mathbf{a}, T)$, where $\varphi_j$ is evaluated in terms of activities $a_i = \gamma_i c_i / c^\circ$.

**Thermodynamic consistency.**
Forward and reverse rate constants are never independent: $k_r(T) = k_f(T) / K(T)$.
Storing $k_r$ separately breaks this identity whenever temperature changes.

**Equilibrium (DAE) mode.**
When a reaction is fast relative to the transport timescale, its ODE is replaced by the algebraic constraint $\ln Q(\mathbf{a}, T) = \ln K(T)$, turning the transport PDE into a PDAE.
This is a timescale-separation argument, not a large-$K$ limit.

**Activity corrections.**
Ionic strength $I = \tfrac{1}{2}\sum_i c_i z_i^2$ is computed once per time step and passed to every activity coefficient module via `PhysicalState`.
The choice of ionic strength model is decoupled from the choice of activity coefficient model; any combination is valid.
```

**Chapters in this part:**

- @implementation-source-term: PDE context and building blocks (`Species`, `Component`, `ReactionModel`).
- @implementation-rate-laws: mass-action rate laws from irreversible to reversible multi-species reactions.
- @implementation-equilibrium: thermodynamic consistency ($k_r = k_f/K(T)$), $K(T)$ via van't Hoff, and the equilibrium (DAE) mode.
- @implementation-activity: activity corrections ($a_i = \gamma_i c_i/c^\circ$), ionic strength models, Debye-Hückel and Davies, and the apparent pKa shift.
- @implementation-kinetics: Arrhenius $k_f(T)$, thermodynamic consistency across temperatures, and saturation kinetics (Michaelis-Menten, Hill).
- @implementation-acid-base: pH, the `pKa` factory, water autoionisation, and Davies corrections.
- @implementation-buffer: buffer capacity, mixed buffers, and ionic strength effects on $\beta$.
- @implementation-interface: the residual/Jacobian contract and planned CADET-Core integration.
