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
This part develops one reaction model, extended one feature at a time:

$$
\varphi(a_i, T) = k^f(T)\prod_i a_i^{e_i^f} - k^r(T)\prod_i a_i^{e_i^r},
\quad k^r(T) = \frac{k^f(T)}{K(T)},
\quad a_i = \frac{\gamma_i c_i}{c^\circ}.
$$

Each chapter adds one feature.
`MassActionReaction` is the existing CADET interface and a special case of `ThermodynamicReaction` with $\gamma_i = 1$, $c^\circ = 1\ \text{M}$, and $k^r$ as a free parameter.
Equilibrium (via $K(T)$) determines the admissible state; kinetics only sets how fast it is reached.

```{admonition} Install
:class: note

Install the library once from the project root:

    pip install -e .
```

```{admonition} Key design decisions
:class: tip
**One object, three knobs.**
`ThermodynamicReaction` is the central object.
Every chapter in this part configures one of three independent arguments:

| Argument               | Controls                                                   | Default                      |
| ---------------------- | ---------------------------------------------------------- | ---------------------------- |
| `equilibrium_constant` | $K(T)$; sets the equilibrium composition                   | required                     |
| `rate_constant`        | $k^f(T)$; sets the relaxation timescale (`mode="kinetic"`) | omit for `mode="equil"`      |
| `activity_coefficient` | $\gamma_i$; corrects for non-ideality                      | `ActivityCoefficientIdeal()` |

`MassActionReaction` is a special case: $\gamma_i = 1$, $c^\circ = 1\ \text{M}$, and $k^r$ treated as a free parameter rather than derived from $K(T)$.

**Kinetic source term.**
The reaction contribution to $\partial c_i / \partial t$ is $r_i = \sum_j \nu_{ij}\, \varphi_j(\mathbf{a}, T)$, where $\varphi_j$ is evaluated in terms of activities $a_i = \gamma_i c_i / c^\circ$.

**Equilibrium (DAE) mode.**
When a reaction is fast relative to the transport timescale, its ODE is replaced by the algebraic constraint $\ln Q(\mathbf{a}, T) = \ln K(T)$, turning the transport PDE into a PDAE.
This is a timescale-separation limit (fast reactions), not a limit in the value of $K$.

**Activity corrections.**
Ionic strength $I = \tfrac{1}{2}\sum_i c_i z_i^2$ is computed once per time step and passed to every activity coefficient module via `PhysicalState`.
The ionic strength model and activity coefficient model are decoupled at the API level and can be combined consistently within their validity ranges.

All reactions are evaluated as $\varphi(\mathbf{a}, T)$; the rest of the API only determines how $\mathbf{a}$, $k^f$, and $K$ are computed.
```

**Chapters in this part:**

- @implementation-source-term: PDE context, building blocks (`Species`, `Component`, `ReactionModel`), and `MassActionReaction` as the existing interface.
- @implementation-equilibrium: `ThermodynamicReaction` in equilibrium mode; $K(T)$ via van't Hoff, Kirchhoff correction, and custom forms.
- @implementation-kinetics: kinetic mode, `RateConstantFixed` and `RateConstantArrhenius`; thermodynamic consistency across temperatures.
- @implementation-energy-balance: coupled energy balance; temperature as a dynamic state; analytic $\partial\varphi/\partial T$ Jacobian. *(optional)*
- @implementation-activity: activity corrections ($a_i = \gamma_i c_i/c^\circ$), ionic strength models, Debye-Hückel and Davies, temperature-dependent $A$ via $\varepsilon_r(T)$, and the apparent pKa shift.
- @implementation-acid-base: pH, the `pKa` factory, water autoionisation, and Davies corrections.
- @implementation-buffer: buffer capacity, mixed buffers, and ionic strength effects on $\beta$.
- @implementation-enzyme: saturation kinetics (`MichaelisMenten`, `HillRate`) as the finite-site constraint expressed as a kinetic rate law.
- @implementation-interface: the residual/Jacobian contract and CADET-Core integration.
