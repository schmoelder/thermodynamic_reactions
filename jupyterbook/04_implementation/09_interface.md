---
kernelspec:
  name: python3
  display_name: Python 3 (ipykernel)
---

(implementation-interface)=
# CADET-Core Interface

The `reactions` library is a self-contained Python prototype.
This chapter describes the mathematical contract it must satisfy to connect to CADET-Core, and identifies the C++ abstractions that a full implementation would require.
CADET already implements the mass action law (MAL); the extensions developed in this book add thermodynamic consistency, activity corrections, and equilibrium mode as new building blocks on top of the existing framework.


## The residual / Jacobian contract

CADET-Core expects each reaction model to satisfy the DAE system

$$
F(t,\, \mathbf{c},\, \dot{\mathbf{c}}) = \mathbf{0}
$$

where the residual vector $F$ has one entry per dynamic species.
The library exposes this interface via two methods on `ReactionModel`:

- `ReactionModel.residual(c, c_dot, T)`: assembles $F$
- `ReactionModel.jacobian(c, c_dot, T)`: returns $\partial F / \partial \mathbf{c}$

Two types of row appear in $F$:

**Kinetic rows** (one per kinetic species $i$):

$$
F_i = \dot{c}_i - \sum_j \nu_{ij}\, \varphi_j(\mathbf{a}, T)
$$

**Algebraic rows** (one per equilibrium reaction, replacing the ODE for the dependent species):

$$
F_\text{dep} = \ln Q_j(\mathbf{a}, T) - \ln K_j(T)
$$

The dependent species is chosen as the product with the largest $|\nu_{ij}|$.
This structure is identical to CADET's existing reaction interface; the only extension is that fluxes $\varphi_j$ are evaluated in terms of activities $a_i = \gamma_i c_i / c^\circ$ rather than concentrations directly.


## C++ building blocks

Three abstractions are needed beyond the existing MAL implementation.

**Equilibrium constant functor $K(T)$.** Both van't Hoff and Arrhenius share the form $A\,\exp(B/T)$; a single C++ functor template with that structure covers `EquilibriumConstantVantHoff`, `RateConstantArrhenius`, and their Kirchhoff-corrected variants.
More general forms (`EquilibriumConstantTabulated`, `EquilibriumConstantCustom`) require a callable interface with the same signature.

**Activity coefficient plugin $\gamma(\mathbf{c}, T)$.** The Davies and Debye-Hückel models are both functions of ionic strength $I(\mathbf{c}) = \frac{1}{2}\sum_i z_i^2 c_i / c^\circ$, itself a linear function of the concentration state.
A plugin interface that takes the concentration vector and charge array and returns $\gamma_i$ covers both models and allows future extensions (Pitzer, custom) without changing the solver.

**Parameter-state dependencies.** Both $K(T)$ and $\gamma(\mathbf{c}, T)$ are state-dependent parameters: $K$ depends on temperature; $\gamma$ depends on composition through $I$.
Recognising this pattern allows the C++ implementation to reuse a common infrastructure for any quantity of the form $f(\text{state})$, rather than treating each model as a special case.
The same infrastructure would support future energy-balance coupling, where $T$ itself becomes a dynamic state variable.


## Jacobian verification

The Jacobian implementation uses analytic derivatives for both kinetic and equilibrium rows.
The following snippet verifies the analytic Jacobian against a finite-difference reference for each model type developed in this part:

```{code-cell} ipython3
import numpy as np
from reactions.api import (
    Species, Component,
    EquilibriumConstant,
    RateConstantFixed,
    MassActionReaction, ThermodynamicReaction,
    ReactionModel, pKa,
)


def check_jacobian(model, c, label, T=298.15, tol=1e-4):
    c_dot = np.zeros(len(c))
    J_ana = model.jacobian(c, c_dot, T=T)
    eps   = 1e-5
    r0    = model.residual(c, c_dot, T=T)
    n     = len(c)
    J_fd  = np.zeros((n, n))
    for k in range(n):
        cp = c.copy(); cp[k] += eps
        J_fd[:, k] = (model.residual(cp, c_dot, T=T) - r0) / eps
    err    = np.max(np.abs(J_ana - J_fd))
    status = "PASSED" if err < tol else "FAILED"
    print(f"  {label:<45}  max|J_ana - J_fd| = {err:.2e}  [{status}]")


comp_ab  = Component("ab", [Species("A"), Species("B")])
model_ma = ReactionModel(
    components=[comp_ab],
    reactions=[MassActionReaction("A <-> B", kf=2.0, kr=0.5)],
)
model_td = ReactionModel(
    components=[comp_ab],
    reactions=[
        ThermodynamicReaction(
            "A <-> B",
            mode="kinetic",
            equilibrium_constant=EquilibriumConstant(4.0),
            rate_constant=RateConstantFixed(2.0),
        )
    ],
)

acetate   = Component("acetate",   [Species("AcOH", charge=0),  Species("AcO-", charge=-1)])
proton    = Component("proton",    [Species("H+", charge=+1)])
hydroxide = Component("hydroxide", [Species("OH-", charge=-1)])
water     = Component("water",     [Species("H2O", charge=0)])
model_ac  = ReactionModel(
    components=[acetate, proton, hydroxide, water],
    reactions=[
        ThermodynamicReaction("AcOH <-> AcO- + H+", mode="equil",
                              equilibrium_constant=pKa(4.756)),
        ThermodynamicReaction("H2O <-> H+ + OH-",   mode="equil",
                              equilibrium_constant=EquilibriumConstant(1e-14)),
    ],
    T=298.15,
)

print("Analytic Jacobian verification:")
check_jacobian(model_ma, np.array([1.0, 0.5]),       "A⇌B MassActionReaction")
check_jacobian(model_td, np.array([500.0, 200.0]),   "A⇌B ThermodynamicReaction kinetic")
check_jacobian(model_ac, np.array([10.0, 90.0, 6.3e-5, 1.6e-4]),
               "Acetic acid (equil)")
```

---

Part 5 applies the full framework to a realistic problem: designing a pH gradient for
ion-exchange chromatography, where buffer composition, ionic strength, and protein
charge state are all coupled (@case-ph-gradient).
