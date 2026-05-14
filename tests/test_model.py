"""
Tests for reactions.model.ReactionModel:
analytic Jacobian (dc and dT) vs finite differences, energy balance Jacobian.
"""

import numpy as np
import pytest
from reactions.equilibrium import (
    EquilibriumConstant,
    EquilibriumConstantVantHoff,
    EquilibriumConstantVantHoffCp,
    pKa,
)
from reactions.model import ReactionModel
from reactions.rate import RateConstantArrhenius, RateConstantFixed
from reactions.reaction import MassActionReaction, ThermodynamicReaction
from reactions.species import Component, Species

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _check_jacobian(model, c, T=298.15, eps=1e-7, tol=1e-4):
    """
    Compare model.jacobian() (analytic) with a relative finite-difference step.

    Uses a step h = eps * max(|c_k|, 1e-10) so that the step scales with
    the concentration magnitude.  Returns the max relative error over all
    nonzero Jacobian entries (absolute error where both analytic and FD are
    near zero).
    """
    c_dot = np.zeros(len(c))
    J_ana = model.jacobian(c, c_dot, T=T)
    r0 = model.residual(c, c_dot, T=T)
    n = len(c)
    J_fd = np.zeros((n, n))
    for k in range(n):
        h = eps * max(abs(float(c[k])), 1e-10)
        cp = c.copy()
        cp[k] += h
        J_fd[:, k] = (model.residual(cp, c_dot, T=T) - r0) / h
    scale = np.maximum(np.abs(J_ana), 1.0)
    return float(np.max(np.abs(J_ana - J_fd) / scale))


def _fd_jacobian_dT(model, c, T, eps=1e-5):
    """Full-model FD approximation of d(residual)/dT."""
    r0 = model.residual(c, np.zeros(len(c)), T)
    return (model.residual(c, np.zeros(len(c)), T + eps) - r0) / eps


def _fd_jacobian_dc(model, c, T, eps=1e-5):
    """Full-model FD approximation of d(residual)/dc."""
    n = len(c)
    r0 = model.residual(c, np.zeros(n), T)
    J = np.zeros((n, n))
    for k in range(n):
        c_p = c.copy()
        c_p[k] += eps
        J[:, k] = (model.residual(c_p, np.zeros(n), T) - r0) / eps
    return J


def _eb_rhs_T(model, c, T, rho_cp):
    """Energy balance scalar rhs: -Q_dot/rho_cp = -Σ_j ΔH_j φ_j / rho_cp."""
    state = model.make_state(c, T)
    aux = model.make_aux(state)
    Q_dot = 0.0
    for j, rxn in enumerate(model.reactions):
        if not model.kinetic_mask[j]:
            continue
        eq = getattr(rxn, "equilibrium_constant", None)
        if eq is None:
            continue
        Q_dot += eq.reaction_enthalpy(T) * rxn.net_rate(state, aux, model.species_index)
    return -Q_dot / rho_cp


def _eb_jac_analytic(model, c, T, rho_cp):
    """Analytic energy balance Jacobian row (∂rhs_T/∂c, ∂rhs_T/∂T)."""
    state = model.make_state(c, T)
    aux = model.make_aux(state)
    n = len(c)
    jac_c = np.zeros(n)
    jac_T = 0.0
    for j, rxn in enumerate(model.reactions):
        if not model.kinetic_mask[j]:
            continue
        eq = getattr(rxn, "equilibrium_constant", None)
        if eq is None:
            continue
        dH = eq.reaction_enthalpy(T)
        dH_dT = eq.d_reaction_enthalpy_dT(T)
        phi = rxn.net_rate(state, aux, model.species_index)
        dphi_dc = rxn.net_rate_jac(state, aux, model.species_index)
        dphi_dT = rxn.net_rate_dT(state, aux, model.species_index)
        jac_c -= dH / rho_cp * dphi_dc
        jac_T -= (dH * dphi_dT + phi * dH_dT) / rho_cp
    return jac_c, jac_T


# ---------------------------------------------------------------------------
# Jacobian dc
# ---------------------------------------------------------------------------


def test_jacobian_mass_action_reaction():
    """Analytic Jacobian of MassActionReaction matches finite differences."""
    comp = Component("ab", [Species("A"), Species("B")])
    model = ReactionModel(
        components=[comp],
        reactions=[MassActionReaction("A <-> B", kf=2.0, kr=0.5)],
    )
    c = np.array([1.0, 0.5])
    err = _check_jacobian(model, c)
    assert err < 1e-4, f"MassActionReaction Jacobian error {err:.2e} exceeds 1e-4"


def test_jacobian_thermodynamic_kinetic():
    """Analytic Jacobian of ThermodynamicReaction (kinetic) matches finite differences."""
    comp = Component("ab", [Species("A"), Species("B")])
    model = ReactionModel(
        components=[comp],
        reactions=[
            ThermodynamicReaction(
                "A <-> B",
                mode="kinetic",
                equilibrium_constant=EquilibriumConstant(4.0),
                rate_constant=RateConstantFixed(kf_value=2.0),
            )
        ],
    )
    c = np.array([500.0, 200.0])
    err = _check_jacobian(model, c)
    assert err < 1e-4, (
        f"ThermodynamicReaction kinetic Jacobian error {err:.2e} exceeds 1e-4"
    )


def test_jacobian_thermodynamic_equil_acetic():
    """Analytic Jacobian of ThermodynamicReaction (equil, acetic acid) matches FD."""
    acetate = Component(
        "acetate", [Species("AcOH", charge=0), Species("AcO-", charge=-1)]
    )
    proton = Component("proton", [Species("H+", charge=+1)])
    hydroxide = Component("hydroxide", [Species("OH-", charge=-1)])
    water = Component("water", [Species("H2O", charge=0)])
    model = ReactionModel(
        components=[acetate, proton, hydroxide, water],
        reactions=[
            ThermodynamicReaction(
                "AcOH <-> AcO- + H+",
                mode="equil",
                equilibrium_constant=pKa(4.756),
            ),
            ThermodynamicReaction(
                "H2O <-> H+ + OH-",
                mode="equil",
                equilibrium_constant=EquilibriumConstant(1e-14),
            ),
        ],
        T=298.15,
    )
    c = np.array([10.0, 90.0, 6.3e-5, 1.6e-4, 1000.0])
    err = _check_jacobian(model, c)
    assert err < 1e-4, (
        f"Acetic acid equil Jacobian relative error {err:.2e} exceeds 1e-4"
    )


@pytest.mark.parametrize("T", [298.15, 310.0, 330.0])
def test_jacobian_dc_vanthoff_arrhenius(T):
    """jacobian (d(residual)/dc) analytic matches FD: VantHoff K + Arrhenius kf."""
    model = ReactionModel(
        components=[Component("A"), Component("B")],
        reactions=[
            ThermodynamicReaction(
                "A <-> B",
                mode="kinetic",
                equilibrium_constant=EquilibriumConstantVantHoff(dH=-20e3, dS=-50.0),
                rate_constant=RateConstantArrhenius(A=1e10, Ea=40e3),
            ),
        ],
    )
    c = np.array([300.0, 700.0])
    ana = model.jacobian(c, np.zeros(2), T)
    fd = _fd_jacobian_dc(model, c, T)
    np.testing.assert_allclose(ana, fd, rtol=1e-4, atol=1e-8)


# ---------------------------------------------------------------------------
# Jacobian dT
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("T", [298.15, 310.0, 330.0])
def test_jacobian_dT_vanthoff_arrhenius(T):
    """jacobian_dT analytic matches FD: VantHoff K + Arrhenius kf."""
    model = ReactionModel(
        components=[Component("A"), Component("B")],
        reactions=[
            ThermodynamicReaction(
                "A <-> B",
                mode="kinetic",
                equilibrium_constant=EquilibriumConstantVantHoff(dH=-20e3, dS=-50.0),
                rate_constant=RateConstantArrhenius(A=1e10, Ea=40e3),
            ),
        ],
    )
    c = np.array([300.0, 700.0])
    ana = model.jacobian_dT(c, np.zeros(2), T)
    fd = _fd_jacobian_dT(model, c, T)
    np.testing.assert_allclose(ana, fd, rtol=1e-4, atol=1e-8)


def test_jacobian_dT_vanthoffcp():
    """jacobian_dT analytic matches FD: VantHoffCp (temperature-dependent dH)."""
    model = ReactionModel(
        components=[Component("A"), Component("B")],
        reactions=[
            ThermodynamicReaction(
                "A <-> B",
                mode="kinetic",
                equilibrium_constant=EquilibriumConstantVantHoffCp(
                    dH=-20e3,
                    dS=-50.0,
                    dCp=100.0,
                ),
                rate_constant=RateConstantArrhenius(A=1e10, Ea=40e3),
            ),
        ],
    )
    c = np.array([300.0, 700.0])
    ana = model.jacobian_dT(c, np.zeros(2), T=310.0)
    fd = _fd_jacobian_dT(model, c, T=310.0)
    np.testing.assert_allclose(ana, fd, rtol=1e-4, atol=1e-8)


def test_jacobian_dT_fixed_is_zero():
    """jacobian_dT is zero for fixed K and fixed kf (no T-dependence)."""
    model = ReactionModel(
        components=[Component("A"), Component("B")],
        reactions=[
            ThermodynamicReaction(
                "A <-> B",
                mode="kinetic",
                equilibrium_constant=EquilibriumConstant(4.0),
                rate_constant=RateConstantFixed(2000.0),
            ),
        ],
    )
    c = np.array([300.0, 700.0])
    ana = model.jacobian_dT(c, np.zeros(2), T=310.0)
    np.testing.assert_allclose(ana, 0.0, atol=1e-12)


def test_jacobian_dT_equil_row():
    """jacobian_dT analytic matches FD for equilibrium (ln Q - ln K) row."""
    model = ReactionModel(
        components=[Component("A"), Component("B")],
        reactions=[
            ThermodynamicReaction(
                "A <-> B",
                mode="equil",
                equilibrium_constant=EquilibriumConstantVantHoff(dH=-20e3, dS=-50.0),
            ),
        ],
    )
    c = np.array([300.0, 700.0])
    ana = model.jacobian_dT(c, np.zeros(2), T=310.0)
    fd = _fd_jacobian_dT(model, c, T=310.0)
    np.testing.assert_allclose(ana, fd, rtol=1e-4, atol=1e-10)


def test_jacobian_dT_mixed_reactions():
    """jacobian_dT analytic matches FD: two reactions, kinetic + equil."""
    model = ReactionModel(
        components=[Component("A"), Component("B"), Component("C")],
        reactions=[
            ThermodynamicReaction(
                "A <-> B",
                mode="kinetic",
                equilibrium_constant=EquilibriumConstantVantHoff(dH=-20e3, dS=-50.0),
                rate_constant=RateConstantArrhenius(A=1e10, Ea=40e3),
            ),
            ThermodynamicReaction(
                "B <-> C",
                mode="equil",
                equilibrium_constant=EquilibriumConstant(2.0),
            ),
        ],
    )
    c = np.array([300.0, 400.0, 300.0])
    ana = model.jacobian_dT(c, np.zeros(3), T=310.0)
    fd = _fd_jacobian_dT(model, c, T=310.0)
    np.testing.assert_allclose(ana, fd, rtol=1e-4, atol=1e-8)


# ---------------------------------------------------------------------------
# Energy balance Jacobian
# ---------------------------------------------------------------------------


def test_energy_balance_jac_analytic_vanthoff_arrhenius():
    """Analytic J[n, :n] and J[n, n] match FD for VantHoff K + Arrhenius kf."""
    water = Species("H2O", molar_mass=0.018, density=1000.0, heat_capacity=75.3)
    model = ReactionModel(
        components=[Component("A"), Component("B"), Component("water", [water])],
        reactions=[
            ThermodynamicReaction(
                "A <-> B",
                mode="kinetic",
                equilibrium_constant=EquilibriumConstantVantHoff(dH=-20e3, dS=-50.0),
                rate_constant=RateConstantArrhenius(A=1e10, Ea=40e3),
            )
        ],
    )
    c = np.array([800.0, 200.0, 1000.0])
    rho_cp = model.volumetric_heat_capacity(c)
    T = 310.0
    eps = 1e-5

    rhs0 = _eb_rhs_T(model, c, T, rho_cp)
    jac_c_fd = np.array(
        [
            (_eb_rhs_T(model, c + np.eye(len(c))[k] * eps, T, rho_cp) - rhs0) / eps
            for k in range(len(c))
        ]
    )
    jac_T_fd = (_eb_rhs_T(model, c, T + eps, rho_cp) - rhs0) / eps

    jac_c, jac_T = _eb_jac_analytic(model, c, T, rho_cp)
    np.testing.assert_allclose(jac_c, jac_c_fd, rtol=1e-5)
    np.testing.assert_allclose(jac_T, jac_T_fd, rtol=1e-5)


def test_energy_balance_jac_analytic_vanthoffcp():
    """With VantHoffCp, dΔH/dT = dCp contributes non-trivially to J[n, n]."""
    dCp = 150.0
    water = Species("H2O", molar_mass=0.018, density=1000.0, heat_capacity=75.3)
    model = ReactionModel(
        components=[Component("A"), Component("B"), Component("water", [water])],
        reactions=[
            ThermodynamicReaction(
                "A <-> B",
                mode="kinetic",
                equilibrium_constant=EquilibriumConstantVantHoffCp(
                    dH=-20e3, dS=-50.0, dCp=dCp
                ),
                rate_constant=RateConstantArrhenius(A=1e10, Ea=40e3),
            )
        ],
    )
    c = np.array([600.0, 400.0, 1000.0])
    rho_cp = model.volumetric_heat_capacity(c)
    T = 305.0
    eps = 1e-5

    rhs0 = _eb_rhs_T(model, c, T, rho_cp)
    jac_T_fd = (_eb_rhs_T(model, c, T + eps, rho_cp) - rhs0) / eps
    _, jac_T = _eb_jac_analytic(model, c, T, rho_cp)


# ---------------------------------------------------------------------------
# check_conservation
# ---------------------------------------------------------------------------


def test_check_conservation_empty_lhs_not_conserved():
    """
    Empty-LHS reaction '<-> H+ + OH-' creates both species from an untracked
    reservoir (water). A component containing both species must be flagged as
    non-conserved.
    """
    from reactions.ionic import IonicStrengthIdeal

    H = Species("H+", charge=1)
    OH = Species("OH-", charge=-1)
    water_pair = Component("water_pair", [H, OH])

    model = ReactionModel(
        components=[water_pair],
        reactions=[
            ThermodynamicReaction(
                "<-> H+ + OH-",
                mode="equil",
                equilibrium_constant=pKa(14),
            )
        ],
        ionic_strength=IonicStrengthIdeal(),
    )
    reports = model.check_conservation()
    assert len(reports) == 1
    assert not reports[0].conserved


def test_check_conservation_closed_stoich_conserved():
    """
    HAc <-> Ac- + H+: total acetate (HAc + Ac-) is conserved because both
    species are on the same side of the closed stoichiometry.
    Contrast with the empty-LHS case where both products appear from nothing.
    """
    from reactions.ionic import IonicStrengthIdeal

    HAc = Species("HAc", charge=0)
    Ac = Species("Ac-", charge=-1)
    H = Species("H+", charge=1)
    acetic = Component("acetic_acid", [HAc, Ac])
    H_comp = Component("H_plus", [H])

    model = ReactionModel(
        components=[acetic, H_comp],
        reactions=[
            ThermodynamicReaction(
                "HAc <-> Ac- + H+",
                mode="equil",
                equilibrium_constant=pKa(4.756),
            )
        ],
        ionic_strength=IonicStrengthIdeal(),
    )
    reports = model.check_conservation()
    assert len(reports) == 1
    assert reports[0].component.name == "acetic_acid"
    assert reports[0].conserved
