"""
Unit tests for equilibrium constant K(T) and rate constant kf(T) models.
"""

import numpy as np
import pytest

from reactions.api import (
    Component,
    EquilibriumConstant,
    EquilibriumConstantCustom,
    EquilibriumConstantPolynomial,
    EquilibriumConstantTabulated,
    EquilibriumConstantVantHoff,
    EquilibriumConstantVantHoffCp,
    R_GAS,
    RateConstantArrhenius,
    RateConstantFixed,
    RateConstantPolynomial,
    ReactionModel,
    Species,
    ThermodynamicReaction,
)
from reactions.solver import simulate

C_REF: float = 1000.0


# ---------------------------------------------------------------------------
# van't Hoff: K(T) temperature sweep
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("T", [298.15, 320.0, 350.0])
def test_vanthoff_converges_at_each_temperature(T):
    """EquilibriumConstantVantHoff: simulation reaches K(T) at each temperature."""
    dH, dS, kf = -20e3, -50.0, 1.0
    K_expected = np.exp(-dH / (8.314462 * T) + dS / 8.314462)

    comp = Component("ab", [Species("A"), Species("B")])
    model = ReactionModel(
        components=[comp],
        reactions=[
            ThermodynamicReaction(
                "A <-> B",
                mode="kinetic",
                equilibrium_constant=EquilibriumConstantVantHoff(dH=dH, dS=dS),
                rate_constant=RateConstantFixed(kf_value=kf),
            )
        ],
    )

    tau = C_REF / (kf * (1.0 + 1.0 / K_expected))
    result = simulate(model, {"A": 1000.0, "B": 0.0}, (0, 10 * tau), T=T, n_points=300)
    assert result.success

    Q = (result["B"][-1] / C_REF) / (result["A"][-1] / C_REF)
    assert abs(Q - K_expected) / K_expected < 1e-3


# ---------------------------------------------------------------------------
# Arrhenius + van't Hoff: kf(T)/kr(T) == K(T)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("T", [250.0, 280.0, 298.15, 320.0, 350.0, 400.0])
def test_arrhenius_vanthoff_kf_kr_ratio(T):
    """kf(T) / kr(T) == K(T) to 1e-12 relative tolerance across temperature."""
    dH = -30e3
    dS = -80.0
    A_arr = 1e8
    Ea = 50e3

    eq_const = EquilibriumConstantVantHoff(dH=dH, dS=dS)
    comp = Component("ab", [Species("A"), Species("B")])
    rxn = ThermodynamicReaction(
        "A <-> B",
        mode="kinetic",
        equilibrium_constant=eq_const,
        rate_constant=RateConstantArrhenius(A=A_arr, Ea=Ea),
    )

    K_val = eq_const.K(T)
    kf_val = rxn.kf(T)
    kr_val = rxn.kr(T)

    ratio = kf_val / kr_val
    rel_err = abs(ratio - K_val) / K_val
    assert rel_err < 1e-12, (
        f"At T={T} K: kf/kr={ratio:.6e}, K={K_val:.6e}, rel_err={rel_err:.2e}"
    )


# ---------------------------------------------------------------------------
# reaction_enthalpy fallbacks
# ---------------------------------------------------------------------------


def test_reaction_enthalpy_fd_fallback_custom():
    """EquilibriumConstantCustom.reaction_enthalpy() falls back to FD via van't Hoff."""
    dH_true = -20e3
    dS = -50.0
    K_fn = lambda T: float(np.exp(-dH_true / (R_GAS * T) + dS / R_GAS))  # noqa: E731
    eq = EquilibriumConstantCustom(K_fn)
    dH_fd = eq.reaction_enthalpy(298.15)
    assert abs(dH_fd - dH_true) / abs(dH_true) < 1e-6


def test_reaction_enthalpy_fd_fallback_tabulated():
    """EquilibriumConstantTabulated.reaction_enthalpy() falls back to FD via van't Hoff."""
    dH_true = -20e3
    dS = -50.0
    T_data = np.linspace(270, 370, 200)
    K_data = np.exp(-dH_true / (R_GAS * T_data) + dS / R_GAS)
    eq = EquilibriumConstantTabulated(T_data, K_data)
    dH_fd = eq.reaction_enthalpy(298.15)
    assert abs(dH_fd - dH_true) / abs(dH_true) < 0.02


def test_reaction_enthalpy_fixed_K_is_zero():
    """EquilibriumConstant (fixed K) implies dlnK/dT = 0, so reaction_enthalpy = 0."""
    eq = EquilibriumConstant(K_eq=4.0)
    assert eq.reaction_enthalpy(298.15) == 0.0


# ---------------------------------------------------------------------------
# Polynomial equilibrium constant
# ---------------------------------------------------------------------------


def test_equilibrium_constant_polynomial_K_and_derivative():
    """EquilibriumConstantPolynomial: K and dlnK/dT match analytic values."""
    a0, a1, a2 = 2.0, -0.01, 3e-5
    eq = EquilibriumConstantPolynomial([a0, a1, a2])
    T = 310.0
    K_expected = np.exp(a0 + a1 * T + a2 * T**2)
    dlnK_expected = a1 + 2 * a2 * T
    dH_expected = R_GAS * T**2 * dlnK_expected
    assert abs(eq.K(T) - K_expected) / K_expected < 1e-12
    assert abs(eq.dlnK_dT(T) - dlnK_expected) < 1e-12
    assert abs(eq.reaction_enthalpy(T) - dH_expected) < 1e-6


def test_equilibrium_constant_polynomial_constant_case():
    """Single-coefficient polynomial is a temperature-independent K."""
    eq = EquilibriumConstantPolynomial([np.log(4.0)])
    assert abs(eq.K(298.15) - 4.0) < 1e-12
    assert eq.dlnK_dT(298.15) == 0.0
    assert eq.reaction_enthalpy(298.15) == 0.0


# ---------------------------------------------------------------------------
# Polynomial rate constant
# ---------------------------------------------------------------------------


def test_rate_constant_polynomial_kf_and_derivative():
    """RateConstantPolynomial: kf and dlnkf/dT match analytic values."""
    b0, b1 = 5.0, -0.005
    rc = RateConstantPolynomial([b0, b1])
    T = 320.0
    kf_expected = np.exp(b0 + b1 * T)
    dlnkf_expected = b1
    assert abs(rc.kf(T) - kf_expected) / kf_expected < 1e-12
    assert abs(rc.dlnkf_dT(T) - dlnkf_expected) < 1e-12


def test_polynomial_jacobian_dT_analytic_vs_fd():
    """jacobian_dT with polynomial K/k is analytic and matches FD."""
    a0, a1 = 1.5, -8e-3
    b0, b1 = 4.0, -5e-3
    model = ReactionModel(
        components=[Component("A"), Component("B")],
        reactions=[
            ThermodynamicReaction(
                "A <-> B",
                mode="kinetic",
                equilibrium_constant=EquilibriumConstantPolynomial([a0, a1]),
                rate_constant=RateConstantPolynomial([b0, b1]),
            ),
        ],
    )
    c = np.array([600.0, 400.0])
    T = 310.0
    analytic = model.jacobian_dT(c, np.zeros(2), T)
    eps = 1e-5
    fd = (model.residual(c, np.zeros(2), T + eps) -
          model.residual(c, np.zeros(2), T - eps)) / (2 * eps)
    np.testing.assert_allclose(analytic, fd, rtol=1e-4)


# ---------------------------------------------------------------------------
# d_reaction_enthalpy_dT
# ---------------------------------------------------------------------------


def test_d_reaction_enthalpy_dT_fixed_k():
    """Fixed K has zero temperature sensitivity — dΔH/dT = 0."""
    eq = EquilibriumConstant(K_eq=10.0)
    assert eq.d_reaction_enthalpy_dT(298.15) == 0.0


def test_d_reaction_enthalpy_dT_vanthoff():
    """VantHoff has constant ΔH — dΔH/dT = 0."""
    eq = EquilibriumConstantVantHoff(dH=-20e3, dS=-50.0)
    assert eq.d_reaction_enthalpy_dT(298.15) == 0.0
    assert eq.d_reaction_enthalpy_dT(350.0) == 0.0


def test_d_reaction_enthalpy_dT_vanthoffcp():
    """VantHoffCp: dΔH/dT = dCp."""
    dCp = 120.0
    eq = EquilibriumConstantVantHoffCp(dH=-20e3, dS=-50.0, dCp=dCp)
    assert eq.d_reaction_enthalpy_dT(298.15) == pytest.approx(dCp)
    assert eq.d_reaction_enthalpy_dT(350.0) == pytest.approx(dCp)


def test_d_reaction_enthalpy_dT_polynomial_analytic_vs_fd():
    """Polynomial d_reaction_enthalpy_dT analytic result matches FD on reaction_enthalpy."""
    coeffs = [5.0, -1e-2, 3e-5]
    eq = EquilibriumConstantPolynomial(coeffs=coeffs)
    for T in [280.0, 298.15, 320.0, 360.0]:
        analytic = eq.d_reaction_enthalpy_dT(T)
        eps = 1e-3
        fd = (eq.reaction_enthalpy(T + eps) - eq.reaction_enthalpy(T - eps)) / (2 * eps)
        assert analytic == pytest.approx(fd, rel=1e-5)


def test_d_reaction_enthalpy_dT_custom_fd_fallback():
    """Custom K uses FD fallback and returns a finite float."""
    eq = EquilibriumConstantCustom(func=lambda T: np.exp(-20e3 / (R_GAS * T)))
    result = eq.d_reaction_enthalpy_dT(298.15)
    assert np.isfinite(result)
