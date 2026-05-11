"""
Integration tests for the reactions library: model assembly, simulation, and solver.

Run with: pytest tests/ -v
Requires: pip install -e .
"""

import warnings

import numpy as np
import pytest

C_REF: float = 1000.0  # mol/m³ — standard-state concentration

from reactions.api import (
    ActivityCoefficientCustom,
    Component,
    EquilibriumConstant,
    EquilibriumConstantCustom,
    EquilibriumConstantVantHoff,
    EquilibriumConstantVantHoffCp,
    MassActionReaction,
    PhysicalState,
    R_GAS,
    RateConstantArrhenius,
    RateConstantFixed,
    ReactionModel,
    Species,
    ThermodynamicReaction,
    pKa,
)
from reactions.solver import simulate, solve_equilibrium


def test_reversible_ab_matches_analytical():
    """A <-> B numerical solution matches closed-form exponential relaxation."""
    kf_val, kr_val = 2.0, 0.5
    K = kf_val / kr_val
    A0, B0 = 1000.0, 0.0

    comp = Component("ab", [Species("A"), Species("B")])
    model = ReactionModel(
        components=[comp],
        reactions=[MassActionReaction("A <-> B", kf=kf_val, kr=kr_val)],
    )

    result = simulate(model, {"A": A0, "B": B0}, (0, 5.0))
    assert result.success

    total = A0 + B0
    A_eq = total / (1 + K)
    B_eq = total * K / (1 + K)
    lam = kf_val + kr_val
    A_ana = A_eq + (A0 - A_eq) * np.exp(-lam * result.t)
    B_ana = B_eq + (B0 - B_eq) * np.exp(-lam * result.t)

    assert np.max(np.abs(result["A"] - A_ana)) < 1e-6
    assert np.max(np.abs(result["B"] - B_ana)) < 1e-6


def test_reversible_ab_conservation():
    """A <-> B: total concentration is conserved throughout."""
    A0, B0 = 1000.0, 0.0
    comp = Component("ab", [Species("A"), Species("B")])
    model = ReactionModel(
        components=[comp],
        reactions=[MassActionReaction("A <-> B", kf=2.0, kr=0.5)],
    )
    result = simulate(model, {"A": A0, "B": B0}, (0, 5.0))
    assert result.success
    assert np.max(np.abs(result["A"] + result["B"] - (A0 + B0))) < 1e-6


def test_thermodynamic_consistency_converges_to_K():
    """ThermodynamicReaction kinetic: simulation reaches Q = K at long times."""
    K_eq, kf, A0, T = 4.0, 2.0, 1000.0, 298.15

    comp = Component("ab", [Species("A"), Species("B")])
    model = ReactionModel(
        components=[comp],
        reactions=[
            ThermodynamicReaction(
                "A <-> B",
                mode="kinetic",
                equilibrium_constant=EquilibriumConstant(K_eq=K_eq),
                rate_constant=RateConstantFixed(kf_value=kf),
            )
        ],
    )

    tau = C_REF / (kf * (1.0 + 1.0 / K_eq))
    result = simulate(model, {"A": A0, "B": 0.0}, (0, 10.0 * tau), T=T, n_points=500)
    assert result.success

    Q_final = result["B"][-1] / result["A"][-1]
    assert abs(Q_final - K_eq) / K_eq < 1e-3


def test_thermodynamic_consistency_kr_derived():
    """ThermodynamicReaction: kr is always kf / K, never stored independently."""
    K_eq, kf, T = 4.0, 2.0, 298.15
    comp = Component("ab", [Species("A"), Species("B")])
    model = ReactionModel(
        components=[comp],
        reactions=[
            ThermodynamicReaction(
                "A <-> B",
                mode="kinetic",
                equilibrium_constant=EquilibriumConstant(K_eq=K_eq),
                rate_constant=RateConstantFixed(kf_value=kf),
            )
        ],
    )
    rxn = model.reactions[0]
    assert abs(rxn.kr(T) - kf / K_eq) < 1e-12


def test_acid_base_henderson_hasselbalch():
    """Acetic acid speciation matches Henderson-Hasselbalch (ideal activity)."""
    pKa_val, c_tot, T = 4.76, 100.0, 298.15

    acetate = Component("acetate", [Species("AcOH", charge=0), Species("AcO-", charge=-1)])
    proton = Component("proton", [Species("H+", charge=+1)])
    hydroxide = Component("hydroxide", [Species("OH-", charge=-1)])
    water = Component("water", [Species("H2O", charge=0)])

    model = ReactionModel(
        components=[acetate, proton, hydroxide, water],
        reactions=[
            ThermodynamicReaction(
                "AcOH <-> AcO- + H+",
                mode="equil",
                equilibrium_constant=pKa(pKa_val),
            ),
            ThermodynamicReaction(
                "H2O <-> H+ + OH-",
                mode="equil",
                equilibrium_constant=EquilibriumConstant(1e-14),
            ),
        ],
        T=T,
    )

    errors = []
    for pH in np.linspace(2.0, 8.0, 13):
        H = 10.0 ** (-pH) * C_REF
        OH = 1e-14 * C_REF**2 / H
        AcO_HH = c_tot / (1.0 + 10.0 ** (pKa_val - pH))
        AcOH_HH = c_tot - AcO_HH
        c0 = {
            "AcOH": max(AcOH_HH, 1e-10),
            "AcO-": max(AcO_HH, 1e-10),
            "H+": H,
            "OH-": OH,
        }
        c_eq = solve_equilibrium(model, c0, T=T, prescribed={"H2O": C_REF})
        errors.append(abs(c_eq["AcOH"] - AcOH_HH))

    assert max(errors) < 1e-6


def test_conservation_multi_reaction():
    """A + B <-> C, C <-> D + E: atom balances conserved throughout."""
    comp = Component("all", [
        Species("A"), Species("B"), Species("C"), Species("D"), Species("E"),
    ])
    model = ReactionModel(
        components=[comp],
        reactions=[
            MassActionReaction("A + B -> C", kf=1.0, kr=0.1),
            MassActionReaction("C -> D + E", kf=0.5, kr=0.2),
        ],
    )

    c0 = {"A": 500.0, "B": 300.0, "C": 0.0, "D": 0.0, "E": 0.0}
    result = simulate(model, c0, (0, 10.0))
    assert result.success

    C1 = result["A"] + result["C"] + result["D"]
    C2 = result["B"] + result["C"] + result["E"]
    assert np.max(np.abs(C1 - C1[0])) < 1e-6
    assert np.max(np.abs(C2 - C2[0])) < 1e-6


# ---------------------------------------------------------------------------
# ThermodynamicReaction reduces to MassActionReaction under ideal conditions
#
# ThermodynamicReaction uses activities a_i = c_i / C_REF (ideal, γ=1):
#   v = kf_thermo * prod(a_i^e_fwd) = (kf_thermo / C_REF^n) * prod(c_i^e_fwd)
# MassActionReaction uses concentrations directly:
#   v = kf_MA * prod(c_i^e_fwd)
# Equivalence requires kf_thermo = kf_MA * C_REF^n.
# For the equilibrium constant: K_thermo = K_conc * C_REF^(n_fwd - n_bwd),
# where n_fwd and n_bwd are the sums of forward and backward exponents.
# ---------------------------------------------------------------------------


def test_thermodynamic_reduces_to_mass_action_n1():
    """n=1 (A⇌B): ThermodynamicReaction with kf*C_REF matches MassActionReaction exactly."""
    kf_MA, kr_MA = 2.0, 0.5
    K_conc = kf_MA / kr_MA

    comp = Component("ab", [Species("A"), Species("B")])
    model_ma = ReactionModel(
        components=[comp],
        reactions=[MassActionReaction("A <-> B", kf=kf_MA, kr=kr_MA)],
    )
    model_td = ReactionModel(
        components=[comp],
        reactions=[
            ThermodynamicReaction(
                "A <-> B",
                mode="kinetic",
                equilibrium_constant=EquilibriumConstant(K_eq=K_conc),
                rate_constant=RateConstantFixed(kf_value=kf_MA * C_REF),
            )
        ],
    )

    c0 = {"A": 800.0, "B": 200.0}
    result_ma = simulate(model_ma, c0, (0, 2.0))
    result_td = simulate(model_td, c0, (0, 2.0))

    assert np.max(np.abs(result_ma["A"] - result_td["A"])) < 1e-6
    assert np.max(np.abs(result_ma["B"] - result_td["B"])) < 1e-6


def test_thermodynamic_reduces_to_mass_action_n2():
    """n=2 (2A⇌B): ThermodynamicReaction with kf*C_REF^2, K*C_REF matches MassActionReaction."""
    kf_MA = 1e-3
    kr_MA = 2e-4
    K_conc = kf_MA / kr_MA

    comp = Component("ab", [Species("A"), Species("B")])
    model_ma = ReactionModel(
        components=[comp],
        reactions=[MassActionReaction("2 A <-> B", kf=kf_MA, kr=kr_MA)],
    )
    model_td = ReactionModel(
        components=[comp],
        reactions=[
            ThermodynamicReaction(
                "2 A <-> B",
                mode="kinetic",
                equilibrium_constant=EquilibriumConstant(K_eq=K_conc * C_REF),
                rate_constant=RateConstantFixed(kf_value=kf_MA * C_REF**2),
            )
        ],
    )

    c0 = {"A": 400.0, "B": 100.0}
    result_ma = simulate(model_ma, c0, (0, 5.0))
    result_td = simulate(model_td, c0, (0, 5.0))

    assert np.max(np.abs(result_ma["A"] - result_td["A"])) < 1e-4
    assert np.max(np.abs(result_ma["B"] - result_td["B"])) < 1e-4


# ---------------------------------------------------------------------------
# Helper
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


# ---------------------------------------------------------------------------
# 1. Jacobian verification
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
    assert err < 1e-4, f"ThermodynamicReaction kinetic Jacobian error {err:.2e} exceeds 1e-4"


def test_jacobian_thermodynamic_equil_acetic():
    """Analytic Jacobian of ThermodynamicReaction (equil, acetic acid) matches FD."""
    acetate = Component("acetate", [Species("AcOH", charge=0), Species("AcO-", charge=-1)])
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
    assert err < 1e-4, f"Acetic acid equil Jacobian relative error {err:.2e} exceeds 1e-4"


# ---------------------------------------------------------------------------
# 2. Phosphate polyprotic: Bjerrum fraction validation at pH 7.2
# ---------------------------------------------------------------------------


def test_phosphate_bjerrum_fractions_ph72():
    """
    Three-step phosphate dissociation: species concentrations at pH 7.2
    match the Bjerrum fraction formula.
    """
    pKa1, pKa2, pKa3 = 2.148, 7.198, 12.350
    Ka1 = 10.0 ** (-pKa1)
    Ka2 = 10.0 ** (-pKa2)
    Ka3 = 10.0 ** (-pKa3)
    c_tot_phos = 100.0
    T = 298.15
    pH = 7.2

    phosphate = Component("phosphate", [
        Species("H3PO4",  charge=0),
        Species("H2PO4-", charge=-1),
        Species("HPO4-2", charge=-2),
        Species("PO4-3",  charge=-3),
    ])
    proton = Component("proton", [Species("H+", charge=+1)])
    hydroxide = Component("hydroxide", [Species("OH-", charge=-1)])
    water = Component("water", [Species("H2O", charge=0)])

    model_phos = ReactionModel(
        components=[phosphate, proton, hydroxide, water],
        reactions=[
            ThermodynamicReaction(
                "H3PO4 <-> H2PO4- + H+",
                mode="equil",
                equilibrium_constant=pKa(pKa1),
            ),
            ThermodynamicReaction(
                "H2PO4- <-> HPO4-2 + H+",
                mode="equil",
                equilibrium_constant=pKa(pKa2),
            ),
            ThermodynamicReaction(
                "HPO4-2 <-> PO4-3 + H+",
                mode="equil",
                equilibrium_constant=pKa(pKa3),
            ),
            ThermodynamicReaction(
                "H2O <-> H+ + OH-",
                mode="equil",
                equilibrium_constant=EquilibriumConstant(1e-14),
            ),
        ],
        T=T,
    )

    h = 10.0 ** (-pH)
    D = h**3 + Ka1 * h**2 + Ka1 * Ka2 * h + Ka1 * Ka2 * Ka3
    alpha0 = h**3 / D
    alpha1 = Ka1 * h**2 / D
    alpha2 = Ka1 * Ka2 * h / D
    alpha3 = Ka1 * Ka2 * Ka3 / D

    H = 10.0 ** (-pH) * C_REF
    OH = 1e-14 * C_REF**2 / H
    c0 = {
        "H3PO4":  max(alpha0 * c_tot_phos, 1e-10),
        "H2PO4-": max(alpha1 * c_tot_phos, 1e-10),
        "HPO4-2": max(alpha2 * c_tot_phos, 1e-10),
        "PO4-3":  max(alpha3 * c_tot_phos, 1e-10),
        "H+": H,
        "OH-": OH,
    }
    c_eq = solve_equilibrium(model_phos, c0, T=T, prescribed={"H2O": C_REF})

    tol = 1e-4
    for sp, alpha in [
        ("H3PO4",  alpha0),
        ("H2PO4-", alpha1),
        ("HPO4-2", alpha2),
        ("PO4-3",  alpha3),
    ]:
        expected = alpha * c_tot_phos
        err = abs(c_eq[sp] - expected)
        assert err < tol, (
            f"{sp}: solver {c_eq[sp]:.6f}, Bjerrum {expected:.6f}, err {err:.2e}"
        )


# ---------------------------------------------------------------------------
# 3. Edge cases
# ---------------------------------------------------------------------------


def test_extreme_ph_acid_ph1():
    """Acetic acid at pH 1: essentially fully protonated; equilibrium solver converges."""
    pKa_val = 4.76
    c_tot = 100.0
    T = 298.15

    acetate = Component("acetate", [Species("AcOH", charge=0), Species("AcO-", charge=-1)])
    proton = Component("proton", [Species("H+", charge=+1)])
    hydroxide = Component("hydroxide", [Species("OH-", charge=-1)])
    water = Component("water", [Species("H2O", charge=0)])

    model = ReactionModel(
        components=[acetate, proton, hydroxide, water],
        reactions=[
            ThermodynamicReaction(
                "AcOH <-> AcO- + H+",
                mode="equil",
                equilibrium_constant=pKa(pKa_val),
            ),
            ThermodynamicReaction(
                "H2O <-> H+ + OH-",
                mode="equil",
                equilibrium_constant=EquilibriumConstant(1e-14),
            ),
        ],
        T=T,
    )

    pH = 1.0
    H = 10.0 ** (-pH) * C_REF
    OH = 1e-14 * C_REF**2 / H
    AcO_hh = c_tot / (1.0 + 10.0 ** (pKa_val - pH))
    AcOH_hh = c_tot - AcO_hh
    c0 = {
        "AcOH": max(AcOH_hh, 1e-10),
        "AcO-": max(AcO_hh, 1e-10),
        "H+": H,
        "OH-": OH,
    }
    c_eq = solve_equilibrium(model, c0, T=T, prescribed={"H2O": C_REF})

    frac_deprotonated = c_eq["AcO-"] / c_tot
    assert frac_deprotonated < 1e-3
    assert abs(c_eq["AcOH"] - AcOH_hh) < 1e-4


def test_extreme_ph_acid_ph13():
    """Acetic acid at pH 13: essentially fully deprotonated; equilibrium solver converges."""
    pKa_val = 4.76
    c_tot = 100.0
    T = 298.15

    acetate = Component("acetate", [Species("AcOH", charge=0), Species("AcO-", charge=-1)])
    proton = Component("proton", [Species("H+", charge=+1)])
    hydroxide = Component("hydroxide", [Species("OH-", charge=-1)])
    water = Component("water", [Species("H2O", charge=0)])

    model = ReactionModel(
        components=[acetate, proton, hydroxide, water],
        reactions=[
            ThermodynamicReaction(
                "AcOH <-> AcO- + H+",
                mode="equil",
                equilibrium_constant=pKa(pKa_val),
            ),
            ThermodynamicReaction(
                "H2O <-> H+ + OH-",
                mode="equil",
                equilibrium_constant=EquilibriumConstant(1e-14),
            ),
        ],
        T=T,
    )

    pH = 13.0
    H = 10.0 ** (-pH) * C_REF
    OH = 1e-14 * C_REF**2 / H
    AcO_hh = c_tot / (1.0 + 10.0 ** (pKa_val - pH))
    AcOH_hh = c_tot - AcO_hh
    c0 = {
        "AcOH": max(AcOH_hh, 1e-10),
        "AcO-": max(AcO_hh, 1e-10),
        "H+": H,
        "OH-": OH,
    }
    c_eq = solve_equilibrium(model, c0, T=T, prescribed={"H2O": C_REF})

    frac_protonated = c_eq["AcOH"] / c_tot
    assert frac_protonated < 1e-7


# ---------------------------------------------------------------------------
# Jacobian vs finite differences
# ---------------------------------------------------------------------------


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
                    dH=-20e3, dS=-50.0, dCp=100.0,
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
# Energy balance (coupled T simulation)
# ---------------------------------------------------------------------------

# Water: rho=997 kg/m³, M=0.018 kg/mol, Cp=75.3 J/(mol·K)
# -> rho_cp = (997/0.018)*75.3 ≈ 4.175e6 J/(m³·K)
_WATER = Species("H2O", charge=0,
                 molar_mass=0.018, density=997.0, heat_capacity=75.3)
_WATER_COMPONENT = Component("water", [_WATER])
_WATER_X = {"H2O": 1.0}
_RHO_CP_WATER = (997.0 / 0.018) * 75.3


def test_volumetric_heat_capacity_single_solvent():
    """volumetric_heat_capacity returns correct value for pure water."""
    model = ReactionModel(
        components=[Component("A"), _WATER_COMPONENT],
        reactions=[
            ThermodynamicReaction(
                "A <-> A",
                mode="kinetic",
                equilibrium_constant=EquilibriumConstant(1.0),
                rate_constant=RateConstantFixed(1.0),
            ),
        ],
    )
    c = np.array([0.0, _WATER.c_ref])
    rho_cp = model.volumetric_heat_capacity(c)
    np.testing.assert_allclose(rho_cp, _RHO_CP_WATER, rtol=1e-10)


def test_volumetric_heat_capacity_independent_of_solute():
    """volumetric_heat_capacity is unaffected by solute concentration."""
    model = ReactionModel(
        components=[Component("A"), _WATER_COMPONENT],
        reactions=[
            ThermodynamicReaction(
                "A <-> A",
                mode="kinetic",
                equilibrium_constant=EquilibriumConstant(1.0),
                rate_constant=RateConstantFixed(1.0),
            ),
        ],
    )
    c_no_solute = np.array([0.0, _WATER.c_ref])
    c_with_solute = np.array([500.0, _WATER.c_ref])
    assert model.volumetric_heat_capacity(c_no_solute) == model.volumetric_heat_capacity(c_with_solute)


def test_volumetric_heat_capacity_mixture():
    """volumetric_heat_capacity weighted sum for two-solvent mixture."""
    water = Species("H2O",  molar_mass=0.018, density=997.0, heat_capacity=75.3)
    mecn  = Species("MeCN", molar_mass=0.041, density=786.0, heat_capacity=91.4)
    model = ReactionModel(
        components=[
            Component("A"),
            Component("water",  [water]),
            Component("MeCN",   [mecn]),
        ],
        reactions=[
            ThermodynamicReaction(
                "A <-> A",
                mode="kinetic",
                equilibrium_constant=EquilibriumConstant(1.0),
                rate_constant=RateConstantFixed(1.0),
            ),
        ],
    )
    c = np.array([0.0, 0.7 * water.c_ref, 0.3 * mecn.c_ref])
    expected = 0.7 * (997.0 / 0.018) * 75.3 + 0.3 * (786.0 / 0.041) * 91.4
    np.testing.assert_allclose(model.volumetric_heat_capacity(c), expected, rtol=1e-10)


def _make_ab_model_with_water(K, kf):
    """ReactionModel for A <-> B in water, used by energy balance tests."""
    return ReactionModel(
        components=[Component("A"), Component("B"), _WATER_COMPONENT],
        reactions=[
            ThermodynamicReaction(
                "A <-> B",
                mode="kinetic",
                equilibrium_constant=K,
                rate_constant=kf,
            ),
        ],
    )


def test_coupled_energy_balance_adiabatic():
    """Adiabatic energy conservation: rho_cp * dT = -dH * dc_B."""
    dH = -20e3
    T0 = 298.15
    model = _make_ab_model_with_water(
        EquilibriumConstantVantHoff(dH=dH, dS=-50.0),
        RateConstantFixed(kf_value=1e5),
    )
    result = simulate(
        model,
        c0={"A": 1000.0, "B": 0.0},
        t_span=(0, 0.5),
        T=T0,
        solvent_composition=_WATER_X,
    )
    assert result.success
    dc_B = result["B"][-1]
    dT_sim = result.T_profile[-1] - T0
    dT_exp = -dH * dc_B / _RHO_CP_WATER
    np.testing.assert_allclose(dT_sim, dT_exp, rtol=1e-6)


def test_coupled_mass_balance():
    """Mass balance holds throughout the coupled energy-balance simulation."""
    c_tot = 1000.0
    model = _make_ab_model_with_water(
        EquilibriumConstantVantHoff(dH=-20e3, dS=-50.0),
        RateConstantArrhenius(A=1e10, Ea=40e3),
    )
    result = simulate(
        model,
        c0={"A": c_tot, "B": 0.0},
        t_span=(0, 10.0),
        T=298.15,
        solvent_composition=_WATER_X,
    )
    total = result["A"] + result["B"]
    np.testing.assert_allclose(total, c_tot, rtol=1e-8)


def test_coupled_T_profile_populated():
    """T_profile is always populated when solvent_composition is given."""
    model = _make_ab_model_with_water(
        EquilibriumConstantVantHoff(dH=-20e3, dS=-50.0),
        RateConstantFixed(kf_value=1e5),
    )
    result = simulate(
        model,
        c0={"A": 1000.0, "B": 0.0},
        t_span=(0, 5.0),
        T=298.15,
        solvent_composition=_WATER_X,
    )
    assert result.T_profile is not None
    assert len(result.T_profile) == len(result.t)


def test_coupled_callable_T_raises():
    """Combining solvent_composition with callable T raises ValueError."""
    model = _make_ab_model_with_water(
        EquilibriumConstant(4.0),
        RateConstantFixed(2000.0),
    )
    with pytest.raises(ValueError, match="solvent_composition cannot be combined"):
        simulate(
            model,
            c0={"A": 1000.0},
            t_span=(0, 1.0),
            T=lambda t: 298.15 + t,
            solvent_composition=_WATER_X,
        )


def test_coupled_callable_solvent_composition():
    """Callable solvent_composition (gradient) is evaluated at each time step."""
    water = Species("H2O",  molar_mass=0.018, density=997.0, heat_capacity=75.3)
    mecn  = Species("MeCN", molar_mass=0.041, density=786.0, heat_capacity=91.4)
    model = ReactionModel(
        components=[
            Component("A"), Component("B"),
            Component("water", [water]),
            Component("MeCN",  [mecn]),
        ],
        reactions=[
            ThermodynamicReaction(
                "A <-> B",
                mode="kinetic",
                equilibrium_constant=EquilibriumConstantVantHoff(dH=-20e3, dS=-50.0),
                rate_constant=RateConstantFixed(kf_value=1e5),
            ),
        ],
    )
    t_end = 5.0
    x_comp = {
        "H2O":  lambda t: 1.0 - 0.3 * t / t_end,
        "MeCN": lambda t: 0.3 * t / t_end,
    }
    result = simulate(
        model,
        c0={"A": 1000.0, "B": 0.0},
        t_span=(0, t_end),
        T=298.15,
        solvent_composition=x_comp,
    )
    assert result.success
    assert result.T_profile is not None


def test_energy_balance_custom_K_matches_vanthoff():
    """Energy balance with EquilibriumConstantCustom matches VantHoff result."""
    dH_true = -20e3
    dS = -50.0
    K_fn = lambda T: float(np.exp(-dH_true / (R_GAS * T) + dS / R_GAS))  # noqa: E731
    model_custom = ReactionModel(
        components=[Component("A"), Component("B"), Component("water", [_WATER])],
        reactions=[
            ThermodynamicReaction(
                "A <-> B",
                mode="kinetic",
                equilibrium_constant=EquilibriumConstantCustom(K_fn),
                rate_constant=RateConstantArrhenius(A=1e10, Ea=40e3),
            ),
        ],
    )
    model_vH = ReactionModel(
        components=[Component("A"), Component("B"), Component("water", [_WATER])],
        reactions=[
            ThermodynamicReaction(
                "A <-> B",
                mode="kinetic",
                equilibrium_constant=EquilibriumConstantVantHoff(dH=dH_true, dS=dS),
                rate_constant=RateConstantArrhenius(A=1e10, Ea=40e3),
            ),
        ],
    )
    kw = dict(c0={"A": 1000.0}, t_span=(0, 2.0), T=298.15, solvent_composition=_WATER_X)
    r_custom = simulate(model_custom, **kw)
    r_vH    = simulate(model_vH,    **kw)
    assert r_custom.success and r_vH.success
    np.testing.assert_allclose(r_custom.T_profile, r_vH.T_profile, rtol=1e-4)


def test_energy_balance_mass_action_warns():
    """MassActionReaction in an energy balance model raises UserWarning."""
    model = ReactionModel(
        components=[Component("A"), Component("B"), Component("water", [_WATER])],
        reactions=[MassActionReaction("A <-> B", kf=1.0, kr=0.25)],
    )
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        simulate(model, c0={"A": 1000.0}, t_span=(0, 0.1), T=298.15,
                 solvent_composition=_WATER_X)
    assert any(issubclass(warning.category, UserWarning) for warning in w)
    assert any("MassActionReaction" in str(warning.message) for warning in w)


# ---------------------------------------------------------------------------
# Bundle A: d_reaction_enthalpy_dT and analytic energy balance Jacobian row
# ---------------------------------------------------------------------------


def _eb_rhs_T(model, c, T, rho_cp):
    """Energy balance scalar rhs: -Q_dot/rho_cp = -Σ_j ΔH_j φ_j / rho_cp."""
    state = model.make_state(c, T)
    Q_dot = 0.0
    for j, rxn in enumerate(model.reactions):
        if not model.kinetic_mask[j]:
            continue
        eq = getattr(rxn, "equilibrium_constant", None)
        if eq is None:
            continue
        Q_dot += eq.reaction_enthalpy(T) * rxn.net_rate(state, model.species_index, model.charges)
    return -Q_dot / rho_cp


def _eb_jac_analytic(model, c, T, rho_cp):
    """Analytic energy balance Jacobian row (∂rhs_T/∂c, ∂rhs_T/∂T)."""
    state = model.make_state(c, T)
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
        phi = rxn.net_rate(state, model.species_index, model.charges)
        dphi_dc = rxn.net_rate_jac(state, model.species_index, model.charges)
        dphi_dT = rxn.net_rate_dT(state, model.species_index, model.charges)
        jac_c -= dH / rho_cp * dphi_dc
        jac_T -= (dH * dphi_dT + phi * dH_dT) / rho_cp
    return jac_c, jac_T


def test_energy_balance_jac_analytic_vanthoff_arrhenius():
    """Analytic J[n, :n] and J[n, n] match FD for VantHoff K + Arrhenius kf."""
    water = Species(
        "H2O", molar_mass=0.018, density=1000.0, heat_capacity=75.3
    )
    model = ReactionModel(
        components=[Component("A"), Component("B"), Component("water", [water])],
        reactions=[ThermodynamicReaction(
            "A <-> B",
            mode="kinetic",
            equilibrium_constant=EquilibriumConstantVantHoff(dH=-20e3, dS=-50.0),
            rate_constant=RateConstantArrhenius(A=1e10, Ea=40e3),
        )],
    )
    c = np.array([800.0, 200.0, 1000.0])
    rho_cp = model.volumetric_heat_capacity(c)
    T = 310.0
    eps = 1e-5

    rhs0 = _eb_rhs_T(model, c, T, rho_cp)
    jac_c_fd = np.array([
        (_eb_rhs_T(model, c + np.eye(len(c))[k] * eps, T, rho_cp) - rhs0) / eps
        for k in range(len(c))
    ])
    jac_T_fd = (_eb_rhs_T(model, c, T + eps, rho_cp) - rhs0) / eps

    jac_c, jac_T = _eb_jac_analytic(model, c, T, rho_cp)
    np.testing.assert_allclose(jac_c, jac_c_fd, rtol=1e-5)
    np.testing.assert_allclose(jac_T, jac_T_fd, rtol=1e-5)


def test_energy_balance_jac_analytic_vanthoffcp():
    """With VantHoffCp, dΔH/dT = dCp contributes non-trivially to J[n, n]."""
    dCp = 150.0
    water = Species(
        "H2O", molar_mass=0.018, density=1000.0, heat_capacity=75.3
    )
    model = ReactionModel(
        components=[Component("A"), Component("B"), Component("water", [water])],
        reactions=[ThermodynamicReaction(
            "A <-> B",
            mode="kinetic",
            equilibrium_constant=EquilibriumConstantVantHoffCp(
                dH=-20e3, dS=-50.0, dCp=dCp
            ),
            rate_constant=RateConstantArrhenius(A=1e10, Ea=40e3),
        )],
    )
    c = np.array([600.0, 400.0, 1000.0])
    rho_cp = model.volumetric_heat_capacity(c)
    T = 305.0
    eps = 1e-5

    rhs0 = _eb_rhs_T(model, c, T, rho_cp)
    jac_T_fd = (_eb_rhs_T(model, c, T + eps, rho_cp) - rhs0) / eps
    _, jac_T = _eb_jac_analytic(model, c, T, rho_cp)
    np.testing.assert_allclose(jac_T, jac_T_fd, rtol=1e-5)


# ---------------------------------------------------------------------------
# Solvent species in state.c (E/P3)
# ---------------------------------------------------------------------------


def test_prescribed_at_call_site_holds_constant():
    """prescribed kwarg in simulate holds a species constant at the given value."""
    model = ReactionModel(
        components=[Component("A"), Component("water", [Species("H2O")])],
        reactions=[MassActionReaction("A <-> A", kf=1.0)],
    )
    result = simulate(model, {"A": 100.0}, (0.0, 1.0), prescribed={"H2O": C_REF})
    np.testing.assert_allclose(result["H2O"], C_REF, rtol=1e-10)


def test_solvent_included_in_state_c():
    """make_state includes all species, including solvents, in state.c."""
    water = Component("water", [Species("H2O")])
    model = ReactionModel(
        components=[Component("A"), water],
        reactions=[MassActionReaction("A <-> A", kf=1.0)],
    )
    c = np.array([100.0, 1000.0])
    state = model.make_state(c, T=300.0)
    assert len(state.c) == 2
    np.testing.assert_array_equal(state.c, c)


def test_custom_ac_sees_solvent_in_state_c():
    """Custom activity coefficient receives solvent concentration in state.c."""
    captured_c_h2o = []
    h2o_idx_ref = [None]

    def custom_ac(state, charges):
        if h2o_idx_ref[0] is not None:
            captured_c_h2o.append(float(state.c[h2o_idx_ref[0]]))
        return np.ones(len(state.c))

    model = ReactionModel(
        components=[Component("A"), Component("B"), Component("water", [_WATER])],
        reactions=[ThermodynamicReaction(
            "A <-> B",
            mode="kinetic",
            equilibrium_constant=EquilibriumConstantVantHoff(dH=-20e3, dS=-50.0),
            rate_constant=RateConstantFixed(1000.0),
            activity_coefficient=ActivityCoefficientCustom(fn=custom_ac),
        )],
    )
    h2o_idx_ref[0] = model.species_index["H2O"]
    simulate(
        model,
        c0={"A": 1000.0},
        t_span=(0, 0.05),
        T=298.15,
        solvent_composition=_WATER_X,
        n_points=10,
    )
    assert len(captured_c_h2o) > 0
    np.testing.assert_allclose(captured_c_h2o, _WATER.c_ref, rtol=1e-6)


def test_gradient_solvent_varies_in_state_c():
    """Callable solvent_composition causes c[H2O] to vary in state.c during simulation."""
    captured_c_h2o = []
    h2o_idx_ref = [None]

    def custom_ac(state, charges):
        if h2o_idx_ref[0] is not None:
            captured_c_h2o.append(float(state.c[h2o_idx_ref[0]]))
        return np.ones(len(state.c))

    def x_water(t):
        return 1.0 - 0.5 * t

    model = ReactionModel(
        components=[Component("A"), Component("B"), Component("water", [_WATER])],
        reactions=[ThermodynamicReaction(
            "A <-> B",
            mode="kinetic",
            equilibrium_constant=EquilibriumConstantVantHoff(dH=-20e3, dS=-50.0),
            rate_constant=RateConstantFixed(1000.0),
            activity_coefficient=ActivityCoefficientCustom(fn=custom_ac),
        )],
    )
    h2o_idx_ref[0] = model.species_index["H2O"]
    simulate(
        model,
        c0={"A": 1000.0},
        t_span=(0, 1.0),
        T=298.15,
        solvent_composition={"H2O": x_water},
        n_points=10,
    )
    assert any(c < 0.99 * _WATER.c_ref for c in captured_c_h2o)


# ---------------------------------------------------------------------------
# Prescribed species (E/P2)
# ---------------------------------------------------------------------------


def test_prescribed_species_constant():
    """Prescribed constant species: value is fixed, ODE species responds correctly."""
    k = 0.1
    A0 = 1000.0

    comp = Component("ab", [Species("A"), Species("B")])
    model = ReactionModel(
        components=[comp],
        reactions=[MassActionReaction("A -> B", kf=k, kr=0.0)],
    )

    result = simulate(
        model,
        {"B": 0.0},
        (0.0, 5.0),
        prescribed={"A": lambda t: A0},
        n_points=200,
    )
    assert result.success

    np.testing.assert_allclose(result["A"], A0, rtol=1e-6)

    B_ana = k * A0 * result.t
    np.testing.assert_allclose(result["B"], B_ana, rtol=1e-4)


def test_prescribed_species_ramp():
    """Prescribed ramp species: A(t) = A0 - r*t drives B accumulation."""
    k = 0.1
    A0, r = 1000.0, 50.0

    comp = Component("ab", [Species("A"), Species("B")])
    model = ReactionModel(
        components=[comp],
        reactions=[MassActionReaction("A -> B", kf=k, kr=0.0)],
    )

    t_end = 4.0
    result = simulate(
        model,
        {"B": 0.0},
        (0.0, t_end),
        prescribed={"A": lambda t: A0 - r * t},
        n_points=200,
    )
    assert result.success

    np.testing.assert_allclose(result["A"], A0 - r * result.t, rtol=1e-6)

    B_ana = k * (A0 * result.t - 0.5 * r * result.t**2)
    np.testing.assert_allclose(result["B"], B_ana, rtol=1e-4)


# ---------------------------------------------------------------------------
# Solution (formulation helper)
# ---------------------------------------------------------------------------


def _water():
    return Component("water", [Species("H2O", charge=0, molar_mass=0.018015, density=1000.0)])


def test_solution_pure_solvent_c0():
    """Pure solvent: c0['H2O'] == water.c_ref."""
    water = _water()
    from reactions.api import Solution
    sol = Solution(water)
    assert sol.c0["H2O"] == pytest.approx(water.c_ref, rel=1e-9)
    assert sol.prescribed["H2O"] == pytest.approx(water.c_ref, rel=1e-9)


def test_solution_solutes_added():
    """Solutes appear in c0 but not in prescribed."""
    water = _water()
    from reactions.api import Solution
    sol = Solution(water, solutes={"H+": 1e-4, "OH-": 1e-7})
    assert sol.c0["H+"] == pytest.approx(1e-4)
    assert sol.c0["OH-"] == pytest.approx(1e-7)
    assert "H+" not in sol.prescribed
    assert "OH-" not in sol.prescribed


def test_solution_mixed_solvents():
    """Two-solvent mix: each concentration is phi * c_ref."""
    ethanol = Component(
        "ethanol",
        [Species("EtOH", charge=0, molar_mass=0.04607, density=789.0)],
    )
    water = _water()
    from reactions.api import Solution
    sol = Solution({water: 0.9, ethanol: 0.1})
    assert sol.c0["H2O"] == pytest.approx(0.9 * water.c_ref, rel=1e-9)
    assert sol.c0["EtOH"] == pytest.approx(0.1 * ethanol.c_ref, rel=1e-9)
    assert set(sol.prescribed.keys()) == {"H2O", "EtOH"}


def test_solution_fractions_must_sum_to_one():
    water = _water()
    from reactions.api import Solution
    with pytest.raises(ValueError, match="sum to 1"):
        Solution({water: 0.5})


def test_solution_missing_density_raises():
    no_density = Component("X", [Species("X")])
    from reactions.api import Solution
    with pytest.raises(ValueError, match="density"):
        Solution(no_density)


def test_solution_c0_and_prescribed_are_copies():
    """Mutating the returned dict does not affect the Solution."""
    water = _water()
    from reactions.api import Solution
    sol = Solution(water)
    c0 = sol.c0
    c0["H2O"] = 0.0
    assert sol.c0["H2O"] != 0.0
