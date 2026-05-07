"""
Validation tests for the reactions library.

Run with: pytest tests/ -v
Requires: pip install -e .
"""

import numpy as np
import pytest

from reactions.api import (
    C_REF,
    ActivityCoefficientDavies,
    Component,
    EquilibriumConstant,
    EquilibriumConstantCustom,
    EquilibriumConstantPolynomial,
    EquilibriumConstantTabulated,
    EquilibriumConstantVantHoff,
    EquilibriumConstantVantHoffCp,
    IonicStrengthBackground,
    IonicStrengthFixed,
    MassActionReaction,
    PhysicalState,
    R_GAS,
    RateConstantArrhenius,
    RateConstantFixed,
    RateConstantPolynomial,
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
    water = Component("water", [Species("H2O", charge=0, is_solvent=True)])

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
        c_eq = solve_equilibrium(model, c0, T=T)
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
    K_conc = kf_MA / kr_MA   # K unchanged for n_fwd=n_bwd=1

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
    kf_MA = 1e-3   # m³/(mol·s)
    kr_MA = 2e-4   # 1/s
    K_conc = kf_MA / kr_MA
    # n_fwd=2, n_bwd=1 → K_thermo = K_conc * C_REF^(2-1) = K_conc * C_REF

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
    # Relative error where |J_ana| > 1; absolute error elsewhere
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
    water = Component("water", [Species("H2O", charge=0, is_solvent=True)])
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
    # Use concentrations from the implementation notebook (pH ~ 5.15, above pKa):
    # [AcOH]=10, [AcO-]=90, [H+]=6.3e-5 mol/m3, [OH-]=1.6e-4 mol/m3
    # These give moderate Jacobian entries suitable for finite-difference comparison.
    c = np.array([10.0, 90.0, 6.3e-5, 1.6e-4])
    err = _check_jacobian(model, c)
    assert err < 1e-4, f"Acetic acid equil Jacobian relative error {err:.2e} exceeds 1e-4"


# ---------------------------------------------------------------------------
# 2. Davies activity corrections
# ---------------------------------------------------------------------------


def test_davies_ph_shift_direction():
    """
    Davies corrections shift the apparent pKa in the correct direction.

    The thermodynamic equilibrium constant is defined in terms of activities:
        Ka = a(AcO-) * a(H+) / a(AcOH)
           = gamma(AcO-) * gamma(H+) * [AcO-] * [H+] / (gamma(AcOH) * [AcOH] * C_REF)

    With gamma(AcOH) = 1 (neutral) and gamma(AcO-) = gamma(H+) = gamma < 1:
        Ka_apparent = [AcO-] * [H+] / ([AcOH] * C_REF) = Ka / gamma^2

    So pKa_app = pKa - log10(gamma^2) = pKa - 2*log10(gamma) = pKa - 2*A*f(I)
    where A = 0.509 and f(I) = sqrt(I_L)/(1+sqrt(I_L)) - 0.3*I_L.

    The test measures the concentration-pH at half-neutralisation ([AcO-]=[AcOH])
    by starting from that half-half initial condition and letting the solver
    adjust the dependent species (H+ and AcO-) to satisfy both equilibrium
    constraints. The non-dependent species (AcOH and OH-) are pinned by the
    initial condition, so the solver cleanly isolates the pKa shift.
    """
    pKa_val = 4.756
    T = 298.15
    A_davies = 0.509

    acetate = Component("acetate", [Species("AcOH", charge=0), Species("AcO-", charge=-1)])
    proton = Component("proton", [Species("H+", charge=+1)])
    hydroxide = Component("hydroxide", [Species("OH-", charge=-1)])
    water = Component("water", [Species("H2O", charge=0, is_solvent=True)])

    def make_model(I_fixed):
        return ReactionModel(
            components=[acetate, proton, hydroxide, water],
            reactions=[
                ThermodynamicReaction(
                    "AcOH <-> AcO- + H+",
                    mode="equil",
                    equilibrium_constant=pKa(pKa_val),
                    activity_coefficient=ActivityCoefficientDavies(),
                ),
                ThermodynamicReaction(
                    "H2O <-> H+ + OH-",
                    mode="equil",
                    equilibrium_constant=EquilibriumConstant(1e-14),
                    activity_coefficient=ActivityCoefficientDavies(),
                ),
            ],
            ionic_strength=IonicStrengthFixed(I=I_fixed),
            T=T,
        )

    Ka = 10.0 ** (-pKa_val)
    # Initial guess for H+ at ideal pKa; OH- from water equilibrium
    H_guess = Ka * C_REF
    OH_guess = 1e-14 * C_REF**2 / H_guess
    # Start at exact half-neutralisation: [AcO-] = [AcOH] = 50 mol/m3
    # The solver pins AcOH and OH- (the non-dependent species) and adjusts
    # H+ and AcO- to satisfy both equilibrium constraints.
    c0 = {
        "AcOH": 50.0,
        "AcO-": 50.0,
        "H+": H_guess,
        "OH-": OH_guess,
    }

    # Ideal model (I=0): H+ at equilibrium should equal Ka * C_REF
    model_ideal = make_model(I_fixed=0.0)
    c_eq_ideal = solve_equilibrium(model_ideal, c0, T=T)
    pH_ideal = -np.log10(c_eq_ideal["H+"] / C_REF)

    # Davies model (I=100 mol/m3 = 0.1 mol/L)
    I_bg = 100.0   # mol/m3
    model_ionic = make_model(I_fixed=I_bg)
    c_eq_davies = solve_equilibrium(model_ionic, c0, T=T)
    pH_davies = -np.log10(c_eq_davies["H+"] / C_REF)

    # Analytical prediction: pKa_app = pKa - 2*A*f(I)
    I_L = I_bg / 1000.0
    sqrt_I = np.sqrt(I_L)
    f_I = sqrt_I / (1.0 + sqrt_I) - 0.3 * I_L
    pKa_app_analytical = pKa_val - 2.0 * A_davies * f_I

    # Ideal case recovers thermodynamic pKa
    assert abs(pH_ideal - pKa_val) < 1e-6, (
        f"Ideal case: pH at half-neutralisation {pH_ideal:.6f} != pKa {pKa_val}"
    )
    # Davies case: apparent pKa is lower (acid appears stronger)
    assert pH_davies < pKa_val, (
        f"Expected pH_davies ({pH_davies:.4f}) < pKa ({pKa_val}) at I > 0"
    )
    # Quantitative agreement with analytical Davies shift (within 0.001 pK units)
    assert abs(pH_davies - pKa_app_analytical) < 0.001, (
        f"Davies pH at half-neutralisation {pH_davies:.4f} vs analytical "
        f"pKa_app {pKa_app_analytical:.4f}, diff {abs(pH_davies-pKa_app_analytical):.4f}"
    )


# ---------------------------------------------------------------------------
# 3. Phosphate polyprotic: Bjerrum fraction validation at pH 7.2
# ---------------------------------------------------------------------------


def test_phosphate_bjerrum_fractions_ph72():
    """
    Three-step phosphate dissociation: species concentrations at pH 7.2
    match the Bjerrum fraction formula.

    Bjerrum fractions for H3PO4/H2PO4-/HPO4-2/PO4-3:
        D = h^3 + Ka1*h^2 + Ka1*Ka2*h + Ka1*Ka2*Ka3
        alpha0 = h^3 / D
        alpha1 = Ka1*h^2 / D
        alpha2 = Ka1*Ka2*h / D
        alpha3 = Ka1*Ka2*Ka3 / D
    """
    pKa1, pKa2, pKa3 = 2.148, 7.198, 12.350
    Ka1 = 10.0 ** (-pKa1)
    Ka2 = 10.0 ** (-pKa2)
    Ka3 = 10.0 ** (-pKa3)
    c_tot_phos = 100.0   # mol/m3
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
    water = Component("water", [Species("H2O", charge=0, is_solvent=True)])

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

    # Bjerrum fractions (dimensionless activities in mol/L, but ratios cancel)
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
    c_eq = solve_equilibrium(model_phos, c0, T=T)

    tol = 1e-4  # mol/m3
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
# 4. Arrhenius + van't Hoff consistency: kf(T)/kr(T) == K(T)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("T", [250.0, 280.0, 298.15, 320.0, 350.0, 400.0])
def test_arrhenius_vanthoff_kf_kr_ratio(T):
    """kf(T) / kr(T) == K(T) to 1e-12 relative tolerance across temperature."""
    dH = -30e3       # J/mol
    dS = -80.0       # J/(mol K)
    A_arr = 1e8      # pre-exponential (same units as kf)
    Ea = 50e3        # J/mol activation energy

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

    # kr is derived as kf/K, so ratio must equal K exactly within FP precision
    ratio = kf_val / kr_val
    rel_err = abs(ratio - K_val) / K_val
    assert rel_err < 1e-12, (
        f"At T={T} K: kf/kr={ratio:.6e}, K={K_val:.6e}, rel_err={rel_err:.2e}"
    )


# ---------------------------------------------------------------------------
# 5. Edge cases
# ---------------------------------------------------------------------------


def test_extreme_ph_acid_ph1():
    """Acetic acid at pH 1: essentially fully protonated; equilibrium solver converges."""
    pKa_val = 4.76
    c_tot = 100.0
    T = 298.15

    acetate = Component("acetate", [Species("AcOH", charge=0), Species("AcO-", charge=-1)])
    proton = Component("proton", [Species("H+", charge=+1)])
    hydroxide = Component("hydroxide", [Species("OH-", charge=-1)])
    water = Component("water", [Species("H2O", charge=0, is_solvent=True)])

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
    c_eq = solve_equilibrium(model, c0, T=T)

    # At pH 1, [AcO-]/c_tot should be tiny (fraction ~ 10^(1-4.76) ~ 1.7e-4)
    frac_deprotonated = c_eq["AcO-"] / c_tot
    assert frac_deprotonated < 1e-3, (
        f"Deprotonated fraction {frac_deprotonated:.2e} unexpectedly large at pH 1"
    )
    # Check against Henderson-Hasselbalch
    assert abs(c_eq["AcOH"] - AcOH_hh) < 1e-4


def test_extreme_ph_acid_ph13():
    """Acetic acid at pH 13: essentially fully deprotonated; equilibrium solver converges."""
    pKa_val = 4.76
    c_tot = 100.0
    T = 298.15

    acetate = Component("acetate", [Species("AcOH", charge=0), Species("AcO-", charge=-1)])
    proton = Component("proton", [Species("H+", charge=+1)])
    hydroxide = Component("hydroxide", [Species("OH-", charge=-1)])
    water = Component("water", [Species("H2O", charge=0, is_solvent=True)])

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
    c_eq = solve_equilibrium(model, c0, T=T)

    # At pH 13, essentially all acetate deprotonated
    frac_protonated = c_eq["AcOH"] / c_tot
    assert frac_protonated < 1e-7, (
        f"Protonated fraction {frac_protonated:.2e} unexpectedly large at pH 13"
    )


def test_near_zero_ionic_strength_gives_ideal_activity():
    """Davies activity coefficients converge monotonically toward 1 as I -> 0.

    Two properties are verified:
    1. The neutral species (z=0) has gamma = 1 exactly at all ionic strengths.
    2. The maximum deviation from 1 across all species decreases monotonically
       as I decreases from 1000 to 0.001 mol/m3, confirming the I -> 0 limit
       approaches ideal (gamma -> 1).
    """
    dav = ActivityCoefficientDavies()
    charges = np.array([0.0, -1.0, 1.0, -2.0])
    c_dummy = np.ones(4) * 100.0

    # Neutral species has gamma = 1 exactly (z=0 => log_gamma = 0)
    for I in [0.001, 1.0, 100.0, 1000.0]:
        state = PhysicalState(c=c_dummy, T=298.15, I=I)
        gamma = dav.activity(state, charges)
        assert abs(gamma[0] - 1.0) < 1e-12, (
            f"At I={I}, neutral species gamma = {gamma[0]:.12f} != 1"
        )

    # Deviations from 1 must decrease monotonically as I decreases in the
    # low-to-moderate range where the Davies equation is designed to apply
    # (I < 500 mol/m3 = 0.5 mol/L).  Above that the -0.3*I_L term causes
    # f(I) to turn over, so non-monotone behaviour at very high I is expected.
    I_vals = [100.0, 10.0, 1.0, 0.1, 0.01, 0.001]
    deviations = []
    for I in I_vals:
        state = PhysicalState(c=c_dummy, T=298.15, I=I)
        gamma = dav.activity(state, charges)
        deviations.append(float(np.max(np.abs(gamma - 1.0))))

    for i in range(len(deviations) - 1):
        assert deviations[i] >= deviations[i + 1], (
            f"Activity deviation not monotone in valid range: "
            f"I={I_vals[i]} dev={deviations[i]:.6f}, "
            f"I={I_vals[i+1]} dev={deviations[i+1]:.6f}"
        )

    # Confirm that the smallest I in the sweep gives a much smaller deviation
    # than the largest I (qualitative convergence toward ideal).
    assert deviations[-1] < deviations[0] / 10.0, (
        f"Deviation at smallest I ({deviations[-1]:.4f}) should be much less than "
        f"deviation at largest I ({deviations[0]:.4f})"
    )


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
_WATER = Species("H2O", charge=0, is_solvent=True,
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
    rho_cp = model.volumetric_heat_capacity(_WATER_X)
    np.testing.assert_allclose(rho_cp, _RHO_CP_WATER, rtol=1e-10)


def test_volumetric_heat_capacity_defaults_single_solvent():
    """volumetric_heat_capacity defaults x_k=1 when only one solvent present."""
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
    assert model.volumetric_heat_capacity() == model.volumetric_heat_capacity(_WATER_X)


def test_volumetric_heat_capacity_mixture():
    """volumetric_heat_capacity weighted sum for two-solvent mixture."""
    water = Species("H2O",  is_solvent=True, molar_mass=0.018, density=997.0, heat_capacity=75.3)
    mecn  = Species("MeCN", is_solvent=True, molar_mass=0.041, density=786.0, heat_capacity=91.4)
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
    x = {"H2O": 0.7, "MeCN": 0.3}
    expected = 0.7 * (997.0 / 0.018) * 75.3 + 0.3 * (786.0 / 0.041) * 91.4
    np.testing.assert_allclose(model.volumetric_heat_capacity(x), expected, rtol=1e-10)


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
    # 30% MeCN ramp over 10 s — ρCp decreases as water is replaced
    water = Species("H2O",  is_solvent=True, molar_mass=0.018, density=997.0, heat_capacity=75.3)
    mecn  = Species("MeCN", is_solvent=True, molar_mass=0.041, density=786.0, heat_capacity=91.4)
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


def test_reaction_enthalpy_fd_fallback_custom():
    """EquilibriumConstantCustom.reaction_enthalpy() falls back to FD via van't Hoff."""
    dH_true = -20e3
    dS = -50.0
    K_fn = lambda T: float(np.exp(-dH_true / (R_GAS * T) + dS / R_GAS))  # noqa: E731
    eq = EquilibriumConstantCustom(K_fn)
    # analytic: ΔH = R T² · dlnK/dT = dH_true
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
    assert abs(dH_fd - dH_true) / abs(dH_true) < 0.02  # linear interpolation limits accuracy


def test_reaction_enthalpy_fixed_K_is_zero():
    """EquilibriumConstant (fixed K) implies dlnK/dT = 0, so reaction_enthalpy = 0."""
    eq = EquilibriumConstant(K_eq=4.0)
    assert eq.reaction_enthalpy(298.15) == 0.0


def test_equilibrium_constant_polynomial_K_and_derivative():
    """EquilibriumConstantPolynomial: K and dlnK/dT match analytic values."""
    # ln K(T) = a0 + a1*T + a2*T^2
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
    import warnings
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
    # coeffs chosen so both dlnK/dT and d²lnK/dT² are non-trivial
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
        "H2O", is_solvent=True, molar_mass=0.018, density=1000.0, heat_capacity=75.3
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
    rho_cp = model.volumetric_heat_capacity({"H2O": 1.0})
    c = np.array([800.0, 200.0])
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
        "H2O", is_solvent=True, molar_mass=0.018, density=1000.0, heat_capacity=75.3
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
    rho_cp = model.volumetric_heat_capacity({"H2O": 1.0})
    c = np.array([600.0, 400.0])
    T = 305.0
    eps = 1e-5

    rhs0 = _eb_rhs_T(model, c, T, rho_cp)
    jac_T_fd = (_eb_rhs_T(model, c, T + eps, rho_cp) - rhs0) / eps
    _, jac_T = _eb_jac_analytic(model, c, T, rho_cp)
    np.testing.assert_allclose(jac_T, jac_T_fd, rtol=1e-5)
