"""
Tests for reactions.solver.simulate():
ODE integration, thermodynamic consistency, energy balance, prescribed species.
"""

import warnings

import numpy as np
import pytest

from reactions.api import (
    ActivityCoefficientCustom,
    Component,
    EquilibriumConstant,
    EquilibriumConstantCustom,
    EquilibriumConstantVantHoff,
    MassActionReaction,
    R_GAS,
    RateConstantArrhenius,
    RateConstantFixed,
    ReactionModel,
    Species,
    ThermodynamicReaction,
)
from reactions.solver import simulate

C_REF: float = 1000.0  # mol/m³ — standard-state concentration

# ---------------------------------------------------------------------------
# Water solvent constants shared by energy balance and solvent-state tests
# ---------------------------------------------------------------------------

_WATER = Species("H2O", charge=0, molar_mass=0.018, density=997.0, heat_capacity=75.3)
_WATER_COMPONENT = Component("water", [_WATER])
_WATER_X = {"H2O": 1.0}
_RHO_CP_WATER = (997.0 / 0.018) * 75.3


# ---------------------------------------------------------------------------
# Basic ODE integration
# ---------------------------------------------------------------------------


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
# Energy balance (coupled T simulation)
# ---------------------------------------------------------------------------


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
# Solvent species in state.c
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
# Prescribed species (callable)
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
