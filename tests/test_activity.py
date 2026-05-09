"""
Unit tests for activity coefficient and ionic strength models.
"""

import warnings

import numpy as np
import pytest

from reactions.api import (
    ActivityCoefficientDavies,
    ActivityCoefficientDebyeHuckel,
    Component,
    IonicStrengthFixed,
    PhysicalState,
    Species,
    ThermodynamicReaction,
    EquilibriumConstant,
    ReactionModel,
    _water_epsilon_r,
    pKa,
)
from reactions.solver import solve_equilibrium

C_REF: float = 1000.0


# ---------------------------------------------------------------------------
# Davies / Debye-Hückel unit tests
# ---------------------------------------------------------------------------


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

    for I in [0.001, 1.0, 100.0, 1000.0]:
        state = PhysicalState(c=c_dummy, T=298.15, I=I)
        gamma = dav.activity(state, charges)
        assert abs(gamma[0] - 1.0) < 1e-12, (
            f"At I={I}, neutral species gamma = {gamma[0]:.12f} != 1"
        )

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

    assert deviations[-1] < deviations[0] / 10.0


def test_davies_ac_ignores_solvent_concentration():
    """Davies AC for ionic species is independent of solvent concentration (z=0 → gamma=1)."""
    dav = ActivityCoefficientDavies()
    charges = np.array([1.0, -1.0, 0.0])
    I = 50.0

    state_low = PhysicalState(c=np.array([1e-4, 1e-4, 100.0]), T=298.15, I=I)
    state_high = PhysicalState(c=np.array([1e-4, 1e-4, 1000.0]), T=298.15, I=I)

    gamma_low = dav.activity(state_low, charges)
    gamma_high = dav.activity(state_high, charges)

    np.testing.assert_array_equal(gamma_low, gamma_high)


def test_davies_ph_shift_direction():
    """
    Davies corrections shift the apparent pKa in the correct direction.

    At half-neutralisation ([AcO-]=[AcOH]), pH = pKa_app.
    With Davies corrections: pKa_app = pKa - 2*A*f(I) < pKa.
    """
    pKa_val = 4.756
    T = 298.15
    A_davies = 0.509

    acetate = Component("acetate", [Species("AcOH", charge=0), Species("AcO-", charge=-1)])
    proton = Component("proton", [Species("H+", charge=+1)])
    hydroxide = Component("hydroxide", [Species("OH-", charge=-1)])
    water = Component("water", [Species("H2O", charge=0)])

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
    H_guess = Ka * C_REF
    OH_guess = 1e-14 * C_REF**2 / H_guess
    c0 = {
        "AcOH": 50.0,
        "AcO-": 50.0,
        "H+": H_guess,
        "OH-": OH_guess,
    }

    model_ideal = make_model(I_fixed=0.0)
    c_eq_ideal = solve_equilibrium(model_ideal, c0, T=T, prescribed={"H2O": C_REF})
    pH_ideal = -np.log10(c_eq_ideal["H+"] / C_REF)

    I_bg = 100.0
    model_ionic = make_model(I_fixed=I_bg)
    c_eq_davies = solve_equilibrium(model_ionic, c0, T=T, prescribed={"H2O": C_REF})
    pH_davies = -np.log10(c_eq_davies["H+"] / C_REF)

    I_L = I_bg / 1000.0
    sqrt_I = np.sqrt(I_L)
    f_I = sqrt_I / (1.0 + sqrt_I) - 0.3 * I_L
    pKa_app_analytical = pKa_val - 2.0 * A_davies * f_I

    assert abs(pH_ideal - pKa_val) < 1e-6
    assert pH_davies < pKa_val
    assert abs(pH_davies - pKa_app_analytical) < 0.001


# ---------------------------------------------------------------------------
# E/P4 — T-dependent A in DH/Davies (epsilon_r parameter)
# ---------------------------------------------------------------------------


def test_davies_A_increases_above_ambient():
    """With T-dependent epsilon_r, Davies A increases above 298 K (εr·T decreases)."""
    charges = np.array([1.0, -1.0])
    I = 100.0

    state_ambient = PhysicalState(c=np.array([1e-4, 1e-4]), T=298.15, I=I)
    state_warm = PhysicalState(c=np.array([1e-4, 1e-4]), T=320.0, I=I)

    dav = ActivityCoefficientDavies(epsilon_r=_water_epsilon_r)
    gamma_ambient = dav.activity(state_ambient, charges)
    gamma_warm = dav.activity(state_warm, charges)

    assert gamma_warm[0] < gamma_ambient[0], (
        f"Expected gamma(320K)={gamma_warm[0]:.4f} < gamma(298K)={gamma_ambient[0]:.4f}"
    )


def test_dh_A_increases_above_ambient():
    """With T-dependent epsilon_r, DH A increases above 298 K."""
    charges = np.array([1.0, -1.0])
    I = 50.0

    state_ambient = PhysicalState(c=np.array([1e-4, 1e-4]), T=298.15, I=I)
    state_warm = PhysicalState(c=np.array([1e-4, 1e-4]), T=320.0, I=I)

    dh = ActivityCoefficientDebyeHuckel(epsilon_r=_water_epsilon_r)
    gamma_ambient = dh.activity(state_ambient, charges)
    gamma_warm = dh.activity(state_warm, charges)

    assert gamma_warm[0] < gamma_ambient[0]


def test_davies_scalar_epsilon_r_overrides_default():
    """A scalar epsilon_r overrides the stored 25 °C A: lower εr → higher A → lower γ."""
    charges = np.array([1.0, -1.0])
    state = PhysicalState(c=np.array([1e-4, 1e-4]), T=298.15, I=100.0)

    dav_default = ActivityCoefficientDavies()
    dav_custom = ActivityCoefficientDavies(epsilon_r=40.0)

    gamma_default = dav_default.activity(state, charges)
    gamma_custom = dav_custom.activity(state, charges)

    assert gamma_custom[0] < gamma_default[0]


def test_davies_warning_fires_at_non_ambient_T():
    """UserWarning emitted when T deviates > 5 K from 298.15 and epsilon_r is None."""
    charges = np.array([1.0, -1.0])
    state = PhysicalState(c=np.array([1e-4, 1e-4]), T=320.0, I=100.0)
    dav = ActivityCoefficientDavies()

    with pytest.warns(UserWarning, match="epsilon_r was not provided"):
        dav.activity(state, charges)


def test_dh_warning_fires_at_non_ambient_T():
    """UserWarning emitted when T deviates > 5 K from 298.15 and epsilon_r is None."""
    charges = np.array([1.0, -1.0])
    state = PhysicalState(c=np.array([1e-4, 1e-4]), T=320.0, I=50.0)
    dh = ActivityCoefficientDebyeHuckel()

    with pytest.warns(UserWarning, match="epsilon_r was not provided"):
        dh.activity(state, charges)


def test_davies_no_warning_at_ambient_T():
    """No warning at T = 298.15 K even without epsilon_r."""
    charges = np.array([1.0, -1.0])
    state = PhysicalState(c=np.array([1e-4, 1e-4]), T=298.15, I=100.0)
    dav = ActivityCoefficientDavies()

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        dav.activity(state, charges)


def test_davies_no_warning_with_epsilon_r():
    """No warning at non-ambient T when epsilon_r is provided."""
    charges = np.array([1.0, -1.0])
    state = PhysicalState(c=np.array([1e-4, 1e-4]), T=320.0, I=100.0)
    dav = ActivityCoefficientDavies(epsilon_r=_water_epsilon_r)

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        dav.activity(state, charges)


def test_davies_backward_compat_at_ambient():
    """Default Davies gives the standard A=0.509 result at 298.15 K."""
    charges = np.array([1.0, -1.0, 0.0])
    I = 100.0
    state = PhysicalState(c=np.array([1e-4, 1e-4, 1000.0]), T=298.15, I=I)

    gamma = ActivityCoefficientDavies().activity(state, charges)

    A, I_L = 0.509, 0.1
    sqrt_I = np.sqrt(I_L)
    expected = 10.0 ** (-A * (sqrt_I / (1.0 + sqrt_I) - 0.3 * I_L))

    np.testing.assert_allclose(gamma[0], expected, rtol=1e-10)
    np.testing.assert_allclose(gamma[2], 1.0, rtol=1e-10)
