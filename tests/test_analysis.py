"""
Tests for reactions.analysis.buffer_capacity().
All expected values are derived analytically; no fitting.

Van Slyke formula: β = ln(10) · c · Σ_{i<j} (z_i − z_j)² · α_i · α_j
At a pKa, adjacent Bjerrum fractions satisfy α_i = α_j = 0.5, so
β_max = ln(10) · c · (Δz)² · 0.25 for each adjacent pair.
"""

import math

import numpy as np
import pytest

from reactions.analysis import buffer_capacity
from reactions.common import H_plus, OH_minus, water
from reactions.equilibrium import pKa
from reactions.model import ReactionModel
from reactions.reaction import ThermodynamicReaction
from reactions.species import Component, Species


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _acetic_model():
    acetic = Component("acetic", [Species("HAc", charge=0), Species("Ac-", charge=-1)])
    return ReactionModel(
        components=[acetic, H_plus, OH_minus, water],
        reactions=[
            ThermodynamicReaction("HAc <-> Ac- + H+", mode="equil", equilibrium_constant=pKa(4.756)),
            ThermodynamicReaction("H2O <-> H+ + OH-", mode="equil", equilibrium_constant=pKa(14.00)),
        ],
    )


def _phosphate_model():
    phosphate = Component(
        "phosphate",
        [
            Species("H3PO4", charge=0),
            Species("H2PO4-", charge=-1),
            Species("HPO4-2", charge=-2),
            Species("PO4-3", charge=-3),
        ],
    )
    return ReactionModel(
        components=[phosphate, H_plus, OH_minus, water],
        reactions=[
            ThermodynamicReaction("H3PO4 <-> H2PO4- + H+", mode="equil", equilibrium_constant=pKa(2.148)),
            ThermodynamicReaction("H2PO4- <-> HPO4-2 + H+", mode="equil", equilibrium_constant=pKa(7.198)),
            ThermodynamicReaction("HPO4-2 <-> PO4-3 + H+", mode="equil", equilibrium_constant=pKa(12.350)),
            ThermodynamicReaction("H2O <-> H+ + OH-", mode="equil", equilibrium_constant=pKa(14.00)),
        ],
    )


# ---------------------------------------------------------------------------
# Monoprotic: acetic acid
# ---------------------------------------------------------------------------


def test_monoprotic_peak_at_pKa():
    """β peaks at pKa with height ln(10)·c·0.25 for Δz=1."""
    c_tot = 100.0
    pKa_val = 4.756
    model = _acetic_model()

    pH = np.array([pKa_val])
    beta = buffer_capacity(model, {"acetic": c_tot}, pH)

    expected = math.log(10) * c_tot * 0.25
    assert pytest.approx(beta["acetic"][0], rel=1e-4) == expected


def test_monoprotic_returns_water_key():
    model = _acetic_model()
    beta = buffer_capacity(model, {"acetic": 100.0}, np.array([7.0]))
    assert "water" in beta
    assert "acetic" in beta


def test_monoprotic_zero_total_omitted():
    """Component with c_total=0 does not appear in the result."""
    model = _acetic_model()
    beta = buffer_capacity(model, {}, np.array([4.756]))
    assert "acetic" not in beta
    assert "water" in beta


def test_monoprotic_symmetry():
    """β is symmetric around pKa: β(pKa − Δ) ≈ β(pKa + Δ)."""
    model = _acetic_model()
    pKa_val = 4.756
    delta = 0.5
    pH = np.array([pKa_val - delta, pKa_val + delta])
    beta = buffer_capacity(model, {"acetic": 100.0}, pH)
    assert pytest.approx(beta["acetic"][0], rel=1e-6) == beta["acetic"][1]


# ---------------------------------------------------------------------------
# Water / autoionisation
# ---------------------------------------------------------------------------


def test_water_beta_at_neutral_pH():
    """β_water at pH 7 = ln(10)·([H+]+[OH-]) ≈ 4.6e-4 mol/(m³·pH)."""
    model = _acetic_model()
    pH = np.array([7.0])
    beta = buffer_capacity(model, {}, pH)
    c_ref = H_plus.species[0].c_ref  # 1000 mol/m³
    c_H = 1e-7 * c_ref
    c_OH = 1e-14 * c_ref**2 / c_H
    expected = math.log(10) * (c_H + c_OH)
    assert pytest.approx(beta["water"][0], rel=1e-6) == expected


def test_water_beta_symmetric_around_pH7():
    """β_water is symmetric in [H+]+[OH-]: β(pH 7+x) = β(pH 7-x)."""
    model = _acetic_model()
    delta = 2.0
    pH = np.array([7.0 - delta, 7.0 + delta])
    beta = buffer_capacity(model, {}, pH)
    assert pytest.approx(beta["water"][0], rel=1e-6) == beta["water"][1]


# ---------------------------------------------------------------------------
# Polyprotic: phosphate (Δz=1 for all adjacent transitions)
# ---------------------------------------------------------------------------


def test_phosphate_three_equal_peaks():
    """
    All three phosphate peaks have equal height: adjacent Δz=1 throughout.
    β_max = ln(10)·c·0.25 at each pKa.
    """
    c_tot = 100.0
    model = _phosphate_model()
    pKas = [2.148, 7.198, 12.350]
    pH = np.array(pKas)

    beta = buffer_capacity(model, {"phosphate": c_tot}, pH)
    expected = math.log(10) * c_tot * 0.25

    for i, pk in enumerate(pKas):
        assert pytest.approx(beta["phosphate"][i], rel=1e-3) == expected, (
            f"Peak at pKa={pk} expected {expected:.4f}, got {beta['phosphate'][i]:.4f}"
        )


def test_phosphate_peaks_all_equal_to_each_other():
    """Peak heights are mutually equal (ratio 1:1:1, not 1:4:9)."""
    c_tot = 100.0
    model = _phosphate_model()
    pH = np.array([2.148, 7.198, 12.350])
    beta = buffer_capacity(model, {"phosphate": c_tot}, pH)["phosphate"]

    assert pytest.approx(beta[1] / beta[0], rel=1e-3) == 1.0
    assert pytest.approx(beta[2] / beta[0], rel=1e-3) == 1.0


def test_phosphate_linearity_in_concentration():
    """β scales linearly with total concentration."""
    model = _phosphate_model()
    pH = np.array([7.198])
    b1 = buffer_capacity(model, {"phosphate": 50.0}, pH)["phosphate"][0]
    b2 = buffer_capacity(model, {"phosphate": 100.0}, pH)["phosphate"][0]
    assert pytest.approx(b2 / b1, rel=1e-6) == 2.0


def test_phosphate_array_input():
    """buffer_capacity accepts array pH and returns array output."""
    model = _phosphate_model()
    pH = np.linspace(4.0, 10.0, 50)
    beta = buffer_capacity(model, {"phosphate": 100.0}, pH)
    assert beta["phosphate"].shape == (50,)
    assert beta["water"].shape == (50,)


# ---------------------------------------------------------------------------
# Autoionisation identification robustness
# ---------------------------------------------------------------------------


def test_kw_not_confused_with_ka1():
    """
    In a phosphate + water model, _find_kw must return the water Kw,
    not Ka1 of H3PO4 (both reactions produce H+ and a charge-(-1) species).
    Incorrect Kw ≈ Ka1·c_ref² ≈ 7e5 mol²/m⁶; correct Kw ≈ 1e-8 mol²/m⁶.
    A wrong Kw gives β_water ≈ 1e8 at pH 7 instead of ~4.6e-4.
    """
    model = _phosphate_model()
    pH = np.array([7.0])
    beta = buffer_capacity(model, {"phosphate": 100.0}, pH)
    assert beta["water"][0] < 1.0, (
        f"β_water={beta['water'][0]:.2e} suggests Kw was misidentified as Ka1"
    )
