"""
Tests for parse_stoichiometry and reaction classes.
"""

import pytest

from reactions.reaction import parse_stoichiometry


def test_basic_reversible():
    nu = parse_stoichiometry("A + 2 B <-> C")
    assert nu == {"A": -1.0, "B": -2.0, "C": 1.0}


def test_charged_species():
    nu = parse_stoichiometry("H2O <-> H+ + OH-")
    assert nu == {"H2O": -1.0, "H+": 1.0, "OH-": 1.0}


def test_phosphate():
    nu = parse_stoichiometry("H2PO4- <-> HPO4-2 + H+")
    assert nu == {"H2PO4-": -1.0, "HPO4-2": 1.0, "H+": 1.0}


def test_empty_lhs():
    """<-> H+ + OH- omits the water term (a_{H2O}=1 absorbed into K_w)."""
    nu = parse_stoichiometry("<-> H+ + OH-")
    assert nu == {"H+": 1.0, "OH-": 1.0}


def test_empty_rhs():
    nu = parse_stoichiometry("H2O <->")
    assert nu == {"H2O": -1.0}


def test_empty_rhs_solubility_product():
    """Ca2+ + CO3-2 <-> represents a solubility product: solid phase has unit activity."""
    nu = parse_stoichiometry("Ca2+ + CO3-2 <->")
    assert nu == {"Ca2+": -1.0, "CO3-2": -1.0}


def test_irreversible():
    nu = parse_stoichiometry("A -> B")
    assert nu == {"A": -1.0, "B": 1.0}


def test_fractional_coefficient():
    nu = parse_stoichiometry("2 A <-> B")
    assert nu["A"] == pytest.approx(-2.0)
    assert nu["B"] == pytest.approx(1.0)
