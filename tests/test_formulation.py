"""
Tests for reactions.formulation.Solution.
"""

import pytest

from reactions.formulation import Solution
from reactions.species import Component, Species


def _water():
    return Component("water", [Species("H2O", charge=0, molar_mass=0.018015, density=1000.0)])


def test_solution_pure_solvent_c0():
    """Pure solvent: c0['H2O'] == water.c_ref."""
    water = _water()
    sol = Solution(water)
    assert sol.c0["H2O"] == pytest.approx(water.c_ref, rel=1e-9)
    assert sol.prescribed["H2O"] == pytest.approx(water.c_ref, rel=1e-9)


def test_solution_solutes_added():
    """Solutes appear in c0 but not in prescribed."""
    water = _water()
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
    sol = Solution({water: 0.9, ethanol: 0.1})
    assert sol.c0["H2O"] == pytest.approx(0.9 * water.c_ref, rel=1e-9)
    assert sol.c0["EtOH"] == pytest.approx(0.1 * ethanol.c_ref, rel=1e-9)
    assert set(sol.prescribed.keys()) == {"H2O", "EtOH"}


def test_solution_fractions_must_sum_to_one():
    water = _water()
    with pytest.raises(ValueError, match="sum to 1"):
        Solution({water: 0.5})


def test_solution_missing_density_raises():
    no_density = Component("X", [Species("X")])
    with pytest.raises(ValueError, match="density"):
        Solution(no_density)


def test_solution_c0_and_prescribed_are_copies():
    """Mutating the returned dict does not affect the Solution."""
    water = _water()
    sol = Solution(water)
    c0 = sol.c0
    c0["H2O"] = 0.0
    assert sol.c0["H2O"] != 0.0
