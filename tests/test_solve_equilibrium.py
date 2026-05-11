"""
Tests for reactions.solver.solve_equilibrium():
acid-base speciation, polyprotic buffers, edge cases, conservation enforcement.
"""

import math

import numpy as np
import pytest

from reactions.api import (
    Component,
    EquilibriumConstant,
    IonicStrengthIdeal,
    ReactionModel,
    Solution,
    Species,
    ThermodynamicReaction,
    acetic_acid,
    acetic_acid_equilibria,
    autoionisation,
    H_plus,
    OH_minus,
    pKa,
    water,
)
from reactions.solver import solve_equilibrium

C_REF: float = 1000.0  # mol/m³ — standard-state concentration


def test_acid_base_henderson_hasselbalch():
    """Acetic acid speciation matches Henderson-Hasselbalch (ideal activity)."""
    pKa_val, c_tot, T = 4.76, 100.0, 298.15

    acetate = Component("acetate", [Species("AcOH", charge=0), Species("AcO-", charge=-1)])
    proton = Component("proton", [Species("H+", charge=+1)])
    hydroxide = Component("hydroxide", [Species("OH-", charge=-1)])
    water_local = Component("water", [Species("H2O", charge=0)])

    model = ReactionModel(
        components=[acetate, proton, hydroxide, water_local],
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
    water_local = Component("water", [Species("H2O", charge=0)])

    model_phos = ReactionModel(
        components=[phosphate, proton, hydroxide, water_local],
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


def test_extreme_ph_acid_ph1():
    """Acetic acid at pH 1: essentially fully protonated; equilibrium solver converges."""
    pKa_val = 4.76
    c_tot = 100.0
    T = 298.15

    acetate = Component("acetate", [Species("AcOH", charge=0), Species("AcO-", charge=-1)])
    proton = Component("proton", [Species("H+", charge=+1)])
    hydroxide = Component("hydroxide", [Species("OH-", charge=-1)])
    water_local = Component("water", [Species("H2O", charge=0)])

    model = ReactionModel(
        components=[acetate, proton, hydroxide, water_local],
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
    water_local = Component("water", [Species("H2O", charge=0)])

    model = ReactionModel(
        components=[acetate, proton, hydroxide, water_local],
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
# Conservation enforcement
# ---------------------------------------------------------------------------


def test_equilibrium_conservation_pure_weak_acid():
    """
    100 mM pure HAc in water: pH ~2.88, total acetate conserved.
    Prior to the conservation fix, pH was determined solely by the initial
    c_OH- guess and total acetate was not conserved.
    """
    model = ReactionModel(
        components=[acetic_acid, H_plus, OH_minus, water],
        reactions=[*acetic_acid_equilibria(), *autoionisation()],
        ionic_strength=IonicStrengthIdeal(),
    )
    sol = Solution(water, solutes={"HAc": 100.0, "Ac-": 1e-6, "H+": 1e-4, "OH-": 1e-7})
    c_eq = solve_equilibrium(model, sol.c0, prescribed=sol.prescribed)

    pH = -math.log10(c_eq["H+"] / 1000)
    total_acetate = c_eq["HAc"] + c_eq["Ac-"]

    assert pytest.approx(pH, abs=0.01) == 2.88
    assert pytest.approx(total_acetate, rel=1e-6) == 100.0


def test_equilibrium_conservation_independent_of_initial_oh():
    """
    solve_equilibrium must give the same pH regardless of c_OH- initial guess.
    """
    model = ReactionModel(
        components=[acetic_acid, H_plus, OH_minus, water],
        reactions=[*acetic_acid_equilibria(), *autoionisation()],
        ionic_strength=IonicStrengthIdeal(),
    )
    c_base = {"HAc": 100.0, "Ac-": 1e-6, "H+": 1e-4, "H2O": water.c_ref}
    pHs = []
    for c_oh in [1e-5, 1e-7, 1e-9]:
        c0 = {**c_base, "OH-": c_oh}
        c_eq = solve_equilibrium(model, c0, prescribed={"H2O": water.c_ref})
        pHs.append(-math.log10(c_eq["H+"] / 1000))

    assert max(pHs) - min(pHs) < 1e-4
