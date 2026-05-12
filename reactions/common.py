"""
reactions.common
================
Pre-defined Component instances and equilibrium reaction factories for
commonly used buffer species and solvents.

All pKa values are at 298.15 K.  Van't Hoff ΔH° values are included where
the temperature sensitivity is practically significant (water, Tris); other
buffers are treated as temperature-independent at this level.

Usage pattern::

    from reactions.common import (
        water, H_plus, OH_minus,
        acetic_acid, autoionisation, acetic_acid_equilibria,
    )
    from reactions.activity import ActivityCoefficientDavies
    from reactions.formulation import Solution
    from reactions.ionic import IonicStrengthIdeal
    from reactions.model import ReactionModel

    davies = ActivityCoefficientDavies()
    model = ReactionModel(
        components=[acetic_acid, H_plus, OH_minus, water],
        reactions=[
            *acetic_acid_equilibria(activity_coefficient=davies),
            *autoionisation(activity_coefficient=davies),
        ],
        ionic_strength=IonicStrengthIdeal(),
    )
    sol = Solution(water, solutes={"HAc": 100.0, "H+": 1e-4, "OH-": 1e-7})
    c_eq = solve_equilibrium(model, sol.c0, prescribed=sol.prescribed)
"""

from __future__ import annotations

from .equilibrium import pKa
from .reaction import ThermodynamicReaction
from .species import Component, Species

__all__ = [
    "water",
    "H_plus",
    "OH_minus",
    "acetic_acid",
    "phosphate",
    "citric_acid",
    "tris",
    "hepes",
    "mops",
    "autoionisation",
    "acetic_acid_equilibria",
    "phosphate_equilibria",
    "citric_acid_equilibria",
    "tris_equilibria",
    "hepes_equilibria",
    "mops_equilibria",
]

# ---------------------------------------------------------------------------
# Solvent
# ---------------------------------------------------------------------------

water = Component(
    "water",
    [Species("H2O", charge=0, molar_mass=0.018015, density=1000.0)],
)

# ---------------------------------------------------------------------------
# Proton carriers  (shared across buffer systems)
# ---------------------------------------------------------------------------

H_plus = Component("H+", [Species("H+", charge=1, molar_mass=1.008e-3)])
OH_minus = Component("OH-", [Species("OH-", charge=-1, molar_mass=17.007e-3)])

# ---------------------------------------------------------------------------
# Buffer components
# Molar masses [kg/mol] from IUPAC atomic weights (H=1.008, C=12.011,
# N=14.007, O=15.999, P=30.974, S=32.06).  Each dissociation step loses
# one H, so molar_mass decreases by _M_H = 1.008e-3 kg/mol per step.
# ---------------------------------------------------------------------------

_M_H = 1.008e-3  # kg/mol

#: Acetic acid / acetate  (pKa = 4.756)
#: C2H4O2 = 60.052 g/mol
acetic_acid = Component(
    "acetic_acid",
    [
        Species("HAc", charge=0,  molar_mass=60.052e-3),
        Species("Ac-", charge=-1, molar_mass=60.052e-3 - _M_H),
    ],
)

#: Phosphoric acid / phosphate  (pKa = 2.148, 7.198, 12.35)
#: H3PO4 = 97.994 g/mol
phosphate = Component(
    "phosphate",
    [
        Species("H3PO4",  charge=0,  molar_mass=97.994e-3),
        Species("H2PO4-", charge=-1, molar_mass=97.994e-3 - _M_H),
        Species("HPO4-2", charge=-2, molar_mass=97.994e-3 - 2 * _M_H),
        Species("PO4-3",  charge=-3, molar_mass=97.994e-3 - 3 * _M_H),
    ],
)

#: Citric acid / citrate  (pKa = 3.128, 4.761, 6.396)
#: C6H8O7 = 192.123 g/mol
citric_acid = Component(
    "citric_acid",
    [
        Species("H3Cit",  charge=0,  molar_mass=192.123e-3),
        Species("H2Cit-", charge=-1, molar_mass=192.123e-3 - _M_H),
        Species("HCit-2", charge=-2, molar_mass=192.123e-3 - 2 * _M_H),
        Species("Cit-3",  charge=-3, molar_mass=192.123e-3 - 3 * _M_H),
    ],
)

#: Tris  (pKa = 8.072, ΔH° = -47.45 kJ/mol; strongly temperature-dependent)
#: C4H11NO3 = 121.136 g/mol (free base); TrisH+ = 122.144 g/mol
tris = Component(
    "tris",
    [
        Species("TrisH+", charge=1,  molar_mass=122.144e-3),
        Species("Tris",   charge=0,  molar_mass=122.144e-3 - _M_H),
    ],
)

#: HEPES  (pKa = 7.55; neutral zwitterion <=> anion)
#: C8H18N2O4S = 238.302 g/mol (protonated form)
hepes = Component(
    "hepes",
    [
        Species("HEPESH", charge=0,  molar_mass=238.302e-3),
        Species("HEPES-", charge=-1, molar_mass=238.302e-3 - _M_H),
    ],
)

#: MOPS  (pKa = 7.20; neutral zwitterion <=> anion)
#: C7H15NO4S = 209.260 g/mol (protonated form)
mops = Component(
    "mops",
    [
        Species("MOPSH", charge=0,  molar_mass=209.260e-3),
        Species("MOPS-", charge=-1, molar_mass=209.260e-3 - _M_H),
    ],
)

# ---------------------------------------------------------------------------
# Equilibrium reaction factories
# ---------------------------------------------------------------------------
# Each factory returns list[ThermodynamicReaction] so the calling convention
# is always *factory(...) regardless of how many reactions a system has.


def autoionisation(activity_coefficient=None) -> list[ThermodynamicReaction]:
    """
    Water autoionisation: ``H2O <-> H+ + OH-``  (pKw = 14.00, ΔH° = 55.8 kJ/mol).

    Requires components: ``water``, ``H_plus``, ``OH_minus``.
    """
    return [
        ThermodynamicReaction(
            "H2O <-> H+ + OH-",
            mode="equil",
            equilibrium_constant=pKa(14.00, dH=55800.0),
            activity_coefficient=activity_coefficient,
        ),
    ]


def acetic_acid_equilibria(activity_coefficient=None) -> list[ThermodynamicReaction]:
    """
    ``HAc <-> Ac- + H+``  (pKa = 4.756; temperature-independent).

    Requires components: ``acetic_acid``, ``H_plus``.
    """
    return [
        ThermodynamicReaction(
            "HAc <-> Ac- + H+",
            mode="equil",
            equilibrium_constant=pKa(4.756),
            activity_coefficient=activity_coefficient,
        ),
    ]


def phosphate_equilibria(activity_coefficient=None) -> list[ThermodynamicReaction]:
    """
    Three dissociation steps of phosphoric acid (pKa = 2.148, 7.198, 12.35).

    Requires components: ``phosphate``, ``H_plus``.
    """
    return [
        ThermodynamicReaction(
            "H3PO4 <-> H2PO4- + H+",
            mode="equil",
            equilibrium_constant=pKa(2.148),
            activity_coefficient=activity_coefficient,
        ),
        ThermodynamicReaction(
            "H2PO4- <-> HPO4-2 + H+",
            mode="equil",
            equilibrium_constant=pKa(7.198),
            activity_coefficient=activity_coefficient,
        ),
        ThermodynamicReaction(
            "HPO4-2 <-> PO4-3 + H+",
            mode="equil",
            equilibrium_constant=pKa(12.35),
            activity_coefficient=activity_coefficient,
        ),
    ]


def citric_acid_equilibria(activity_coefficient=None) -> list[ThermodynamicReaction]:
    """
    Three dissociation steps of citric acid (pKa = 3.128, 4.761, 6.396).

    Requires components: ``citric_acid``, ``H_plus``.
    """
    return [
        ThermodynamicReaction(
            "H3Cit <-> H2Cit- + H+",
            mode="equil",
            equilibrium_constant=pKa(3.128),
            activity_coefficient=activity_coefficient,
        ),
        ThermodynamicReaction(
            "H2Cit- <-> HCit-2 + H+",
            mode="equil",
            equilibrium_constant=pKa(4.761),
            activity_coefficient=activity_coefficient,
        ),
        ThermodynamicReaction(
            "HCit-2 <-> Cit-3 + H+",
            mode="equil",
            equilibrium_constant=pKa(6.396),
            activity_coefficient=activity_coefficient,
        ),
    ]


def tris_equilibria(activity_coefficient=None) -> list[ThermodynamicReaction]:
    """
    ``TrisH+ <-> Tris + H+``  (pKa = 8.072, ΔH° = −47.45 kJ/mol).

    Temperature sensitivity: pKa drops ~0.03 per °C — use van't Hoff for
    any simulation where temperature deviates from 298.15 K.

    Requires components: ``tris``, ``H_plus``.
    """
    return [
        ThermodynamicReaction(
            "TrisH+ <-> Tris + H+",
            mode="equil",
            equilibrium_constant=pKa(8.072, dH=-47450.0),
            activity_coefficient=activity_coefficient,
        ),
    ]


def hepes_equilibria(activity_coefficient=None) -> list[ThermodynamicReaction]:
    """
    ``HEPESH <-> HEPES- + H+``  (pKa = 7.55; temperature-independent).

    Requires components: ``hepes``, ``H_plus``.
    """
    return [
        ThermodynamicReaction(
            "HEPESH <-> HEPES- + H+",
            mode="equil",
            equilibrium_constant=pKa(7.55),
            activity_coefficient=activity_coefficient,
        ),
    ]


def mops_equilibria(activity_coefficient=None) -> list[ThermodynamicReaction]:
    """
    ``MOPSH <-> MOPS- + H+``  (pKa = 7.20; temperature-independent).

    Requires components: ``mops``, ``H_plus``.
    """
    return [
        ThermodynamicReaction(
            "MOPSH <-> MOPS- + H+",
            mode="equil",
            equilibrium_constant=pKa(7.20),
            activity_coefficient=activity_coefficient,
        ),
    ]
