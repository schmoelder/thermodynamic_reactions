"""
reactions.common.
================
Pre-defined Component instances and equilibrium reaction factories for
commonly used buffer species and solvents.

All pKa values are at 298.15 K.  Van't Hoff ΔH° values are included where
the temperature sensitivity is practically significant (water, Tris); other
buffers are treated as temperature-independent at this level.

Usage pattern::

    from reactions.activity import ActivityCoefficientDavies
    from reactions.common import (
        water,
        H_plus,
        OH_minus,
        acetic_acid,
        autoionization,
        acetic_acid_equilibria,
    )
    from reactions.formulation import Solution
    from reactions.ionic import IonicStrengthIdeal
    from reactions.model import ReactionModel

    model = ReactionModel(
        components=[acetic_acid, H_plus, OH_minus, water],
        reactions=[
            *acetic_acid_equilibria(),
            *autoionization(),
        ],
        ionic_strength=IonicStrengthIdeal(),
        activity_coefficient=ActivityCoefficientDavies(),
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
    "bis_tris",
    "citric_acid",
    "dap",
    "dea",
    "hepes",
    "imidazole",
    "lactic_acid",
    "mops",
    "phosphate",
    "piperazine",
    "tris",
    "autoionization",
    "acetic_acid_equilibria",
    "bis_tris_equilibria",
    "citric_acid_equilibria",
    "dap_equilibria",
    "dea_equilibria",
    "hepes_equilibria",
    "imidazole_equilibria",
    "lactic_acid_equilibria",
    "mops_equilibria",
    "phosphate_equilibria",
    "piperazine_equilibria",
    "tris_equilibria",
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
# Buffer components  (alphabetical)
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
        Species("HAc", charge=0, molar_mass=60.052e-3),
        Species("Ac-", charge=-1, molar_mass=60.052e-3 - _M_H),
    ],
)

#: Bis-Tris  (pKa = 6.484; Δz² = 0, ionic-strength-insensitive)
#: C8H19NO5 = 209.242 g/mol (free base)
bis_tris = Component(
    "bis_tris",
    [
        Species("BisH+", charge=1, molar_mass=209.242e-3 + _M_H),
        Species("Bis", charge=0, molar_mass=209.242e-3),
    ],
)

#: Citric acid / citrate  (pKa = 3.128, 4.761, 6.396)
#: C6H8O7 = 192.123 g/mol
citric_acid = Component(
    "citric_acid",
    [
        Species("H3Cit", charge=0, molar_mass=192.123e-3),
        Species("H2Cit-", charge=-1, molar_mass=192.123e-3 - _M_H),
        Species("HCit-2", charge=-2, molar_mass=192.123e-3 - 2 * _M_H),
        Species("Cit-3", charge=-3, molar_mass=192.123e-3 - 3 * _M_H),
    ],
)

#: 1,3-Diaminopropane (DAP)  (pKa = 8.640, 10.470)
#: C3H10N2 = 74.127 g/mol (free base)
#: Step 1 (H2DAP2+ → HDAP+): Δz² = −2 (pKa_app rises with I).
#: Step 2 (HDAP+ → DAP): Δz² = 0.
dap = Component(
    "dap",
    [
        Species("H2DAP2+", charge=2, molar_mass=74.127e-3 + 2 * _M_H),
        Species("HDAP+", charge=1, molar_mass=74.127e-3 + _M_H),
        Species("DAP", charge=0, molar_mass=74.127e-3),
    ],
)

#: Diethanolamine (DEA)  (pKa = 8.883; Δz² = 0, ionic-strength-insensitive)
#: C4H11NO2 = 105.137 g/mol (free base)
dea = Component(
    "dea",
    [
        Species("DEAH+", charge=1, molar_mass=105.137e-3 + _M_H),
        Species("DEA", charge=0, molar_mass=105.137e-3),
    ],
)

#: HEPES  (pKa = 7.55; neutral zwitterion <=> anion)
#: C8H18N2O4S = 238.302 g/mol (protonated form)
hepes = Component(
    "hepes",
    [
        Species("HEPESH", charge=0, molar_mass=238.302e-3),
        Species("HEPES-", charge=-1, molar_mass=238.302e-3 - _M_H),
    ],
)

#: Imidazole  (pKa = 6.993; Δz² = 0, ionic-strength-insensitive)
#: C3H4N2 = 68.079 g/mol (free base)
imidazole = Component(
    "imidazole",
    [
        Species("ImH+", charge=1, molar_mass=68.079e-3 + _M_H),
        Species("Im", charge=0, molar_mass=68.079e-3),
    ],
)

#: Lactic acid / lactate  (pKa = 3.860)
#: C3H6O3 = 90.078 g/mol
lactic_acid = Component(
    "lactic_acid",
    [
        Species("LacH", charge=0, molar_mass=90.078e-3),
        Species("Lac-", charge=-1, molar_mass=90.078e-3 - _M_H),
    ],
)

#: MOPS  (pKa = 7.20; neutral zwitterion <=> anion)
#: C7H15NO4S = 209.260 g/mol (protonated form)
mops = Component(
    "mops",
    [
        Species("MOPSH", charge=0, molar_mass=209.260e-3),
        Species("MOPS-", charge=-1, molar_mass=209.260e-3 - _M_H),
    ],
)

#: Phosphoric acid / phosphate  (pKa = 2.148, 7.198, 12.35)
#: H3PO4 = 97.994 g/mol
phosphate = Component(
    "phosphate",
    [
        Species("H3PO4", charge=0, molar_mass=97.994e-3),
        Species("H2PO4-", charge=-1, molar_mass=97.994e-3 - _M_H),
        Species("HPO4-2", charge=-2, molar_mass=97.994e-3 - 2 * _M_H),
        Species("PO4-3", charge=-3, molar_mass=97.994e-3 - 3 * _M_H),
    ],
)

#: Piperazine  (pKa = 5.333, 9.731)
#: C4H10N2 = 86.138 g/mol (free base)
#: Step 1 (H2Pip2+ → HPip+): Δz² = −2 (pKa_app rises with I).
#: Step 2 (HPip+ → Pip): Δz² = 0.
piperazine = Component(
    "piperazine",
    [
        Species("H2Pip2+", charge=2, molar_mass=86.138e-3 + 2 * _M_H),
        Species("HPip+", charge=1, molar_mass=86.138e-3 + _M_H),
        Species("Pip", charge=0, molar_mass=86.138e-3),
    ],
)

#: Tris  (pKa = 8.072, ΔH° = -47.45 kJ/mol; strongly temperature-dependent)
#: C4H11NO3 = 121.136 g/mol (free base); TrisH+ = 122.144 g/mol
tris = Component(
    "tris",
    [
        Species("TrisH+", charge=1, molar_mass=122.144e-3),
        Species("Tris", charge=0, molar_mass=122.144e-3 - _M_H),
    ],
)

# ---------------------------------------------------------------------------
# Equilibrium reaction factories  (alphabetical; autoionization first)
# ---------------------------------------------------------------------------
# Each factory returns list[ThermodynamicReaction] so the calling convention
# is always *factory(...) regardless of how many reactions a system has.


def autoionization() -> list[ThermodynamicReaction]:
    r"""
    Water autoionization constraint for dilute aqueous systems (pKw = 14.00, ΔH° = 55.8 kJ/mol).

    Uses empty-LHS stoichiometry ``<-> H+ + OH-``, which absorbs $a_{H_2O} = 1$
    into Kw.  Water is therefore **not** required as a component.

    Use this factory when water is the background solvent and is not tracked
    explicitly.  For mixed-solvent systems where $a_{H_2O} \neq 1$, include
    water as a component and write the full reaction ``H2O <-> H+ + OH-``
    with a composition-dependent equilibrium constant instead.

    Requires components: ``H_plus``, ``OH_minus``.
    """
    return [
        ThermodynamicReaction(
            "<-> H+ + OH-",
            mode="equil",
            equilibrium_constant=pKa(14.00, dH=55800.0),
        ),
    ]


def acetic_acid_equilibria() -> list[ThermodynamicReaction]:
    """
    ``HAc <-> Ac- + H+``  (pKa = 4.756; temperature-independent).

    Requires components: ``acetic_acid``, ``H_plus``.
    """
    return [
        ThermodynamicReaction(
            "HAc <-> Ac- + H+",
            mode="equil",
            equilibrium_constant=pKa(4.756),
        ),
    ]


def bis_tris_equilibria() -> list[ThermodynamicReaction]:
    """
    ``BisH+ <-> Bis + H+``  (pKa = 6.484; Δz² = 0).

    Requires components: ``bis_tris``, ``H_plus``.
    """
    return [
        ThermodynamicReaction(
            "BisH+ <-> Bis + H+",
            mode="equil",
            equilibrium_constant=pKa(6.484),
        ),
    ]


def citric_acid_equilibria() -> list[ThermodynamicReaction]:
    """
    Three dissociation steps of citric acid (pKa = 3.128, 4.761, 6.396).

    Requires components: ``citric_acid``, ``H_plus``.
    """
    return [
        ThermodynamicReaction(
            "H3Cit <-> H2Cit- + H+",
            mode="equil",
            equilibrium_constant=pKa(3.128),
        ),
        ThermodynamicReaction(
            "H2Cit- <-> HCit-2 + H+",
            mode="equil",
            equilibrium_constant=pKa(4.761),
        ),
        ThermodynamicReaction(
            "HCit-2 <-> Cit-3 + H+",
            mode="equil",
            equilibrium_constant=pKa(6.396),
        ),
    ]


def dap_equilibria() -> list[ThermodynamicReaction]:
    """
    Two-step DAP (1,3-diaminopropane) dissociation (pKa = 8.640, 10.470).

    Step 1 (Δz² = −2): pKa_app rises slightly with ionic strength.
    Step 2 (Δz² = 0): ionic-strength-insensitive.

    Requires components: ``dap``, ``H_plus``.
    """
    return [
        ThermodynamicReaction(
            "H2DAP2+ <-> HDAP+ + H+",
            mode="equil",
            equilibrium_constant=pKa(8.640),
        ),
        ThermodynamicReaction(
            "HDAP+ <-> DAP + H+",
            mode="equil",
            equilibrium_constant=pKa(10.470),
        ),
    ]


def dea_equilibria() -> list[ThermodynamicReaction]:
    """
    ``DEAH+ <-> DEA + H+``  (pKa = 8.883; Δz² = 0).

    Requires components: ``dea``, ``H_plus``.
    """
    return [
        ThermodynamicReaction(
            "DEAH+ <-> DEA + H+",
            mode="equil",
            equilibrium_constant=pKa(8.883),
        ),
    ]


def hepes_equilibria() -> list[ThermodynamicReaction]:
    """
    ``HEPESH <-> HEPES- + H+``  (pKa = 7.55; temperature-independent).

    Requires components: ``hepes``, ``H_plus``.
    """
    return [
        ThermodynamicReaction(
            "HEPESH <-> HEPES- + H+",
            mode="equil",
            equilibrium_constant=pKa(7.55),
        ),
    ]


def imidazole_equilibria() -> list[ThermodynamicReaction]:
    """
    ``ImH+ <-> Im + H+``  (pKa = 6.993; Δz² = 0).

    Requires components: ``imidazole``, ``H_plus``.
    """
    return [
        ThermodynamicReaction(
            "ImH+ <-> Im + H+",
            mode="equil",
            equilibrium_constant=pKa(6.993),
        ),
    ]


def lactic_acid_equilibria() -> list[ThermodynamicReaction]:
    """
    ``LacH <-> Lac- + H+``  (pKa = 3.860; temperature-independent).

    Requires components: ``lactic_acid``, ``H_plus``.
    """
    return [
        ThermodynamicReaction(
            "LacH <-> Lac- + H+",
            mode="equil",
            equilibrium_constant=pKa(3.860),
        ),
    ]


def mops_equilibria() -> list[ThermodynamicReaction]:
    """
    ``MOPSH <-> MOPS- + H+``  (pKa = 7.20; temperature-independent).

    Requires components: ``mops``, ``H_plus``.
    """
    return [
        ThermodynamicReaction(
            "MOPSH <-> MOPS- + H+",
            mode="equil",
            equilibrium_constant=pKa(7.20),
        ),
    ]


def phosphate_equilibria() -> list[ThermodynamicReaction]:
    """
    Three dissociation steps of phosphoric acid (pKa = 2.148, 7.198, 12.35).

    Requires components: ``phosphate``, ``H_plus``.
    """
    return [
        ThermodynamicReaction(
            "H3PO4 <-> H2PO4- + H+",
            mode="equil",
            equilibrium_constant=pKa(2.148),
        ),
        ThermodynamicReaction(
            "H2PO4- <-> HPO4-2 + H+",
            mode="equil",
            equilibrium_constant=pKa(7.198),
        ),
        ThermodynamicReaction(
            "HPO4-2 <-> PO4-3 + H+",
            mode="equil",
            equilibrium_constant=pKa(12.35),
        ),
    ]


def piperazine_equilibria() -> list[ThermodynamicReaction]:
    """
    Two-step piperazine dissociation (pKa = 5.333, 9.731).

    Step 1 (Δz² = −2): pKa_app rises slightly with ionic strength.
    Step 2 (Δz² = 0): ionic-strength-insensitive.

    Requires components: ``piperazine``, ``H_plus``.
    """
    return [
        ThermodynamicReaction(
            "H2Pip2+ <-> HPip+ + H+",
            mode="equil",
            equilibrium_constant=pKa(5.333),
        ),
        ThermodynamicReaction(
            "HPip+ <-> Pip + H+",
            mode="equil",
            equilibrium_constant=pKa(9.731),
        ),
    ]


def tris_equilibria() -> list[ThermodynamicReaction]:
    """
    ``TrisH+ <-> Tris + H+``  (pKa = 8.072, ΔH° = +47.45 kJ/mol).

    Temperature sensitivity: pKa drops ~0.028 per °C (endothermic ionisation;
    Goldberg et al. 2002, J. Phys. Chem. Ref. Data 31, 231).
    Note: dH is the *ionisation* (dissociation) enthalpy, positive for endothermic
    reactions.  ITC literature reports the opposite sign (protonation enthalpy).
    Use van't Hoff for any simulation where temperature deviates from 298.15 K.

    Requires components: ``tris``, ``H_plus``.
    """
    return [
        ThermodynamicReaction(
            "TrisH+ <-> Tris + H+",
            mode="equil",
            equilibrium_constant=pKa(8.072, dH=47450.0),
        ),
    ]
