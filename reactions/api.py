"""
ReactionModel API
=================
Public-facing re-export facade.

All symbols are defined in their own submodules:
    reactions.species    — constants, PhysicalState, Species, Component
    reactions.ionic      — IonicStrengthIdeal, IonicStrengthBackground, IonicStrengthFixed
    reactions.activity   — ActivityCoefficient*, _water_epsilon_r
    reactions.equilibrium — EquilibriumConstant*, pKa
    reactions.rate       — RateConstant*, RateBase, MichaelisMenten, HillRate, CustomRate
    reactions.reaction   — parse_stoichiometry, ReactionBase, Thermodynamic/MassAction/...
    reactions.model      — ConservationReport, ReactionModel
    reactions.analysis   — buffer_capacity, speciation_fractions, solve_equilibrium_sweep

Units (SI throughout)
---------------------
Concentration  : mol / m³
Temperature    : K
Energy         : J / mol
Molar mass     : kg / mol
Density        : kg / m³

Equilibrium constants are dimensionless (activities, not concentrations).
Activities are a_i = γ_i · c_i / c_ref_i where c_ref_i is per-species (Species.c_ref,
default 1000 mol/m³).
"""

from .activity import (
    ActivityCoefficientBase,
    ActivityCoefficientCustom,
    ActivityCoefficientDavies,
    ActivityCoefficientDebyeHuckel,
    ActivityCoefficientIdeal,
    _DH_A_CONST,
    _DH_B_L_CONST,
    _water_epsilon_r,
)
from .equilibrium import (
    EquilibriumConstant,
    EquilibriumConstantBase,
    EquilibriumConstantCustom,
    EquilibriumConstantPolynomial,
    EquilibriumConstantTabulated,
    EquilibriumConstantVantHoff,
    EquilibriumConstantVantHoffCp,
    pKa,
)
from .ionic import (
    IonicStrengthBackground,
    IonicStrengthBase,
    IonicStrengthFixed,
    IonicStrengthIdeal,
)
from .model import ConservationReport, ReactionModel
from .rate import (
    CustomRate,
    HillRate,
    MichaelisMenten,
    RateBase,
    RateConstantArrhenius,
    RateConstantBase,
    RateConstantFixed,
    RateConstantPolynomial,
    RateConstantTabulated,
)
from .reaction import (
    CustomReaction,
    EnzymaticReaction,
    MassActionReaction,
    ReactionBase,
    ThermodynamicReaction,
    parse_stoichiometry,
)
from .common import (
    H_plus,
    OH_minus,
    acetic_acid,
    acetic_acid_equilibria,
    autoionisation,
    citric_acid,
    citric_acid_equilibria,
    hepes,
    hepes_equilibria,
    mops,
    mops_equilibria,
    phosphate,
    phosphate_equilibria,
    tris,
    tris_equilibria,
    water,
)
from .analysis import buffer_capacity, speciation_fractions, solve_equilibrium_sweep
from .formulation import Solution
from .species import (
    KB,
    R_GAS,
    H_PLANCK,
    Component,
    PhysicalState,
    Species,
)

__all__ = [
    # constants
    "R_GAS",
    "KB",
    "H_PLANCK",
    # state
    "PhysicalState",
    # species / components
    "Species",
    "Component",
    # ionic strength
    "IonicStrengthBase",
    "IonicStrengthIdeal",
    "IonicStrengthBackground",
    "IonicStrengthFixed",
    # activity coefficients
    "ActivityCoefficientBase",
    "ActivityCoefficientIdeal",
    "ActivityCoefficientDebyeHuckel",
    "ActivityCoefficientDavies",
    "ActivityCoefficientCustom",
    "_DH_A_CONST",
    "_DH_B_L_CONST",
    "_water_epsilon_r",
    # equilibrium constants
    "EquilibriumConstantBase",
    "EquilibriumConstant",
    "EquilibriumConstantVantHoff",
    "EquilibriumConstantVantHoffCp",
    "EquilibriumConstantCustom",
    "EquilibriumConstantTabulated",
    "EquilibriumConstantPolynomial",
    "pKa",
    # rate constants
    "RateConstantBase",
    "RateConstantFixed",
    "RateConstantArrhenius",
    "RateConstantPolynomial",
    "RateConstantTabulated",
    # enzymatic rate laws
    "RateBase",
    "MichaelisMenten",
    "HillRate",
    "CustomRate",
    # reactions
    "parse_stoichiometry",
    "ReactionBase",
    "ThermodynamicReaction",
    "MassActionReaction",
    "EnzymaticReaction",
    "CustomReaction",
    # model
    "ConservationReport",
    "ReactionModel",
    # analysis
    "buffer_capacity",
    "speciation_fractions",
    "solve_equilibrium_sweep",
    # formulation
    "Solution",
    # common components and reaction factories
    "water",
    "H_plus",
    "OH_minus",
    "acetic_acid",
    "acetic_acid_equilibria",
    "autoionisation",
    "citric_acid",
    "citric_acid_equilibria",
    "hepes",
    "hepes_equilibria",
    "mops",
    "mops_equilibria",
    "phosphate",
    "phosphate_equilibria",
    "tris",
    "tris_equilibria",
]


# ---------------------------------------------------------------------------
# Backwards-compatibility guard for removed names
# ---------------------------------------------------------------------------


def __getattr__(name: str):
    if name == "C_REF":
        raise AttributeError(
            "'C_REF' has been removed from reactions.api. "
            "Use Species.c_ref (default 1000.0 mol/m³) for per-species "
            "standard-state concentrations, or write 1000.0 directly."
        )
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
