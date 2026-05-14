"""
reactions — thermodynamic and kinetic reaction models for CADET.

Primary imports
---------------
from reactions.species import R_GAS, KB, PhysicalState, Species, Component
from reactions.activity import ActivityCoefficientDebyeHuckel, ActivityCoefficientDavies
from reactions.ionic import IonicStrengthIdeal, IonicStrengthBackground
from reactions.equilibrium import (
    EquilibriumConstant,
    EquilibriumConstantVantHoff,
    EquilibriumConstantVantHoffCp,
    pKa,
)
from reactions.rate import RateConstantFixed, RateConstantArrhenius, MichaelisMenten, HillRate
from reactions.reaction import MassActionReaction, ThermodynamicReaction, EnzymaticReaction
from reactions.model import ReactionModel, ConservationReport
from reactions.formulation import Solution
from reactions.common import water, H_plus, OH_minus, autoionisation, acetic_acid_equilibria
from reactions.solver import simulate, solve_equilibrium, SimulationResult
"""
