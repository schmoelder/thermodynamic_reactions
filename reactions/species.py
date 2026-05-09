"""
Core types: constants, state, species, components.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

R_GAS: float = 8.314462        # J / (mol K)
KB: float = 1.380649e-23       # J / K
H_PLANCK: float = 6.626070e-34 # J s


# ---------------------------------------------------------------------------
# Physical state
# ---------------------------------------------------------------------------


@dataclass
class PhysicalState:
    """
    State at a single spatial point and instant in time.

    Built by ReactionModel.make_state(). Reaction modules receive this
    and do not need to know where T or I came from.

    Attributes
    ----------
    c : np.ndarray
        Species concentrations [mol/m³], shape (n_species,).
    T : float
        Temperature [K].
    I : float
        Ionic strength [mol/m³], derived by IonicStrengthModule.
    c_ref : np.ndarray
        Per-species standard-state concentrations [mol/m³], shape (n_species,).
    """

    c: np.ndarray
    T: float
    I: float = 0.0
    c_ref: np.ndarray = field(default_factory=lambda: np.array([]))


# ---------------------------------------------------------------------------
# Species and Component
# ---------------------------------------------------------------------------


@dataclass
class Species:
    """
    A specific chemical form of a component.

    Parameters
    ----------
    name : str
    charge : int
        Ionic charge. Default 0.
    c_ref : float
        Standard-state concentration [mol/m³]. Default 1000.
        Set to ``rho/M`` for solvents (e.g. 55500 mol/m³ for water)
        so that activity ``a = gamma * c / c_ref`` approaches 1 in
        dilute solution.
    molar_mass : float, optional
        Molar mass [kg/mol].
    density : float, optional
        Pure-component density [kg/m³].
    heat_capacity : float, optional
        Molar heat capacity [J/(mol·K)].
        Species with all three physical fields set contribute to
        ``volumetric_heat_capacity``; others are ignored.
    """

    name: str
    charge: int = 0
    c_ref: float = 1000.0          # mol/m³ standard-state concentration
    molar_mass: Optional[float] = None
    density: Optional[float] = None
    heat_capacity: Optional[float] = None  # molar Cp [J/(mol·K)]


class Component:
    """
    A conserved chemical entity that may exist as multiple species.

    Parameters
    ----------
    name : str
    species : list[Species]

    Examples
    --------
    Component("A")               # shorthand: auto-creates Species("A", charge=0)
    Component("H+", charge=+1)   # kwargs forwarded to Species

    Component("sodium", [Species("Na+", charge=+1)])

    Component("phosphate", [
        Species("H3PO4",  charge=0),
        Species("H2PO4-", charge=-1),
        Species("HPO4-2", charge=-2),
        Species("PO4-3",  charge=-3),
    ])
    """

    def __init__(self, name: str, species: list[Species] | None = None, **kwargs) -> None:
        self.name = name
        self.species = species if species is not None else [Species(name, **kwargs)]
        self._validate()

    def _validate(self) -> None:
        if not self.species:
            raise ValueError(f"Component '{self.name}' needs at least one species.")
        names = [s.name for s in self.species]
        if len(names) != len(set(names)):
            raise ValueError(f"Component '{self.name}' has duplicate species names.")

    def __repr__(self) -> str:
        return f"Component('{self.name}', species={[s.name for s in self.species]})"
