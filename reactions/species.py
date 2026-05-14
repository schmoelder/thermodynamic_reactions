"""Species and components: constants, Species, Component."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

__all__ = [
    "R_GAS",
    "KB",
    "H_PLANCK",
    "Species",
    "Component",
]

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

R_GAS: float = 8.314462  # J / (mol K)
KB: float = 1.380649e-23  # J / K
H_PLANCK: float = 6.626070e-34  # J s


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
    c_ref : float, optional
        Standard-state concentration [mol/m³].
        If not given and both ``density`` and ``molar_mass`` are set,
        defaults to ``density / molar_mass`` (the pure-component molar
        concentration, e.g. ~55 556 mol/m³ for water).
        Otherwise defaults to 1000 mol/m³ (1 mol/L, the conventional
        aqueous standard state for solutes).
    molar_mass : float, optional
        Molar mass [kg/mol].
        Together with ``density``, marks this species as a volume-accounting
        solvent: its concentration is derived from the remaining volume rather
        than specified by the user.
    density : float, optional
        Pure-component density [kg/m³].
    heat_capacity : float, optional
        Molar heat capacity [J/(mol·K)].
        Species with all three physical fields set contribute to
        ``volumetric_heat_capacity``; others are ignored.
    """

    name: str
    charge: int = 0
    c_ref: Optional[float] = None
    molar_mass: Optional[float] = None
    density: Optional[float] = None
    heat_capacity: Optional[float] = None  # molar Cp [J/(mol·K)]

    def __post_init__(self) -> None:
        if self.c_ref is None:
            if self.density is not None and self.molar_mass is not None:
                self.c_ref = self.density / self.molar_mass
            else:
                self.c_ref = 1000.0


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

    def __init__(
        self, name: str, species: list[Species] | None = None, **kwargs
    ) -> None:
        self.name = name
        self.species = species if species is not None else [Species(name, **kwargs)]
        self._validate()

    def _validate(self) -> None:
        if not self.species:
            raise ValueError(f"Component '{self.name}' needs at least one species.")
        names = [s.name for s in self.species]
        if len(names) != len(set(names)):
            raise ValueError(f"Component '{self.name}' has duplicate species names.")

    def _single_attr(self, attr: str):
        if len(self.species) != 1:
            raise AttributeError(
                f"Component '{self.name}' has {len(self.species)} species; "
                f"use .species[i].{attr} explicitly."
            )
        return getattr(self.species[0], attr)

    @property
    def charge(self) -> int:
        return self._single_attr("charge")

    @property
    def c_ref(self) -> float:
        return self._single_attr("c_ref")

    @property
    def molar_mass(self) -> Optional[float]:
        return self._single_attr("molar_mass")

    @property
    def density(self) -> Optional[float]:
        return self._single_attr("density")

    @property
    def heat_capacity(self) -> Optional[float]:
        return self._single_attr("heat_capacity")

    def __repr__(self) -> str:
        return f"Component('{self.name}', species={[s.name for s in self.species]})"
