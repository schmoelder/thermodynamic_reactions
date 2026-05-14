"""
formulation.py.
==============
Solution: converts a recipe (solutes + solvent fractions) to c0 and prescribed.
"""

from __future__ import annotations

from typing import Union

from .species import Component

__all__ = ["Solution"]


class Solution:
    """
    Convert a buffer recipe to concentration dicts for the solver.

    Parameters
    ----------
    solvents : Component or dict[Component, float]
        Solvent(s). A single ``Component`` is treated as pure solvent
        (volume fraction 1.0). A dict maps each ``Component`` to its volume
        fraction; fractions must sum to 1.
        Each solvent species must carry ``density`` and ``molar_mass`` so that
        the pure-component concentration ``c_ref = density / molar_mass`` is
        defined.
    solutes : dict[str, float], optional
        Solute concentrations {species_name: value [mol/m³]}.

    Attributes
    ----------
    c0 : dict[str, float]
        All species concentrations [mol/m³]. Solvent species are set to
        ``phi * c_ref``; solute species are taken directly from the recipe.
    prescribed : dict[str, float]
        Solvent concentrations for passing to the solver as ``prescribed``.
        Identical to the solvent entries in ``c0``.

    Examples
    --------
    Pure-solvent shorthand::

        sol = Solution(water, solutes={"HAc": 100.0, "H+": 1e-4, "OH-": 1e-7})
        c_eq = solve_equilibrium(model, sol.c0, T=298.15, prescribed=sol.prescribed)

    Mixed solvents::

        sol = Solution({water: 0.9, ethanol: 0.1}, solutes={"HAc": 100.0})
    """

    def __init__(
        self,
        solvents: Union[Component, dict],
        solutes: dict[str, float] | None = None,
    ) -> None:
        if isinstance(solvents, Component):
            solvents = {solvents: 1.0}

        total = sum(solvents.values())
        if abs(total - 1.0) > 1e-9:
            raise ValueError(
                f"Solvent volume fractions must sum to 1.0; got {total:.6f}."
            )

        for comp in solvents:
            for sp in comp.species:
                if sp.density is None or sp.molar_mass is None:
                    raise ValueError(
                        f"Solvent species '{sp.name}' requires density and "
                        "molar_mass to determine its pure-component concentration."
                    )

        self._c0: dict[str, float] = {}
        self._prescribed: dict[str, float] = {}

        for comp, phi in solvents.items():
            for sp in comp.species:
                c = phi * sp.c_ref
                self._c0[sp.name] = c
                self._prescribed[sp.name] = c

        for name, conc in (solutes or {}).items():
            self._c0[name] = conc

    @property
    def c0(self) -> dict[str, float]:
        return dict(self._c0)

    @property
    def prescribed(self) -> dict[str, float]:
        return dict(self._prescribed)
