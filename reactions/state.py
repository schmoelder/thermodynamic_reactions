"""State and AuxiliaryState: runtime simulation state containers."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import xarray as xr

__all__ = [
    "State",
    "AuxiliaryState",
]


class State:
    """
    True ODE state at a single point in time.

    Holds the quantities the solver integrates (concentrations, and
    temperature when an energy balance is active). The backing array
    ``_s`` is the single source of truth and is never replaced after
    construction — only mutated in-place — so DataArray views remain
    valid across solver steps.

    Parameters
    ----------
    name : str
        Identifier (e.g. ``"bulk"``).
    entries : dict[str, list[str] | int]
        Entry name → list of coordinate labels (strings) or a plain
        count (int, converted to ``range(n)``).
        Example: ``{"c": ["HAc", "Ac-", "H+", "OH-"]}``.
    T : float
        Temperature [K].  Stored separately from ``_s``; becomes part
        of the solver buffer when energy balance is added later.
    """

    def __init__(
        self,
        name: str,
        entries: dict[str, list[str] | int],
        T: float = 298.15,
    ) -> None:
        self.name = name
        self.entries: dict[str, list] = {
            k: v if isinstance(v, list) else list(range(v)) for k, v in entries.items()
        }
        self.T = T

        n_total = sum(len(coords) for coords in self.entries.values())
        self._s = np.zeros(n_total)

        self._entry_slices: dict[str, slice] = {}
        i = 0
        for entry, coords in self.entries.items():
            n = len(coords)
            self._entry_slices[entry] = slice(i, i + n)
            i += n

        self._field: dict[str, xr.DataArray] = self._build_field()

    def _build_field(self) -> dict[str, xr.DataArray]:
        result = {}
        for entry, coords in self.entries.items():
            data = self._s[self._entry_slices[entry]]
            dim = "species" if coords and isinstance(coords[0], str) else "component"
            result[entry] = xr.DataArray(
                data,
                dims=[dim],
                coords={dim: coords},
                name=entry,
            )
        return result

    def wire_to_buffer(self, buf: np.ndarray, offset: int) -> None:
        """Replace ``_s`` with a view into ``buf[offset:offset+n_dof]``.

        Called once before the first solver step. Rebuilds DataArray views
        so they point at the new memory. Any previously held DataArray
        references become stale after this call.
        """
        self._s = buf[offset : offset + self.n_dof]
        self._field = self._build_field()

    @property
    def n_dof(self) -> int:
        return len(self._s)

    @property
    def s_flat(self) -> np.ndarray:
        return self._s

    @s_flat.setter
    def s_flat(self, value: np.ndarray) -> None:
        self._s[...] = value

    @property
    def c(self) -> np.ndarray:
        """Concentration array [mol/m³], shape (n_species,). View into ``_s``."""
        return self._s[self._entry_slices["c"]]

    @c.setter
    def c(self, value: np.ndarray) -> None:
        self._s[self._entry_slices["c"]] = value

    @property
    def field(self) -> dict[str, xr.DataArray]:
        """Zero-copy DataArray views per entry. Always reflects current ``_s``."""
        return self._field

    def __getitem__(self, entry: str) -> xr.DataArray:
        try:
            return self._field[entry]
        except KeyError:
            raise KeyError(f"'{entry}' is not a valid state entry.")

    def __setitem__(self, entry: str, value: np.ndarray) -> None:
        if entry not in self._entry_slices:
            raise KeyError(f"'{entry}' is not a valid state entry.")
        self._s[self._entry_slices[entry]] = value

    def copy(self) -> "State":
        """Return a new State with the same entries, T, and a copy of ``_s``."""
        new = State(self.name, {k: list(v) for k, v in self.entries.items()}, T=self.T)
        new._s[:] = self._s
        return new

    def with_T(self, T: float) -> "State":
        """Return a copy of this state with a different temperature."""
        new = self.copy()
        new.T = T
        return new

    def __repr__(self) -> str:
        return (
            f"State(name='{self.name}', "
            f"entries={list(self.entries.keys())}, "
            f"n_dof={self.n_dof}, T={self.T})"
        )


@dataclass
class AuxiliaryState:
    """
    Quantities derived from State; not integrated by the solver.

    Built by ReactionModel.make_aux() at each residual evaluation.

    Attributes
    ----------
    I : float
        Ionic strength [mol/m³].
    c_ref : np.ndarray
        Per-species standard-state concentrations [mol/m³], shape (n_species,).
    """

    I: float
    c_ref: np.ndarray
