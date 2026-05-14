"""State-dependent parameter functions (modulators) with optional analytic gradients."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from .species import R_GAS
from .state import AuxiliaryState, State

if TYPE_CHECKING:
    from .rate import RateBase

__all__ = [
    "StateModulatorBase",
    "ConstantModulator",
    "LinearModulator",
    "PowerLawModulator",
    "ExponentialModulator",
    "SigmoidalModulator",
    "ArrheniusModulator",
    "BellModulator",
    "ModulatedRate",
]


class StateModulatorBase(ABC):
    """
    Abstract base for state-dependent scalar multipliers.

    Concrete classes implement value(); analytic gradient methods are optional.
    Returning None signals that the caller must fall back to finite differences.
    """

    @abstractmethod
    def value(
        self,
        state: State,
        aux: AuxiliaryState,
        species_index: dict[str, int],
    ) -> float:
        """Return scalar multiplier evaluated at the current state."""

    def dvalue_dc(
        self,
        state: State,
        aux: AuxiliaryState,
        species_index: dict[str, int],
    ) -> np.ndarray | None:
        """Return gradient of value w.r.t. concentration vector, or None (FD fallback)."""
        return None

    def dvalue_dT(
        self,
        state: State,
        aux: AuxiliaryState,
        species_index: dict[str, int],
    ) -> float | None:
        """Return derivative of value w.r.t. temperature, or None (FD fallback)."""
        return None


@dataclass
class ConstantModulator(StateModulatorBase):
    """Constant multiplier independent of state."""

    c: float

    def value(
        self, state: State, aux: AuxiliaryState, species_index: dict[str, int]
    ) -> float:
        return self.c

    def dvalue_dc(
        self,
        state: State,
        aux: AuxiliaryState,
        species_index: dict[str, int],
    ) -> np.ndarray:
        return np.zeros(len(state.c))

    def dvalue_dT(
        self,
        state: State,
        aux: AuxiliaryState,
        species_index: dict[str, int],
    ) -> float:
        return 0.0


@dataclass
class LinearModulator(StateModulatorBase):
    """slope * c[species] + intercept."""

    species: str
    slope: float
    intercept: float = 0.0

    def value(
        self, state: State, aux: AuxiliaryState, species_index: dict[str, int]
    ) -> float:
        return self.slope * state.c[species_index[self.species]] + self.intercept

    def dvalue_dc(
        self,
        state: State,
        aux: AuxiliaryState,
        species_index: dict[str, int],
    ) -> np.ndarray:
        g = np.zeros(len(state.c))
        g[species_index[self.species]] = self.slope
        return g

    def dvalue_dT(
        self,
        state: State,
        aux: AuxiliaryState,
        species_index: dict[str, int],
    ) -> float:
        return 0.0


@dataclass
class PowerLawModulator(StateModulatorBase):
    """scale * c[species]^exponent."""

    species: str
    exponent: float
    scale: float = 1.0

    def value(
        self, state: State, aux: AuxiliaryState, species_index: dict[str, int]
    ) -> float:
        c = state.c[species_index[self.species]]
        return self.scale * c**self.exponent

    def dvalue_dc(
        self,
        state: State,
        aux: AuxiliaryState,
        species_index: dict[str, int],
    ) -> np.ndarray:
        c = state.c[species_index[self.species]]
        g = np.zeros(len(state.c))
        g[species_index[self.species]] = (
            self.scale * self.exponent * c ** (self.exponent - 1)
        )
        return g

    def dvalue_dT(
        self,
        state: State,
        aux: AuxiliaryState,
        species_index: dict[str, int],
    ) -> float:
        return 0.0


@dataclass
class ExponentialModulator(StateModulatorBase):
    """scale * exp(rate_const * c[species])."""

    species: str
    rate_const: float
    scale: float = 1.0

    def value(
        self, state: State, aux: AuxiliaryState, species_index: dict[str, int]
    ) -> float:
        c = state.c[species_index[self.species]]
        return self.scale * np.exp(self.rate_const * c)

    def dvalue_dc(
        self,
        state: State,
        aux: AuxiliaryState,
        species_index: dict[str, int],
    ) -> np.ndarray:
        g = np.zeros(len(state.c))
        g[species_index[self.species]] = self.rate_const * self.value(
            state, aux, species_index
        )
        return g

    def dvalue_dT(
        self,
        state: State,
        aux: AuxiliaryState,
        species_index: dict[str, int],
    ) -> float:
        return 0.0


@dataclass
class SigmoidalModulator(StateModulatorBase):
    """1 / (1 + exp(-k * (c[species] - c_half)))."""

    species: str
    k: float
    c_half: float

    def value(
        self, state: State, aux: AuxiliaryState, species_index: dict[str, int]
    ) -> float:
        c = state.c[species_index[self.species]]
        return 1.0 / (1.0 + np.exp(-self.k * (c - self.c_half)))

    def dvalue_dc(
        self,
        state: State,
        aux: AuxiliaryState,
        species_index: dict[str, int],
    ) -> np.ndarray:
        v = self.value(state, aux, species_index)
        g = np.zeros(len(state.c))
        g[species_index[self.species]] = self.k * v * (1.0 - v)
        return g

    def dvalue_dT(
        self,
        state: State,
        aux: AuxiliaryState,
        species_index: dict[str, int],
    ) -> float:
        return 0.0


@dataclass
class ArrheniusModulator(StateModulatorBase):
    """A * exp(-Ea / (R * T)). Analytic dvalue_dT."""

    A: float
    Ea: float

    def value(
        self, state: State, aux: AuxiliaryState, species_index: dict[str, int]
    ) -> float:
        return self.A * np.exp(-self.Ea / (R_GAS * state.T))

    def dvalue_dc(
        self,
        state: State,
        aux: AuxiliaryState,
        species_index: dict[str, int],
    ) -> np.ndarray:
        return np.zeros(len(state.c))

    def dvalue_dT(
        self,
        state: State,
        aux: AuxiliaryState,
        species_index: dict[str, int],
    ) -> float:
        return self.value(state, aux, species_index) * self.Ea / (R_GAS * state.T**2)


@dataclass
class BellModulator(StateModulatorBase):
    """
    Three-state protonation activity fraction.

    Models pH-dependent enzyme activity via the middle Bjerrum fraction:
    f_active = Ka1 * a(H+) / (a(H+)^2 + Ka1 * a(H+) + Ka1 * Ka2)

    where Ka1 = 10^(-pKa1) and Ka2 = 10^(-pKa2) are thermodynamic constants
    and a(H+) = gamma * c(H+) / c_ref is the proton activity.

    f_active is maximised at a(H+) = sqrt(Ka1 * Ka2), i.e. pH_opt = (pKa1 + pKa2) / 2.
    The peak value approaches 1 as pKa2 - pKa1 → ∞ but is always less than 1 for finite
    pKa separation.
    """

    pKa1: float
    pKa2: float
    proton_species: str

    def _activity_H(
        self,
        state: State,
        aux: AuxiliaryState,
        species_index: dict[str, int],
    ) -> float:
        i = species_index[self.proton_species]
        return aux.gamma[i] * state.c[i] / aux.c_ref[i]

    def value(
        self, state: State, aux: AuxiliaryState, species_index: dict[str, int]
    ) -> float:
        a_H = self._activity_H(state, aux, species_index)
        Ka1 = 10.0 ** (-self.pKa1)
        Ka2 = 10.0 ** (-self.pKa2)
        denom = a_H**2 + Ka1 * a_H + Ka1 * Ka2
        return Ka1 * a_H / denom

    def dvalue_dc(
        self,
        state: State,
        aux: AuxiliaryState,
        species_index: dict[str, int],
    ) -> np.ndarray:
        i = species_index[self.proton_species]
        a_H = self._activity_H(state, aux, species_index)
        Ka1 = 10.0 ** (-self.pKa1)
        Ka2 = 10.0 ** (-self.pKa2)
        denom = a_H**2 + Ka1 * a_H + Ka1 * Ka2
        # df/da_H = Ka1 * (Ka1*Ka2 - a_H^2) / denom^2
        df_daH = Ka1 * (Ka1 * Ka2 - a_H**2) / denom**2
        da_H_dci = aux.gamma[i] / aux.c_ref[i]
        g = np.zeros(len(state.c))
        g[i] = df_daH * da_H_dci
        return g

    def dvalue_dT(
        self,
        state: State,
        aux: AuxiliaryState,
        species_index: dict[str, int],
    ) -> float:
        return 0.0


@dataclass
class ModulatedRate:
    """
    Product of a base rate and a state-dependent modulator.

    v(state, aux) = modulator.value(state, aux) * base_rate(state, aux)

    Acts as a RateBase: pass it as the ``rate`` argument of EnzymaticReaction.
    """

    base_rate: "RateBase"
    modulator: StateModulatorBase

    def __call__(
        self,
        state: State,
        aux: AuxiliaryState,
        species_index: dict[str, int],
    ) -> float:
        return self.modulator.value(state, aux, species_index) * self.base_rate(
            state,
            aux,
            species_index,
        )
