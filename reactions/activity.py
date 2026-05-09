"""
Activity coefficient models: γ(state, charges).
"""

from __future__ import annotations

import warnings
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Callable, Optional, Union

import numpy as np

from .species import PhysicalState

# ---------------------------------------------------------------------------
# Debye-Hückel physical constants
# A_L(T) = _DH_A_CONST / (εr·T)^(3/2)   [(mol/L)^(-1/2)]
# B_SI(T) = _DH_B_SI_CONST / sqrt(εr·T)  [m^(-1)·(mol/m³)^(-1/2)]
# ---------------------------------------------------------------------------

_DH_A_CONST: float = 1.8246e6
_DH_B_SI_CONST: float = 1.5908e10


def _water_epsilon_r(T: float) -> float:
    """Malmberg-Maryott (1956) relative permittivity of liquid water, valid 0–60 °C."""
    t = T - 273.15
    return 87.740 - 0.4008 * t + 9.398e-4 * t**2 - 1.410e-6 * t**3


# ---------------------------------------------------------------------------
# Activity coefficient base and models
# ---------------------------------------------------------------------------


class ActivityCoefficientBase(ABC):
    """
    Abstract base for activity coefficient models.

    Returns γᵢ for each dynamic species.
    Thermodynamic activity: aᵢ = γᵢ · cᵢ / c_ref_i.
    """

    @abstractmethod
    def activity(self, state: PhysicalState, charges: np.ndarray) -> np.ndarray:
        """
        Return activity coefficients γᵢ, shape (n_species,).

        Parameters
        ----------
        state : PhysicalState
        charges : np.ndarray
            Ionic charges for each dynamic species.
        """


@dataclass
class ActivityCoefficientIdeal(ActivityCoefficientBase):
    """γᵢ = 1 for all species."""

    def activity(self, state: PhysicalState, charges: np.ndarray) -> np.ndarray:
        return np.ones(len(state.c))


@dataclass
class ActivityCoefficientDebyeHuckel(ActivityCoefficientBase):
    """
    Extended Debye-Hückel:
        log10(γᵢ) = -A · zᵢ² · √I / (1 + B · a_ion · √I)

    Valid up to I ~ 100 mol/m³ (0.1 mol/L).

    Parameters
    ----------
    A : float
        [(mol/m³)^(-1/2)]. Default: 25 °C water value (0.509/√1000).
        Ignored when epsilon_r is provided.
    B : float
        [m^(-1)·(mol/m³)^(-1/2)]. Default: 25 °C water value (3.28e9/√1000).
        Ignored when epsilon_r is provided.
    a_ion : float
        Mean ion-size parameter [m]. Typical: 3e-10 m.
    epsilon_r : float or callable, optional
        Relative permittivity εr of the solvent.  If callable, called as
        epsilon_r(T) at each evaluation.  When provided, A and B are computed
        from εr and T via the Debye-Hückel formula (A ∝ (εr·T)^(-3/2),
        B ∝ (εr·T)^(-1/2)) and the stored A and B fields are ignored.
        When None (default), A and B are used as-is and a UserWarning is
        emitted when state.T deviates > 5 K from 298.15 K.
    """

    A: float = 0.509 / (1000.0**0.5)
    B: float = 3.28e9 / (1000.0**0.5)
    a_ion: float = 3e-10
    epsilon_r: Optional[Union[float, Callable[[float], float]]] = None

    def activity(self, state: PhysicalState, charges: np.ndarray) -> np.ndarray:
        T = state.T
        if self.epsilon_r is not None:
            er = float(self.epsilon_r(T)) if callable(self.epsilon_r) else float(self.epsilon_r)
            A = _DH_A_CONST / ((er * T) ** 1.5 * 1000.0 ** 0.5)
            B = _DH_B_SI_CONST / (er * T) ** 0.5
        else:
            if abs(T - 298.15) > 5.0:
                warnings.warn(
                    f"ActivityCoefficientDebyeHuckel: T={T:.1f} K deviates more than "
                    "5 K from 298.15 K and epsilon_r was not provided. "
                    "A and B are fixed at 25 °C water values. "
                    "Pass epsilon_r=_water_epsilon_r for T-dependent behaviour.",
                    UserWarning,
                    stacklevel=2,
                )
            A = self.A
            B = self.B
        sqrt_I = np.sqrt(state.I)
        log_gamma = (
            -A * charges**2 * sqrt_I
            / (1.0 + B * self.a_ion * sqrt_I)
        )
        return 10.0**log_gamma


@dataclass
class ActivityCoefficientDavies(ActivityCoefficientBase):
    """
    Davies equation:
        log10(γᵢ) = -A · zᵢ² · (√I_L / (1 + √I_L) - 0.3 · I_L)

    where I_L = I / 1000 is ionic strength in mol/L.
    Valid up to I ~ 500 mol/m³ (0.5 mol/L).

    Parameters
    ----------
    A : float
        At 25 °C in water: 0.509 (L/mol)^0.5. Ignored when epsilon_r is provided.
    epsilon_r : float or callable, optional
        Relative permittivity εr of the solvent.  If callable, called as
        epsilon_r(T) at each evaluation.  When provided, A is computed from
        εr and T via A ∝ (εr·T)^(-3/2) and the stored A field is ignored.
        When None (default), A is used as-is and a UserWarning is emitted
        when state.T deviates > 5 K from 298.15 K.
    """

    A: float = 0.509
    epsilon_r: Optional[Union[float, Callable[[float], float]]] = None

    def activity(self, state: PhysicalState, charges: np.ndarray) -> np.ndarray:
        T = state.T
        if self.epsilon_r is not None:
            er = float(self.epsilon_r(T)) if callable(self.epsilon_r) else float(self.epsilon_r)
            A = _DH_A_CONST / (er * T) ** 1.5
        else:
            if abs(T - 298.15) > 5.0:
                warnings.warn(
                    f"ActivityCoefficientDavies: T={T:.1f} K deviates more than "
                    "5 K from 298.15 K and epsilon_r was not provided. "
                    "A is fixed at the 25 °C water value. "
                    "Pass epsilon_r=_water_epsilon_r for T-dependent behaviour.",
                    UserWarning,
                    stacklevel=2,
                )
            A = self.A
        I_L = state.I / 1000.0
        sqrt_I = np.sqrt(I_L)
        log_gamma = -A * charges**2 * (sqrt_I / (1.0 + sqrt_I) - 0.3 * I_L)
        return 10.0**log_gamma


@dataclass
class ActivityCoefficientCustom(ActivityCoefficientBase):
    """
    User-supplied activity coefficient function.

    Parameters
    ----------
    fn : callable(state, charges) -> np.ndarray
    """

    fn: Callable[[PhysicalState, np.ndarray], np.ndarray]

    def activity(self, state: PhysicalState, charges: np.ndarray) -> np.ndarray:
        return self.fn(state, charges)
