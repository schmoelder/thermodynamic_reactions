"""Rate constant models: kf(T). Enzymatic rate laws: v(state, species_index)."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np

from .species import R_GAS
from .state import State

__all__ = [
    "RateConstantBase",
    "RateConstantFixed",
    "RateConstantArrhenius",
    "RateConstantPolynomial",
    "RateConstantTabulated",
    "RateBase",
    "MichaelisMenten",
    "HillRate",
    "CustomRate",
]


# ---------------------------------------------------------------------------
# Forward rate constant  kf(T)
# kr is never stored — always derived as kf(T) / K(T) by ThermodynamicReaction
# ---------------------------------------------------------------------------


class RateConstantBase(ABC):
    """
    Abstract base for forward rate constant models.

    Units convention
    ----------------
    Rate constants in this framework operate on dimensionless activities
    (a_i = gamma_i * c_i / C_REF) when used in ThermodynamicReaction,
    and on concentrations [mol/m³] when used in MassActionReaction.

    This means kf has different effective units depending on reaction type
    and order (n = sum of forward stoichiometric exponents):

        MassActionReaction:     kf  in  (m³/mol)^(n-1) / s
        ThermodynamicReaction:  kf  in  mol/(m³·s)  for any order

    For any reaction order n, the numerical value of kf differs between
    the two classes by a factor of C_REF^n:

        kf_thermo = kf_massaction * C_REF^n

    where C_REF = 1000 mol/m³.  For first-order reactions (n=1) this
    means kf_thermo = kf_massaction * 1000: the same numerical value
    gives dynamics ~1000× slower in ThermodynamicReaction.
    When fitting kf to experimental data, be consistent: use the same
    reaction class for fitting and prediction.
    """

    @abstractmethod
    def kf(self, T: float) -> float:
        """Forward rate constant at temperature T [K]."""

    def dlnkf_dT(self, T: float) -> Optional[float]:
        """
        Return d(ln kf)/dT [1/K], or None if no analytic form is available.

        None signals that the caller must fall back to finite differences.
        """
        return None


@dataclass
class RateConstantFixed(RateConstantBase):
    """
    Temperature-independent forward rate constant.

    Parameters
    ----------
    kf_value : float
        Forward rate constant. Units depend on which reaction class this
        is used with and the reaction order — see RateConstantBase docstring.
    """

    kf_value: float

    def kf(self, T: float) -> float:
        return self.kf_value

    def dlnkf_dT(self, T: float) -> float:
        return 0.0


@dataclass
class RateConstantArrhenius(RateConstantBase):
    """
    k(T) = A · exp(-Ea / (R·T)).

    Parameters
    ----------
    A : float
        Pre-exponential factor.  Same unit convention as RateConstantFixed —
        see RateConstantBase docstring.
    Ea : float
        Activation energy [J/mol].
    """

    A: float
    Ea: float

    def kf(self, T: float) -> float:
        return self.A * float(np.exp(-self.Ea / (R_GAS * T)))

    def dlnkf_dT(self, T: float) -> float:
        return self.Ea / (R_GAS * T**2)


@dataclass
class RateConstantPolynomial(RateConstantBase):
    """
    kf(T) = exp(b₀ + b₁T + b₂T² + ...).

    coeffs = [b₀, b₁, b₂, ...] in ascending power order.

    Parameters
    ----------
    coeffs : array-like
        Polynomial coefficients [b₀, b₁, ...], ascending power order.
    """

    coeffs: np.ndarray

    def __post_init__(self) -> None:
        self.coeffs = np.asarray(self.coeffs, dtype=float)
        self._powers = np.arange(len(self.coeffs), dtype=float)
        self._deriv_coeffs = self._powers * self.coeffs

    def kf(self, T: float) -> float:
        return float(np.exp(np.dot(self.coeffs, T**self._powers)))

    def dlnkf_dT(self, T: float) -> float:
        powers_m1 = np.where(self._powers > 0, self._powers - 1, 0.0)
        return float(np.dot(self._deriv_coeffs, T**powers_m1))


@dataclass
class RateConstantTabulated(RateConstantBase):
    """
    kf(T) from linearly interpolated tabulated data.

    Parameters
    ----------
    T_data : array-like   Temperature values [K], monotonically increasing.
    kf_data : array-like  Forward rate constants at each temperature.
    """

    T_data: np.ndarray
    kf_data: np.ndarray

    def __post_init__(self) -> None:
        self.T_data = np.asarray(self.T_data, dtype=float)
        self.kf_data = np.asarray(self.kf_data, dtype=float)

    def kf(self, T: float) -> float:
        return float(np.interp(T, self.T_data, self.kf_data))


# ---------------------------------------------------------------------------
# Enzymatic rate laws  v(state, species_index)
# ---------------------------------------------------------------------------


class RateBase(ABC):
    """Abstract base for enzymatic / custom rate laws."""

    @abstractmethod
    def __call__(self, state: State, species_index: dict[str, int]) -> float:
        """Net reaction rate [mol/(m³·s)]."""


@dataclass
class MichaelisMenten(RateBase):
    """
    v = Vmax · [S] / (Km + [S]).

    Parameters
    ----------
    Vmax : float      Maximum rate [mol/(m³·s)].
    Km : float        Michaelis constant [mol/m³].
    substrate : str   Name of the substrate species.
    """

    Vmax: float
    Km: float
    substrate: str

    def __call__(self, state: State, species_index: dict[str, int]) -> float:
        S = state.c[species_index[self.substrate]]
        return self.Vmax * S / (self.Km + S)


@dataclass
class HillRate(RateBase):
    """
    v = Vmax · [S]^n / (Km^n + [S]^n).

    Parameters
    ----------
    Vmax : float      Maximum rate [mol/(m³·s)].
    Km : float        Half-saturation constant [mol/m³].
    n : float         Hill coefficient [-].
    substrate : str   Name of the substrate species.
    """

    Vmax: float
    Km: float
    n: float
    substrate: str

    def __call__(self, state: State, species_index: dict[str, int]) -> float:
        S = state.c[species_index[self.substrate]]
        return self.Vmax * S**self.n / (self.Km**self.n + S**self.n)


@dataclass
class CustomRate(RateBase):
    """
    User-supplied rate function.

    Parameters
    ----------
    fn : callable(state, species_index) -> float
    """

    fn: Callable[[State, dict[str, int]], float]

    def __call__(self, state: State, species_index: dict[str, int]) -> float:
        return self.fn(state, species_index)
