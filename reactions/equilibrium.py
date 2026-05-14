"""Equilibrium constant models: K(T)."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np

from .species import R_GAS

__all__ = [
    "EquilibriumConstantBase",
    "EquilibriumConstant",
    "EquilibriumConstantVantHoff",
    "EquilibriumConstantVantHoffCp",
    "EquilibriumConstantCustom",
    "EquilibriumConstantTabulated",
    "EquilibriumConstantPolynomial",
    "pKa",
]


class EquilibriumConstantBase(ABC):
    """Abstract base for equilibrium constant models."""

    @abstractmethod
    def K(self, T: float) -> float:  # noqa: N802
        """Return dimensionless equilibrium constant at temperature T [K]."""

    def dlnK_dT(self, T: float) -> Optional[float]:  # noqa: N802
        """
        Return d(ln K)/dT [1/K], or None if no analytic form is available.

        None signals that the caller must fall back to finite differences.
        """
        return None

    def reaction_enthalpy(self, T: float, eps: float = 1e-4) -> float:
        """
        Return ΔrH°(T) [J/mol].

        Uses the analytic form if dlnK_dT() returns a value; otherwise falls
        back to central-difference FD via the van't Hoff relation:
            ΔrH°(T) = R T² · d(ln K)/dT
        """
        dlnK = self.dlnK_dT(T)
        if dlnK is None:
            dlnK = (np.log(self.K(T + eps)) - np.log(self.K(T - eps))) / (2 * eps)
        return R_GAS * T**2 * dlnK

    def d_reaction_enthalpy_dT(self, T: float, eps: float = 1e-4) -> float:
        """
        Return d(ΔrH°)/dT [J/(mol·K)].

        FD fallback via central difference on reaction_enthalpy(T).
        Subclasses override with analytic forms where available.
        """
        return (self.reaction_enthalpy(T + eps) - self.reaction_enthalpy(T - eps)) / (
            2 * eps
        )


@dataclass
class EquilibriumConstant(EquilibriumConstantBase):
    """
    Temperature-independent equilibrium constant.

    Parameters
    ----------
    K_eq : float
        Dimensionless equilibrium constant.

    Notes
    -----
    For acid dissociation, use the pKa() factory:
        pKa(4.76)  ->  EquilibriumConstant(K_eq=10**-4.76)
    """

    K_eq: float

    def K(self, T: float) -> float:  # noqa: N802
        return self.K_eq

    def dlnK_dT(self, T: float) -> float:  # noqa: N802
        return 0.0

    def d_reaction_enthalpy_dT(self, T: float) -> float:
        return 0.0


@dataclass
class EquilibriumConstantVantHoff(EquilibriumConstantBase):
    """
    ln K(T) = -dH / (R·T) + dS / R.

    Parameters
    ----------
    dH : float    Standard enthalpy [J/mol].
    dS : float    Standard entropy [J/(mol·K)].
    T_ref : float Reference temperature [K] (documentation only).

    Notes
    -----
    For acid dissociation with known pKa and dH, use:
        pKa(value, dH=...)  ->  EquilibriumConstantVantHoff(...)
    """

    dH: float
    dS: float
    T_ref: float = 298.15

    def K(self, T: float) -> float:  # noqa: N802
        return float(np.exp(-self.dH / (R_GAS * T) + self.dS / R_GAS))

    def dlnK_dT(self, T: float) -> float:  # noqa: N802
        return self.dH / (R_GAS * T**2)

    def reaction_enthalpy(self, T: float) -> float:
        return self.dH

    def d_reaction_enthalpy_dT(self, T: float) -> float:
        return 0.0


@dataclass
class EquilibriumConstantVantHoffCp(EquilibriumConstantBase):
    """
    Van't Hoff with heat capacity correction (Kirchhoff's law):
        dH(T) = dH_ref + dCp · (T - T_ref)
        dS(T) = dS_ref + dCp · ln(T / T_ref).

    Parameters
    ----------
    dH : float    Standard enthalpy at T_ref [J/mol].
    dS : float    Standard entropy at T_ref [J/(mol·K)].
    dCp : float   Heat capacity difference [J/(mol·K)].
    T_ref : float Reference temperature [K].
    """

    dH: float
    dS: float
    dCp: float
    T_ref: float = 298.15

    def K(self, T: float) -> float:  # noqa: N802
        dH_T = self.dH + self.dCp * (T - self.T_ref)
        dS_T = self.dS + self.dCp * np.log(T / self.T_ref)
        return float(np.exp(-dH_T / (R_GAS * T) + dS_T / R_GAS))

    def dlnK_dT(self, T: float) -> float:  # noqa: N802
        dH_T = self.dH + self.dCp * (T - self.T_ref)
        return dH_T / (R_GAS * T**2)

    def reaction_enthalpy(self, T: float) -> float:
        return self.dH + self.dCp * (T - self.T_ref)

    def d_reaction_enthalpy_dT(self, T: float) -> float:
        return self.dCp


@dataclass
class EquilibriumConstantCustom(EquilibriumConstantBase):
    """K(T) from any callable — use for fitted polynomials, exponentials, or lookup tables.

    Parameters
    ----------
    func : callable
        Any callable ``(T: float) -> float`` returning the dimensionless K at T [K].
    """

    func: Callable[[float], float]

    def K(self, T: float) -> float:  # noqa: N802
        return float(self.func(T))


@dataclass
class EquilibriumConstantTabulated(EquilibriumConstantBase):
    """
    K(T) from linearly interpolated tabulated data.

    Parameters
    ----------
    T_data : array-like   Temperature values [K], monotonically increasing.
    K_data : array-like   Equilibrium constants at each temperature.
    """

    T_data: np.ndarray
    K_data: np.ndarray

    def __post_init__(self) -> None:
        self.T_data = np.asarray(self.T_data, dtype=float)
        self.K_data = np.asarray(self.K_data, dtype=float)

    def K(self, T: float) -> float:  # noqa: N802
        return float(np.interp(T, self.T_data, self.K_data))


@dataclass
class EquilibriumConstantPolynomial(EquilibriumConstantBase):
    """
    K(T) = exp(a₀ + a₁T + a₂T² + ...).

    coeffs = [a₀, a₁, a₂, ...] in ascending power order.

    Analytic derivatives:
        d(ln K)/dT = a₁ + 2a₂T + 3a₃T² + ...
        ΔrH°(T)    = R T² · d(ln K)/dT

    Parameters
    ----------
    coeffs : array-like
        Polynomial coefficients [a₀, a₁, ...], ascending power order.
    """

    coeffs: np.ndarray

    def __post_init__(self) -> None:
        self.coeffs = np.asarray(self.coeffs, dtype=float)
        self._powers = np.arange(len(self.coeffs), dtype=float)
        self._deriv_coeffs = self._powers * self.coeffs
        self._deriv2_coeffs = self._powers * (self._powers - 1) * self.coeffs

    def K(self, T: float) -> float:  # noqa: N802
        return float(np.exp(np.dot(self.coeffs, T**self._powers)))

    def dlnK_dT(self, T: float) -> float:  # noqa: N802
        powers_m1 = np.where(self._powers > 0, self._powers - 1, 0.0)
        return float(np.dot(self._deriv_coeffs, T**powers_m1))

    def reaction_enthalpy(self, T: float) -> float:
        return R_GAS * T**2 * self.dlnK_dT(T)

    def d_reaction_enthalpy_dT(self, T: float) -> float:
        dlnK = self.dlnK_dT(T)
        powers_m2 = np.where(self._powers > 1, self._powers - 2, 0.0)
        d2lnK = float(np.dot(self._deriv2_coeffs, T**powers_m2))
        return R_GAS * (2.0 * T * dlnK + T**2 * d2lnK)


def pKa(
    value: float,
    dH: Optional[float] = None,
    T_ref: float = 298.15,
) -> EquilibriumConstantBase:
    """
    Construct an equilibrium constant from a pKa value.

    Without dH: returns EquilibriumConstant(K_eq=10**-value).
    With dH:    returns EquilibriumConstantVantHoff with dS back-calculated
                from pKa and dH at T_ref.

    Parameters
    ----------
    value : float   pKa at T_ref (dimensionless).
    dH : float, optional
        Standard enthalpy of dissociation [J/mol].
    T_ref : float   Reference temperature [K].
    """
    Ka_ref = 10.0 ** (-value)
    if dH is None:
        return EquilibriumConstant(K_eq=Ka_ref)
    dS = (dH + R_GAS * T_ref * np.log(Ka_ref)) / T_ref
    return EquilibriumConstantVantHoff(dH=dH, dS=dS, T_ref=T_ref)
