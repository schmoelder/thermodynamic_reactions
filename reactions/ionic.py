"""
Ionic strength models: I(c, charges).
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np

from .species import PhysicalState

__all__ = [
    "IonicStrengthBase",
    "IonicStrengthIdeal",
    "IonicStrengthBackground",
    "IonicStrengthFixed",
]


class IonicStrengthBase(ABC):
    """Abstract base for ionic strength models."""

    @abstractmethod
    def evaluate(self, c: np.ndarray, charges: np.ndarray) -> float:
        """
        Return ionic strength [mol/m³].

        Parameters
        ----------
        c : np.ndarray
            Dynamic species concentrations [mol/m³].
        charges : np.ndarray
            Ionic charges, extracted from components by ReactionModel.
        """


@dataclass
class IonicStrengthIdeal(IonicStrengthBase):
    """I = 0.5 · sum(cᵢ · zᵢ²)  [mol/m³]"""

    def evaluate(self, c: np.ndarray, charges: np.ndarray) -> float:
        return 0.5 * float(np.sum(c * charges**2))


@dataclass
class IonicStrengthBackground(IonicStrengthBase):
    """
    I = I_bg + 0.5 · sum(cᵢ · zᵢ²)  [mol/m³]

    Parameters
    ----------
    I_bg : float
        Background ionic strength [mol/m³].
    """

    I_bg: float

    def evaluate(self, c: np.ndarray, charges: np.ndarray) -> float:
        return self.I_bg + 0.5 * float(np.sum(c * charges**2))


@dataclass
class IonicStrengthFixed(IonicStrengthBase):
    """
    I = constant  [mol/m³].

    Parameters
    ----------
    I : float
        Fixed ionic strength [mol/m³].
    """

    I: float

    def evaluate(self, c: np.ndarray, charges: np.ndarray) -> float:
        return self.I
