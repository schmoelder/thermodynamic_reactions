"""
Reaction classes: parse_stoichiometry, ReactionBase and subclasses.
"""

from __future__ import annotations

import re
from abc import ABC, abstractmethod
from typing import Optional, Union

import numpy as np

from .activity import ActivityCoefficientBase, ActivityCoefficientIdeal
from .equilibrium import EquilibriumConstantBase
from .rate import RateBase, RateConstantBase, RateConstantFixed
from .species import PhysicalState

__all__ = [
    "parse_stoichiometry",
    "ReactionBase",
    "ThermodynamicReaction",
    "MassActionReaction",
    "EnzymaticReaction",
    "CustomReaction",
]


# ---------------------------------------------------------------------------
# Stoichiometry parser
# ---------------------------------------------------------------------------


def parse_stoichiometry(stoichiometry_str: str) -> dict[str, float]:
    """
    Parse a reaction string into stoichiometric coefficients.

    Reactants are negative, products are positive.
    Splits on ' + ' (with surrounding spaces) to preserve charged
    species names like H+, OH-, HPO4-2.

    Parameters
    ----------
    stoichiometry_str : str
        E.g. "H2PO4- <-> HPO4-2 + H+"  or  "2 A + B -> C"

    Returns
    -------
    dict[str, float]
        {species_name: stoichiometric_coefficient}

    Examples
    --------
    >>> parse_stoichiometry("A + 2 B <-> C")
    {'A': -1.0, 'B': -2.0, 'C': 1.0}
    >>> parse_stoichiometry("H2O <-> H+ + OH-")
    {'H2O': -1.0, 'H+': 1.0, 'OH-': 1.0}
    >>> parse_stoichiometry("H2PO4- <-> HPO4-2 + H+")
    {'H2PO4-': -1.0, 'HPO4-2': 1.0, 'H+': 1.0}
    """
    reversible = "<->" in stoichiometry_str
    sep = "<->" if reversible else "->"
    parts = stoichiometry_str.split(sep)
    if len(parts) != 2:
        raise ValueError(
            f"Expected exactly one '{sep}' in '{stoichiometry_str}'."
        )

    def _parse_side(s: str) -> dict[str, float]:
        result: dict[str, float] = {}
        for term in re.split(r"\s+\+\s+", s.strip()):
            term = term.strip()
            if not term:
                continue
            m = re.match(
                r"^(\d+(?:\.\d+)?)?\s*([A-Za-z][A-Za-z0-9_]*[\-\+]?\d*)$",
                term,
            )
            if not m:
                raise ValueError(
                    f"Cannot parse stoichiometry term '{term}'."
                )
            coeff = float(m.group(1)) if m.group(1) else 1.0
            name = m.group(2)
            result[name] = result.get(name, 0.0) + coeff
        return result

    lhs = _parse_side(parts[0])
    rhs = _parse_side(parts[1])

    nu: dict[str, float] = {}
    for name, coeff in lhs.items():
        nu[name] = nu.get(name, 0.0) - coeff
    for name, coeff in rhs.items():
        nu[name] = nu.get(name, 0.0) + coeff
    return nu


# ---------------------------------------------------------------------------
# Reaction base
# ---------------------------------------------------------------------------


class ReactionBase(ABC):
    """
    Abstract base for all reaction types.

    Subclasses implement residual() for their specific rate law.
    jacobian() falls back to finite differences unless overridden.
    """

    mode: str  # "kinetic" | "equil", set by subclass

    def __init__(self, stoichiometry: str) -> None:
        self.stoichiometry_str = stoichiometry
        self.nu = parse_stoichiometry(stoichiometry)

    @abstractmethod
    def residual(
        self,
        state: PhysicalState,
        species_index: dict[str, int],
    ) -> np.ndarray:
        """
        Contribution to the residual vector at this state.

        Returns
        -------
        np.ndarray, shape (n_species,)
        """

    def jacobian(
        self,
        state: PhysicalState,
        species_index: dict[str, int],
    ) -> np.ndarray:
        """Finite-difference Jacobian. Override analytically if needed."""
        eps = 1e-6
        n = len(state.c)
        res0 = self.residual(state, species_index)
        J = np.zeros((len(res0), n))
        for i in range(n):
            c_pert = state.c.copy()
            c_pert[i] += eps
            s_pert = PhysicalState(c=c_pert, T=state.T, I=state.I, c_ref=state.c_ref)
            J[:, i] = (self.residual(s_pert, species_index) - res0) / eps
        return J

    def species_names(self) -> list[str]:
        return list(self.nu.keys())

    def _build_exponent_arrays(
        self,
        species_index: dict[str, int],
        n_species: int,
        exponent_fwd: Optional[np.ndarray] = None,
        exponent_bwd: Optional[np.ndarray] = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Build forward and backward exponent arrays from stoichiometry,
        or use the provided overrides.

        Default (stoichiometric convention, CADET eq. 25):
            e_fwd[i] = max(0, -nu[i])   reactant exponents
            e_bwd[i] = max(0, +nu[i])   product exponents
        """
        if exponent_fwd is not None and exponent_bwd is not None:
            return (
                np.asarray(exponent_fwd, dtype=float),
                np.asarray(exponent_bwd, dtype=float),
            )
        e_fwd = np.zeros(n_species)
        e_bwd = np.zeros(n_species)
        for name, coeff in self.nu.items():
            if name not in species_index:
                continue
            i = species_index[name]
            e_fwd[i] = max(0.0, -coeff)
            e_bwd[i] = max(0.0, +coeff)
        return e_fwd, e_bwd

    @staticmethod
    def _mass_action_rate(
        c: np.ndarray,
        kf: float,
        kr: float,
        e_fwd: np.ndarray,
        e_bwd: np.ndarray,
    ) -> float:
        """Net mass-action rate: v = kf * prod(c^e_fwd) - kr * prod(c^e_bwd)."""
        c_safe = np.maximum(c, 0.0)
        vf = kf * float(np.prod(c_safe ** e_fwd))
        vr = kr * float(np.prod(c_safe ** e_bwd))
        return vf - vr

    @staticmethod
    def _mass_action_rate_jac(
        c: np.ndarray,
        kf: float,
        kr: float,
        e_fwd: np.ndarray,
        e_bwd: np.ndarray,
    ) -> np.ndarray:
        """Analytic dv/dc_k for mass-action rate. Returns shape (n_species,)."""
        c_safe = np.maximum(c, 1e-300)
        prod_fwd = kf * float(np.prod(c_safe ** e_fwd))
        prod_bwd = kr * float(np.prod(c_safe ** e_bwd))
        dv = np.zeros(len(c))
        for k in range(len(c)):
            if c_safe[k] > 0:
                dv[k] = (
                    prod_fwd * e_fwd[k] / c_safe[k]
                    - prod_bwd * e_bwd[k] / c_safe[k]
                )
        return dv

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}('{self.stoichiometry_str}')"


# ---------------------------------------------------------------------------
# ThermodynamicReaction
# ---------------------------------------------------------------------------


class ThermodynamicReaction(ReactionBase):
    """
    Thermodynamically consistent reaction.

    Always requires an EquilibriumConstantBase (K(T)).
    For kinetic mode, also requires a RateConstantBase (kf(T)).
    kr is always derived as kf(T) / K(T) — never set independently.

    Rate law and units
    ------------------
    The net rate uses dimensionless activities a_i = gamma_i * c_i / c_ref_i:

        v = kf(T) · prod(a_i ^ e_fwd_i) - kr(T) · prod(a_i ^ e_bwd_i)

    Because activities are dimensionless, kf must have units of
    [mol/(m³·s)] for the rate to be in [mol/(m³·s)], regardless of
    reaction order.

    Parameters
    ----------
    stoichiometry : str
    mode : str
        "kinetic" — rate law drives dynamics (rate_constant required).
        "equil"   — algebraic K constraint, no explicit rate law.
    equilibrium_constant : EquilibriumConstantBase
        Provides K(T).
    rate_constant : RateConstantBase, optional
        Provides kf(T) in [mol/(m³·s)]. Required for mode="kinetic".
    activity_coefficient : ActivityCoefficientBase, optional
        Default: ActivityCoefficientIdeal().
    """

    def __init__(
        self,
        stoichiometry: str,
        mode: str,
        equilibrium_constant: EquilibriumConstantBase,
        rate_constant: Optional[RateConstantBase] = None,
        activity_coefficient: Optional[ActivityCoefficientBase] = None,
    ) -> None:
        super().__init__(stoichiometry)
        if mode not in ("kinetic", "equil"):
            raise ValueError("mode must be 'kinetic' or 'equil'.")
        if mode == "kinetic" and rate_constant is None:
            raise ValueError("rate_constant required for mode='kinetic'.")
        self.mode = mode
        self.equilibrium_constant = equilibrium_constant
        self.rate_constant = rate_constant
        self.activity_coefficient = activity_coefficient or ActivityCoefficientIdeal()

    def kf(self, T: float) -> float:
        if self.rate_constant is None:
            raise ValueError("No rate_constant — reaction is equil only.")
        return self.rate_constant.kf(T)

    def kr(self, T: float) -> float:
        """Reverse rate constant — always derived from kf / K(T)."""
        return self.kf(T) / self.equilibrium_constant.K(T)

    def net_rate(
        self,
        state: PhysicalState,
        species_index: dict[str, int],
        charges: np.ndarray,
    ) -> float:
        """Net reaction rate [mol/(m³·s)] for kinetic mode."""
        gamma = self.activity_coefficient.activity(state, charges)
        n = len(state.c)
        a = gamma * state.c / state.c_ref

        e_fwd, e_bwd = self._build_exponent_arrays(species_index, n)
        return self._mass_action_rate(
            a,
            self.kf(state.T),
            self.kr(state.T),
            e_fwd,
            e_bwd,
        )

    def net_rate_jac(
        self,
        state: PhysicalState,
        species_index: dict[str, int],
        charges: np.ndarray,
    ) -> np.ndarray:
        """Analytic dv/dc, shape (n_species,)."""
        gamma = self.activity_coefficient.activity(state, charges)
        n = len(state.c)
        a = gamma * state.c / state.c_ref
        e_fwd, e_bwd = self._build_exponent_arrays(species_index, n)
        dv_da = self._mass_action_rate_jac(
            a,
            self.kf(state.T),
            self.kr(state.T),
            e_fwd,
            e_bwd,
        )
        return dv_da * gamma / state.c_ref

    def net_rate_dT(
        self,
        state: PhysicalState,
        species_index: dict[str, int],
        charges: np.ndarray,
        eps: float = 1e-6,
    ) -> float:
        """
        d(phi_j)/dT using the analytic formula:

            dphi/dT = (d ln kf / dT) * phi  +  kr * (d ln K / dT) * P_bwd

        Falls back to finite differences when either derivative is unavailable.
        """
        dlnkf = self.rate_constant.dlnkf_dT(state.T)
        dlnK = self.equilibrium_constant.dlnK_dT(state.T)

        if dlnkf is None or dlnK is None:
            phi0 = self.net_rate(state, species_index, charges)
            s_pert = PhysicalState(c=state.c, T=state.T + eps, I=state.I, c_ref=state.c_ref)
            return (self.net_rate(s_pert, species_index, charges) - phi0) / eps

        gamma = self.activity_coefficient.activity(state, charges)
        n = len(state.c)
        a = gamma * state.c / state.c_ref
        e_fwd, e_bwd = self._build_exponent_arrays(species_index, n)
        a_safe = np.maximum(a, 0.0)
        P_bwd = float(np.prod(a_safe ** e_bwd))
        phi = self._mass_action_rate(a, self.kf(state.T), self.kr(state.T), e_fwd, e_bwd)
        return dlnkf * phi + self.kr(state.T) * dlnK * P_bwd

    def log_K_residual(
        self,
        state: PhysicalState,
        species_index: dict[str, int],
        charges: np.ndarray,
    ) -> float:
        """Algebraic equilibrium residual: ln(Q) - ln(K)."""
        gamma = self.activity_coefficient.activity(state, charges)

        def _a(name: str) -> float:
            if name not in species_index:
                return 1.0
            i = species_index[name]
            return float(gamma[i] * state.c[i] / state.c_ref[i])

        ln_Q = sum(
            coeff * np.log(max(_a(name), 1e-300))
            for name, coeff in self.nu.items()
        )
        return float(ln_Q - np.log(self.equilibrium_constant.K(state.T)))

    def residual(
        self, state: PhysicalState, species_index: dict[str, int]
    ) -> np.ndarray:
        raise NotImplementedError(
            "ReactionModel assembles residuals — do not call residual() directly."
        )


# ---------------------------------------------------------------------------
# MassActionReaction
# ---------------------------------------------------------------------------


class MassActionReaction(ReactionBase):
    """
    Mass-action kinetics. kf and kr are independent — no thermodynamic
    consistency enforced.

    The net rate is:
        v = kf · prod(c_i ^ e_fwd_i) - kr · prod(c_i ^ e_bwd_i)

    Use ThermodynamicReaction when K(T) has physical meaning.

    Parameters
    ----------
    stoichiometry : str
    kf : float or RateConstantBase
    kr : float or RateConstantBase
        Use 0.0 for irreversible reactions.
    exponent_fwd : np.ndarray, optional
        Forward exponents, shape (n_species,). Defaults to stoichiometry.
    exponent_bwd : np.ndarray, optional
        Backward exponents, shape (n_species,). Defaults to stoichiometry.
    """

    mode = "kinetic"

    def __init__(
        self,
        stoichiometry: str,
        kf: Union[float, RateConstantBase] = 1.0,
        kr: Union[float, RateConstantBase] = 0.0,
        exponent_fwd: Optional[np.ndarray] = None,
        exponent_bwd: Optional[np.ndarray] = None,
    ) -> None:
        super().__init__(stoichiometry)
        self._kf: RateConstantBase = (
            RateConstantFixed(kf) if isinstance(kf, (int, float)) else kf
        )
        self._kr: RateConstantBase = (
            RateConstantFixed(kr) if isinstance(kr, (int, float)) else kr
        )
        self._exponent_fwd_override = (
            np.asarray(exponent_fwd, dtype=float)
            if exponent_fwd is not None else None
        )
        self._exponent_bwd_override = (
            np.asarray(exponent_bwd, dtype=float)
            if exponent_bwd is not None else None
        )

    def kf(self, T: float) -> float:
        return self._kf.kf(T)

    def kr(self, T: float) -> float:
        return self._kr.kf(T)

    def net_rate(
        self,
        state: PhysicalState,
        species_index: dict[str, int],
        charges: np.ndarray,
    ) -> float:
        """Net mass-action rate [mol/(m³·s)] using concentrations directly."""
        n = len(state.c)
        e_fwd, e_bwd = self._build_exponent_arrays(
            species_index, n,
            self._exponent_fwd_override,
            self._exponent_bwd_override,
        )
        return self._mass_action_rate(
            state.c,
            self.kf(state.T),
            self.kr(state.T),
            e_fwd,
            e_bwd,
        )

    def net_rate_jac(
        self,
        state: PhysicalState,
        species_index: dict[str, int],
        charges: np.ndarray,
    ) -> np.ndarray:
        """Analytic dv/dc, shape (n_species,)."""
        n = len(state.c)
        e_fwd, e_bwd = self._build_exponent_arrays(
            species_index, n,
            self._exponent_fwd_override,
            self._exponent_bwd_override,
        )
        return self._mass_action_rate_jac(
            state.c,
            self.kf(state.T),
            self.kr(state.T),
            e_fwd,
            e_bwd,
        )

    def residual(
        self, state: PhysicalState, species_index: dict[str, int]
    ) -> np.ndarray:
        raise NotImplementedError(
            "ReactionModel assembles residuals — do not call residual() directly."
        )


# ---------------------------------------------------------------------------
# EnzymaticReaction
# ---------------------------------------------------------------------------


class EnzymaticReaction(ReactionBase):
    """
    Saturation kinetics. Concrete rate law supplied via a RateBase instance.

    Parameters
    ----------
    stoichiometry : str
    rate : RateBase
        E.g. MichaelisMenten(...) or HillRate(...).
    """

    mode = "kinetic"

    def __init__(self, stoichiometry: str, rate: RateBase) -> None:
        super().__init__(stoichiometry)
        self.rate = rate

    def net_rate(
        self,
        state: PhysicalState,
        species_index: dict[str, int],
        charges: np.ndarray,
    ) -> float:
        """Net rate from the supplied RateBase instance."""
        return self.rate(state, species_index)

    def residual(
        self, state: PhysicalState, species_index: dict[str, int]
    ) -> np.ndarray:
        raise NotImplementedError(
            "ReactionModel assembles residuals — do not call residual() directly."
        )


# ---------------------------------------------------------------------------
# CustomReaction
# ---------------------------------------------------------------------------


class CustomReaction(ReactionBase):
    """
    User-supplied rate wrapped in the reaction interface.

    Parameters
    ----------
    stoichiometry : str
    rate : RateBase
    mode : str   "kinetic" | "equil"
    """

    def __init__(
        self,
        stoichiometry: str,
        rate: RateBase,
        mode: str = "kinetic",
    ) -> None:
        super().__init__(stoichiometry)
        if mode not in ("kinetic", "equil"):
            raise ValueError("mode must be 'kinetic' or 'equil'.")
        self.mode = mode
        self.rate = rate

    def residual(
        self, state: PhysicalState, species_index: dict[str, int]
    ) -> np.ndarray:
        raise NotImplementedError(
            "CustomReaction.residual() — to be implemented."
        )
