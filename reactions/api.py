"""
ReactionModel API
=================
Public-facing classes for defining thermodynamic and kinetic reaction models.

Units (SI throughout)
---------------------
Concentration  : mol / m³
Temperature    : K
Energy         : J / mol
Molar mass     : kg / mol
Density        : kg / m³

Equilibrium constants are dimensionless (activities, not concentrations).
Concentrations are divided by C_REF = 1000 mol/m³ when computing activities.

Formatting: ruff-compatible (line length 88, double quotes).
"""

from __future__ import annotations

import re
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Callable, Optional, Union

import numpy as np

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

R_GAS: float = 8.314462        # J / (mol K)
KB: float = 1.380649e-23       # J / K
H_PLANCK: float = 6.626070e-34 # J s
C_REF: float = 1000.0          # mol/m³  (standard state, 1 mol/L)


# ---------------------------------------------------------------------------
# Physical state — pure container, no logic
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
        Dynamic species concentrations [mol/m³], shape (n_species,).
        Solvent species (is_solvent=True) are excluded.
    T : float
        Temperature [K].
    I : float
        Ionic strength [mol/m³], derived by IonicStrengthModule.
    """

    c: np.ndarray
    T: float
    I: float = 0.0


# ---------------------------------------------------------------------------
# Conservation report — returned by ReactionModel.check_conservation()
# ---------------------------------------------------------------------------


@dataclass
class ConservationReport:
    """
    Result for one named Component from ReactionModel.check_conservation().

    Attributes
    ----------
    component : Component
        The component that was checked.
    conserved : bool
        True if the sum of the component's dynamic species is a conserved
        moiety of the reaction network (within tolerance).
    residual : float
        ||v - Q Q^T v|| / ||v||, where v is the component sum vector and
        Q spans the left null space of N.  Zero means exactly conserved.
    moiety_vector : np.ndarray
        The null-space vector most aligned with v, shape (n_dynamic_species,).
        Species ordering matches ReactionModel.species.  None if the null
        space is empty.
    """

    component: "Component"
    conserved: bool
    residual: float
    moiety_vector: Optional[np.ndarray]


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
    is_solvent : bool
        If True, activity is fixed at 1 and the species is excluded
        from the dynamic state vector. Default False.
    molar_mass : float, optional
        Molar mass [kg/mol].
    density : float, optional
        Pure-component density [kg/m³].
    """

    name: str
    charge: int = 0
    is_solvent: bool = False
    molar_mass: Optional[float] = None
    density: Optional[float] = None


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
    Component("H2O", is_solvent=True)

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


# ---------------------------------------------------------------------------
# Ionic strength  I(c, charges)
# Charges are extracted from components by ReactionModel and passed in.
# ---------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------
# Activity coefficient  γ(state, charges)
# ---------------------------------------------------------------------------


class ActivityCoefficientBase(ABC):
    """
    Abstract base for activity coefficient models.

    Returns γᵢ for each dynamic species.
    Thermodynamic activity: aᵢ = γᵢ · cᵢ / C_REF.
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
        [(m³/mol)^0.5]. At 25 °C: 0.509 (L/mol)^0.5 → 0.509/√1000.
    B : float
        [(m³/mol)^0.5 / m]. At 25 °C: 3.28e9/√1000.
    a_ion : float
        Mean ion-size parameter [m]. Typical: 3e-10 m.
    """

    A: float = 0.509 / (1000.0**0.5)
    B: float = 3.28e9 / (1000.0**0.5)
    a_ion: float = 3e-10

    def activity(self, state: PhysicalState, charges: np.ndarray) -> np.ndarray:
        sqrt_I = np.sqrt(state.I)
        log_gamma = (
            -self.A * charges**2 * sqrt_I
            / (1.0 + self.B * self.a_ion * sqrt_I)
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
        At 25 °C in water: 0.509 (L/mol)^0.5.
    """

    A: float = 0.509

    def activity(self, state: PhysicalState, charges: np.ndarray) -> np.ndarray:
        I_L = state.I / 1000.0
        sqrt_I = np.sqrt(I_L)
        log_gamma = -self.A * charges**2 * (sqrt_I / (1.0 + sqrt_I) - 0.3 * I_L)
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


# ---------------------------------------------------------------------------
# Equilibrium constant  K(T)
# ---------------------------------------------------------------------------


class EquilibriumConstantBase(ABC):
    """Abstract base for equilibrium constant models."""

    @abstractmethod
    def K(self, T: float) -> float:  # noqa: N802
        """Return dimensionless equilibrium constant at temperature T [K]."""


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


@dataclass
class EquilibriumConstantVantHoff(EquilibriumConstantBase):
    """
    ln K(T) = -dH / (R·T) + dS / R

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


@dataclass
class EquilibriumConstantVantHoffCp(EquilibriumConstantBase):
    """
    Van't Hoff with heat capacity correction (Kirchhoff's law):
        dH(T) = dH_ref + dCp · (T - T_ref)
        dS(T) = dS_ref + dCp · ln(T / T_ref)

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


# ---------------------------------------------------------------------------
# pKa factory — convenience constructor, not a class
# ---------------------------------------------------------------------------


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

    Examples
    --------
    pKa(4.76)                -> EquilibriumConstant(K_eq=1.738e-5)
    pKa(7.2, dH=+4000)       -> EquilibriumConstantVantHoff(...)
    """
    Ka_ref = 10.0 ** (-value)
    if dH is None:
        return EquilibriumConstant(K_eq=Ka_ref)
    dS = (dH + R_GAS * T_ref * np.log(Ka_ref)) / T_ref
    return EquilibriumConstantVantHoff(dH=dH, dS=dS, T_ref=T_ref)


# ---------------------------------------------------------------------------
# Rate constant  kf(T)
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

    For first-order reactions (n=1) the two conventions coincide.
    For higher-order reactions, the numerical value of kf differs by
    powers of C_REF = 1000 mol/m³:

        kf_thermo = kf_massaction * C_REF^(n-1)

    When fitting kf to experimental data, be consistent: use the same
    reaction class for fitting and prediction.  If you have literature
    values in mol/L units, multiply by 1000^(n-1) to convert to mol/m³.
    """

    @abstractmethod
    def kf(self, T: float) -> float:
        """Forward rate constant at temperature T [K]."""


@dataclass
class RateConstantFixed(RateConstantBase):
    """
    Temperature-independent forward rate constant.

    Parameters
    ----------
    kf_value : float
        Forward rate constant.
        Units depend on which reaction class this is used with and the
        reaction order — see RateConstantBase docstring for the convention.

    Examples
    --------
    First-order A -> B, kf = 1 /s (same in both MassAction and Thermodynamic):
        RateConstantFixed(kf_value=1.0)

    Second-order A + B -> C in MassActionReaction, kf = 1e-3 m³/(mol·s):
        RateConstantFixed(kf_value=1e-3)

    Same reaction in ThermodynamicReaction, kf = 1e-3 * C_REF = 1.0 mol/(m³·s):
        RateConstantFixed(kf_value=1.0)
    """

    kf_value: float

    def kf(self, T: float) -> float:
        return self.kf_value


@dataclass
class RateConstantArrhenius(RateConstantBase):
    """
    k(T) = A · exp(-Ea / (R·T))

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
    def __call__(
        self, state: PhysicalState, species_index: dict[str, int]
    ) -> float:
        """Net reaction rate [mol/(m³·s)]."""


@dataclass
class MichaelisMenten(RateBase):
    """
    v = Vmax · [S] / (Km + [S])

    Parameters
    ----------
    Vmax : float      Maximum rate [mol/(m³·s)].
    Km : float        Michaelis constant [mol/m³].
    substrate : str   Name of the substrate species.
    """

    Vmax: float
    Km: float
    substrate: str

    def __call__(
        self, state: PhysicalState, species_index: dict[str, int]
    ) -> float:
        S = state.c[species_index[self.substrate]]
        return self.Vmax * S / (self.Km + S)


@dataclass
class HillRate(RateBase):
    """
    v = Vmax · [S]^n / (Km^n + [S]^n)

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

    def __call__(
        self, state: PhysicalState, species_index: dict[str, int]
    ) -> float:
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

    fn: Callable[[PhysicalState, dict[str, int]], float]

    def __call__(
        self, state: PhysicalState, species_index: dict[str, int]
    ) -> float:
        return self.fn(state, species_index)


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
            s_pert = PhysicalState(c=c_pert, T=state.T, I=state.I)
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

        Parameters
        ----------
        species_index : dict[str, int]
        n_species : int
        exponent_fwd : np.ndarray, optional
            Override forward exponents, shape (n_species,).
        exponent_bwd : np.ndarray, optional
            Override backward exponents, shape (n_species,).

        Returns
        -------
        e_fwd, e_bwd : np.ndarray, shape (n_species,)
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
                continue   # solvent — skip
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
        """
        Net mass-action rate: v = kf * prod(c^e_fwd) - kr * prod(c^e_bwd).
        Zero or negative concentrations are clamped to avoid domain errors.
        """
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
        """
        Analytic derivative dv/dc_k for mass-action rate.
        dv/dc_k = kf * e_fwd[k]/c[k] * prod(c^e_fwd)
                - kr * e_bwd[k]/c[k] * prod(c^e_bwd)
        Returns shape (n_species,).
        """
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
    The net rate uses dimensionless activities a_i = gamma_i * c_i / C_REF:

        v = kf(T) · prod(a_i ^ e_fwd_i) - kr(T) · prod(a_i ^ e_bwd_i)

    Because activities are dimensionless, kf must have units of
    [mol/(m³·s)] for the rate to be in [mol/(m³·s)], regardless of
    reaction order.  This differs from MassActionReaction where kf
    operates on concentrations and its units depend on reaction order.

    For a reaction of order n (= sum of forward exponents):

        kf_ThermodynamicReaction = kf_MassActionReaction * C_REF^(n-1)

    where C_REF = 1000 mol/m³.  For first-order reactions (n=1) the
    two are identical.

    Relaxation timescale
    --------------------
    Because kf operates on activities, the relaxation timescale in
    concentration space is:

        tau = C_REF / (kf * (1 + 1/K))   [s]  for first-order A <-> B

    compared to tau = 1 / (kf + kr) for MassActionReaction.  If you
    specify kf = 1.0 mol/(m³·s), the system relaxes ~1000x slower in
    concentration units than a MassActionReaction with kf = 1.0 1/s.
    This is physically correct — it reflects the activity normalisation
    — but can surprise users coming from concentration-based frameworks.

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

    Examples
    --------
    Rapid equilibrium, pKa-specified, Davies activity:
        ThermodynamicReaction(
            "H2PO4- <-> HPO4-2 + H+",
            mode="equil",
            equilibrium_constant=pKa(7.2),
            activity_coefficient=ActivityCoefficientDavies(),
        )

    Kinetic, van't Hoff K(T), Arrhenius kf in [mol/(m³·s)]:
        ThermodynamicReaction(
            "A <-> B",
            mode="kinetic",
            equilibrium_constant=EquilibriumConstantVantHoff(dH=-20e3, dS=-50.0),
            rate_constant=RateConstantArrhenius(A=1e10, Ea=40e3),
        )

    Converting from a MassActionReaction kf (first-order, n=1 — no change needed):
        # MassActionReaction("A <-> B", kf=2.0, kr=0.5)
        # kf_thermo = kf_mass * C_REF^(1-1) = 2.0 * 1 = 2.0  mol/(m³·s)
        ThermodynamicReaction(
            "A <-> B",
            mode="kinetic",
            equilibrium_constant=EquilibriumConstant(K_eq=4.0),
            rate_constant=RateConstantFixed(kf_value=2.0),
        )
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
        """
        Net reaction rate [mol/(m³·s)] for kinetic mode.

        Uses activities a_i = gamma_i * c_i / C_REF instead of
        concentrations. Solvent species have activity = 1.
        Exponents default to stoichiometric convention (CADET eq. 25).
        """
        gamma = self.activity_coefficient.activity(state, charges)
        # Build activity array: a_i = gamma_i * c_i / C_REF
        # Solvent species (not in species_index) get activity 1.0
        n = len(state.c)
        a = np.ones(n)
        for i in range(n):
            a[i] = gamma[i] * state.c[i] / C_REF

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
        """
        Analytic dv/dc, shape (n_species,).
        Chain rule: dv/dc_k = dv/da_k * da_k/dc_k = dv/da_k * gamma_k/C_REF.
        """
        gamma = self.activity_coefficient.activity(state, charges)
        n = len(state.c)
        a = gamma * state.c / C_REF
        e_fwd, e_bwd = self._build_exponent_arrays(species_index, n)
        dv_da = self._mass_action_rate_jac(
            a,
            self.kf(state.T),
            self.kr(state.T),
            e_fwd,
            e_bwd,
        )
        # chain rule: dv/dc_k = dv/da_k * gamma_k / C_REF
        return dv_da * gamma / C_REF

    def log_K_residual(
        self,
        state: PhysicalState,
        species_index: dict[str, int],
        charges: np.ndarray,
    ) -> float:
        """
        Algebraic equilibrium residual: ln(Q) - ln(K).
        Q = prod(a_i^nu_i) over all species in this reaction.
        """
        gamma = self.activity_coefficient.activity(state, charges)

        def _a(name: str) -> float:
            if name not in species_index:
                return 1.0
            return float(gamma[species_index[name]] * state.c[species_index[name]] / C_REF)

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

    By default the exponent arrays are derived from stoichiometry
    (CADET convention, eq. 25):
        e_fwd[i] = max(0, -nu[i])   reactant exponents
        e_bwd[i] = max(0, +nu[i])   product exponents

    These can be overridden for non-elementary reactions or empirical
    power laws where reaction order differs from stoichiometry.

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

    Examples
    --------
    Elementary reversible:
        MassActionReaction("A <-> B", kf=1.0, kr=0.1)

    Irreversible second-order:
        MassActionReaction("2 A -> B", kf=0.5)

    Non-elementary (override exponents):
        MassActionReaction("A + B -> C", kf=1.0,
                           exponent_fwd=np.array([0.5, 1.0, 0.0]))
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
        # Store overrides — resolved to arrays in net_rate() once
        # species_index is known (at ReactionModel init time via cache)
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
        """
        Net mass-action rate [mol/(m³·s)] using concentrations directly.
        """
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

    Examples
    --------
    EnzymaticReaction("S -> P",
                      rate=MichaelisMenten(Vmax=1e-3, Km=50.0, substrate="S"))
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


# ---------------------------------------------------------------------------
# ReactionModel
# ---------------------------------------------------------------------------


class ReactionModel:
    """
    Self-contained reaction model for a single unit operation.

    Assembles PhysicalState from (c, T), dispatches to reaction modules,
    and exposes residuals / Jacobians for the DAE solver (CADET interface).

    Charges are extracted from components here and passed to the ionic
    strength module — modules themselves stay ignorant of component structure.

    Parameters
    ----------
    components : list[Component]
        System-level components. Referenced, not owned.
    reactions : list[ReactionBase]
    ionic_strength : IonicStrengthBase, optional
        Default: IonicStrengthIdeal().
    T : float
        Default temperature [K].

    Attributes
    ----------
    species : list[Species]
        Dynamic species only (is_solvent=False), in component order.
    species_index : dict[str, int]
        Maps species name to index in state.c.
    charges : np.ndarray
        Ionic charges for dynamic species, extracted from components.
    """

    def __init__(
        self,
        components: list[Component],
        reactions: list[ReactionBase],
        ionic_strength: Optional[IonicStrengthBase] = None,
        T: float = 298.15,
    ) -> None:
        self.components = components
        self.reactions = reactions
        self.ionic_strength = ionic_strength or IonicStrengthIdeal()
        self.T = T

        # All species — needed for stoichiometry validation
        self._all_species: list[Species] = [
            sp for comp in components for sp in comp.species
        ]
        # Dynamic species only — these form state.c
        self.species: list[Species] = [
            sp for sp in self._all_species if not sp.is_solvent
        ]
        self.species_index: dict[str, int] = {
            sp.name: i for i, sp in enumerate(self.species)
        }
        # Charges extracted here — ionic strength module stays dumb
        self.charges: np.ndarray = np.array(
            [sp.charge for sp in self.species], dtype=float
        )
        self._validate()
        # Cache stoichiometric matrix — built once, not per residual call
        self.nu, self.kinetic_mask, self.equil_dep = (
            self._build_stoichiometric_matrix()
        )

    def _validate(self) -> None:
        all_names = [sp.name for sp in self._all_species]
        if len(all_names) != len(set(all_names)):
            seen, dupes = set(), set()
            for n in all_names:
                (dupes if n in seen else seen).add(n)
            raise ValueError(f"Duplicate species names across components: {dupes}")
        known = set(all_names)
        for rxn in self.reactions:
            for name in rxn.species_names():
                if name not in known:
                    raise ValueError(
                        f"Species '{name}' in {rxn!r} not found in components."
                    )

    def _build_stoichiometric_matrix(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Build stoichiometric matrix and identify dependent species for
        equilibrium reactions.

        Returns
        -------
        nu : np.ndarray, shape (n_species, n_reactions)
            Stoichiometric coefficients. Rows = dynamic species,
            cols = reactions. Solvent species are excluded.
        kinetic_mask : np.ndarray of bool, shape (n_reactions,)
            True for kinetic reactions.
        equil_dep : np.ndarray of int, shape (n_equil_reactions,)
            Index of the dependent species for each equilibrium reaction
            (species with largest |nu_ij| — its ODE is replaced by the
            algebraic constraint).
        """
        n_s = len(self.species)
        n_r = len(self.reactions)
        nu = np.zeros((n_s, n_r))

        for j, rxn in enumerate(self.reactions):
            for name, coeff in rxn.nu.items():
                if name in self.species_index:
                    nu[self.species_index[name], j] = coeff
                # solvent species: excluded from state, skip

        kinetic_mask = np.array(
            [rxn.mode == "kinetic" for rxn in self.reactions]
        )

        # For each equilibrium reaction, pick the dependent species:
        # largest |nu_ij|; ties broken by preferring products (nu > 0).
        def _pick_dep(col: np.ndarray) -> int:
            abs_nu = np.abs(col)
            max_val = abs_nu.max()
            candidates = np.where(abs_nu == max_val)[0]
            # prefer product among candidates
            products = [i for i in candidates if col[i] > 0]
            return int(products[0] if products else candidates[0])

        equil_dep = np.array([
            _pick_dep(nu[:, j])
            for j in range(n_r)
            if not kinetic_mask[j]
        ], dtype=int)

        return nu, kinetic_mask, equil_dep

    def make_state(
        self, c: np.ndarray, T: Optional[float] = None
    ) -> PhysicalState:
        """
        Build a PhysicalState. Ionic strength computed from charges
        extracted from components — modules never see components directly.
        """
        T_val = T if T is not None else self.T
        I = self.ionic_strength.evaluate(c, self.charges)
        return PhysicalState(c=c, T=T_val, I=I)

    def residual(
        self,
        c: np.ndarray,
        c_dot: np.ndarray,
        T: Optional[float] = None,
    ) -> np.ndarray:
        """
        DAE residual vector, shape (n_dynamic_species,).

        For kinetic reactions (ODE rows):
            r_i = c_dot_i - sum_j nu_ij * v_j(c, T)

        For equilibrium reactions (algebraic rows, one per equil reaction):
            r_dep = ln(Q_j) - ln(K_j(T))

        where dep is the dependent species of reaction j — its ODE row
        is replaced by the algebraic constraint.

        Parameters
        ----------
        c : np.ndarray
            Dynamic species concentrations [mol/m³], shape (n_species,).
            This is the slice already provided by the system.
        c_dot : np.ndarray
            Time derivative dc/dt [mol/(m³·s)], shape (n_species,).
        T : float, optional
            Temperature [K]. Uses model default if not provided.
        """
        state = self.make_state(c, T)
        nu = self.nu
        kinetic_mask = self.kinetic_mask
        equil_dep = self.equil_dep
        n_r = len(self.reactions)

        # Start from ODE residual: r_i = c_dot_i - sum_j nu_ij * v_j
        # (will be overwritten for dependent species of equil reactions)
        r = c_dot.copy().astype(float)

        equil_counter = 0
        for j, rxn in enumerate(self.reactions):
            if kinetic_mask[j]:
                # kinetic: subtract nu_ij * v_j from all species
                v = rxn.net_rate(state, self.species_index, self.charges)
                r -= nu[:, j] * v
            else:
                # equil: replace ODE of dependent species with ln(Q)-ln(K)
                dep = equil_dep[equil_counter]
                equil_counter += 1
                r[dep] = rxn.log_K_residual(
                    state, self.species_index, self.charges
                )

        return r

    def jacobian(
        self,
        c: np.ndarray,
        c_dot: np.ndarray,
        T: Optional[float] = None,
        eps: float = 1e-6,
    ) -> np.ndarray:
        """
        Jacobian d(residual)/d(c), shape (n_species, n_species).

        Analytic for MassActionReaction and ThermodynamicReaction (kinetic).
        Analytic for equilibrium algebraic rows (ln Q derivative).
        Falls back to finite differences for EnzymaticReaction and CustomReaction.

        Analytic structure
        ------------------
        Kinetic row i, column k:
            dr_i/dc_k = -sum_j nu_ij * dv_j/dc_k

        For mass-action v_j = kf * prod(c^|nu|, reactants) - kr * prod(c^nu, products):
            dv_j/dc_k = kf * |nu_kj| / c_k * prod_reactants
                      - kr *  nu_kj  / c_k * prod_products
                      (only when nu_kj != 0)

        Equilibrium row dep_j, column k:
            dr_dep/dc_k = nu_kj / c_k   (derivative of ln Q w.r.t. c_k)

        The alpha * I term (d(residual)/d(c_dot)) is identity for kinetic rows,
        zero for equilibrium rows — handled by the DAE solver externally.
        """
        state = self.make_state(c, T)
        T_val = state.T
        nu = self.nu
        n_s = len(self.species)
        n_r = len(self.reactions)
        J = np.zeros((n_s, n_s))

        equil_counter = 0
        for j, rxn in enumerate(self.reactions):
            if self.kinetic_mask[j]:
                # --- analytic dv/dc for mass-action and thermodynamic ---
                if hasattr(rxn, "net_rate_jac"):
                    # analytic dv/dc available — use it
                    dv_dc = rxn.net_rate_jac(state, self.species_index, self.charges)
                    J -= np.outer(nu[:, j], dv_dc)
                else:
                    # fallback: finite differences for this reaction
                    v0 = rxn.net_rate(state, self.species_index, self.charges)
                    for k in range(n_s):
                        c_pert = c.copy()
                        c_pert[k] += eps
                        s_pert = PhysicalState(c=c_pert, T=state.T, I=state.I)
                        v_pert = rxn.net_rate(s_pert, self.species_index, self.charges)
                        J[:, k] -= nu[:, j] * (v_pert - v0) / eps

            else:
                # --- equilibrium row: dr_dep/dc_k = nu_kj / c_k ---
                dep = self.equil_dep[equil_counter]
                equil_counter += 1
                # zero out the ODE row for dep species (replaced by algebraic)
                J[dep, :] = 0.0
                for k in range(n_s):
                    sp_name = self.species[k].name
                    if sp_name not in rxn.nu:
                        continue
                    ck = float(c[k])
                    if ck == 0.0:
                        continue
                    # d(ln Q)/dc_k = nu_kj / (gamma_k * c_k / C_REF) * gamma_k/C_REF
                    #              = nu_kj / c_k  (gamma cancels)
                    J[dep, k] = rxn.nu[sp_name] / ck

        return J

    def check_conservation(
        self,
        tol: Optional[float] = None,
        report_all: bool = False,
    ) -> "list[ConservationReport] | tuple[list[ConservationReport], list[np.ndarray]]":
        """
        Check whether each multi-species Component corresponds to a conserved
        moiety of the reaction network.

        The left null space of the stoichiometric matrix N is computed via SVD.
        A vector v (1 for species in the component, 0 elsewhere) is conserved
        iff it lies in that null space.

        Parameters
        ----------
        tol : float, optional
            Singular-value threshold for null-space detection.
            None, 0, or negative uses machine precision:
            eps * max(singular_values).
        report_all : bool
            If True, also return all left null-space vectors as a second
            element.  Each vector has shape (n_dynamic_species,); species
            ordering matches self.species.

        Returns
        -------
        reports : list[ConservationReport]
            One entry per Component that has more than one dynamic species.
            Single-species components are skipped (trivially one variable).
        moieties : list[np.ndarray]
            Returned only when report_all=True.  All left null-space vectors
            of N.
        """
        N = self.nu  # shape (n_species, n_reactions)
        n_s = N.shape[0]

        U, s, _ = np.linalg.svd(N, full_matrices=True)

        if tol is None or tol <= 0:
            sv_tol = np.finfo(float).eps * (float(s.max()) if len(s) > 0 else 1.0)
        else:
            sv_tol = float(tol)

        # Columns of U with near-zero singular values span the left null space.
        # U has shape (n_s, n_s); s has length min(n_s, n_reactions).
        null_mask = np.concatenate(
            [s < sv_tol, np.ones(n_s - len(s), dtype=bool)]
        )
        Q = U[:, null_mask]  # shape (n_s, n_null)

        proj_tol = sv_tol ** 0.5  # projection residuals accumulate rounding error

        reports: list[ConservationReport] = []
        for comp in self.components:
            dyn = [
                sp for sp in comp.species
                if not sp.is_solvent and sp.name in self.species_index
            ]
            if len(dyn) <= 1:
                continue

            v = np.zeros(n_s)
            for sp in dyn:
                v[self.species_index[sp.name]] = 1.0
            v_norm = float(np.linalg.norm(v))

            if v_norm == 0.0 or Q.shape[1] == 0:
                reports.append(ConservationReport(
                    component=comp,
                    conserved=False,
                    residual=1.0,
                    moiety_vector=None,
                ))
                continue

            v_unit = v / v_norm
            proj = Q @ (Q.T @ v_unit)
            res = float(np.linalg.norm(v_unit - proj))
            dots = np.abs(Q.T @ v_unit)
            closest = Q[:, int(np.argmax(dots))] if dots.size > 0 else None
            reports.append(ConservationReport(
                component=comp,
                conserved=res < proj_tol,
                residual=res,
                moiety_vector=closest,
            ))

        if report_all:
            moieties = [Q[:, i] for i in range(Q.shape[1])]
            return reports, moieties
        return reports

    @property
    def parameters(self) -> dict[str, float]:
        """
        Flat dict of all fittable scalar float parameters.
        Suitable for scipy.optimize.
        """
        params: dict[str, float] = {}
        for i, rxn in enumerate(self.reactions):
            prefix = f"reactions[{i}]"
            for attr in (
                "equilibrium_constant",
                "rate_constant",
                "activity_coefficient",
                "_kf",
                "_kr",
                "rate",
            ):
                module = getattr(rxn, attr, None)
                if module is None:
                    continue
                for k, v in vars(module).items():
                    if isinstance(v, float):
                        params[f"{prefix}.{attr}.{k}"] = v
        return params

    def __repr__(self) -> str:
        return (
            f"ReactionModel("
            f"{len(self.components)} components, "
            f"{len(self.reactions)} reactions, "
            f"{len(self.species)} dynamic species)"
        )

