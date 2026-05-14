"""Analytical post-processing of ReactionModel equilibria."""

from __future__ import annotations

import numpy as np

__all__ = ["buffer_capacity", "speciation_fractions", "solve_equilibrium_sweep"]


def buffer_capacity(
    model,
    c_total: dict[str, float],
    pH,
    T: float = 298.15,
) -> dict[str, np.ndarray]:
    """
    Analytical buffer capacity β(pH) for a ReactionModel.

    Uses ideal Bjerrum fractions derived from the model's equilibrium
    constants; activity corrections are not applied (ideal β).
    Proton-transfer reactions are identified by stoichiometry: any
    reaction where H⁺ appears as a product (ν > 0).

    Parameters
    ----------
    model : ReactionModel
    c_total : dict[str, float]
        Total concentration of each component [mol/m³].
        Components absent from the dict contribute zero.
    pH : array-like
        pH values at which to evaluate β.
    T : float
        Temperature [K], used to evaluate K(T).

    Returns
    -------
    dict[str, np.ndarray]
        β [mol/(m³·pH)] for each multi-species component that has a
        complete proton-transfer ladder, plus a ``"water"`` key for the
        water/autoionization contribution.

    Notes
    -----
    Fixes two issues present in the CADET-Process implementation:

    * Species are sorted by charge (descending) rather than relying on
      reaction iteration order, so the pKa ladder is always assembled
      correctly regardless of how reactions were defined.
    * Charge differences use actual ``Species.charge`` values rather than
      index differences, so non-unit charge steps are handled correctly.
    """
    pH = np.asarray(pH, dtype=float)

    h_name, h_cref = _find_proton(model)
    c_H = 10.0 ** (-pH) * h_cref  # mol/m³

    result: dict[str, np.ndarray] = {}
    result["water"] = _beta_water(c_H, model, h_name, h_cref, T)

    for component in model.components:
        if len(component.species) < 2:
            continue

        species = sorted(component.species, key=lambda s: -s.charge)
        charges = np.array([s.charge for s in species], dtype=float)

        K_ladder = _k_ladder(model, species, h_name, T)
        if K_ladder is None:
            continue

        c_tot = c_total.get(component.name, 0.0)
        if c_tot == 0.0:
            continue

        alpha = _bjerrum(K_ladder, c_H, h_cref)

        n = len(species)
        beta = np.zeros_like(c_H)
        for j in range(n):
            for i in range(j):
                dz = charges[i] - charges[j]
                beta += dz**2 * alpha[i] * alpha[j]
        result[component.name] = np.log(10) * c_tot * beta

    return result


def speciation_fractions(
    pH,
    pKa_list: list[float],
) -> np.ndarray:
    """
    Analytical Bjerrum speciation fractions for a polyprotic acid.

    Parameters
    ----------
    pH : array-like
        pH values at which to evaluate the fractions.
    pKa_list : list[float]
        pKa values in dissociation order (Ka1, Ka2, …).
        Length n gives n+1 species: fully protonated (index 0)
        through fully deprotonated (index n).

    Returns
    -------
    np.ndarray, shape (n_species, len(pH))
        ``fractions[k]`` is the mole fraction of the species that has
        donated k protons relative to the fully protonated form.

    Examples
    --------
    Monoprotic (acetic acid, pKa = 4.756):

    >>> f = speciation_fractions([4.756], [4.756])
    >>> f[:, 0]  # at pH = pKa both fractions are 0.5
    array([0.5, 0.5])

    Phosphate (three pKa values):

    >>> f = speciation_fractions(np.linspace(0, 14, 500), [2.148, 7.198, 12.350])
    >>> f.shape
    (4, 500)
    """
    pH = np.atleast_1d(np.asarray(pH, dtype=float))
    n = len(pKa_list)
    h = 10.0 ** (-pH)
    Ka = [10.0 ** (-pk) for pk in pKa_list]

    weights = np.zeros((n + 1, len(h)))
    K_prod = 1.0
    for k in range(n + 1):
        weights[k] = K_prod * h ** (n - k)
        if k < n:
            K_prod *= Ka[k]

    D = weights.sum(axis=0)
    return weights / D


def solve_equilibrium_sweep(
    model,
    pH,
    c0: dict[str, float],
    prescribed: dict | None = None,
    T: float = 298.15,
) -> dict[str, np.ndarray]:
    """
    Solve equilibrium across a pH range by prescribing [H⁺] at each point.

    At each step [H⁺] = 10^(−pH) · c_ref is added to ``prescribed`` and
    passed to ``solve_equilibrium``.  The result of each solve is used as
    the warm-start for the next point, so only the first pH point requires
    a reasonable ``c0``; the rest follow automatically.

    Parameters
    ----------
    model : ReactionModel
    pH : array-like
        pH values to sweep over.
    c0 : dict[str, float]
        Initial concentration guess for the first pH point [mol/m³].
    prescribed : dict, optional
        Additional species to hold fixed throughout the sweep (e.g.
        ``{"H2O": water.c_ref}``).  H⁺ is always added automatically.
    T : float
        Temperature [K].

    Returns
    -------
    dict[str, np.ndarray]
        Equilibrium concentrations [mol/m³] for every species, each an
        array of length ``len(pH)``.  Points where the solver did not
        converge are filled with ``np.nan``.

    Notes
    -----
    Prescribing [H⁺] is equivalent to the pH-stat pattern
    (@implementation-practical): speciation adjusts freely while pH is
    held fixed.  This is more robust than providing Henderson–Hasselbalch
    initial guesses at every pH because the warm-start propagates the
    previous equilibrium state across the sweep.
    """
    from .solver import solve_equilibrium

    pH = np.asarray(pH, dtype=float)
    h_name, h_cref = _find_proton(model)
    species_names = [sp.name for sp in model.species]

    result = {name: np.full(len(pH), np.nan) for name in species_names}
    c_curr = dict(c0)

    for i, ph in enumerate(pH):
        c_H = 10.0 ** (-ph) * h_cref
        c_curr[h_name] = c_H
        step_prescribed = {h_name: c_H, **(prescribed or {})}
        try:
            c_eq = solve_equilibrium(model, c_curr, T=T, prescribed=step_prescribed)
            for name in species_names:
                result[name][i] = c_eq.get(name, np.nan)
            c_curr = dict(c_eq)
        except RuntimeError:
            pass  # NaN already filled; warm-start from last successful result

    return result


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _find_proton(model) -> tuple[str, float]:
    """Return (species_name, c_ref) for H⁺ (charge +1, single-species component)."""
    for comp in model.components:
        if len(comp.species) == 1 and comp.species[0].charge == 1:
            sp = comp.species[0]
            return sp.name, sp.c_ref
    raise ValueError(
        "No H⁺ species found. Add a single-species component with charge=+1."
    )


def _beta_water(
    c_H: np.ndarray,
    model,
    h_name: str,
    h_cref: float,
    T: float,
) -> np.ndarray:
    """
    β_water = ln(10) · ([H⁺] + [OH⁻]).

    Kw is taken from the autoionization reaction if present (enabling
    van't Hoff temperature dependence); otherwise falls back to
    Kw = 10⁻¹⁴ at 298.15 K, converted to mol²/m⁶.
    """
    kw = _find_kw(model, h_name, h_cref, T)
    c_OH = kw / c_H
    return np.log(10) * (c_H + c_OH)


def _find_kw(model, h_name: str, h_cref: float, T: float) -> float:
    """
    Return Kw in mol²/m⁶ (= c_H · c_OH at equilibrium).

    Checks two patterns (tried in order):

    1. **Empty-LHS** (preferred): no reactants, products are H⁺ plus exactly
       one charge-(-1) species.  This is the ``<-> H+ + OH-`` stoichiometry
       produced by ``autoionization()``.

    2. **Explicit-water** (mixed-solvent): sole reactant is a solvent species
       (density set), same product pattern.  Discriminating on density
       distinguishes H₂O → H⁺ + OH⁻ from acid dissociation reactions such as
       H₃PO₄ → H₂PO₄⁻ + H⁺, where H₂PO₄⁻ also carries charge −1.

    Falls back to 1e-14 · c_ref² (25 °C value in pure water).
    """
    for rxn in model.reactions:
        reactants = {k: v for k, v in rxn.nu.items() if v < 0}
        products = {k: v for k, v in rxn.nu.items() if v > 0}
        if h_name not in products:
            continue
        other_products = {k: v for k, v in products.items() if k != h_name}
        if len(other_products) != 1:
            continue
        other_name = next(iter(other_products))
        if _species_charge(model, other_name) != -1:
            continue
        # Pattern 1: empty-LHS (autoionization() default)
        if len(reactants) == 0:
            return rxn.equilibrium_constant.K(T) * h_cref**2
        # Pattern 2: explicit water reactant (mixed-solvent models)
        if len(reactants) == 1:
            reactant_name = next(iter(reactants))
            if _species_is_solvent(model, reactant_name):
                return rxn.equilibrium_constant.K(T) * h_cref**2
    return 1e-14 * h_cref**2


def _k_ladder(model, species_sorted: list, h_name: str, T: float) -> list | None:
    """
    Build the K ladder for adjacent species pairs.

    For each pair (species[i], species[i+1]), find the reaction
    species[i] <-> species[i+1] + H⁺ and return K(T).
    Returns None if any step is missing from the model.
    """
    K_ladder = []
    for i in range(len(species_sorted) - 1):
        sp_from = species_sorted[i].name
        sp_to = species_sorted[i + 1].name
        K = None
        for rxn in model.reactions:
            nu = rxn.nu
            if (
                nu.get(sp_from, 0) < 0
                and nu.get(sp_to, 0) > 0
                and nu.get(h_name, 0) > 0
            ):
                K = rxn.equilibrium_constant.K(T)
                break
        if K is None:
            return None
        K_ladder.append(K)
    return K_ladder


def _bjerrum(K_ladder: list[float], c_H: np.ndarray, h_cref: float) -> np.ndarray:
    """
    Ideal Bjerrum fractions α[k, :] for n = len(K_ladder)+1 species.

    K_ladder[i] is the dimensionless equilibrium constant for
    species[i] <-> species[i+1] + H⁺, so K_eff[i] = K_ladder[i] * h_cref
    has units mol/m³ and gives the concentration-space Ka.
    """
    n = len(K_ladder) + 1
    K_eff = np.array(K_ladder) * h_cref  # mol/m³

    weights = np.zeros((n, len(c_H)))
    K_prod = 1.0
    for k in range(n):
        weights[k] = K_prod * c_H ** (n - 1 - k)
        if k < n - 1:
            K_prod *= K_eff[k]

    D = weights.sum(axis=0)
    return weights / D


def _species_charge(model, name: str) -> int | None:
    for comp in model.components:
        for sp in comp.species:
            if sp.name == name:
                return sp.charge
    return None


def _species_is_solvent(model, name: str) -> bool:
    """Return True if the named species has density set (solvent, not solute)."""
    for comp in model.components:
        for sp in comp.species:
            if sp.name == name:
                return sp.density is not None
    return False
