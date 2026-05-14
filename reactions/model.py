"""ReactionModel: assembles state, dispatches to reaction modules, exposes DAE interface."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from .ionic import IonicStrengthBase, IonicStrengthIdeal
from .reaction import ReactionBase
from .species import Component, Species
from .state import AuxiliaryState, State

__all__ = [
    "ConservationReport",
    "ReactionModel",
]


# ---------------------------------------------------------------------------
# Conservation report
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
        True if the sum of the component's species is a conserved
        moiety of the reaction network (within tolerance).
    residual : float
        ||v - Q Q^T v|| / ||v||, where v is the component sum vector and
        Q spans the left null space of N.  Zero means exactly conserved.
    moiety_vector : np.ndarray
        The null-space vector most aligned with v, shape (n_dynamic_species,).
        Species ordering matches ReactionModel.species.  None if the null
        space is empty.
    """

    component: Component
    conserved: bool
    residual: float
    moiety_vector: Optional[np.ndarray]


# ---------------------------------------------------------------------------
# ReactionModel
# ---------------------------------------------------------------------------


class ReactionModel:
    """
    Self-contained reaction model for a single unit operation.

    Assembles State and AuxiliaryState from (c, T), dispatches to reaction modules,
    and exposes residuals / Jacobians for the DAE solver (CADET interface).

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
        All species, in component order.
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

        self._all_species: list[Species] = [
            sp for comp in components for sp in comp.species
        ]
        self.species: list[Species] = self._all_species
        self.species_index: dict[str, int] = {
            sp.name: i for i, sp in enumerate(self.species)
        }
        self.charges: np.ndarray = np.array(
            [sp.charge for sp in self.species], dtype=float
        )
        self._validate()
        self.nu, self.kinetic_mask, self.equil_dep = self._build_stoichiometric_matrix()

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
        kinetic_mask : np.ndarray of bool, shape (n_reactions,)
        equil_dep : np.ndarray of int, shape (n_equil_reactions,)
            Index of the dependent species (largest |nu_ij|) for each
            equilibrium reaction.
        """
        n_s = len(self.species)
        n_r = len(self.reactions)
        nu = np.zeros((n_s, n_r))

        for j, rxn in enumerate(self.reactions):
            for name, coeff in rxn.nu.items():
                if name in self.species_index:
                    nu[self.species_index[name], j] = coeff

        kinetic_mask = np.array([rxn.mode == "kinetic" for rxn in self.reactions])

        def _pick_dep(col: np.ndarray) -> int:
            abs_nu = np.abs(col)
            eligible = [i for i in range(len(col)) if abs_nu[i] > 0]
            max_val = max(abs_nu[i] for i in eligible)
            candidates = [i for i in eligible if abs_nu[i] == max_val]
            products = [i for i in candidates if col[i] > 0]
            return int(products[0] if products else candidates[0])

        equil_dep = np.array(
            [_pick_dep(nu[:, j]) for j in range(n_r) if not kinetic_mask[j]], dtype=int
        )

        return nu, kinetic_mask, equil_dep

    def make_state(
        self,
        c: np.ndarray,
        T: Optional[float] = None,
    ) -> State:
        """Build a State from a concentration array."""
        T_val = T if T is not None else self.T
        species_names = [sp.name for sp in self.species]
        state = State("bulk", {"c": species_names}, T=T_val)
        state.c = c
        return state

    def make_aux(self, state: State) -> AuxiliaryState:
        """Compute AuxiliaryState (I, c_ref) from the current State."""
        I = self.ionic_strength.evaluate(state.c, self.charges)
        c_ref = np.array([sp.c_ref for sp in self.species])
        return AuxiliaryState(I=I, c_ref=c_ref)

    def residual(
        self,
        c: np.ndarray,
        c_dot: np.ndarray,
        T: Optional[float] = None,
    ) -> np.ndarray:
        """
        DAE residual vector, shape (n_species,).

        For kinetic reactions (ODE rows):
            r_i = c_dot_i - sum_j nu_ij * v_j(c, T)

        For equilibrium reactions (algebraic rows, one per equil reaction):
            r_dep = ln(Q_j) - ln(K_j(T))
        """
        state = self.make_state(c, T)
        aux = self.make_aux(state)
        nu = self.nu
        kinetic_mask = self.kinetic_mask
        equil_dep = self.equil_dep

        r = c_dot.copy().astype(float)

        equil_counter = 0
        for j, rxn in enumerate(self.reactions):
            if kinetic_mask[j]:
                v = rxn.net_rate(state, aux, self.species_index, self.charges)
                r -= nu[:, j] * v
            else:
                dep = equil_dep[equil_counter]
                equil_counter += 1
                r[dep] = rxn.log_K_residual(
                    state, aux, self.species_index, self.charges
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
        """
        state = self.make_state(c, T)
        aux = self.make_aux(state)
        nu = self.nu
        n_s = len(self.species)
        J = np.zeros((n_s, n_s))

        equil_counter = 0
        for j, rxn in enumerate(self.reactions):
            if self.kinetic_mask[j]:
                if hasattr(rxn, "net_rate_jac"):
                    dv_dc = rxn.net_rate_jac(
                        state, aux, self.species_index, self.charges
                    )
                    J -= np.outer(nu[:, j], dv_dc)
                else:
                    v0 = rxn.net_rate(state, aux, self.species_index, self.charges)
                    for k in range(n_s):
                        s_pert = state.copy()
                        s_pert.c[k] += eps
                        aux_pert = self.make_aux(s_pert)
                        v_pert = rxn.net_rate(
                            s_pert, aux_pert, self.species_index, self.charges
                        )
                        J[:, k] -= nu[:, j] * (v_pert - v0) / eps

            else:
                dep = self.equil_dep[equil_counter]
                equil_counter += 1
                J[dep, :] = 0.0
                for k in range(n_s):
                    sp_name = self.species[k].name
                    if sp_name not in rxn.nu:
                        continue
                    ck = float(c[k])
                    if ck == 0.0:
                        continue
                    J[dep, k] = rxn.nu[sp_name] / ck

        return J

    def volumetric_heat_capacity(
        self,
        c: np.ndarray,
    ) -> float:
        """
        Volumetric heat capacity ρCp [J/(m³·K)] from species with heat_capacity set.

        ρCp = Σ_k c_k · Cp_k
        where c_k [mol/m³] and Cp_k [J/(mol·K)].
        """
        thermal = [
            (i, sp) for i, sp in enumerate(self.species) if sp.heat_capacity is not None
        ]
        if not thermal:
            raise ValueError(
                "No species with molar_mass, density, and heat_capacity set. "
                "Energy balance requires at least one such species."
            )
        return sum(c[i] * sp.heat_capacity for i, sp in thermal)

    def jacobian_dT(
        self,
        c: np.ndarray,
        c_dot: np.ndarray,
        T: Optional[float] = None,
        eps: float = 1e-6,
    ) -> np.ndarray:
        """
        d(residual)/dT, shape (n_species,).

        Kinetic row i:   dr_i/dT = -sum_j nu_ij * dphi_j/dT
        Equilibrium dep: dr_dep/dT = -d(ln K_j)/dT

        Analytic where available; falls back to finite differences per reaction.
        """
        state = self.make_state(c, T)
        aux = self.make_aux(state)
        n_s = len(self.species)
        drdT = np.zeros(n_s)

        equil_counter = 0
        for j, rxn in enumerate(self.reactions):
            if self.kinetic_mask[j]:
                if hasattr(rxn, "net_rate_dT"):
                    dvdT = rxn.net_rate_dT(
                        state, aux, self.species_index, self.charges, eps
                    )
                else:
                    v0 = rxn.net_rate(state, aux, self.species_index, self.charges)
                    s_pert = state.with_T(state.T + eps)
                    dvdT = (
                        rxn.net_rate(s_pert, aux, self.species_index, self.charges) - v0
                    ) / eps
                drdT -= self.nu[:, j] * dvdT
            else:
                dep = self.equil_dep[equil_counter]
                equil_counter += 1
                dlnK = rxn.equilibrium_constant.dlnK_dT(state.T)
                if dlnK is None:
                    r0 = rxn.log_K_residual(
                        state, aux, self.species_index, self.charges
                    )
                    s_pert = state.with_T(state.T + eps)
                    r1 = rxn.log_K_residual(
                        s_pert, aux, self.species_index, self.charges
                    )
                    drdT[dep] = (r1 - r0) / eps
                else:
                    drdT[dep] = -dlnK

        return drdT

    def check_conservation(
        self,
        tol: Optional[float] = None,
        report_all: bool = False,
    ) -> "list[ConservationReport] | tuple[list[ConservationReport], list[np.ndarray]]":
        """
        Check whether each multi-species Component corresponds to a conserved
        moiety of the reaction network.

        Parameters
        ----------
        tol : float, optional
            Singular-value threshold for null-space detection.
            None, 0, or negative uses machine precision.
        report_all : bool
            If True, also return all left null-space vectors as a second element.

        Returns
        -------
        reports : list[ConservationReport]
            One entry per Component with more than one species.
        moieties : list[np.ndarray]
            Returned only when report_all=True.
        """
        N = self.nu
        n_s = N.shape[0]

        U, s, _ = np.linalg.svd(N, full_matrices=True)

        if tol is None or tol <= 0:
            sv_tol = np.finfo(float).eps * (float(s.max()) if len(s) > 0 else 1.0)
        else:
            sv_tol = float(tol)

        null_mask = np.concatenate([s < sv_tol, np.ones(n_s - len(s), dtype=bool)])
        Q = U[:, null_mask]

        proj_tol = sv_tol**0.5

        reports: list[ConservationReport] = []
        for comp in self.components:
            dyn = [sp for sp in comp.species if sp.name in self.species_index]
            if len(dyn) <= 1:
                continue

            v = np.zeros(n_s)
            for sp in dyn:
                v[self.species_index[sp.name]] = 1.0
            v_norm = float(np.linalg.norm(v))

            if v_norm == 0.0 or Q.shape[1] == 0:
                reports.append(
                    ConservationReport(
                        component=comp,
                        conserved=False,
                        residual=1.0,
                        moiety_vector=None,
                    )
                )
                continue

            v_unit = v / v_norm
            proj = Q @ (Q.T @ v_unit)
            res = float(np.linalg.norm(v_unit - proj))
            dots = np.abs(Q.T @ v_unit)
            closest = Q[:, int(np.argmax(dots))] if dots.size > 0 else None
            reports.append(
                ConservationReport(
                    component=comp,
                    conserved=res < proj_tol,
                    residual=res,
                    moiety_vector=closest,
                )
            )

        if report_all:
            moieties = [Q[:, i] for i in range(Q.shape[1])]
            return reports, moieties
        return reports

    @property
    def parameters(self) -> dict[str, float]:
        """Flat dict of all fittable scalar float parameters."""
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
