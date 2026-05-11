"""
solver.py
=========
ODE/DAE integration and equilibrium solving for ReactionModel.

Two entry points:
    simulate()          — time integration via scipy Radau (stiff ODE)
    solve_equilibrium() — steady-state via Newton in log space

Units: SI throughout (mol/m³, K, s).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional, Union

import numpy as np
from scipy.integrate import solve_ivp

from .model import ReactionModel

__all__ = [
    "SimulationResult",
    "simulate",
    "solve_equilibrium",
]


# ---------------------------------------------------------------------------
# Simulation result
# ---------------------------------------------------------------------------


@dataclass
class SimulationResult:
    """
    Result of a ReactionModel simulation.

    Attributes
    ----------
    t : np.ndarray
        Time points [s].
    c : np.ndarray
        Concentrations [mol/m³], shape (n_t, n_species).
    species : list[str]
        Species names, matching columns of c.
    success : bool
    message : str
    """

    t: np.ndarray
    c: np.ndarray
    species: list[str]
    success: bool
    message: str
    T_profile: Optional[np.ndarray] = None  # None when T was scalar/default

    def __getitem__(self, name: str) -> np.ndarray:
        """result['A'] returns the concentration array of species A."""
        return self.c[:, self.species.index(name)]

    def plot(self, ax=None, species: Optional[list[str]] = None, **kwargs):
        """
        Plot concentrations vs time.

        Parameters
        ----------
        ax : matplotlib Axes, optional
        species : list[str], optional
            Subset of species to plot. Default: all.
        """
        import matplotlib.pyplot as plt
        if ax is None:
            _, ax = plt.subplots()
        names = species or self.species
        for name in names:
            ax.plot(self.t, self[name], label=name, **kwargs)
        ax.set_xlabel("time [s]")
        ax.set_ylabel("concentration [mol/m³]")
        ax.legend()
        return ax


# ---------------------------------------------------------------------------
# simulate
# ---------------------------------------------------------------------------


def simulate(
    model: ReactionModel,
    c0: dict[str, float],
    t_span: tuple[float, float],
    T: Union[float, Callable[[float], float], None] = None,
    solvent_composition: Optional[dict] = None,
    prescribed: Optional[dict] = None,
    n_points: int = 500,
    rtol: float = 1e-8,
    atol: float = 1e-10,
    method: str = "radau",
) -> SimulationResult:
    """
    Integrate a ReactionModel forward in time.

    Parameters
    ----------
    model : ReactionModel
    c0 : dict[str, float]
        Initial concentrations {species_name: value [mol/m³]}.
        Missing species default to 0.
    t_span : tuple[float, float]
        (t_start, t_end) [s].
    T : float or callable or None
        Temperature [K].  Pass a float for isothermal simulations, or a
        callable T(t) -> float for a prescribed temperature programme.
        Uses model default if not provided.
        Must be a float (initial condition) when solvent_composition is given.
    solvent_composition : dict, optional
        Solvent mole fractions {species_name: value}, where value is either
        a float (fixed) or a callable x_k(t) -> float (gradient programme).
        When provided, T is treated as a dynamic state solved from the
        energy balance:

            ρCp(t) · dT/dt = -∑_j ΔH_j(T) · φ_j(c, T)

        ρCp is computed at each step via model.volumetric_heat_capacity().
        Cannot be combined with a callable T.
    prescribed : dict, optional
        Species held at a fixed or time-varying concentration
        {species_name: float | callable(t) -> float}.
        Prescribed species are excluded from the ODE and held at their
        specified value throughout the integration.
        Intended for parametric sweeps (e.g. fixing c_H+ to scan pH).
        Solvent concentrations are handled by the ``Solution`` class,
        which passes the appropriate entries here.
    n_points : int
        Number of output time points.
    rtol, atol : float
        Relative and absolute tolerances.
    method : {'radau', 'ida'}
        ODE/DAE solver backend.  'radau' uses scipy Radau (default);
        'ida' uses SUNDIALS IDA via scikits.odes.  The IDA backend treats
        prescribed species and equilibrium constraints as true algebraic
        rows rather than using injection.  Requires scikits.odes to be
        installed.  Not yet supported with solvent_composition.

    Returns
    -------
    SimulationResult
        result.T_profile is populated when T is callable or
        solvent_composition is given; None otherwise.
    """
    species_names = [sp.name for sp in model.species]
    c_init = np.array([c0.get(name, 0.0) for name in species_names])
    n = len(c_init)
    t_eval = np.linspace(t_span[0], t_span[1], n_points)

    # --- prescribed species: inject c[i] = fn(t) before every RHS call ---
    _prescribed = prescribed or {}
    _prescribed_fns: dict[str, Callable[[float], float]] = {}
    for _pname, _pval in _prescribed.items():
        if callable(_pval):
            _prescribed_fns[_pname] = _pval
        else:
            _v = float(_pval)
            _prescribed_fns[_pname] = lambda t, v=_v: v
    prescribed_mask = np.array([sp.name in _prescribed_fns for sp in model.species])
    has_prescribed = bool(prescribed_mask.any())
    dynamic_idx = np.where(~prescribed_mask)[0]
    prescribed_idx = np.where(prescribed_mask)[0]
    prescribed_fns = [_prescribed_fns[model.species[i].name] for i in prescribed_idx]
    if has_prescribed:
        # Initialise prescribed species from their function at t0
        for i, fn in zip(prescribed_idx, prescribed_fns):
            c_init[i] = fn(t_span[0])

    # --- IDA path (scikits.odes) ---
    if method == "ida":
        try:
            from scikits.odes.dae import dae as ida_dae
        except ImportError:
            raise ImportError(
                "scikits.odes is required for method='ida'. "
                "Install with: pip install scikits.odes"
            )
        if solvent_composition is not None:
            raise NotImplementedError(
                "method='ida' is not yet supported with solvent_composition."
            )

        T_func_ida: Callable[[float], Optional[float]]
        if callable(T):
            T_func_ida = T
        else:
            T_func_ida = lambda t: T  # noqa: E731

        # Algebraic variable indices: equilibrium-dependent + prescribed (deduped)
        algebraic_idx = list(dict.fromkeys(
            [int(i) for i in model.equil_dep] + [int(i) for i in prescribed_idx]
        ))
        ode_rows = [i for i in range(n) if i not in set(algebraic_idx)]

        # Consistent initial ydot: rhs for ODE rows, 0 for algebraic
        r0 = model.residual(c_init, np.zeros(n), T_func_ida(t_span[0]))
        ydot0 = np.zeros(n)
        for i in ode_rows:
            ydot0[i] = -r0[i]

        def ida_residual(t, y, ydot, result):
            c = y.copy()
            for idx, fn in zip(prescribed_idx, prescribed_fns):
                c[idx] = fn(t)
            result[:] = model.residual(c, ydot, T_func_ida(t))
            for idx, fn in zip(prescribed_idx, prescribed_fns):
                result[idx] = y[idx] - fn(t)

        def ida_jac(t, y, ydot, residual, cj, J):
            c = y.copy()
            for idx, fn in zip(prescribed_idx, prescribed_fns):
                c[idx] = fn(t)
            J[:, :] = model.jacobian(c, np.zeros(n), T_func_ida(t))
            # Add cj·I for ODE rows (dF/dydot = I for those rows)
            for i in ode_rows:
                J[i, i] += cj
            # Prescribed rows: F_i = y_i - fn(t), independent of ydot
            for idx in prescribed_idx:
                J[idx, :] = 0.0
                J[idx, idx] = 1.0

        solver_kwargs: dict = dict(
            rtol=rtol,
            atol=atol,
            jacfn=ida_jac,
            old_api=False,
        )
        if algebraic_idx:
            solver_kwargs["algebraic_vars_idx"] = algebraic_idx

        ida_solver = ida_dae("ida", ida_residual, **solver_kwargs)
        output = ida_solver.solve(t_eval, c_init, ydot0)

        c_out = output.values.y.copy()
        for idx, fn in zip(prescribed_idx, prescribed_fns):
            c_out[:, idx] = np.array([fn(t) for t in output.values.t])

        return SimulationResult(
            t=output.values.t,
            c=c_out,
            species=species_names,
            success=output.flag >= 0,
            message=str(output.message),
            T_profile=(
                np.array([T_func_ida(t) for t in output.values.t])
                if callable(T) else None
            ),
        )

    # --- coupled energy balance ---
    if solvent_composition is not None:
        if callable(T):
            raise ValueError(
                "solvent_composition cannot be combined with a callable T. "
                "Pass T as a float (initial temperature)."
            )
        T0 = float(T if T is not None else model.T)

        # Solvent indices in the full state vector (solvents are now in c)
        solvent_idx = {
            name: model.species_index[name]
            for name in solvent_composition
            if name in model.species_index
        }

        def _inject_solvents(c: np.ndarray, t: float) -> np.ndarray:
            """Return a copy of c with solvent concentrations set from solvent_composition."""
            c = c.copy()
            for name, val in solvent_composition.items():
                if name in solvent_idx:
                    x = val(t) if callable(val) else val
                    c[solvent_idx[name]] = x * model.species[solvent_idx[name]].c_ref
            return c

        # Initialise solvent concentrations in c_init
        for name, val in solvent_composition.items():
            if name in solvent_idx:
                x0 = val(t_span[0]) if callable(val) else val
                c_init[solvent_idx[name]] = x0 * model.species[solvent_idx[name]].c_ref

        y_init = np.append(c_init, T0)

        # Warn once for any reactions that cannot contribute to the energy balance.
        import warnings
        for rxn in model.reactions:
            if getattr(rxn, "equilibrium_constant", None) is None:
                warnings.warn(
                    f"{type(rxn).__name__} has no equilibrium constant and cannot "
                    "contribute to the energy balance. Its heat of reaction will be "
                    "treated as zero. Use ThermodynamicReaction with an "
                    "EquilibriumConstantVantHoff (or similar) to include it.",
                    UserWarning,
                    stacklevel=3,
                )

        def rhs_coupled(t: float, y: np.ndarray) -> np.ndarray:
            c = _inject_solvents(y[:n], t)
            T_cur = float(y[n])
            dc_dt = -model.residual(c, np.zeros(n), T_cur)
            # Solvents are prescribed — zero their ODE rows
            for i in solvent_idx.values():
                dc_dt[i] = 0.0
            rho_cp = model.volumetric_heat_capacity(c)
            state = model.make_state(c, T_cur)
            Q_dot = 0.0
            for j, rxn in enumerate(model.reactions):
                if not model.kinetic_mask[j]:
                    continue
                eq = getattr(rxn, "equilibrium_constant", None)
                if eq is None:
                    continue
                dH = eq.reaction_enthalpy(T_cur)
                Q_dot += dH * rxn.net_rate(state, model.species_index, model.charges)
            return np.append(dc_dt, -Q_dot / rho_cp)

        def jac_coupled(t: float, y: np.ndarray) -> np.ndarray:
            c = _inject_solvents(y[:n], t)
            T_cur = float(y[n])
            state = model.make_state(c, T_cur)
            rho_cp = model.volumetric_heat_capacity(c)
            J = np.zeros((n + 1, n + 1))
            J[:n, :n] = -model.jacobian(c, np.zeros(n), T_cur)
            J[:n, n] = -model.jacobian_dT(c, np.zeros(n), T_cur)
            # Zero Jacobian rows for solvent species (prescribed, no ODE)
            for i in solvent_idx.values():
                J[i, :] = 0.0
            # Analytic energy-balance row:
            #   ∂rhs_T/∂c_k = -Σ_j [ΔH_j / ρCp] · ∂φ_j/∂c_k
            #   ∂rhs_T/∂T   = -Σ_j [ΔH_j · ∂φ_j/∂T + φ_j · d(ΔH_j)/dT] / ρCp
            for j, rxn in enumerate(model.reactions):
                if not model.kinetic_mask[j]:
                    continue
                eq = getattr(rxn, "equilibrium_constant", None)
                if eq is None:
                    continue
                dH = eq.reaction_enthalpy(T_cur)
                dH_dT = eq.d_reaction_enthalpy_dT(T_cur)
                phi = rxn.net_rate(state, model.species_index, model.charges)
                dphi_dc = rxn.net_rate_jac(state, model.species_index, model.charges)
                dphi_dT = rxn.net_rate_dT(state, model.species_index, model.charges)
                J[n, :n] -= dH / rho_cp * dphi_dc
                J[n, n] -= (dH * dphi_dT + phi * dH_dT) / rho_cp
            return J

        sol = solve_ivp(
            rhs_coupled,
            t_span,
            y_init,
            method="Radau",
            t_eval=t_eval,
            jac=jac_coupled,
            rtol=rtol,
            atol=atol,
            dense_output=False,
        )
        return SimulationResult(
            t=sol.t,
            c=sol.y[:n, :].T,
            species=species_names,
            success=sol.success,
            message=sol.message,
            T_profile=sol.y[n, :],
        )

    # --- isothermal or prescribed T(t) ---
    T_func: Callable[[float], Optional[float]]
    if callable(T):
        T_func = T
    else:
        T_func = lambda t: T  # noqa: E731  (scalar or None — constant)

    if has_prescribed:
        # Solve only for dynamic species; inject prescribed values before each call.
        n_dyn = len(dynamic_idx)

        def _full_c(c_dyn: np.ndarray, t: float) -> np.ndarray:
            c_full = np.empty(n)
            c_full[dynamic_idx] = c_dyn
            for i, fn in zip(prescribed_idx, prescribed_fns):
                c_full[i] = fn(t)
            return c_full

        def rhs(t: float, c_dyn: np.ndarray) -> np.ndarray:
            c_full = _full_c(c_dyn, t)
            return (-model.residual(c_full, np.zeros(n), T_func(t)))[dynamic_idx]

        def jac(t: float, c_dyn: np.ndarray) -> np.ndarray:
            c_full = _full_c(c_dyn, t)
            J_full = -model.jacobian(c_full, np.zeros(n), T_func(t))
            return J_full[np.ix_(dynamic_idx, dynamic_idx)]

        sol = solve_ivp(
            rhs,
            t_span,
            c_init[dynamic_idx],
            method="Radau",
            t_eval=t_eval,
            jac=jac,
            rtol=rtol,
            atol=atol,
            dense_output=False,
        )

        c_out = np.zeros((len(sol.t), n))
        c_out[:, dynamic_idx] = sol.y.T
        for i, fn in zip(prescribed_idx, prescribed_fns):
            c_out[:, i] = np.array([fn(t) for t in sol.t])

        return SimulationResult(
            t=sol.t,
            c=c_out,
            species=species_names,
            success=sol.success,
            message=sol.message,
            T_profile=np.array([T_func(t) for t in sol.t]) if callable(T) else None,
        )

    def rhs(t: float, c: np.ndarray) -> np.ndarray:
        # residual(c, c_dot=0) = c_dot - f(c)  =>  f(c) = -residual(c, 0)
        return -model.residual(c, np.zeros_like(c), T_func(t))

    def jac(t: float, c: np.ndarray) -> np.ndarray:
        return -model.jacobian(c, np.zeros_like(c), T_func(t))

    sol = solve_ivp(
        rhs,
        t_span,
        c_init,
        method="Radau",
        t_eval=t_eval,
        jac=jac,
        rtol=rtol,
        atol=atol,
        dense_output=False,
    )

    T_profile = np.array([T_func(t) for t in sol.t]) if callable(T) else None

    return SimulationResult(
        t=sol.t,
        c=sol.y.T,
        species=species_names,
        success=sol.success,
        message=sol.message,
        T_profile=T_profile,
    )


# ---------------------------------------------------------------------------
# solve_equilibrium
# ---------------------------------------------------------------------------


def solve_equilibrium(
    model: ReactionModel,
    c0: dict[str, float],
    T: Optional[float] = None,
    prescribed: Optional[dict] = None,
    max_iter: int = 200,
    tol: float = 1e-12,
) -> dict[str, float]:
    """
    Solve for equilibrium concentrations using Newton's method in log space.

    Working in y = ln(c) keeps concentrations strictly positive throughout
    the iteration.  The full DAE residual with c_dot = 0 is used:

        kinetic rows:   -sum_j nu_ij * v_j(c) = 0   (zero net flux)
        equil rows:     ln(Q_j) - ln(K_j)     = 0   (K satisfied)

    The Jacobian is transformed via the chain rule:
        dr/dy_k = dr/dc_k * c_k

    A backtracking line search prevents divergence.

    Parameters
    ----------
    model : ReactionModel
    c0 : dict[str, float]
        Initial guess {species_name: concentration [mol/m³]}.
        Missing species default to 1e-10.
    T : float, optional
        Temperature [K].  Uses model default if not provided.
    prescribed : dict, optional
        Species held fixed {species_name: float | callable(t) -> float}.
        These are excluded from the Newton iteration and held at their
        specified value (evaluated at t=0 for callables).
        Intended for parametric sweeps (e.g. fixing c_H+ to scan pH).
        Solvent concentrations are handled by the ``Solution`` class.
    max_iter : int
        Maximum Newton iterations.
    tol : float
        Convergence tolerance on residual norm.

    Returns
    -------
    dict[str, float]
        Equilibrium concentrations [mol/m³].

    Raises
    ------
    RuntimeError
        If the solver does not converge within max_iter iterations.
    """
    species_names = [sp.name for sp in model.species]
    n = len(species_names)

    # Prescribed species: keep fixed at the specified concentration.
    _prescribed = prescribed or {}
    _prescribed_fns: dict[str, Callable[[float], float]] = {}
    for _pname, _pval in _prescribed.items():
        if callable(_pval):
            _prescribed_fns[_pname] = _pval
        else:
            _v = float(_pval)
            _prescribed_fns[_pname] = lambda t, v=_v: v
    prescribed_mask = np.array([sp.name in _prescribed_fns for sp in model.species])
    dynamic_idx = np.where(~prescribed_mask)[0]
    prescribed_idx = np.where(prescribed_mask)[0]
    n_dyn = len(dynamic_idx)

    c_init = np.array([max(c0.get(name, 1e-10), 1e-10) for name in species_names])
    for i in prescribed_idx:
        c_init[i] = _prescribed_fns[model.species[i].name](0.0)

    # -----------------------------------------------------------------------
    # Conservation enforcement
    #
    # For an equilibrium-only (DAE) system, some ODE rows are trivially zero
    # (species that are neither dep in any equil reaction nor involved in any
    # kinetic reaction).  These rows carry no gradient information and the
    # Newton step leaves the corresponding species at their initial values,
    # violating stoichiometric conservation.
    #
    # Fix: identify those trivially unconstrained dynamic species, compute
    # the left null space of the stoichiometric matrix (the conservation
    # moieties), and replace their zero residual rows with explicit
    # conservation equations  w^T c = w^T c0.
    # -----------------------------------------------------------------------
    equil_dep_set = set(int(d) for d in model.equil_dep)
    kinetic_col_idx = [j for j, is_k in enumerate(model.kinetic_mask) if is_k]

    # Local indices (within dynamic_idx) of trivially unconstrained species.
    trivial_local: list[int] = []
    for local_i, global_i in enumerate(dynamic_idx):
        if global_i in equil_dep_set:
            continue  # algebraic row — already constrained
        if kinetic_col_idx and np.any(model.nu[global_i, kinetic_col_idx] != 0):
            continue  # has a kinetic ODE contribution — not trivial
        trivial_local.append(local_i)

    # Left null space of N_dyn: conservation moiety vectors, shape (n_dyn, n_moiety)
    if trivial_local:
        N_dyn = model.nu[np.ix_(dynamic_idx, list(range(len(model.reactions))))]
        U_s, sv, _ = np.linalg.svd(N_dyn, full_matrices=True)
        sv_tol = max(N_dyn.shape) * np.finfo(float).eps * (sv.max() if len(sv) > 0 else 1.0)
        n_sv = min(N_dyn.shape)
        is_null = np.concatenate([sv < sv_tol, np.ones(n_dyn - n_sv, dtype=bool)])
        moiety_vecs = U_s[:, is_null]               # (n_dyn, n_moiety)
        moiety_targets = moiety_vecs.T @ c_init[dynamic_idx]  # conserved quantities
    else:
        moiety_vecs = np.zeros((n_dyn, 0))
        moiety_targets = np.zeros(0)
    n_moiety = moiety_vecs.shape[1]

    def _apply_conservation(r: np.ndarray, c_dyn: np.ndarray) -> np.ndarray:
        """Replace trivial ODE rows with conservation equations."""
        r = r.copy()
        for slot, local_i in enumerate(trivial_local):
            if slot >= n_moiety:
                break
            w = moiety_vecs[:, slot]
            r[local_i] = np.dot(w, c_dyn) - moiety_targets[slot]
        return r

    def _conservation_jac_rows(J_y: np.ndarray, c_dyn: np.ndarray) -> np.ndarray:
        """Set trivial rows of J_y to conservation equation Jacobians."""
        J_y = J_y.copy()
        for slot, local_i in enumerate(trivial_local):
            if slot >= n_moiety:
                break
            w = moiety_vecs[:, slot]
            J_y[local_i, :] = w * c_dyn   # ∂(w^T c)/∂y_k = w_k * c_k
        return J_y

    # Newton in log space, reduced to dynamic species only.
    y = np.log(c_init[dynamic_idx])
    c_dot_zero = np.zeros(n)

    for iteration in range(max_iter):
        c = c_init.copy()
        c[dynamic_idx] = np.exp(y)
        c_dyn = np.exp(y)

        r_full = model.residual(c, c_dot_zero, T)
        r = _apply_conservation(r_full[dynamic_idx], c_dyn)
        r_norm = float(np.linalg.norm(r))

        if r_norm < tol:
            break

        # Jacobian in log space: dr/dy_k = dr/dc_k * c_k
        J_c_full = model.jacobian(c, c_dot_zero, T)
        J_c = J_c_full[np.ix_(dynamic_idx, dynamic_idx)]
        J_y = _conservation_jac_rows(J_c * c_dyn[np.newaxis, :], c_dyn)

        # Tikhonov regularisation for near-singular systems
        J_y += np.eye(n_dyn) * 1e-12

        try:
            dy = np.linalg.solve(J_y, -r)
        except np.linalg.LinAlgError:
            dy = np.linalg.lstsq(J_y, -r, rcond=None)[0]

        # Backtracking line search
        alpha = 1.0
        for _ in range(20):
            y_new = y + alpha * dy
            c_new = c_init.copy()
            c_new[dynamic_idx] = np.exp(y_new)
            c_dyn_new = np.exp(y_new)
            r_new_full = model.residual(c_new, c_dot_zero, T)[dynamic_idx]
            r_new = _apply_conservation(r_new_full, c_dyn_new)
            if np.linalg.norm(r_new) < r_norm * (1.0 - 1e-4 * alpha):
                break
            alpha *= 0.5

        y = y + alpha * dy

    c_sol = c_init.copy()
    c_sol[dynamic_idx] = np.exp(y)
    c_dyn_sol = c_sol[dynamic_idx]
    r_sol_full = model.residual(c_sol, c_dot_zero, T)[dynamic_idx]
    r_final_norm = float(np.linalg.norm(
        _apply_conservation(r_sol_full, c_dyn_sol)
    ))

    if r_final_norm > 1e-6:
        raise RuntimeError(
            f"solve_equilibrium did not converge after {max_iter} iterations. "
            f"Final residual norm = {r_final_norm:.2e}.  "
            f"Try a better initial guess or check that the system is well-posed."
        )

    return dict(zip(species_names, c_sol))
