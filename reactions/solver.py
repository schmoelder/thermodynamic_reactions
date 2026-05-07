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

from reactions.api import ReactionModel


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
    n_points: int = 500,
    rtol: float = 1e-8,
    atol: float = 1e-10,
) -> SimulationResult:
    """
    Integrate a ReactionModel forward in time using scipy Radau.

    Radau is an implicit solver suited for stiff systems.  The analytic
    Jacobian from ReactionModel.jacobian() is used automatically.

    For systems with equilibrium reactions, the algebraic constraints are
    satisfied at each step via the residual structure — equil rows replace
    ODE rows for the dependent species, so the system remains determined.
    A proper DAE solver (SUNDIALS IDA) would handle the index more
    rigorously; this is a good approximation for prototyping.

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
    n_points : int
        Number of output time points.
    rtol, atol : float
        Relative and absolute tolerances.

    Returns
    -------
    SimulationResult
        result.T_profile is populated when T is callable or heat_capacity
        is given; None otherwise.
    """
    species_names = [sp.name for sp in model.species]
    c_init = np.array([c0.get(name, 0.0) for name in species_names])
    n = len(c_init)
    t_eval = np.linspace(t_span[0], t_span[1], n_points)

    # --- coupled energy balance ---
    if solvent_composition is not None:
        if callable(T):
            raise ValueError(
                "solvent_composition cannot be combined with a callable T. "
                "Pass T as a float (initial temperature)."
            )
        T0 = float(T if T is not None else model.T)
        y_init = np.append(c_init, T0)

        def _x_k(t: float) -> dict:
            return {
                name: (val(t) if callable(val) else val)
                for name, val in solvent_composition.items()
            }

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
            c = y[:n]
            T_cur = float(y[n])
            dc_dt = -model.residual(c, np.zeros(n), T_cur)
            rho_cp = model.volumetric_heat_capacity(_x_k(t))
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
            c = y[:n]
            T_cur = float(y[n])
            state = model.make_state(c, T_cur)
            rho_cp = model.volumetric_heat_capacity(_x_k(t))
            J = np.zeros((n + 1, n + 1))
            J[:n, :n] = -model.jacobian(c, np.zeros(n), T_cur)
            J[:n, n] = -model.jacobian_dT(c, np.zeros(n), T_cur)
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
    c_init = np.array([max(c0.get(name, 1e-10), 1e-10) for name in species_names])
    y = np.log(c_init)
    c_dot_zero = np.zeros(len(y))

    for iteration in range(max_iter):
        c = np.exp(y)
        r = model.residual(c, c_dot_zero, T)
        r_norm = float(np.linalg.norm(r))

        if r_norm < tol:
            break

        # Jacobian in log space: dr/dy_k = dr/dc_k * c_k
        J_c = model.jacobian(c, c_dot_zero, T)
        J_y = J_c * c[np.newaxis, :]

        # Tikhonov regularisation for near-singular systems
        J_y += np.eye(len(y)) * 1e-12

        try:
            dy = np.linalg.solve(J_y, -r)
        except np.linalg.LinAlgError:
            dy = np.linalg.lstsq(J_y, -r, rcond=None)[0]

        # Backtracking line search
        alpha = 1.0
        for _ in range(20):
            y_new = y + alpha * dy
            c_new = np.exp(y_new)
            r_new = model.residual(c_new, c_dot_zero, T)
            if np.linalg.norm(r_new) < r_norm * (1.0 - 1e-4 * alpha):
                break
            alpha *= 0.5

        y = y + alpha * dy

    c_sol = np.exp(y)
    r_final = model.residual(c_sol, c_dot_zero, T)
    r_final_norm = float(np.linalg.norm(r_final))

    if r_final_norm > 1e-6:
        raise RuntimeError(
            f"solve_equilibrium did not converge after {max_iter} iterations. "
            f"Final residual norm = {r_final_norm:.2e}.  "
            f"Try a better initial guess or check that the system is well-posed."
        )

    return dict(zip(species_names, c_sol))
