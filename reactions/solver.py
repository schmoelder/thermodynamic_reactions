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
from typing import Optional

import numpy as np
import matplotlib.pyplot as plt
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
    T: Optional[float] = None,
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
    T : float, optional
        Temperature [K].  Uses model default if not provided.
    n_points : int
        Number of output time points.
    rtol, atol : float
        Relative and absolute tolerances.

    Returns
    -------
    SimulationResult
    """
    species_names = [sp.name for sp in model.species]
    c_init = np.array([c0.get(name, 0.0) for name in species_names])

    def rhs(t: float, c: np.ndarray) -> np.ndarray:
        # residual(c, c_dot=0) = c_dot - f(c)  =>  f(c) = -residual(c, 0)
        return -model.residual(c, np.zeros_like(c), T)

    def jac(t: float, c: np.ndarray) -> np.ndarray:
        return -model.jacobian(c, np.zeros_like(c), T)

    t_eval = np.linspace(t_span[0], t_span[1], n_points)

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

    return SimulationResult(
        t=sol.t,
        c=sol.y.T,
        species=species_names,
        success=sol.success,
        message=sol.message,
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
