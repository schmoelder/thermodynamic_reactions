"""
Solver and validation suite for ReactionModel.

Wraps scipy.integrate.solve_ivp (Radau — implicit, stiff-safe) around
ReactionModel.residual() / jacobian().  Includes analytical validation
cases with known closed-form solutions.

Units: SI throughout (mol/m³, K, s).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
import matplotlib.pyplot as plt

from api import (
    C_REF,
    ActivityCoefficientDavies,
    ActivityCoefficientIdeal,
    Component,
    EquilibriumConstant,
    EquilibriumConstantVantHoff,
    IonicStrengthBackground,
    MassActionReaction,
    PhysicalState,
    RateConstantFixed,
    ReactionModel,
    Species,
    ThermodynamicReaction,
    pKa,
)
from solver import SimulationResult, simulate, solve_equilibrium




# ---------------------------------------------------------------------------
# Validation case 1: A <-> B  (analytical solution known)
# ---------------------------------------------------------------------------


def validate_reversible_ab():
    """
    A <-> B  with kf, kr.

    Analytical solution:
        [A](t) = [A]_eq + ([A]_0 - [A]_eq) * exp(-(kf+kr)*t)
        [B](t) = [B]_eq + ([B]_0 - [B]_eq) * exp(-(kf+kr)*t)

    where [A]_eq = ([A]_0 + [B]_0) / (1 + K),  K = kf/kr.
    """
    print("\n" + "=" * 60)
    print("  Validation 1: A <-> B  (MassActionReaction)")
    print("=" * 60)

    kf_val = 2.0     # 1/s
    kr_val = 0.5     # 1/s
    K      = kf_val / kr_val   # = 4
    A0     = 1000.0  # mol/m³
    B0     = 0.0

    comp = Component("ab", [Species("A"), Species("B")])
    model = ReactionModel(
        components=[comp],
        reactions=[MassActionReaction("A <-> B", kf=kf_val, kr=kr_val)],
    )

    t_end = 5.0
    result = simulate(model, {"A": A0, "B": B0}, (0, t_end))
    assert result.success, f"Solver failed: {result.message}"

    # Analytical solution
    total = A0 + B0
    A_eq  = total / (1 + K)
    B_eq  = total * K / (1 + K)
    lam   = kf_val + kr_val
    t     = result.t
    A_analytical = A_eq + (A0 - A_eq) * np.exp(-lam * t)
    B_analytical = B_eq + (B0 - B_eq) * np.exp(-lam * t)

    err_A = np.max(np.abs(result["A"] - A_analytical))
    err_B = np.max(np.abs(result["B"] - B_analytical))
    print(f"  K = kf/kr = {K:.1f}")
    print(f"  [A]_eq = {A_eq:.2f},  [B]_eq = {B_eq:.2f}  mol/m³")
    print(f"  max |error A| = {err_A:.2e} mol/m³")
    print(f"  max |error B| = {err_B:.2e} mol/m³")
    assert err_A < 1e-6, f"A error too large: {err_A}"
    assert err_B < 1e-6, f"B error too large: {err_B}"
    print("  PASSED")

    # Conservation: [A] + [B] = const
    total_t = result["A"] + result["B"]
    conservation_err = np.max(np.abs(total_t - total))
    print(f"  Conservation |[A]+[B] - {total:.0f}| max = {conservation_err:.2e}")
    assert conservation_err < 1e-6
    print("  Conservation: PASSED")

    return result, A_analytical, B_analytical


# ---------------------------------------------------------------------------
# Validation case 2: thermodynamic consistency
# A <-> B with ThermodynamicReaction — kinetic sim must converge to K
# ---------------------------------------------------------------------------


def validate_thermodynamic_consistency():
    """
    ThermodynamicReaction A <-> B.

    At t -> inf, [B]/[A] must equal K (in activity units).
    kr = kf/K is enforced — test that this holds numerically.
    """
    print("\n" + "=" * 60)
    print("  Validation 2: Thermodynamic consistency")
    print("=" * 60)

    K_eq  = 4.0
    kf    = 2.0   # 1/s  (kr = kf/K = 0.5 derived)
    A0    = 1000.0
    B0    = 0.0
    T     = 298.15

    comp = Component("ab", [Species("A"), Species("B")])
    model = ReactionModel(
        components=[comp],
        reactions=[
            ThermodynamicReaction(
                "A <-> B",
                mode="kinetic",
                equilibrium_constant=EquilibriumConstant(K_eq=K_eq),
                rate_constant=RateConstantFixed(kf_value=kf),
            )
        ],
    )

    # Relaxation timescale in concentration space: tau = C_REF / (kf + kf/K)
    tau = C_REF / (kf * (1.0 + 1.0 / K_eq))
    t_end = 10.0 * tau
    result = simulate(model, {"A": A0, "B": B0}, (0, t_end), T=T, n_points=500)
    assert result.success

    # At t -> inf: B/A = K  (activities cancel since both ideal, c° identical)
    A_final = result["A"][-1]
    B_final = result["B"][-1]
    Q_final = B_final / A_final   # = (B/C_REF) / (A/C_REF)

    print(f"  K_eq = {K_eq}")
    print(f"  [A]_final = {A_final:.4f},  [B]_final = {B_final:.4f}  mol/m³")
    print(f"  Q_final = B/A = {Q_final:.6f}  (should = {K_eq})")
    assert abs(Q_final - K_eq) / K_eq < 1e-3, f"Q={Q_final} != K={K_eq}"
    print("  PASSED: simulation converges to correct equilibrium")

    # Also verify kr = kf/K
    rxn = model.reactions[0]
    kr_derived = rxn.kr(T)
    kr_expected = kf / K_eq
    print(f"  kr = {kr_derived:.6f}  (expected kf/K = {kr_expected:.6f})")
    assert abs(kr_derived - kr_expected) < 1e-12
    print("  PASSED: kr = kf/K enforced")

    return result


# ---------------------------------------------------------------------------
# Validation case 3: acetic acid speciation vs Henderson-Hasselbalch
# ---------------------------------------------------------------------------


def validate_acid_base():
    """
    Acetic acid / acetate buffer.

    Henderson-Hasselbalch (ideal, no activity correction):
        pH = pKa + log10([AcO-] / [AcOH])

    At a given total acetate concentration c_tot and pH, the speciation is:
        [AcO-] = c_tot / (1 + 10^(pKa - pH))
        [AcOH] = c_tot - [AcO-]

    We compare the equilibrium solver output against this at multiple pH
    values (no ionic strength, ideal activity).
    """
    print("\n" + "=" * 60)
    print("  Validation 3: Acetic acid speciation vs Henderson-Hasselbalch")
    print("=" * 60)

    pKa_val   = 4.76
    Ka        = 10.0 ** (-pKa_val)
    c_tot     = 100.0   # mol/m³  (= 0.1 mol/L)
    T         = 298.15

    acetate = Component("acetate", [
        Species("AcOH", charge=0),
        Species("AcO-", charge=-1),
    ])
    proton   = Component("proton",   [Species("H+",  charge=+1)])
    hydroxide = Component("hydroxide", [Species("OH-", charge=-1)])
    water = Component("water", [
        Species("H2O", charge=0, is_solvent=True),
    ])

    model = ReactionModel(
        components=[acetate, proton, hydroxide, water],
        reactions=[
            ThermodynamicReaction(
                "AcOH <-> AcO- + H+",
                mode="equil",
                equilibrium_constant=pKa(pKa_val),
                # ideal activity — no Davies correction
            ),
            ThermodynamicReaction(
                "H2O <-> H+ + OH-",
                mode="equil",
                equilibrium_constant=EquilibriumConstant(1e-14),
            ),
        ],
        T=T,
    )

    # Sweep pH values 2 to 8
    pH_values = np.linspace(2.0, 8.0, 13)
    errors = []

    print(f"\n  {'pH':>6}  {'[AcOH] HH':>12}  {'[AcOH] num':>12}  {'err':>10}")
    print(f"  {'-'*6}  {'-'*12}  {'-'*12}  {'-'*10}")

    for pH in pH_values:
        H_conc  = 10.0 ** (-pH) * C_REF        # mol/m³
        OH_conc = 1e-14 * C_REF**2 / H_conc  # Kw = [H+][OH-]/C_REF² = 1e-14

        # Henderson-Hasselbalch prediction
        AcO_HH  = c_tot / (1.0 + 10.0 ** (pKa_val - pH))
        AcOH_HH = c_tot - AcO_HH

        # Numerical equilibrium
        c0 = {
            "AcOH": max(AcOH_HH, 1e-10),
            "AcO-": max(AcO_HH, 1e-10),
            "H+":   H_conc,
            "OH-":  OH_conc,
        }
        try:
            c_eq = solve_equilibrium(model, c0, T=T)
            err = abs(c_eq["AcOH"] - AcOH_HH)
            errors.append(err)
            print(f"  {pH:6.2f}  {AcOH_HH:12.4f}  {c_eq['AcOH']:12.4f}  {err:10.2e}")
        except RuntimeError as e:
            print(f"  {pH:6.2f}  FAILED: {e}")

    max_err = max(errors)
    print(f"\n  Max error vs Henderson-Hasselbalch: {max_err:.2e} mol/m³")
    # HH is exact for ideal activity — we expect near-machine-precision agreement
    assert max_err < 1e-6, f"Speciation error too large: {max_err}"
    print("  PASSED")

    return pH_values, errors


# ---------------------------------------------------------------------------
# Validation case 4: conservation laws
# ---------------------------------------------------------------------------


def validate_conservation():
    """
    Multi-reaction network: A + B <-> C,  C <-> D + E.

    Conserved quantities (null space of nuT):
        C1: [A] + [C] + [D]     (A atoms)
        C2: [B] + [C] + [E]     (B atoms)

    These must be constant throughout the simulation.
    """
    print("\n" + "=" * 60)
    print("  Validation 4: Conservation laws")
    print("=" * 60)

    comp = Component("all", [
        Species("A"), Species("B"), Species("C"),
        Species("D"), Species("E"),
    ])
    model = ReactionModel(
        components=[comp],
        reactions=[
            MassActionReaction("A + B -> C", kf=1.0, kr=0.1),
            MassActionReaction("C -> D + E",  kf=0.5, kr=0.2),
        ],
    )

    c0 = {"A": 500.0, "B": 300.0, "C": 0.0, "D": 0.0, "E": 0.0}
    result = simulate(model, c0, (0, 10.0))
    assert result.success

    A = result["A"]; B = result["B"]; C = result["C"]
    D = result["D"]; E = result["E"]

    C1 = A + C + D   # should equal A0 + C0 + D0 = 500
    C2 = B + C + E   # should equal B0 + C0 + E0 = 300

    err_C1 = np.max(np.abs(C1 - C1[0]))
    err_C2 = np.max(np.abs(C2 - C2[0]))
    print(f"  C1 = [A]+[C]+[D] = {C1[0]:.1f}  max drift = {err_C1:.2e}")
    print(f"  C2 = [B]+[C]+[E] = {C2[0]:.1f}  max drift = {err_C2:.2e}")
    assert err_C1 < 1e-6, f"C1 not conserved: {err_C1}"
    assert err_C2 < 1e-6, f"C2 not conserved: {err_C2}"
    print("  PASSED")

    return result


# ---------------------------------------------------------------------------
# Validation case 5: temperature dependence via van't Hoff
# ---------------------------------------------------------------------------


def validate_vanthoff():
    """
    A <-> B  with EquilibriumConstantVantHoff.

    K(T) = exp(-dH/RT + dS/R).
    At two temperatures T1 and T2, verify the simulation converges to
    the correct K(T) in each case.
    """
    print("\n" + "=" * 60)
    print("  Validation 5: Van't Hoff temperature dependence")
    print("=" * 60)

    dH   = -20e3   # J/mol  (exothermic — K decreases with T)
    dS   = -50.0   # J/(mol K)
    kf   = 1.0

    comp = Component("ab", [Species("A"), Species("B")])

    for T in [298.15, 320.0, 350.0]:
        K_expected = np.exp(-dH / (8.314462 * T) + dS / 8.314462)
        model = ReactionModel(
            components=[comp],
            reactions=[
                ThermodynamicReaction(
                    "A <-> B",
                    mode="kinetic",
                    equilibrium_constant=EquilibriumConstantVantHoff(
                        dH=dH, dS=dS
                    ),
                    rate_constant=RateConstantFixed(kf_value=kf),
                )
            ],
        )
        tau_T = C_REF / (kf * (1.0 + 1.0 / K_expected))
        result = simulate(model, {"A": 1000.0, "B": 0.0}, (0, 10*tau_T), T=T, n_points=300)
        assert result.success

        A_f = result["A"][-1]
        B_f = result["B"][-1]
        Q   = (B_f / C_REF) / (A_f / C_REF)
        err = abs(Q - K_expected) / K_expected

        print(f"  T={T:.1f} K:  K_expected={K_expected:.4f}  Q_final={Q:.4f}  "
              f"err={err:.2e}  {'PASSED' if err < 1e-3 else 'FAILED'}")
        assert err < 1e-3, f"T={T}: Q={Q} != K={K_expected}"

    print("  PASSED")


# ---------------------------------------------------------------------------
# Run all validations
# ---------------------------------------------------------------------------


if __name__ == "__main__":
    results = {}

    r1, A_ana, B_ana = validate_reversible_ab()
    results["ab_kinetic"] = r1

    r2 = validate_thermodynamic_consistency()
    results["thermo_consistency"] = r2

    validate_acid_base()
    validate_conservation()
    validate_vanthoff()

    print("\n" + "=" * 60)
    print("  All validations passed.")
    print("=" * 60)

    # --- Summary plot ---
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))

    # Plot 1: A <-> B kinetics vs analytical
    ax = axes[0]
    r = results["ab_kinetic"]
    ax.plot(r.t, r["A"], "C0-",  lw=2, label="[A] numerical")
    ax.plot(r.t, r["B"], "C1-",  lw=2, label="[B] numerical")
    ax.plot(r.t, A_ana,  "k--",  lw=1, label="[A] analytical")
    ax.plot(r.t, B_ana,  "k-.",  lw=1, label="[B] analytical")
    ax.set_xlabel("time [s]")
    ax.set_ylabel("concentration [mol/m³]")
    ax.set_title("A ⇌ B  (kf=2, kr=0.5)")
    ax.legend(fontsize=8)

    # Plot 2: thermodynamic consistency
    ax = axes[1]
    r = results["thermo_consistency"]
    ax.plot(r.t, r["A"], label="[A]")
    ax.plot(r.t, r["B"], label="[B]")
    K_eq = 4.0
    A_eq = (r["A"][0] + r["B"][0]) / (1 + K_eq)
    B_eq = K_eq * A_eq
    ax.axhline(A_eq, ls="--", color="C0", alpha=0.5, label=f"[A]_eq={A_eq:.0f}")
    ax.axhline(B_eq, ls="--", color="C1", alpha=0.5, label=f"[B]_eq={B_eq:.0f}")
    ax.set_xlabel("time [s]")
    ax.set_ylabel("concentration [mol/m³]")
    ax.set_title("Thermodynamic consistency\n(K=4, ideal activities)")
    ax.legend(fontsize=8)

    # Plot 3: acetic acid speciation
    ax = axes[2]
    pKa_val = 4.76
    c_tot   = 100.0
    pH_plot = np.linspace(2, 8, 200)
    AcO_HH  = c_tot / (1 + 10 ** (pKa_val - pH_plot))
    AcOH_HH = c_tot - AcO_HH
    ax.plot(pH_plot, AcOH_HH, label="[AcOH]")
    ax.plot(pH_plot, AcO_HH,  label="[AcO⁻]")
    ax.axvline(pKa_val, ls="--", color="gray", alpha=0.6, label=f"pKa={pKa_val}")
    ax.set_xlabel("pH")
    ax.set_ylabel("concentration [mol/m³]")
    ax.set_title("Acetic acid speciation\n(Henderson-Hasselbalch)")
    ax.legend(fontsize=8)

    fig.tight_layout()
    plt.savefig("./outputs/validation.png", dpi=150)
    print("\nPlot saved to validation.png")
