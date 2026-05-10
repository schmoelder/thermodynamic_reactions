---
kernelspec:
  name: python3
  display_name: Python 3 (ipykernel)
---

(nonidealities)=
# Pure-Fluid Non-Ideality: Equations of State and Fugacity

Non-ideality has two conceptually distinct origins.
Even a pure fluid can deviate from ideal behaviour because intermolecular interactions alter the equation of state; this is captured by **fugacity**, which replaces pressure in the chemical potential.
Mixtures introduce a second, independent layer: the interactions between unlike species can differ from those in the pure components, producing excess mixing properties and activity coefficients.
Fugacity handles non-ideal state behaviour; activity coefficients handle non-ideal mixing behaviour.
The two corrections are independent: a mixture of ideal gases has fugacity coefficients of unity but can still have non-unity activity coefficients, and a pure non-ideal fluid has a fugacity correction but no mixing non-ideality.

This chapter addresses the first source: pure-fluid departures from the ideal equation of state.
Mixture non-idealities and excess Gibbs energy are developed in @mixing.


## Real gases: the van der Waals equation

An **equation of state** is a relation between the macroscopic state variables $P$, $V$, $T$, and $n$ of a substance.
The ideal gas law $PV = nRT$ is the simplest; it works well at low pressure where molecules are far apart and interactions are negligible.
It fails at higher pressures because it ignores two effects: molecules have finite size (they cannot overlap), and they attract each other at short range.
The simplest equation of state that accounts for both is the **van der Waals equation** (1873):

$$
\left(P + \frac{an^2}{V^2}\right)\left(V - nb\right) = nRT
$$

The correction $nb$ subtracts the volume actually occupied by the molecules themselves, so $V - nb$ is the free volume available for motion.
The term $an^2/V^2$ corrects for attractive interactions: molecules pulling on each other reduce the force they exert on the walls, lowering the effective pressure.
The constants $a$ and $b$ are specific to each gas and must be measured; when $a = b = 0$ the ideal gas law is recovered.

An **isotherm** traces the $P$–$V_m$ relationship at a single fixed temperature; plotting several isotherms together reveals how the phase structure of a fluid changes with temperature.

Below $T_c$, the van der Waals equation predicts an S-shaped isotherm with multiple mathematical roots at a given pressure: a small-volume root (**liquid**) and a large-volume root (**vapour**), separated by an intermediate region.
This intermediate branch is not physically realised: the mechanical stability condition $(\partial P/\partial V)_T < 0$ is violated there, making the system unstable against infinitesimal fluctuations in volume.
Instead of following this branch, the system lowers its Gibbs free energy by separating into two phases at a single pressure, replacing the unstable region by a horizontal coexistence plateau — the macroscopic manifestation of global Gibbs minimisation.

The liquid is a dense, condensed phase in which attractive interactions hold molecules close together; the vapour is a dilute, disordered phase in which they are negligible.
At the **critical point** $(T_c, P_c, V_c)$ the distinction vanishes; above $T_c$ only a single fluid phase exists.
Every substance with attractive intermolecular interactions has a critical point.
The van der Waals equation locates it from $a$ and $b$ alone:

$$
T_c = \frac{8a}{27Rb}, \qquad P_c = \frac{a}{27b^2}, \qquad V_c = 3b.
$$


```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-vdw-isotherms

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.optimize import brentq
from scipy.integrate import quad

R = 0.08206   # L·atm/(mol·K)
a = 3.640     # CO₂, L²·atm/mol²
b = 0.04267   # CO₂, L/mol

Tc = 8*a / (27*R*b)
Pc = a / (27*b**2)
Vc = 3*b

def P_vdw(Vm, T):
    return R*T / (Vm - b) - a / Vm**2

def coexistence(T):
    def dPdV(Vm):
        return -R*T / (Vm - b)**2 + 2*a / Vm**3
    Vm_s1 = brentq(dPdV, b*1.001, Vc*0.9999)
    Vm_s2 = brentq(dPdV, Vc*1.0001, 20.0)
    P_min = P_vdw(Vm_s1, T)   # local pressure minimum (liquid spinodal)
    P_max = P_vdw(Vm_s2, T)   # local pressure maximum (vapor spinodal)
    V_hi  = 500.0
    # P_coex lies in (P_min, P_max); P_min can be negative so clamp from below
    # by P_vdw(V_hi) so that the vapor-root search always has a sign change
    P_lo = max(P_min, P_vdw(V_hi, T)) + 1e-6
    P_hi = P_max - 1e-6
    def area(P_star):
        Vl = brentq(lambda V: P_vdw(V, T) - P_star, b*1.001, Vm_s1)
        Vv = brentq(lambda V: P_vdw(V, T) - P_star, Vm_s2, V_hi)
        return quad(lambda V: P_vdw(V, T) - P_star, Vl, Vv)[0]
    Pc_co = brentq(area, P_lo, P_hi)
    Vl = brentq(lambda V: P_vdw(V, T) - Pc_co, b*1.001, Vm_s1)
    Vv = brentq(lambda V: P_vdw(V, T) - Pc_co, Vm_s2, V_hi)
    return Pc_co, Vl, Vv

def isotherm_start(T):
    """Return Vm start, skipping the unphysical negative-P liquid branch when present."""
    def dPdV(Vm): return -R*T / (Vm - b)**2 + 2*a / Vm**3
    try:
        Vm_s1 = brentq(dPdV, b*1.001, Vc*0.9999)
    except ValueError:
        return b * 1.005   # supercritical: no spinodal
    if P_vdw(Vm_s1, T) < 0:
        Vm_s2 = brentq(dPdV, Vc*1.0001, 20.0)
        return brentq(lambda V: P_vdw(V, T), Vm_s1, Vm_s2)
    return b * 1.005

XLIM1 = 0.8
# Restrict saturation curve to T where Vv ≤ XLIM1 so dome outline closes at plot edge
T_cross = brentq(lambda T: coexistence(T)[2] - XLIM1, 0.70*Tc, 0.80*Tc)
T_coex1 = np.linspace(T_cross * 0.999, 0.998*Tc, 80)
coex1   = [coexistence(T) for T in T_coex1]
Pc_1    = [c[0] for c in coex1]
Vl_1    = [c[1] for c in coex1]
Vv_1    = [c[2] for c in coex1]

def draw_fills(ax, Pc_arr, Vl_arr, Vv_arr, xlim):
    # Gas: right of vapor branch down to dome bottom only (avoids linear artifact)
    ax.fill_betweenx(
        [Pc] + list(reversed(Pc_arr)),
        [Vc] + list(reversed(Vv_arr)),
        xlim,
        color="#fdae61", alpha=0.18,
    )
    # Two-phase dome
    ax.fill(
        Vl_arr + [Vc] + list(reversed(Vv_arr)) + [Vv_arr[0], Vl_arr[0]],
        Pc_arr + [Pc] + list(reversed(Pc_arr)) + [0.0, 0.0],
        color="#888888", alpha=0.12,
    )
    # Liquid: left of liquid branch, extend to axis
    ax.fill_betweenx(
        [0.0] + Pc_arr + [Pc],
        0,
        [Vl_arr[0]] + Vl_arr + [Vc],
        color="#4393c3", alpha=0.18,
    )
    # Saturation curve outline
    ax.plot(
        Vl_arr + [Vc] + list(reversed(Vv_arr)),
        Pc_arr + [Pc] + list(reversed(Pc_arr)),
        color="#888888", lw=1.4, ls="--", label="saturation curve",
    )

T_show = [0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15]
# Subcritical: dark→light blue as T→Tc; Tc: black; supercritical: light→dark red
_n_sub = sum(1 for Tr in T_show if Tr < 1.0)
_n_sup = sum(1 for Tr in T_show if Tr > 1.0)
_sub_v = np.linspace(0.85, 0.45, _n_sub)
_sup_v = np.linspace(0.45, 0.85, _n_sup)
_si = _ri = 0
colors = []
for Tr in T_show:
    if   Tr < 1.0: colors.append(cm.Blues(_sub_v[_si])); _si += 1
    elif Tr > 1.0: colors.append(cm.Reds(_sup_v[_ri]));  _ri += 1
    else:          colors.append("black")

fig, ax = plt.subplots(figsize=(5.5, 4.5))
draw_fills(ax, Pc_1, Vl_1, Vv_1, XLIM1)

Vm_base = np.linspace(b*1.005, XLIM1 + 0.2, 6000)
for Tr, col in zip(T_show, colors):
    Vm_s = isotherm_start(Tr * Tc)
    Vm   = Vm_base[Vm_base >= Vm_s]
    P    = P_vdw(Vm, Tr * Tc)
    lbl  = r"$T = T_c$" if Tr == 1.00 else rf"$T = {Tr:.2f}\,T_c$"
    ax.plot(Vm, np.where(P < 0, np.nan, P), color=col, lw=2.2 if Tr == 1.00 else 1.8, label=lbl)
    if Tr < 1.0:
        Pc_co, Vl_co, Vv_co = coexistence(Tr * Tc)
        ax.plot([Vl_co, Vv_co], [Pc_co, Pc_co], color=col, ls="--", lw=1.0, zorder=4)
        ax.plot([Vl_co, Vv_co], [Pc_co, Pc_co], "o", color=col, ms=4, zorder=5)

ax.plot(Vc, Pc, "ko", ms=6, zorder=5)
ax.annotate(
    "critical\npoint", xy=(Vc, Pc), xytext=(Vc + 0.02, Pc + 2),
    fontsize=8, ha="center",
)

ax.text(0.025, 42, "liquid", fontsize=8, color="#2166ac",
        va="center", ha="center", rotation=90)
ax.text(0.3, 10, "two-phase\nregion", fontsize=9, color="#555555",
        va="center", ha="center")
ax.text(0.50, 60, "gas", fontsize=9, color="#d6604d", va="center")

ax.axhline(Pc, color="#888888", lw=0.7, ls=":", alpha=0.7)
ax.text(0.4, 90, "supercritical fluid",
        fontsize=8, color="#888888", va="top", ha="center")

ax.set_xlim(0, XLIM1)
ax.set_ylim(0, 100)
ax.set_xlabel(r"$V_m$ / (L mol$^{-1}$)")
ax.set_ylabel(r"$P$ / atm")
ax.legend(fontsize=8, loc="upper right")
fig.tight_layout()
```

```{figure} #cell-vdw-isotherms
:name: fig-vdw-isotherms

Van der Waals isotherms for $\ce{CO2}$ ($a = 3.640\ \mathrm{L^2\,atm\,mol^{-2}}$, $b = 0.04267\ \mathrm{L\,mol^{-1}}$) spanning subcritical, critical, and supercritical temperatures.
The saturation curve (dashed grey) is the region of liquid–vapour coexistence and separates the liquid, gas, and two-phase regions.
Below $T_c$, each van der Waals isotherm develops a pressure loop; the physical coexistence pressure is found by the Maxwell equal-area construction, and the dashed tie-line connects the saturated liquid ($V_l$) and saturated vapour ($V_v$), with phase fractions determined by the lever rule.
As $T \to T_c$, the coexistence region contracts and the loop collapses to an inflection point at the critical point $(T_c, P_c, V_c)$.
Above $P_c$ (dotted line) the liquid–gas distinction disappears and the substance becomes a supercritical fluid.
Each isotherm is shown only where $P \geq 0$; the negative-pressure branch corresponds to metastable tension states.
```


## Fugacity: the effective pressure

The ideal chemical potential for a gas, $\mu = \mu^\circ + RT\ln(P/P^\circ)$, was derived assuming $PV = nRT$.
For a real gas this expression is no longer exact.
The **fugacity** $f$ is defined as the quantity that makes the chemical potential expression exact for any gas:

$$
\mu = \mu^\circ(T) + RT\ln\frac{f}{P^\circ}
$$

It has units of pressure and can be thought of as an effective pressure that accounts for intermolecular interactions.
The ratio $\varphi = f/P$ is the **fugacity coefficient**: it measures how far the gas deviates from ideal behaviour.
By construction, $\varphi \to 1$ as $P \to 0$ (all gases become ideal at low pressure), so $f \to P$ in the dilute limit.

For the van der Waals gas, the fugacity coefficient can be derived analytically.
At moderate pressures, $\ln\varphi \approx (b - a/RT) \cdot P/RT$: repulsion ($b > 0$) increases $f$ above $P$; attraction ($a > 0$) decreases it.

```{admonition} Intuition
:class: tip

Attractions reduce a molecule's tendency to escape or react: the gas behaves as if at lower pressure than it actually is, so $f < P$ and $\varphi < 1$.
At very high pressure, hard-core repulsion dominates and $f > P$.
The limit $\varphi \to 1$ as $P \to 0$ reflects that all gases become ideal at low density.
```

---

The equation-of-state correction developed here applies to pure fluids and to mixture components individually.
Mixing these non-ideal components introduces an additional source of departure: the Gibbs energy of mixing deviates from the ideal baseline when unlike interactions differ from pure-component ones.
The next chapter develops that second layer, derives activity coefficients as partial molar derivatives of the excess Gibbs energy, and introduces the electrolyte models used throughout Parts 3 and 4.
