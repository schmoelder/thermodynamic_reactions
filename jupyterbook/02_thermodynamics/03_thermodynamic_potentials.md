---
kernelspec:
  name: python3
  display_name: Python 3 (ipykernel)
---

(thermodynamic-potentials)=
# Thermodynamic Potentials

The laws established in @laws-of-thermodynamics tell us what is conserved and what direction change takes.
This chapter builds the mathematical framework: the fundamental relation combines the first and second laws into a single equation, Legendre transforms then produce the thermodynamic potentials suited to different experimental conditions, and $G$ emerges as the right one for chemistry.

## The fundamental relation

Substituting the reversible equalities $\delta Q = T\,dS$ and $\delta W = P\,dV$ into the first law gives the **fundamental relation**:

$$
dU = T\,dS - P\,dV + \sum_i \mu_i\,dn_i
$$

This is the first and second laws combined into a single equation.
Because $U$, $S$, $V$, and $n_i$ are all state functions, this relation holds generally, not only for reversible processes, and is the foundation for everything that follows.

The three terms each pair an intensive variable ($T$, $P$, $\mu_i$) with the differential of an extensive one ($S$, $V$, $n_i$).
These are the **conjugate pairs** of thermodynamics.
Each product has units of energy: $T\,dS$ is heat, $P\,dV$ is mechanical work, and $\mu_i\,dn_i$ is chemical work.


## The integrated fundamental relation

Because $U$, $S$, $V$, and $n_i$ are all extensive (they scale linearly with system size), Euler's homogeneity theorem for first-order homogeneous functions gives:

$$
U = TS - PV + \sum_i \mu_i n_i
$$

This expresses total energy as a sum of conjugate products.
Having both the differential and the integrated form makes it possible to extract a third relation, as the next section shows.


## The Gibbs-Duhem equation

Differentiating the integrated form:

$$
dU = T\,dS + S\,dT - P\,dV - V\,dP + \sum_i \mu_i\,dn_i + \sum_i n_i\,d\mu_i
$$

Subtracting the fundamental relation cancels the $T\,dS$, $P\,dV$, and $\mu_i\,dn_i$ terms, leaving only the intensive differentials:

$$
S\,dT - V\,dP + \sum_i n_i\,d\mu_i = 0
$$

This is the **Gibbs-Duhem equation**: $T$, $P$, and $\{\mu_i\}$ are not all independent; changing one forces changes in the others.


## The four potentials

The natural variables of $U$ are $S$ and $V$, making it the appropriate potential for an isolated system where entropy and volume are fixed.
In experiments, however, entropy cannot be directly controlled or imposed; it is determined by the system's microscopic state.
Instead, temperature and pressure are the variables that can be fixed by external reservoirs.
To obtain a potential expressed in these experimentally controlled variables, $U$ is transformed via **Legendre transforms**, which replace an extensive variable with its conjugate intensive one while preserving the underlying state information.
Subtracting $TS$ replaces $S$ with $T$, and adding $PV$ replaces $V$ with $P$.
This produces the thermodynamic potentials $H$, $A$, and $G$, shown in @fig-potentials.

```{figure} figures/thermodynamic_potentials.png
:name: fig-potentials

The four thermodynamic potentials as Legendre transforms of $U$.
Each arrow swaps one extensive natural variable for its conjugate intensive one.
$G$ is the potential for chemistry: it is minimised at the constant $T$ and $P$ conditions of a laboratory.
```

The definitions and natural variables of each potential are summarised in the table below.

| Fixed conditions | Minimised potential   | Definition    | Natural variables |
| ---------------- | --------------------- | ------------- | ----------------- |
| $S$, $V$         | Internal energy $U$   |               | $S$, $V$, $n_i$   |
| $S$, $P$         | Enthalpy $H$          | $U + PV$      | $S$, $P$, $n_i$   |
| $T$, $V$         | Helmholtz energy $A$  | $U - TS$      | $T$, $V$, $n_i$   |
| $T$, $P$         | Gibbs free energy $G$ | $U - TS + PV$ | $T$, $P$, $n_i$   |

The Legendre transform machinery is general: it produces a valid thermodynamic potential for any choice of constraints.
Which potential is useful, however, depends entirely on what is held fixed experimentally.
Enthalpy $H$ arises when volume is replaced by pressure as a control variable, $H = U + PV$.
At constant pressure, changes in enthalpy correspond directly to heat exchange, $dH = \delta Q_P$.
Reaction enthalpies $\Delta_r H$ and the Kirchhoff relations are therefore expressed in terms of $H$.
The Helmholtz energy $A$ is the appropriate potential when temperature and volume are controlled, corresponding to a rigid system in thermal contact with a reservoir.
This is the natural setting for molecular simulation, where $N$, $V$, and $T$ are directly specified, and for certain condensed-matter models.
In typical solution chemistry, reactions occur in open vessels at approximately constant temperature and pressure, so the Gibbs free energy $G$ is the relevant potential; $A$ does not appear further in this book.

## Differentials of the potentials

The differential of each potential is fixed by its natural variables and the fundamental relation:

The differentials of each potential follow from the fundamental relation and their definitions:

$$
dH = T\,dS + V\,dP + \sum_i \mu_i\,dn_i
$$

$$
dA = -S\,dT - P\,dV + \sum_i \mu_i\,dn_i
$$

$$
dG = -S\,dT + V\,dP + \sum_i \mu_i\,dn_i
$$

Note that $\mu_i$ appears with the same coefficient $dn_i$ across all four potentials: chemical potential does not depend on which potential describes the system.

The differentials make abstract potentials measurable: every partial derivative with respect to a natural variable defines a thermodynamic response function, and all thermodynamic response quantities used in this book arise from this structure.

A canonical example of this derivative structure is the **heat capacity at constant pressure**:

$$
C_p = \left.\frac{\partial H}{\partial T}\right|_{P,n}
$$

At constant pressure, changes in enthalpy correspond to heat exchange, so $C_p$ measures the energetic cost of temperature change per kelvin.
Tabulated $C_p$ values allow extrapolation of enthalpies and entropies across temperature via the Kirchhoff relations, discussed in the equilibrium chapter.

```{admonition} Intuition
:class: tip
A large $C_p$ means the system resists temperature change: the same heat input produces a smaller rise.
This buffering capacity and its temperature dependence are what the Kirchhoff relations in @equilibrium-temperature quantify.
```

Other examples that recur frequently are:

- **partial molar volume**  
  $\bar{V}_i = (\partial V/\partial n_i)_{T,P,n_{j\neq i}}$, governing how concentration replaces pressure as a composition variable (@chemical-potential)
- **partial molar enthalpy**  
  $\bar{H}_i = (\partial H/\partial n_i)_{T,P,n_{j\neq i}}$, entering energy balances in reacting systems (@implementation-energy-balance)

In chemistry, temperature and pressure are specified as inputs: a thermostat fixes $T$ and a pressure controller or the surrounding atmosphere fixes $P$, while volume follows as a derived response through the equation of state.
The natural potential is therefore the Gibbs free energy $G(T,P,\{n_i\})$.

At fixed $T$ and $P$, composition is the only variable that changes during reaction, diffusion, or phase change, and all driving forces reduce to differences in chemical potential.

The partial molar quantities above reflect a deeper structural property: thermodynamic potentials are first-order homogeneous functions of their extensive variables; for $G(T,P,\{n_i\})$, this refers to the particle numbers $\{n_i\}$ at fixed $T$ and $P$.
This extensivity implies a global scaling relation that connects local derivatives to integrated forms via Euler's theorem.

Applying Euler's theorem to $G$ gives:

$$
G = \sum_i n_i \left(\frac{\partial G}{\partial n_i}\right)_{T,P} = \sum_i \mu_i n_i
$$

This identity holds generally and is not restricted to equilibrium or to constant $T$ and $P$.
The Gibbs free energy is therefore the composition-weighted sum of chemical potentials, a structure used throughout the chapters on equilibrium and reactions.

```{admonition} Why does $dG$ contain $-S\,dT + V\,dP$ if $G = \sum_i \mu_i n_i$?
:class: note

$G = \sum_i \mu_i n_i$ appears to depend only on composition, but each $\mu_i$ depends on $T$, $P$, and composition.

Differentiating gives:

$$
dG = \sum_i \mu_i\,dn_i + \sum_i n_i\,d\mu_i
$$

The Gibbs–Duhem equation gives $\sum_i n_i\,d\mu_i = -S\,dT + V\,dP$; substituting:

$$
dG = -S\,dT + V\,dP + \sum_i \mu_i\,dn_i
$$

consistent with the direct differential form.
```

## Thermodynamic equilibrium

Each thermodynamic potential is minimised under its own set of natural constraints: $dU = 0$ at fixed $S$ and $V$, $dH = 0$ at fixed $S$ and $P$, $dA = 0$ at fixed $T$ and $V$, and $dG = 0$ at fixed $T$ and $P$.
A system is in **thermodynamic equilibrium** when the relevant potential reaches this minimum.
In a laboratory, $T$ and $P$ are controlled, so $dG = 0$ is the operative condition; the choice of $G$ follows from the experimental constraints, not from any special property of $G$ itself.
Three conditions must hold simultaneously:

- **Thermal equilibrium**: temperature is uniform; no net heat flow occurs.
- **Mechanical equilibrium**: pressure is uniform; no net volume change occurs.
- **Chemical equilibrium**: chemical potentials are equal; no net matter transfer occurs.

The statistical basis for all three was established in @entropy: equilibrium corresponds to maximisation of total entropy.
For two systems exchanging energy, maximising entropy yields equal $\partial S/\partial U$, i.e. equal temperature.
Allowing volume exchange gives equal pressure, and allowing particle exchange gives equal chemical potentials.

```{figure} figures/thermodynamic_equilibrium.png
:name: fig-equilibrium

Two subsystems in contact.
Depending on the type of contact, different quantities equalise at equilibrium: a diathermal wall allows heat exchange and equalises $T$; a movable wall equalises $P$; a permeable membrane equalises $\mu_i$ for the species that can cross.
A rigid, sealed, insulating boundary fixes all three and leaves the system in whatever state it starts in.
```


## Le Chatelier's principle

Equilibrium conditions also determine how a system responds to external perturbations.

> **Le Chatelier's principle:** a perturbation in an intensive variable induces a shift in the system that reduces the resulting change in the corresponding conjugate extensive variable.

This behaviour follows directly from the conjugate structure of the fundamental relation:

- **Temperature** ($T$, conjugate to $S$): increasing temperature shifts equilibrium toward states of higher entropy, typically corresponding to the endothermic direction ($\Delta_r H > 0$), while decreasing temperature favours the exothermic direction.
- **Pressure** ($P$, conjugate to $V$): increasing pressure shifts equilibrium toward states of lower volume.
For gas-phase reactions, this often corresponds to a reduction in the number of moles of gas.
- **Chemical potential** ($\mu_i$, conjugate to $n_i$): increasing $\mu_i$ (e.g. by adding species $i$) drives the system to reduce its chemical potential by consuming it.
In reacting systems, this appears as a nonzero reaction Gibbs energy, $\Delta_r G \neq 0$, which drives the system toward re-equilibration (@reaction-gibbs-energy).

All equilibrium shifts in thermodynamic systems can be interpreted as responses to changes in the conjugate intensive variables of the relevant thermodynamic potential.

Before the operational use of $G$ is developed in @chemical-potential, one concrete application of $dG = 0$ illustrates the framework: phase equilibrium is the first place where $G$ acts as a global potential compared across competing thermodynamic states, and it makes $C_p$ and $\Delta H_\text{trans}$ concrete along the way.


(phase-equilibrium)=
## Phase and phase equilibrium

This section extends the Gibbs framework to the comparison of multiple thermodynamic states.
Phase equilibrium is a global minimisation problem: different phases correspond to distinct thermodynamic manifolds, each with its own Gibbs energy surface, and the stable state is the one with the lowest $G$ at given $(T,P,\{n_i\})$.
Each phase corresponds to a distinct branch of $G$ as a function of composition; the stable state at any given composition is the phase, or mixture of phases, that minimises $G$ globally.
The goal here is not a full treatment of phase diagrams but the minimal variational structure needed to interpret coexistence as equality of chemical potentials.

A **phase** is a homogeneous region of matter with uniform composition and physical state, separated from other regions by a phase boundary.
At equilibrium between two phases $\alpha$ and $\beta$, matter can transfer across the boundary; setting $dG = 0$ yields the phase equilibrium criterion directly:

$$\mu_i(\alpha) = \mu_i(\beta) \quad \text{for all } i$$

When $\mu_i(\alpha) > \mu_i(\beta)$, species $i$ has a thermodynamic driving force to move from phase $\alpha$ to phase $\beta$; the system evolves until the potentials equalise.

The **common tangent construction** gives the criterion a geometric meaning: a straight line simultaneously tangent to $G(x)$ at two compositions $x^\alpha$ and $x^\beta$ identifies the coexistence compositions (@fig-common-tangent).
Tangency at both points is the geometric statement of $\mu_i(\alpha) = \mu_i(\beta)$.
Any overall composition lying between $x^\alpha$ and $x^\beta$ separates into the two coexisting phases; a phase is stable where $G$ is locally convex, and loss of convexity signals spontaneous phase separation.


todo: update figure.
```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-common-tangent

import numpy as np
import matplotlib.pyplot as plt
from reactions.plots import setup_figure, COLORS
from scipy.optimize import brentq

x = np.linspace(1e-4, 1 - 1e-4, 2000)

A = 2.8  # interaction parameter; > 2 for phase separation
B = 0.25  # asymmetry


def G_f(xi):
    return xi * np.log(xi) + (1 - xi) * np.log(1 - xi) + A * xi * (1 - xi) + B * xi


def dG_f(xi):
    return np.log(xi) - np.log(1 - xi) + A * (1 - 2 * xi) + B


def residual(xa):
    m = dG_f(xa)
    try:
        xb = brentq(lambda t: dG_f(t) - m, xa + 0.15, 1 - 1e-3)
    except ValueError:
        return 1.0
    return (G_f(xb) - G_f(xa)) / (xb - xa) - m


xa = brentq(residual, 0.04, 0.44)
m_ct = dG_f(xa)
xb = brentq(lambda t: dG_f(t) - m_ct, xa + 0.15, 1 - 1e-3)

x_line = np.linspace(0, 1, 200)
y_line = G_f(xa) + m_ct * (x_line - xa)

x_mid = x[(x >= xa) & (x <= xb)]
G_mid = G_f(x_mid)
tan_mid = G_f(xa) + m_ct * (x_mid - xa)

fig, ax = setup_figure()
ax.plot(x, G_f(x), color="C0", lw=2.2)
ax.plot(x_line, y_line, color="C2", lw=1.5, ls="--")
ax.fill_between(x_mid, tan_mid, G_mid, alpha=0.12, color="C2")
ax.plot([xa, xb], [G_f(xa), G_f(xb)], "o", color="C2", ms=7, zorder=5)

ax.text(
    2 * (xa + xb) / 3,
    2 * (G_f(xa) + G_f(xb)) / 3,
    "two phases",
    ha="center",
    fontsize=9,
    color="C2",
    style="italic",
)

x_ex = 0.75
G_ex = G_f(x_ex)
tan_ex = G_f(xa) + m_ct * (x_ex - xa)
ax.plot(x_ex, G_ex, "o", color="C3", ms=7, zorder=6)
ax.plot(x_ex, tan_ex, "o", color="C3", ms=7, zorder=6, mfc="white", mew=1.5)
ax.annotate(
    "",
    xy=(x_ex, tan_ex),
    xytext=(x_ex, G_ex),
    arrowprops=dict(arrowstyle="<->", color="C3", lw=1.4),
)
ax.text(
    x_ex + 0.03,
    (G_ex + tan_ex) / 2,
    r"$\Delta G$",
    color="C3",
    fontsize=10,
    va="center",
)
ax.axvline(x_ex, color="C3", lw=0.8, ls=":", ymax=0.92)

# lever rule: two segments at tangent level, split at x_ex
f_beta = (x_ex - xa) / (xb - xa)  # fraction of phase beta
f_alpha = (xb - x_ex) / (xb - xa)  # fraction of phase alpha
ax.annotate(
    "",
    xy=(x_ex, tan_ex),
    xytext=(xa, tan_ex),
    arrowprops=dict(arrowstyle="<->", color="C3", lw=1.0, linestyle="dashed"),
)
ax.annotate(
    "",
    xy=(xb, tan_ex),
    xytext=(x_ex, tan_ex),
    arrowprops=dict(arrowstyle="<->", color="C3", lw=1.0, linestyle="dashed"),
)
ax.text(
    (xa + x_ex) / 2,
    tan_ex + 0.05,
    rf"$f^\beta={f_beta:.2f}$",
    ha="center",
    va="top",
    color="C3",
    fontsize=8.5,
)
ax.text(
    (x_ex + xb) / 2,
    tan_ex + 0.05,
    rf"$f^\alpha={f_alpha:.2f}$",
    ha="center",
    va="top",
    color="C3",
    fontsize=8.5,
)
ax.axvline(xa, color="C3", lw=0.8, ls=":", ymax=0.92)
ax.axvline(xb, color="C3", lw=0.8, ls=":", ymax=0.92)

ax.set_xlabel("Composition $x$")
ax.set_ylabel("Molar Gibbs energy $G$")
ax.set_xlim(-0.02, 1.02)
ax.set_xticks([0, xa, xb, 1])
ax.set_xticklabels(["0", r"$x^\alpha$", r"$x^\beta$", "1"])
ax.set_yticks([])
fig.tight_layout()
```

```{figure} #cell-common-tangent
:name: fig-common-tangent

Common tangent construction for a binary mixture with a tendency to phase-separate.
The molar Gibbs energy $G(x)$ has a non-convex region between $x^\alpha$ and $x^\beta$.
The common tangent (dashed) simultaneously touches both coexistence compositions; tangency at both points is the geometric statement of $\mu_i(\alpha) = \mu_i(\beta)$.
Any overall composition in the shaded region splits into the two coexisting phases.
```

**Clausius-Clapeyron equation.**
Differentiating $\mu^\alpha = \mu^\beta$ along the coexistence curve, and using $d\mu = -(S/n)\,dT + (V/n)\,dP$ for each phase (from $dG = -S\,dT + V\,dP$ divided by $n$), equating across the boundary and rearranging yields:

$$
\frac{dP}{dT} = \frac{\Delta S}{\Delta V}
$$

Using $\Delta S = \Delta H_\text{trans}/T$ at the reversible coexistence condition gives the Clausius–Clapeyron form:

$$\frac{dP}{dT} = \frac{\Delta H_\text{trans}}{T \cdot \Delta V}$$

where $\Delta H_\text{trans}$ is the enthalpy absorbed at constant pressure when crossing the phase boundary (the latent heat) and $\Delta V$ is the volume change.

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-latent-heat

import numpy as np
import matplotlib.pyplot as plt
from reactions.plots import setup_figure, COLORS

Cp_liq = 4.18  # kJ/(kg·K)
Cp_vap = 2.01  # kJ/(kg·K)
dH_trans = 2257.0  # kJ/kg
T_start = 20.0  # °C
T_boil = 100.0  # °C

Q_heat = Cp_liq * (T_boil - T_start)

Q1 = np.linspace(0, Q_heat, 200)
T1 = T_start + Q1 / Cp_liq

Q2 = np.linspace(Q_heat, Q_heat + dH_trans, 200)
T2 = np.full_like(Q2, T_boil)

Q3 = np.linspace(Q_heat + dH_trans, Q_heat + dH_trans + 250, 80)
T3 = T_boil + (Q3 - (Q_heat + dH_trans)) / Cp_vap

Q_all = np.concatenate([Q1, Q2, Q3])
T_all = np.concatenate([T1, T2, T3])

fig, ax = setup_figure()
ax.plot(Q_all, T_all, color="C0", lw=2.2)

ax.annotate(
    "",
    xy=(Q_heat, T_boil + 4),
    xytext=(0, T_boil + 4),
    arrowprops=dict(arrowstyle="<->", color="C1", lw=1.4),
)
ax.text(
    Q_heat / 2,
    T_boil + 9,
    r"$\int C_p\,dT \approx 335\ \mathrm{kJ}$",
    ha="center",
    va="bottom",
    color="C1",
    fontsize=9.5,
)

ax.annotate(
    "",
    xy=(Q_heat + dH_trans, T_boil + 4),
    xytext=(Q_heat, T_boil + 4),
    arrowprops=dict(arrowstyle="<->", color="C2", lw=1.4),
)
ax.text(
    Q_heat + dH_trans / 2,
    T_boil + 9,
    r"$\Delta H_\text{vap} \approx 2260\ \mathrm{kJ}$",
    ha="center",
    va="bottom",
    color="C2",
    fontsize=9.5,
)

ax.axhline(T_boil, color="gray", lw=0.8, linestyle=":")
ax.set_xlabel(r"Heat added $Q$ (kJ per kg)")
ax.set_ylabel(r"Temperature ($^\circ$C)")
ax.set_xlim(-80, Q_all[-1] + 80)
ax.set_ylim(0, 145)
fig.tight_layout()
```

```{figure} #cell-latent-heat
:name: fig-latent-heat

Heat added vs temperature for 1 kg of water.
During liquid heating the slope is $C_p$; at 100 °C the temperature plateaus while the enthalpy of vaporisation $\Delta H_\text{vap} \approx 2260\ \mathrm{kJ/kg}$ is absorbed (roughly seven times the energy of the heating segment).
The Clausius-Clapeyron equation converts this large latent heat into a steep coexistence slope: water boils near 70 °C at altitude and above 100 °C in a pressure cooker.
```

---

The next chapter develops what $G$ tells us in practice: the criterion for spontaneity and the per-species structure through chemical potential.
