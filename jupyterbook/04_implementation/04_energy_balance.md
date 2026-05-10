---
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

(implementation-energy-balance)=
# Coupled Energy Balance

```{admonition} Optional chapter
:class: note

This chapter extends the solver to treat temperature as a dynamic state variable.
It is not required for @implementation-activity or the acid-base chapters that follow.
```

The kinetic and equilibrium frameworks in @implementation-kinetics model temperature as either a fixed parameter or a prescribed function $T(t)$ supplied externally.
Both require the temperature programme to be known before the simulation begins.
In practice, reactions release or absorb heat that immediately changes the fluid temperature, which in turn shifts $K(T)$ and $k_f(T)$ and changes the rate.
Modelling this thermal feedback requires treating $T$ as a dynamic state variable governed by the energy balance.


## Energy balance for a reactive fluid

Consider a closed, well-mixed liquid-phase system at constant pressure with no shaft work.
For liquid-phase systems, density and pressure variations are small enough to neglect: $\rho$ and $C_p$ are treated as constant over the integration, and expansion work $P\,dV$ is negligible.
Under these conditions, the first law reduces to (@laws-of-thermodynamics, @thermodynamic-potentials):

$$
dH = \delta Q_\text{ext}.
$$ (eq-first-law-H)

The total enthalpy $H$ depends on both temperature and composition.
Its total differential is

$$
dH = C_{p,\text{tot}}\,dT + \sum_i \bar{H}_i\,dn_i,
$$

where $C_{p,\text{tot}} = \sum_i n_i C_{p,i}$ is the total heat capacity and $\bar{H}_i$ is the partial molar enthalpy of species $i$.
The composition changes are constrained by the reaction stoichiometry: $dn_i = \sum_j \nu_{ij}\,d\xi_j$ (@reaction-coordinates).
Substituting and collecting over reactions:

$$
dH = C_{p,\text{tot}}\,dT + \sum_j \underbrace{\left(\sum_i \nu_{ij}\,\bar{H}_i\right)}_{\Delta_r H_j}\,d\xi_j,
$$

where $\Delta_r H_j = \sum_i \nu_{ij} \bar{H}_i$ is the reaction enthalpy under current conditions.
Setting equal to {eq}`eq-first-law-H`, dividing by the fluid volume $V$ and by $dt$, and using the reaction flux $\varphi_j = (1/V)\,(d\xi_j/dt)$ (@kinetics):

$$
\rho C_p\,\dot{T}
= -\sum_j \Delta_r H_j\,\varphi_j(\mathbf{c}, T) + \dot{Q}_\text{ext},
$$ (eq-energy-balance)

where $\rho C_p = C_{p,\text{tot}}/V$ is the volumetric heat capacity [J/(m³·K)] and $\dot{Q}_\text{ext} = \delta Q_\text{ext}/(V\,dt)$ is the external heat input per unit volume [W/m³].
For an adiabatic system, $\dot{Q}_\text{ext} = 0$.

In the current implementation, $\Delta_r H_j$ is parameterised using the standard-state quantity $\Delta_r H^\circ_j(T)$ from the equilibrium constant model (@equilibrium-temperature): for `EquilibriumConstantVantHoff`, this is the constant $\Delta H^\circ$; for `EquilibriumConstantVantHoffCp`, it includes the Kirchhoff correction $\Delta H^\circ + \Delta C_p(T - T_\text{ref})$.
This approximation is accurate for dilute liquid-phase systems where activity corrections to the reaction enthalpy are small.

An exothermic reaction ($\Delta_r H^\circ < 0$) acts as a heat source and raises the fluid temperature.
The rising temperature then decreases $K(T)$, reducing the forward driving force — Le Chatelier's principle expressed as a coupled ODE.

```{admonition} Intuition: thermal response scale
:class: tip

The energy balance says reaction enthalpy acts as a source term for temperature.
A rough estimate of the temperature change when a concentration $\Delta c$ of product forms:

$$
\Delta T \sim \frac{|\Delta_r H^\circ|\,\Delta c}{\rho C_p}.
$$

A large $|\Delta_r H^\circ|$ or small $\rho C_p$ means a large thermal response.
Water has an unusually high $\rho C_p \approx 4.2 \times 10^6\ \mathrm{J/(m^3 \cdot K)}$, which is why aqueous reactions change temperature slowly.
Organic solvents and gas-phase systems respond much more strongly to the same reaction heat.
```

For a dilute solution, $\rho C_p$ is determined by the solvent.
Under the dilute approximation (solute molar volumes negligible), it is computed from the solvent species:

$$
\rho C_p = \sum_k x_k\,\frac{\rho_k}{M_k}\,C_{p,k},
$$ (eq-rhocp)

where $x_k$ is the mole fraction of solvent $k$, $\rho_k$ [kg/m³] and $M_k$ [kg/mol] are its density and molar mass, and $C_{p,k}$ [J/(mol·K)] is its molar heat capacity.
The ratio $\rho_k / M_k$ is the molar concentration of pure solvent [mol/m³]; multiplied by $x_k C_{p,k}$ it gives the volumetric contribution.
For pure water ($\rho = 1000\ \mathrm{kg/m^3}$, $M = 0.018\ \mathrm{kg/mol}$, $C_p = 75.3\ \mathrm{J/(mol\cdot K)}$), {eq}`eq-rhocp` gives $\rho C_p \approx 4.18 \times 10^6\ \mathrm{J/(m^3\cdot K)}$.


## Temperature dependence of the reaction flux

Appending $T$ to the concentration state vector requires the partial derivative $\partial\varphi_j/\partial T$ for each reaction.
Splitting the flux using $k_r = k_f/K$, with the activity products $P_\text{fwd}$ and $P_\text{bwd}$ independent of $T$:

$$
\varphi_j = k_f(T)\,P_\text{fwd} - k_r(T)\,P_\text{bwd}.
$$

Differentiating with respect to $T$ and collecting terms:

$$
\frac{\partial\varphi_j}{\partial T}
= \frac{d\ln k_f}{dT}\,\varphi_j
+ k_r\,\frac{d\ln K}{dT}\,P_\text{bwd}.
$$ (eq-dvarphi-dT)

The first term scales the current flux by the logarithmic rate sensitivity; the second is a correction from the shifting equilibrium.
For `RateConstantArrhenius` and `EquilibriumConstantVantHoff`, both $d\ln k_f/dT = E_a/(RT^2)$ and $d\ln K/dT = \Delta_r H^\circ/(RT^2)$ are analytic (@kinetics-temperature, @equilibrium-temperature).
`EquilibriumConstantCustom` and `EquilibriumConstantTabulated` return `None` for these derivatives; the library substitutes a finite-difference estimate automatically.


## Extended Jacobian structure

With the state vector extended to $\mathbf{y} = [\mathbf{c},\,T]^\intercal$, the full residual is

$$
\mathbf{F}(\mathbf{y},\dot{\mathbf{y}}) =
\begin{pmatrix}
\dot{\mathbf{c}} - \mathbf{S}\boldsymbol{\varphi}(\mathbf{c},T) \\
\dot{T} + \bigl(\sum_j \Delta_r H^\circ_j\,\varphi_j\bigr)\!/\rho C_p
\end{pmatrix}.
$$

The Jacobian $\partial\mathbf{F}/\partial\mathbf{y}$ has block structure:

$$
J =
\begin{pmatrix}
-\mathbf{S}\,\partial\boldsymbol{\varphi}/\partial\mathbf{c}
  & -\mathbf{S}\,\partial\boldsymbol{\varphi}/\partial T \\
\partial F_T/\partial\mathbf{c}
  & \partial F_T/\partial T
\end{pmatrix}.
$$

The top-left block and the top-right column are analytic, computed by `ReactionModel.jacobian()` and `ReactionModel.jacobian_dT()`.
The bottom row is also analytic for structured models: $\partial F_T/\partial c_k = \sum_j (\Delta_r H_j/\rho C_p)\,\partial\varphi_j/\partial c_k$ from `net_rate_jac`, and $\partial F_T/\partial T = -\sum_j (\Delta_r H_j\,\partial\varphi_j/\partial T + \varphi_j\,\mathrm{d}\Delta_r H_j/\mathrm{d}T)/\rho C_p$ from `net_rate_dT` and `d_reaction_enthalpy_dT`.
`EquilibriumConstantCustom` and `EquilibriumConstantTabulated` substitute finite-difference estimates for $\Delta_r H_j$ and $\mathrm{d}\Delta_r H_j/\mathrm{d}T$; the remaining blocks stay analytic.


## Specifying solvent properties

The energy balance requires $\rho C_p$ computed from the solvent species.
A species participates in the heat-capacity calculation when all three physical fields are set: `molar_mass` [kg/mol], `density` [kg/m³], and `heat_capacity` [J/(mol·K)].
All three must be provided on at least one species; the library raises `ValueError` if none are.

```{code-cell} ipython3
import numpy as np
from reactions.api import (
    Component,
    Species,
    EquilibriumConstantVantHoff,
    RateConstantArrhenius,
    ThermodynamicReaction,
    ReactionModel,
)
from reactions.solver import simulate

WATER = Species(
    name="water",
    molar_mass=0.018,    # kg/mol
    density=1000.0,      # kg/m³
    heat_capacity=75.3,  # J/(mol·K)
)

dH_rxn = -20e3   # J/mol  (exothermic)
dS_rxn = -50.0   # J/(mol·K)
Ea     = 40e3    # J/mol
A_arr  = 1e10    # mol/(m³·s)
T0     = 298.15  # K
c_tot  = 1000.0  # mol/m³

K_vH  = EquilibriumConstantVantHoff(dH=dH_rxn, dS=dS_rxn)
kf_arr = RateConstantArrhenius(A=A_arr, Ea=Ea)

model = ReactionModel(
    components=[
        Component("A"),
        Component("B"),
        Component("water", species=[WATER]),
    ],
    reactions=[
        ThermodynamicReaction(
            "A <-> B",
            mode="kinetic",
            equilibrium_constant=K_vH,
            rate_constant=kf_arr,
        ),
    ],
)
```

`volumetric_heat_capacity` computes $\rho C_p$ from {eq}`eq-rhocp`:

```{code-cell} ipython3
c_pure = np.zeros(len(model.species))
c_pure[model.species_index["water"]] = 1000.0
rho_cp = model.volumetric_heat_capacity(c_pure)
print(f"ρCp = {rho_cp:.3e} J/(m³·K)")
```

Passing `solvent_composition` to `simulate` triggers the energy balance: $T$ is appended to the state vector and the result gains a `T_profile` array alongside the concentration trajectories.

```{code-cell} ipython3
result = simulate(
    model,
    c0={"A": c_tot, "B": 0.0},
    t_span=(0, 5.0),
    T=T0,
    solvent_composition={"water": 1.0},
)

print(f"T(0)   = {result.T_profile[0]:.3f} K")
print(f"T(end) = {result.T_profile[-1]:.3f} K")
print(f"c_B(end) = {result['B'][-1]:.2f} mol/m³")
```

The temperature rises by a few kelvin as the exothermic reaction produces B, then plateaus once equilibrium is reached.


## Adiabatic simulation and conservation check

For a single reaction $\ce{A <=> B}$ with $\nu_B = +1$ and $\dot{Q}_\text{ext} = 0$, {eq}`eq-energy-balance` integrates to the adiabatic invariant

$$
\rho C_p\,T + \Delta_r H^\circ\,c_B = \text{const}.
$$ (eq-adiabatic-invariant)

An exothermic reaction ($\Delta_r H^\circ < 0$) decreases the second term as $c_B$ grows, so $T$ must rise to maintain the sum.
The solver should preserve {eq}`eq-adiabatic-invariant` to within the specified tolerance:

```{code-cell} ipython3
invariant = rho_cp * result.T_profile + dH_rxn * result["B"]
rel_var = (invariant.max() - invariant.min()) / abs(invariant[0])
print(f"Invariant range:      {invariant.min():.4f} to {invariant.max():.4f} J/m³")
print(f"Relative variation:   {rel_var:.2e}")
```

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-adiabatic

import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(8, 3.8))

axes[0].plot(result.t, result["A"], color="C0", label="A")
axes[0].plot(result.t, result["B"], color="C1", label="B")
axes[0].set_xlabel("time [s]")
axes[0].set_ylabel(r"concentration [mol/m³]")
axes[0].legend()

axes[1].plot(result.t, result.T_profile - 273.15, color="C3")
axes[1].set_xlabel("time [s]")
axes[1].set_ylabel(r"$T$ [°C]")

fig.tight_layout()
```

```{figure} #cell-adiabatic
:name: fig-adiabatic

Adiabatic simulation of $\ce{A <=> B}$ with $\Delta_r H^\circ = -20\ \mathrm{kJ/mol}$ in water ($\rho C_p \approx 4.18 \times 10^6\ \mathrm{J/(m^3\cdot K)}$, $E_a = 40\ \mathrm{kJ/mol}$, $T_0 = 298.15\ \mathrm{K}$).
Left: concentrations; the exothermic reaction produces B while heating the fluid.
Right: temperature rises approximately $4\ \mathrm{K}$ as the reaction proceeds and levels off at equilibrium.
```

A prescribed-temperature simulation at $T_0$ converges to a different equilibrium: without thermal feedback, the solver uses the fixed $K(T_0)$ throughout and overshoots the true adiabatic endpoint.

```{code-cell} ipython3
result_iso = simulate(
    model,
    c0={"A": c_tot, "B": 0.0},
    t_span=(0, 5.0),
    T=T0,
)
```

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-comparison

fig, axes = plt.subplots(1, 2, figsize=(8, 3.8))

K_T0    = K_vH.K(T0)
T_final = float(result.T_profile[-1])
K_Tf    = K_vH.K(T_final)
c_B_iso_eq = c_tot * K_T0 / (1 + K_T0)
c_B_adi_eq = c_tot * K_Tf / (1 + K_Tf)

axes[0].plot(result_iso.t, result_iso["B"], color="C1", ls="--", label="isothermal $T_0$")
axes[0].plot(result.t,     result["B"],     color="C1",           label="adiabatic")
axes[0].axhline(c_B_iso_eq, color="C1", lw=0.8, ls=":")
axes[0].axhline(c_B_adi_eq, color="C1", lw=0.8, ls=":")
axes[0].set_xlabel("time [s]")
axes[0].set_ylabel(r"$c_B$ [mol/m³]")
axes[0].legend()

axes[1].plot(result_iso.t, np.full_like(result_iso.t, T0) - 273.15,
             color="C3", ls="--", label="isothermal $T_0$")
axes[1].plot(result.t, result.T_profile - 273.15, color="C3", label="adiabatic")
axes[1].set_xlabel("time [s]")
axes[1].set_ylabel(r"$T$ [°C]")
axes[1].legend()

fig.tight_layout()
```

```{figure} #cell-comparison
:name: fig-comparison

Isothermal (dashed) vs adiabatic (solid) simulation of the same reaction.
Left: the isothermal model overestimates the long-time $c_B$ because it uses $K(T_0)$ throughout; the adiabatic model reaches a lower equilibrium as the rising temperature shifts $K$ via Le Chatelier.
Dotted horizontal lines mark the respective equilibrium concentrations.
Right: the isothermal model prescribes a flat temperature; the adiabatic model resolves the $\approx 4\ \mathrm{K}$ rise driven by the reaction heat.
```

The gap between the two endpoints is the error introduced by prescribing temperature in a system that self-heats.
For reactions with large $|\Delta_r H^\circ|$ or low $\rho C_p$ (organic solvents, gas-phase), this error grows substantially.


## Gradient solvent composition

`solvent_composition` values may be callables $x_k(t) \to \mathbb{R}$ to represent a time-varying mixture, for example a water/organic gradient in preparative liquid chromatography.
Each callable is evaluated at the current integration time $t$:

```{code-cell} ipython3
MECN = Species(
    name="MeCN",
    molar_mass=0.041,    # kg/mol
    density=786.0,       # kg/m³
    heat_capacity=91.5,  # J/(mol·K)
)

model_mix = ReactionModel(
    components=[
        Component("A"),
        Component("B"),
        Component("water", species=[WATER]),
        Component("MeCN", species=[MECN]),
    ],
    reactions=[
        ThermodynamicReaction(
            "A <-> B",
            mode="kinetic",
            equilibrium_constant=K_vH,
            rate_constant=kf_arr,
        ),
    ],
)

t_end = 5.0
result_grad = simulate(
    model_mix,
    c0={"A": c_tot, "B": 0.0},
    t_span=(0, t_end),
    T=T0,
    solvent_composition={
        "water": lambda t: 1.0 - 0.4 * t / t_end,
        "MeCN":  lambda t: 0.4 * t / t_end,
    },
)
```

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-gradient

x_MeCN = np.array([0.4 * t / t_end for t in result_grad.t])

def _rho_cp_mix(x_mecn: float) -> float:
    c = np.zeros(len(model_mix.species))
    c[model_mix.species_index["water"]] = 1.0 - x_mecn
    c[model_mix.species_index["MeCN"]] = x_mecn
    return model_mix.volumetric_heat_capacity(c)

rho_cp_t = np.array([_rho_cp_mix(x) for x in x_MeCN])

fig, axes = plt.subplots(1, 2, figsize=(9, 3.8))

axes[0].plot(result_grad.t, result_grad["A"], color="C0", label="A")
axes[0].plot(result_grad.t, result_grad["B"], color="C1", label="B")
axes[0].set_xlabel("time [s]")
axes[0].set_ylabel(r"concentration [mol/m³]")
axes[0].legend()
axes[0].set_title("Concentrations (gradient solvent)")

ax2 = axes[1].twinx()
axes[1].plot(result_grad.t, result_grad.T_profile - 273.15, color="C3", label=r"$T$")
ax2.plot(result_grad.t, rho_cp_t / 1e6, color="gray", ls="--", label=r"$\rho C_p$")
axes[1].set_xlabel("time [s]")
axes[1].set_ylabel(r"$T$ [°C]", color="C3")
ax2.set_ylabel(r"$\rho C_p$ [$\mathrm{MJ/(m^3 \cdot K)}$]", color="gray")
axes[1].set_title(r"Temperature and $\rho C_p$")

lines1, labels1 = axes[1].get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
axes[1].legend(lines1 + lines2, labels1 + labels2, fontsize=8)

fig.tight_layout()
```

```{figure} #cell-gradient
:name: fig-gradient

Adiabatic simulation with a linear water/MeCN ramp ($x_\text{MeCN}$: 0 → 0.4 over 5 s).
Left: concentrations; the reaction profile is similar to the pure-water case.
Right: temperature evolution and the time-varying $\rho C_p$ (dashed); as MeCN fraction increases, the lower heat capacity of acetonitrile reduces $\rho C_p$, so the same heat release produces a progressively larger temperature rise per unit volume.
```

`solvent_composition` cannot be combined with a callable $T(t)$: the two approaches prescribe the same degree of freedom.
Passing both raises `ValueError`.

```{admonition} Scope: aqueous activity coefficients
:class: warning

`volumetric_heat_capacity` correctly handles solvent mixtures via {eq}`eq-rhocp`.
The activity coefficient models (`ActivityCoefficientDebyeHuckel`, `ActivityCoefficientDavies`) are currently calibrated for pure water at 25 °C; the Debye-Hückel parameter $A$ depends on the solvent dielectric constant $\varepsilon$, which changes with composition.
Mixed-solvent activity corrections require a dielectric mixing model $\varepsilon(x_k)$ that is not yet implemented.
See the discussion in `LIBRARY.md` under "Known limitations".
```

---

Temperature as a dynamic state requires only solvent species with complete physical fields (`molar_mass`, `density`, `heat_capacity`) and reactions that expose $\Delta_r H^\circ(T)$ via `reaction_enthalpy`.
The solver interface is otherwise unchanged: `simulate` returns the same `SimulationResult` with an additional `T_profile` array.
The following chapter extends the activity term $a_i = \gamma_i c_i / c^\circ$ to non-ideal solutions, where ionic-strength corrections shift the apparent equilibrium composition (@implementation-activity).
