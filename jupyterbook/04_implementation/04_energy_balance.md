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
In practice, reactions release or absorb heat that immediately changes the fluid temperature, which in turn shifts $K(T)$ and $k^f(T)$ and changes the rate.
Modelling this thermal feedback requires treating $T$ as a dynamic state variable governed by the energy balance.


## Energy balance for a reactive fluid

Consider a closed, well-mixed liquid-phase system at constant pressure with no shaft work.
For liquid-phase systems, density variations and expansion work are negligible, so $\rho$ and $C_p$ are treated as constant.
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
Equating with {eq}`eq-first-law-H`, dividing by $V\,dt$, and using the reaction flux $\varphi_j = (1/V)\,(d\xi_j/dt)$ (@kinetics):

$$
\rho C_p\,\dot{T}
= -\sum_j \Delta_r H_j\,\varphi_j(\mathbf{c}, T) + \dot{Q}_\text{ext},
$$ (eq-energy-balance)

where $\rho C_p = C_{p,\text{tot}}/V$ is the volumetric heat capacity [J/(m³·K)] and $\dot{Q}_\text{ext} = \delta Q_\text{ext}/(V\,dt)$ is the external heat input per unit volume [W/m³].
A system is **adiabatic** when it exchanges no heat with its surroundings; for such a system, $\dot{Q}_\text{ext} = 0$.

```{admonition} Implementation note
:class: note

$\Delta_r H_j$ is approximated by the standard-state value $\Delta_r H^\circ_j(T)$ from the equilibrium constant model (@equilibrium-temperature): for `EquilibriumConstantVantHoff`, this is the constant $\Delta H^\circ$; for `EquilibriumConstantVantHoffCp`, it includes the Kirchhoff correction $\Delta H^\circ + \Delta C_p(T - T_\text{ref})$.
This approximation is accurate for dilute liquid-phase systems where activity corrections to the reaction enthalpy are small.
```

An exothermic reaction ($\Delta_r H^\circ < 0$) raises the temperature, which shifts $K(T)$ and feeds back into the reaction rate: Le Chatelier's principle expressed dynamically.

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
    molar_mass=0.018,  # kg/mol
    density=1000.0,  # kg/m³
    heat_capacity=75.3,  # J/(mol·K)
)

dH_rxn = -20e3  # J/mol  (exothermic)
dS_rxn = -50.0  # J/(mol·K)
Ea = 40e3  # J/mol
A_arr = 1e10  # mol/(m³·s)
T0 = 298.15  # K
c_tot = 1000.0  # mol/m³

K_vH = EquilibriumConstantVantHoff(dH=dH_rxn, dS=dS_rxn)
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

Passing `solvent_composition` to `simulate` triggers the energy balance: $T$ is appended to the state vector and the result gains a `result["T"]` array alongside the concentration trajectories.

```{code-cell} ipython3
result = simulate(
    model,
    c0={"A": c_tot, "B": 0.0},
    t_span=(0, 5.0),
    T=T0,
    solvent_composition={"water": 1.0},
)

print(f"T(0)   = {result['T'].values[0]:.3f} K")
print(f"T(end) = {result['T'].values[-1]:.3f} K")
print(f"c_B(end) = {result['c'].sel(species='B').values[-1]:.2f} mol/m³")
```

The temperature rises by a few kelvin as the exothermic reaction produces B, then plateaus once equilibrium is reached (@fig-adiabatic).


## Adiabatic simulation and conservation check

For a single reaction $\ce{A <=> B}$ with $\nu_B = +1$ and $\dot{Q}_\text{ext} = 0$, {eq}`eq-energy-balance` integrates to the adiabatic invariant

$$
\rho C_p\,T + \Delta_r H^\circ\,c_B = \text{const}.
$$ (eq-adiabatic-invariant)

An exothermic reaction ($\Delta_r H^\circ < 0$) decreases the second term as $c_B$ grows, so $T$ must rise to maintain the sum.
The solver should preserve {eq}`eq-adiabatic-invariant` to within the specified tolerance:

```{code-cell} ipython3
invariant = rho_cp * result["T"].values + dH_rxn * result["c"].sel(species="B").values
rel_var = (invariant.max() - invariant.min()) / abs(invariant[0])
print(f"Invariant range:      {invariant.min():.4f} to {invariant.max():.4f} J/m³")
print(f"Relative variation:   {rel_var:.2e}")
```

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-adiabatic

import matplotlib.pyplot as plt
from reactions.plots import setup_figure, COLORS

fig, axes = setup_figure(1, 2)

axes[0].plot(result.coords["time"], result["c"].sel(species="A"), color="C0", label="A")
axes[0].plot(result.coords["time"], result["c"].sel(species="B"), color="C1", label="B")
axes[0].set_xlabel("time [s]")
axes[0].set_ylabel(r"concentration [mol/m³]")
axes[0].legend()

axes[1].plot(result.coords["time"], result["T"] - 273.15, color="C3")
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

To make the thermal effect visible, the comparison uses acetonitrile ($\rho C_p \approx 1.75 \times 10^6\ \mathrm{J/(m^3\cdot K)}$), whose lower heat capacity amplifies the adiabatic temperature rise to $\approx 10\ \mathrm{K}$, compared to $\approx 4\ \mathrm{K}$ in water.
A prescribed-temperature simulation at $T_0$ then overshoots the true adiabatic endpoint by a clearly visible margin (@fig-comparison).

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-comparison

_MECN = Species(
    name="MeCN",
    molar_mass=0.041,
    density=786.0,
    heat_capacity=91.5,
)
_model_cmp = ReactionModel(
    components=[
        Component("A"),
        Component("B"),
        Component("MeCN", species=[_MECN]),
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

result_adi_cmp = simulate(
    _model_cmp,
    c0={"A": c_tot, "B": 0.0},
    t_span=(0, 5.0),
    T=T0,
    solvent_composition={"MeCN": 1.0},
)
result_iso_cmp = simulate(
    _model_cmp,
    c0={"A": c_tot, "B": 0.0},
    t_span=(0, 5.0),
    T=T0,
)

K_T0 = K_vH.K(T0)
T_final_cmp = float(result_adi_cmp["T"].values[-1])
c_B_iso_eq = c_tot * K_T0 / (1 + K_T0)
c_B_adi_eq = c_tot * K_vH.K(T_final_cmp) / (1 + K_vH.K(T_final_cmp))

C_ADI = "C0"
C_ISO = "C1"

fig, axes = setup_figure(1, 2)

axes[0].plot(result_iso_cmp.coords["time"], result_iso_cmp["c"].sel(species="B"), color=C_ISO, ls="--", label="isothermal $T_0$")
axes[0].plot(result_adi_cmp.coords["time"], result_adi_cmp["c"].sel(species="B"), color=C_ADI, label="adiabatic")
axes[0].axhline(c_B_iso_eq, color=C_ISO, lw=0.8, ls="--", alpha=0.5)
axes[0].axhline(c_B_adi_eq, color=C_ADI, lw=0.8, ls=":", alpha=0.5)
axes[0].text(
    0.98, c_B_iso_eq, r"$c_{B,\mathrm{eq}}^{\mathrm{iso}}$",
    transform=axes[0].get_yaxis_transform(),
    ha="right", va="bottom", fontsize=11, color=C_ISO,
)
axes[0].text(
    0.98, c_B_adi_eq, r"$c_{B,\mathrm{eq}}^{\mathrm{adi}}$",
    transform=axes[0].get_yaxis_transform(),
    ha="right", va="top", fontsize=11, color=C_ADI,
)
axes[0].set_xlabel("time [s]")
axes[0].set_ylabel(r"$c_B$ [mol/m³]")
axes[0].legend()

axes[1].plot(
    result_iso_cmp.coords["time"],
    np.full_like(result_iso_cmp.coords["time"].values, T0) - 273.15,
    color=C_ISO,
    ls="--",
    label="isothermal $T_0$",
)
axes[1].plot(result_adi_cmp.coords["time"], result_adi_cmp["T"] - 273.15, color=C_ADI, label="adiabatic")
axes[1].set_xlabel("time [s]")
axes[1].set_ylabel(r"$T$ [°C]")
axes[1].legend()

fig.tight_layout()
```

```{figure} #cell-comparison
:name: fig-comparison

Isothermal (orange, dashed) vs adiabatic (blue, solid) simulation in acetonitrile ($\rho C_p \approx 1.75 \times 10^6\ \mathrm{J/(m^3\cdot K)}$); colors are consistent across both panels.
Left: the isothermal model overestimates the long-time $c_B$ because it uses $K(T_0)$ throughout; the adiabatic model reaches a lower equilibrium as the rising temperature shifts $K$ via Le Chatelier.
Faint horizontal lines mark the respective equilibrium concentrations; the gap between them is the error from prescribing temperature in a self-heating system.
Right: the isothermal model prescribes a flat temperature; the adiabatic model resolves the $\approx 10\ \mathrm{K}$ rise, larger than in water because of acetonitrile's lower heat capacity.
```

The gap between the two endpoints is the error introduced by prescribing temperature in a system that self-heats.
For reactions with large $|\Delta_r H^\circ|$ or low $\rho C_p$ (organic solvents, gas-phase), this error grows substantially.


## Gradient solvent composition

`solvent_composition` values may be callables $x_k(t) \to \mathbb{R}$ to represent a time-varying mixture, for example a water/organic gradient in preparative liquid chromatography.
Each callable is evaluated at the current integration time $t$:

```{code-cell} ipython3
MECN = Species(
    name="MeCN",
    molar_mass=0.041,  # kg/mol
    density=786.0,  # kg/m³
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
        "MeCN": lambda t: 0.4 * t / t_end,
    },
)
```

The resulting concentration and temperature profiles are shown in @fig-gradient.

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-gradient

x_MeCN = np.array([0.4 * t / t_end for t in result_grad.coords["time"].values])


def _rho_cp_mix(x_mecn: float) -> float:
    c = np.zeros(len(model_mix.species))
    c[model_mix.species_index["water"]] = 1.0 - x_mecn
    c[model_mix.species_index["MeCN"]] = x_mecn
    return model_mix.volumetric_heat_capacity(c)


rho_cp_t = np.array([_rho_cp_mix(x) for x in x_MeCN])

fig, axes = setup_figure(1, 2)

axes[0].plot(result_grad.coords["time"], result_grad["c"].sel(species="A"), color="C0", label="A")
axes[0].plot(result_grad.coords["time"], result_grad["c"].sel(species="B"), color="C1", label="B")
axes[0].set_xlabel("time [s]")
axes[0].set_ylabel(r"concentration [mol/m³]")
axes[0].legend()
axes[0].set_title("Concentrations (gradient solvent)")

ax2 = axes[1].twinx()
axes[1].plot(result_grad.coords["time"], result_grad["T"] - 273.15, color="C3", label=r"$T$")
ax2.plot(result_grad.coords["time"], rho_cp_t / 1e6, color="gray", ls="--", label=r"$\rho C_p$")
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

`solvent_composition` cannot be combined with a callable $T(t)$ because both prescribe the thermal state externally.
Passing both raises `ValueError`.


## Temperature dependence of the reaction flux

Appending $T$ to the concentration state vector requires the partial derivative $\partial\varphi_j/\partial T$ for each reaction:

$$
\frac{\partial\varphi_j}{\partial T}
= \frac{d\ln k^f}{dT}\,\varphi_j
+ k^r\,\frac{d\ln K}{dT}\,P^r_j.
$$ (eq-dvarphi-dT)

The first term scales the current flux by the logarithmic rate sensitivity; the second captures the shift in equilibrium.
For `RateConstantArrhenius` and `EquilibriumConstantVantHoff`, both derivatives are analytic (@kinetics-temperature, @equilibrium-temperature).


## Extended Jacobian structure

With the state vector extended to $\mathbf{y} = [\mathbf{c},\,T]^\intercal$, and under the standard-state approximation ($\Delta_r H_j \approx \Delta_r H^\circ_j$) and adiabatic assumption ($\dot{Q}_\text{ext} = 0$), the full residual is

$$
\mathbf{F}(\mathbf{y},\dot{\mathbf{y}}) =
\begin{pmatrix}
\dot{\mathbf{c}} - \mathbf{S}\boldsymbol{\varphi}(\mathbf{c},T) \\
\dot{T} + \bigl(\sum_j \Delta_r H^\circ_j\,\varphi_j\bigr)\!/\rho C_p
\end{pmatrix}.
$$

The Jacobian retains block structure:

$$
J =
\begin{pmatrix}
-\mathbf{S}\,\partial\boldsymbol{\varphi}/\partial\mathbf{c}
  & -\mathbf{S}\,\partial\boldsymbol{\varphi}/\partial T \\
\partial F_T/\partial\mathbf{c}
  & \partial F_T/\partial T
\end{pmatrix}.
$$

All blocks are analytic for structured models and are exposed through the corresponding `ReactionModel` derivative methods.
`EquilibriumConstantCustom` and `EquilibriumConstantTabulated` substitute finite-difference estimates where needed; the remaining blocks stay analytic.


```{admonition} Scope: aqueous activity coefficients
:class: warning

`volumetric_heat_capacity` correctly handles solvent mixtures via {eq}`eq-rhocp`.
The activity coefficient models (`ActivityCoefficientDebyeHuckel`, `ActivityCoefficientDavies`) are currently calibrated for pure water at 25 °C; the Debye-Hückel parameter $A$ depends on the solvent dielectric constant $\varepsilon_r$, which changes with composition.
Mixed-solvent activity corrections require a dielectric mixing model $\varepsilon(x_k)$ that is not yet implemented.
See the discussion in `LIBRARY.md` under "Known limitations".
```

---

Temperature as a dynamic state requires only solvent species with complete physical fields (`molar_mass`, `density`, `heat_capacity`) and reactions that expose $\Delta_r H^\circ(T)$ via `reaction_enthalpy`.
The solver interface is otherwise unchanged; `simulate` simply adds a `result["T"]` array to the returned `xr.Dataset`.
The following chapter extends the activity term $a_i = \gamma_i c_i / c^\circ$ to non-ideal solutions, where ionic-strength corrections shift the apparent equilibrium composition (@implementation-activity).
