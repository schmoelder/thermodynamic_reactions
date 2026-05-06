---
kernelspec:
  name: python3
  display_name: Python 3 (ipykernel)
---

(chemical-potential)=
# Gibbs Energy and Chemical Potential

@thermodynamic-potentials established that the Gibbs free energy $G$ is minimised at equilibrium under constant temperature and pressure.
This chapter develops what $G$ actually tells us: first as a criterion for spontaneity, then through its structure as a sum of per-species contributions.

A **spontaneous process** is one that proceeds on its own without any external driving force.
Heat flows from hot to cold; gases expand into vacuum; salts dissolve in water.
The direction of spontaneous change is not determined by energy alone: dissolving ammonium nitrate in water absorbs heat yet happens spontaneously, while freezing water releases heat yet requires cooling below $0\,\degree\text{C}$.
What decides is the competition between energy and entropy, captured by $G$.

Spontaneity is not limited to chemical reactions.
Any process occurring at constant $T$ and $P$ (phase transitions, mixing, protein folding, membrane transport, electrochemistry) is spontaneous if and only if it lowers $G$.
This universality is what makes $G$ the central quantity of chemical thermodynamics.

One important limitation: $\Delta G < 0$ tells us that a process is thermodynamically possible, not that it will actually occur on any observable timescale.
Rate is a separate question, governed by kinetics and the activation energy barrier (see the note below).

## Spontaneity

From @thermodynamic-potentials, $G = H - TS$, so for a process at constant $T$ and $P$:

$$
\Delta G = \Delta H - T\Delta S
$$

A spontaneous process lowers $G$, so the criterion is $\Delta G < 0$.
This expression is linear in $T$: the intercept at $T = 0$ is $\Delta H$, and the slope is $-\Delta S$.
The sign of $\Delta G$ therefore depends on both the signs of $\Delta H$ and $\Delta S$ and on temperature, as shown in @fig-gibbs-helmholtz.

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-gibbs-helmholtz

import numpy as np
import matplotlib.pyplot as plt

T = np.linspace(0, 2.5, 300)

line_color = "#2c7bb6"
ds_color   = "#2ca25f"
dh_color   = "#d73027"
alpha_curv = 0.04   # small heat-capacity curvature term

fig, axes = plt.subplots(1, 2, figsize=(8, 3.8))

for ax, dH, subtitle in zip(
    axes,
    [1.0, -1.0],
    ["a) endothermic  ($\\Delta H > 0$)", "b) exothermic  ($\\Delta H < 0$)"],
):
    for dS, ds_label, va in [
        ( 0.4, r"$\Delta S > 0$", "top"),
        (-0.4, r"$\Delta S < 0$", "bottom"),
    ]:
        # Lines curve away from each other (heat-capacity effect)
        dG = dH - T * dS - np.sign(dS) * alpha_curv * T**2
        ax.plot(T, dG, color=line_color, lw=2.2)

        # Place label at T=1.2, well away from zero crossing
        t_lab = 1.2
        dG_lab = dH - t_lab * dS - np.sign(dS) * alpha_curv * t_lab**2
        offset = -0.18 if va == "top" else 0.18
        ax.text(t_lab, dG_lab + offset, ds_label,
                color=ds_color, fontsize=9.5, ha="center", va=va)

    sign = ">" if dH > 0 else "<"
    ax.annotate(rf"$\Delta H {sign} 0$",
                xy=(0, dH), xytext=(-0.35, dH),
                color=dh_color, fontsize=9.5, va="center", ha="right",
                annotation_clip=False)
    ax.plot(0, dH, "o", color=dh_color, ms=5, zorder=5)

    ax.spines["left"].set_position("zero")
    ax.spines["bottom"].set_position("zero")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim(-0.5, 2.8)
    ax.set_ylim(-2.2, 2.2)

    ax.text(2.72, 0.12, r"$T$", fontsize=11)
    ax.text(0.08, 2.1, r"$\Delta G$", fontsize=11)
    ax.text(0.5, -0.1, subtitle, transform=ax.transAxes,
            ha="center", fontstyle="italic", fontsize=10)

fig.tight_layout()
```

```{figure} #cell-gibbs-helmholtz
:name: fig-gibbs-helmholtz

$\Delta G = \Delta H - T\Delta S$ as a function of temperature for endothermic (left) and exothermic (right) reactions, each shown for positive and negative $\Delta S$.
The intercept at $T = 0$ is $\Delta H$; the slope is $-\Delta S$.
A reaction becomes spontaneous ($\Delta G < 0$) where its line crosses zero, at the crossover temperature $T^* = \Delta H / \Delta S$.
```

The four sign combinations identified in @fig-gibbs-helmholtz are summarised in @fig-spontaneity.

```{figure} figures/spontaneity_quadrant.png
:name: fig-spontaneity

The four combinations of $\Delta H$ and $\Delta S$ signs and their consequences for spontaneity.
```

```{admonition} Intuition
:class: tip

Four reactions illustrate the four sign combinations:
- $\Delta H < 0$, $\Delta S > 0$: decomposition of $\ce{N2O}$; spontaneous at all temperatures
- $\Delta H < 0$, $\Delta S < 0$: Haber process ($\ce{N2 + 3H2 -> 2NH3}$); spontaneous only below $T^*$
- $\Delta H > 0$, $\Delta S > 0$: calcination ($\ce{CaCO3 -> CaO + CO2}$); spontaneous only above $T^*$ (~840 °C)
- $\Delta H > 0$, $\Delta S < 0$: the reverse of the first case; no driving force exists at any temperature

The same logic applies to phase transitions (melting, freezing) and dissolution: $G = H - TS$ governs them all.
```

```{admonition} Thermodynamics and kinetics
:class: note

A negative $\Delta G$ establishes that a process is thermodynamically spontaneous, but says nothing about how fast it proceeds.
The rate is governed by the activation energy $E_a$ introduced in @maxwell-boltzmann: even a highly favourable reaction may be immeasurably slow if the energy barrier is large.
Thermodynamics predicts whether a process can occur; kinetics determines whether it does so on a relevant timescale.
This distinction is taken up in detail in the chapter on reaction kinetics.
```

## Chemical potential

From $G = \sum_i \mu_i n_i$ (@thermodynamic-potentials), the chemical potential of species $i$ is:

$$
\mu_i = \left.\frac{\partial G}{\partial n_i}\right|_{T,P,n_{j \neq i}}
$$

It is the free energy cost of adding one mole of species $i$ while holding temperature, pressure, and all other amounts fixed.
A species with high $\mu_i$ has a strong thermodynamic driving force to leave (react, diffuse, or precipitate); a species with low $\mu_i$ is stable where it is.

```{admonition} Intuition
:class: tip
$\mu_i$ acts as a chemical pressure: species flow from regions of high $\mu_i$ to regions of low $\mu_i$, just as fluids flow from high mechanical pressure to low.
This single driving force underlies diffusion across a membrane, partitioning between phases, and chemical reaction: all are matter moving downhill in $\mu$.
```

At constant $T$ and $P$, changes in composition change $G$ as:

$$
dG\big|_{T,P} = \sum_i \mu_i\,dn_i
$$

For this to be negative (spontaneous), matter must flow in the direction that lowers $G$.
If species $i$ has higher $\mu_i$ on one side of a membrane than the other, it diffuses toward the lower side.
If the chemical potentials of reactants exceed those of products, the reaction proceeds forward.
Equilibrium is reached when $dG = 0$: chemical potentials are balanced and no net driving force remains.

## The ideal chemical potential

The pressure dependence of $\mu$ follows directly from the differential of $G$.
At constant $T$ and fixed composition, $dG = V\,dP$, and since $G = n\mu$ for a pure substance:

$$
\left.\frac{\partial \mu}{\partial P}\right|_T = V_m
$$

where $V_m = V/n$ is the molar volume.
For an ideal gas, $V_m = RT/P$, giving:

$$
d\mu = \frac{RT}{P}\,dP
$$

Integrating from a reference pressure $P^\circ$ to $P$:

$$
\mu(T, P) = \mu^\circ(T) + RT\ln\frac{P}{P^\circ}
$$

The logarithm is a direct consequence of integrating $1/P$, which in turn comes from the ideal gas equation of state.
Higher pressure raises $\mu$: molecules at higher pressure have a stronger driving force to expand or react.
For real gases, $V_m \neq RT/P$ and the integral gives the fugacity $f$ rather than $P$, taken up in @nonidealities.


## Partial molar volume

The same pattern extends to mixtures.
The **partial molar volume** of species $i$ is:

$$
\bar{V}_i = \left.\frac{\partial V}{\partial n_i}\right|_{T,P,n_{j\neq i}}
$$

It is the volume change when an infinitesimal amount of $i$ is added at constant $T$, $P$, and composition.
For an ideal gas, $\bar{V}_i = V_m = RT/P$ for all components.
For liquids, $\bar{V}_i \approx 10^{-5}\,\text{m}^3/\text{mol}$, roughly three orders of magnitude smaller and nearly independent of composition.


## From gas to solution

The difference in $\bar{V}_i$ between gases and liquids explains why concentration replaces pressure as the composition variable in solution.
For an ideal gas mixture, $\mu_i = \mu_i^\circ + RT\ln(P_i/P^\circ)$ where $P_i = x_i P$ from Dalton's law (@ideal-gas).
In dilute solution, the mole fraction $x_i \approx c_i / c^\circ$, and replacing partial pressure with concentration gives:

$$
\mu_i = \mu_i^\circ(T) + RT \ln \frac{c_i}{c^\circ}
$$

The pressure correction $\bar{V}_i \Delta P$ for a liquid is negligible compared to $RT\ln(c_i/c^\circ)$ at any realistic pressure:

- **Gas**: $\bar{V}_i = RT/P \approx 0.025\,\text{m}^3/\text{mol}$; pressure strongly affects $\mu$
- **Liquid**: $\bar{V}_i \approx 10^{-5}\,\text{m}^3/\text{mol}$; pressure dependence absorbed into $\mu_i^\circ$

This expression holds in the ideal dilute solution limit ($\gamma_i = 1$); deviations are captured through the activity $a_i = \gamma_i c_i / c^\circ$ introduced in @nonidealities.

```{admonition} Intuition
:class: tip

The term $RT \ln(c_i/c^\circ)$ sets the **energy scale** of concentration effects.

At 25 °C, $RT \approx 2.5\ \text{kJ/mol}$:
- a tenfold concentration difference gives $RT \ln 10 \approx 5.7\ \text{kJ/mol}$
- a hundredfold difference gives $\approx 11.4\ \text{kJ/mol}$

These are not small corrections: 5–10 kJ/mol is comparable to weak intermolecular interactions such as hydrogen bonds, and large enough to drive diffusion, dissolution, and chemical reactions.
$\mu_i$ is therefore not an abstract derivative: it is a **free energy per mole**, and concentration differences translate into quantitatively significant driving forces.
```

**Standard states.**
The reference values $P^\circ$ and $c^\circ$ serve the same role: the argument of a logarithm must be dimensionless.
Any choice works; changing them shifts $\mu^\circ$ by a constant, leaving all physical predictions unchanged since those depend on differences $\Delta\mu$.
By convention $P^\circ = 1\,\text{bar}$ and $c^\circ = 1\,\text{mol/L}$, so values expressed in these units enter the logarithm directly.

```{admonition} Connection to statistical mechanics
:class: note

For readers of Part 1: the logarithm and the $RT$ prefactor both trace directly to statistical mechanics.
The entropy of mixing that produces $\ln(c_i/c^\circ)$ is the Gibbs form $\mathcal{H} = -\sum_i p_i \ln p_i$ from @entropy, applied to the distribution of molecules over positions.
The prefactor $RT = k_B T \cdot N_A$ is the Boltzmann energy scale $k_BT$ from @maxwell-boltzmann, converted from per-molecule to per-mole units by $N_A$.
The thermodynamic derivation above and the statistical one are the same calculation viewed at different scales.
```


```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-chempot

import numpy as np
import matplotlib.pyplot as plt

c_ref = 1.0
c = np.logspace(-3, 1, 300) * c_ref
R = 8.314

fig, ax = plt.subplots(figsize=(6, 4))
for T, color in zip([280, 310, 340], ["C0", "C1", "C2"]):
    mu = R * T * np.log(c / c_ref)
    ax.plot(c / c_ref, mu / 1000, color=color, label=f"$T = {T}$ K")

ax.axvline(1, color="gray", linestyle=":", lw=1)
ax.axhline(0, color="gray", linestyle=":", lw=1)
ax.annotate(r"$c = c^\circ$, $\mu = \mu^\circ$", xy=(1, 0), xytext=(2, -3),
            fontsize=9, color="gray",
            arrowprops=dict(arrowstyle="->", color="gray", lw=0.8))
ax.set_xscale("log")
ax.set_xlabel(r"$c / c^\circ$")
ax.set_ylabel(r"$(\mu - \mu^\circ)\ /\ \mathrm{kJ\,mol^{-1}}$")
ax.legend()
fig.tight_layout()
```

```{figure} #cell-chempot
:name: fig-chempot

Ideal chemical potential $\mu = \mu^\circ + RT\ln(c/c^\circ)$ as a function of concentration for three temperatures.
At $c = c^\circ$ the chemical potential equals $\mu^\circ$ by definition.
The slope $RT$ grows with temperature: the same concentration difference drives a larger free energy difference at higher $T$.
As $c \to 0$, $\mu \to -\infty$: an infinitely dilute species has an infinitely strong tendency to draw in more of itself, which is why even trace concentrations always carry a thermodynamic driving force.
```

(phase-equilibrium)=
## Phase and phase equilibrium

A **phase** is a homogeneous region of matter that is uniform in chemical composition and physical state, separated from other regions by a phase boundary.
Examples include liquid water, water vapor, ice, or a crystal of a pure substance; mixtures can also form distinct phases (e.g., oil and water).

At equilibrium between two phases $\alpha$ and $\beta$ with no chemical reaction, matter can transfer across the boundary.
The equilibrium condition follows from $dG = 0$:

$$\mu_i(\alpha) = \mu_i(\beta) \quad \text{for all } i$$

This is the fundamental **phase equilibrium** criterion: the chemical potential of each species must be equal in all phases at equilibrium.
When $\mu_i(\alpha) > \mu_i(\beta)$, species $i$ has a thermodynamic driving force to move from phase $\alpha$ to phase $\beta$; the system evolves until the potentials equalise.

**Clausius-Clapeyron equation.**
Differentiating the equality $\mu_i(\alpha) = \mu_i(\beta)$ along the coexistence curve gives:

$$\frac{dP}{dT} = \frac{\Delta H_\text{trans}}{T \cdot \Delta V}$$

where $\Delta H_\text{trans}$ is the enthalpy change of the transition and $\Delta V$ is the volume change.
This equation predicts how the pressure of phase coexistence changes with temperature.
For liquid-vapor transitions, $\Delta V \approx V_\text{gas} = RT/P$ (since $V_\text{liquid} \ll V_\text{gas}$), and integrating gives the familiar vapor-pressure dependence on temperature.

```{admonition} Intuition
:class: tip
$\Delta H_\text{trans}$ links a phase transition to its energy cost.
Vaporising 1 L of water at 100 °C requires about 2260 kJ, roughly five times more than heating the same water from 0 °C to 100 °C.
The Clausius-Clapeyron equation converts this enthalpy into a slope along the coexistence curve: a large $\Delta H_\text{trans}$ means the boiling point shifts strongly with pressure.
Water boils near 70 °C at high altitude (low pressure) for exactly this reason; a pressure cooker raises the boiling point above 100 °C for the same reason.
Phase transitions are also exploited in refrigeration and heat pumps: evaporation absorbs a large amount of heat at constant temperature, condensation releases it, and the asymmetry between the two sides drives efficient heat transfer.
```


---

The next chapter introduces corrections for non-ideal behaviour: at high concentration, molecules interact and the simple $RT\ln(c/c^\circ)$ expression no longer holds.
