---
kernelspec:
  name: python3
  display_name: Python 3 (ipykernel)
---

(chemical-potential)=
# Gibbs Energy and Chemical Potential

@thermodynamic-potentials established that the Gibbs free energy $G$ is minimised at equilibrium under constant temperature and pressure.
There $G$ was a structural object, a Legendre transform adapted to constant $T$ and $P$, but its operational meaning was left open.
This chapter develops what $G$ actually tells us: first as a criterion for spontaneity, then through its structure as a sum of per-species contributions.


## Spontaneity

A **spontaneous process** is one that proceeds on its own without any external driving force.
Heat flows from hot to cold; gases expand into vacuum; salts dissolve in water.
The direction of spontaneous change is not determined by energy alone: dissolving ammonium nitrate in water absorbs heat yet happens spontaneously, while freezing water releases heat yet requires cooling below $0\,\degree\text{C}$.
What decides is the competition between energy and entropy, captured by $G$.

Spontaneity is not limited to chemical reactions.
Any process occurring at constant $T$ and $P$ (phase transitions, mixing, protein folding, membrane transport, electrochemistry) is spontaneous if and only if it lowers $G$.
This universality is what makes $G$ the central quantity of chemical thermodynamics.

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

Up to this point, $G$ has been used as a criterion for whether processes proceed spontaneously at fixed $T$ and $P$.
We now shift perspective and ask how $G$ is structured as a function of composition, i.e. how it changes when matter is redistributed between species.

From the definition
$$
\mu_i = \left.\frac{\partial G}{\partial n_i}\right|_{T,P,n_{j \neq i}}
$$
the chemical potential is the free energy cost of adding one mole of species $i$ while holding temperature, pressure, and all other amounts fixed.
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

How strongly $\mu_i$ responds to pressure depends on the volume contributed by each mole of species $i$.
In a mixture, this is the **partial molar volume**:

$$
\bar{V}_i = \left.\frac{\partial V}{\partial n_i}\right|_{T,P,n_{j\neq i}}
$$

It is the volume change when an infinitesimal amount of $i$ is added at constant $T$, $P$, and composition.
For an ideal gas, $\bar{V}_i = V_m = RT/P$ for all components.
For liquids, $\bar{V}_i \approx 10^{-5}\,\text{m}^3/\text{mol}$, roughly three orders of magnitude smaller and nearly independent of composition.
The same construction applies to any extensive property: the **partial molar enthalpy** $\bar{H}_i = (\partial H/\partial n_i)_{T,P,n_{j\neq i}}$ is the enthalpy increment when an infinitesimal amount of $i$ is added at fixed $T$, $P$, and remaining composition.
It appears in the energy balance for reacting systems (@implementation-energy-balance).


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

At 25 °C, $RT \approx 2.5\ \mathrm{kJ/mol}$:
- a tenfold concentration difference gives $RT \ln 10 \approx 5.7\ \mathrm{kJ/mol}$
- a hundredfold difference gives $\approx 11.4\ \mathrm{kJ/mol}$

These are not small corrections: 5–10 kJ/mol is comparable to weak intermolecular interactions such as hydrogen bonds, and large enough to drive diffusion, dissolution, and chemical reactions.
$\mu_i$ is therefore not an abstract derivative: it is a **free energy per mole**, and concentration differences translate into quantitatively significant driving forces.
```

**Standard states.**
The reference values $P^\circ$ and $c^\circ$ serve the same role: the argument of a logarithm must be dimensionless.
Any choice works; changing them shifts $\mu^\circ$ by a constant, leaving all physical predictions unchanged since those depend on differences $\Delta\mu$.
By convention $P^\circ = 1\,\mathrm{bar}$ and $c^\circ = 1\,\mathrm{mol/L}$, so values expressed in these units enter the logarithm directly.

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

---

The next chapter introduces the first source of non-ideal behaviour: departures from the ideal equation of state and the fugacity correction that makes the chemical potential exact for real gases.
