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

Because $U$, $S$, $V$, and $n_i$ are all extensive (they scale linearly with system size), integrating the fundamental relation gives:

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
Its practical application to activity coefficients becomes concrete in @nonidealities.

## The four potentials

The additional potentials are obtained from $U$ by **Legendre transforms**: replacing an extensive natural variable with its conjugate intensive one.
Subtracting $TS$ swaps $S$ for $T$; adding $PV$ swaps $V$ for $P$.
The transform relationships are shown in @fig-potentials.

```{figure} figures/thermodynamic_potentials.png
:name: fig-potentials

The four thermodynamic potentials as Legendre transforms of $U$.
Each arrow swaps one extensive natural variable for its conjugate intensive one.
$G$ is the potential for chemistry: it is minimised at the constant $T$ and $P$ conditions of a laboratory.
```

The definitions and natural variables of each potential are summarised in the table below.

| Fixed conditions | Minimised potential   | Definition    | Natural variables |
| ---------------- | --------------------- | ------------- | ----------------- |
| $S$, $V$         | Internal energy $U$   | --            | $S$, $V$, $n_i$   |
| $S$, $P$         | Enthalpy $H$          | $U + PV$      | $S$, $P$, $n_i$   |
| $T$, $V$         | Helmholtz energy $A$  | $U - TS$      | $T$, $V$, $n_i$   |
| $T$, $P$         | Gibbs free energy $G$ | $U - TS + PV$ | $T$, $P$, $n_i$   |

The Legendre transform machinery is general: it produces a valid potential for any combination of constraints.
Whether a given potential is *useful* depends on what the experiment actually controls.
The Helmholtz energy $A$ is minimised when temperature and volume are both fixed: a rigid, thermostatted container.
This is the natural setting for molecular simulation (where particle number, volume, and temperature are controlled directly) and certain condensed-matter problems, but it is rarely the right choice in solution chemistry.
In a laboratory, reactions proceed in vessels open to the atmosphere at controlled temperature, so volume is not fixed and pressure is.
That makes $G$ the relevant potential, and $A$ does not appear further in this book.

## Thermodynamic equilibrium

The table above identifies which potential is minimised under each set of constraints.
A system has reached **thermodynamic equilibrium** when the relevant potential is at its minimum and its differential is zero: no further spontaneous change is possible.
Three conditions must hold simultaneously:

- **Thermal equilibrium**: temperature is uniform throughout; no net heat flows.
- **Mechanical equilibrium**: pressure is uniform throughout; no net volume changes occur.
- **Chemical equilibrium**: the chemical potential of each species is equal in all regions; no net matter transfer occurs.

The statistical basis for all three was established in @entropy: at equilibrium, total entropy is maximised.
For thermal equilibrium, @entropy showed directly that maximising entropy over two systems sharing energy requires equal $\partial S/\partial U$, i.e. equal temperature.
The same maximum-entropy principle extends to pressure (if boundaries can move, volume redistributes until $P$ equalises) and chemical potential (if species can cross a boundary, they redistribute until $\mu_i$ equalises).
Spontaneous change drives the system toward greater total entropy; equilibrium is the state where no such drive remains.

```{figure} figures/thermodynamic_equilibrium.png
:name: fig-equilibrium

Two subsystems in contact.
Depending on the type of contact, different quantities equalise at equilibrium: a diathermal wall allows heat exchange and equalises $T$; a movable wall equalises $P$; a permeable membrane equalises $\mu_i$ for the species that can cross.
A rigid, sealed, insulating boundary fixes all three and leaves the system in whatever state it starts in.
```

## Le Chatelier's principle

The equilibrium conditions also determine how a system responds to external perturbations.

**Le Chatelier's principle**: a perturbation in an intensive variable drives a shift in the system that counteracts this change.

This behaviour follows directly from the conjugate pairs in the fundamental relation:

- **Temperature** ($T$, conjugate to $S$): increasing temperature shifts equilibrium in the endothermic direction ($\Delta_r H > 0$), while decreasing temperature favours the exothermic direction.
- **Pressure** ($P$, conjugate to $V$): increasing pressure shifts equilibrium toward states of lower volume.
For gas-phase reactions, this typically means fewer moles of gas.
- **Chemical potential** ($\mu_i$, conjugate to $n_i$): increasing the chemical potential of species $i$ (e.g. by adding it) drives the system to consume it.
In reactions, this appears quantitatively as a nonzero reaction Gibbs energy, $\Delta_r G \neq 0$, which drives the system back toward equilibrium (@reaction-gibbs-energy).

The last case is the most general: all equilibrium shifts in reacting systems can be expressed in terms of changes in chemical potentials.


## Differentials and the role of $G$

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

Note that $\mu_i$ appears with the same coefficient $dn_i$ in all four differentials: chemical potential does not depend on which potential describes the system.

The differentials make abstract potentials measurable.
At constant pressure and fixed composition, $dH = \delta Q$: enthalpy equals the heat exchanged.
The **heat capacity at constant pressure**,

$$
C_p = \left.\frac{\partial H}{\partial T}\right|_{P,n}
$$

measures this per kelvin.
Tabulated $C_p$ values make it possible to extrapolate $\Delta H$ and $\Delta S$ from one temperature to another via the Kirchhoff relations, which are taken up in the chapter on equilibrium.

```{admonition} Intuition
:class: tip
$C_p$ links a temperature change to an energy cost.
Heating 1 L of water from room temperature to near boiling requires about 335 kJ, which corresponds to roughly four minutes with a 1500 W kettle.
The same integral, $\int C_p\,dT$, determines how reaction enthalpies change with temperature.
This is exactly what the Kirchhoff relations in @equilibrium-temperature compute.
```


(phase-equilibrium)=
## Phase and phase equilibrium

A **phase** is a homogeneous region of matter that is uniform in chemical composition and physical state, separated from other regions by a phase boundary.
Examples include liquid water, water vapour, ice, or a crystal of a pure substance; mixtures can also form distinct phases (e.g., oil and water).

At equilibrium between two phases $\alpha$ and $\beta$ with no chemical reaction, matter can transfer across the boundary.
The equilibrium condition follows from $dG = 0$:

$$\mu_i(\alpha) = \mu_i(\beta) \quad \text{for all } i$$

This is the fundamental **phase equilibrium** criterion: the chemical potential of each species must be equal in all phases at equilibrium.
When $\mu_i(\alpha) > \mu_i(\beta)$, species $i$ has a thermodynamic driving force to move from phase $\alpha$ to phase $\beta$; the system evolves until the potentials equalise.

**Clausius-Clapeyron equation.**
Differentiating the equality $\mu_i(\alpha) = \mu_i(\beta)$ along the coexistence curve, using $dG = -S\,dT + V\,dP$ from the preceding section, gives:

$$\frac{dP}{dT} = \frac{\Delta H_\text{trans}}{T \cdot \Delta V}$$

where $\Delta H_\text{trans}$ is the enthalpy absorbed at constant pressure when crossing the phase boundary (the latent heat) and $\Delta V$ is the volume change.
This equation predicts how the pressure of phase coexistence changes with temperature.
For liquid-vapour transitions, $\Delta V \approx V_\text{gas} = RT/P$ (since $V_\text{liquid} \ll V_\text{gas}$), and integrating gives the familiar vapour-pressure dependence on temperature.

```{admonition} Intuition
:class: tip
$\Delta H_\text{trans}$ links a phase transition to its energy cost.
Vaporising 1 L of water at 100 °C requires about 2260 kJ, roughly five times more than heating the same water from 0 °C to 100 °C (about 335 kJ from the $C_p$ integral).
The Clausius-Clapeyron equation converts this enthalpy into a slope along the coexistence curve: a large $\Delta H_\text{trans}$ means the boiling point shifts strongly with pressure.
Water boils near 70 °C at high altitude (low pressure) for exactly this reason; a pressure cooker raises the boiling point above 100 °C for the same reason.
```


The integrated form also applies to each potential.
For the Gibbs free energy the $TS$ and $PV$ terms cancel entirely:

$$
G = \sum_i \mu_i n_i
$$

This follows from substituting $G = U - TS + PV$ into the integrated fundamental relation $U = TS - PV + \sum_i\mu_i n_i$, so it is an identity that holds in general, not a condition specific to equilibrium or to constant $T$ and $P$.
The Gibbs free energy is the sum of chemical potentials weighted by amount, a result used throughout the chapters on equilibrium and reactions.

```{admonition} Why does $dG$ contain $-S\,dT + V\,dP$ if $G = \sum_i \mu_i n_i$?
:class: note

$G = \sum_i \mu_i n_i$ might suggest that $G$ depends only on composition, but this is misleading: $\mu_i$ itself depends on $T$, $P$, and composition.
So $G$ is implicitly a function of $T$ and $P$ through the $\mu_i$, and a change in temperature or pressure shifts all chemical potentials and therefore changes $G$ even if no matter is added or removed.
Differentiating $G = \sum_i \mu_i n_i$ gives:

$$
dG = \sum_i \mu_i\,dn_i + \sum_i n_i\,d\mu_i
$$

The Gibbs-Duhem equation (derived above) identifies the second sum: $\sum_i n_i\,d\mu_i = -S\,dT + V\,dP$.
Substituting recovers $dG = -S\,dT + V\,dP + \sum_i \mu_i\,dn_i$, consistent with the definition route.
```

---

The next chapter develops what $G$ tells us in practice: the criterion for spontaneity and the per-species structure through chemical potential.
