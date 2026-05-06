---
kernelspec:
  name: python3
  display_name: Python 3 (ipykernel)
---

(nonidealities)=
# Non-idealities

The ideal gas law and the ideal chemical potential are useful starting points, but they rest on a single assumption: molecules do not interact.
For gases at low pressure and solutions at low concentration this is a reasonable approximation.
At high pressure or high concentration it breaks down, and the thermodynamic framework needs corrections.
This chapter introduces those corrections: **fugacity** for gases and **activity** for solutions, then covers the models used in practice.
Readers primarily interested in aqueous electrolyte systems can skip directly to the activity and Debye-Hückel sections.


## Equations of state

An **equation of state** is a relation between the macroscopic state variables $P$, $V$, $T$, and $n$ of a substance.
The ideal gas law $PV = nRT$ is the simplest; it works well at low pressure where molecules are far apart and interactions are negligible.
When those interactions matter, a more realistic equation of state is needed.


## Real gases: the van der Waals equation

The ideal gas law fails because it ignores two effects: molecules have finite size (they cannot overlap), and they attract each other at short range.
The simplest equation of state that accounts for both is the **van der Waals equation** (1873):

$$
\left(P + \frac{an^2}{V^2}\right)\left(V - nb\right) = nRT
$$

The correction $nb$ subtracts the volume actually occupied by the molecules themselves, so $V - nb$ is the free volume available for motion.
The term $an^2/V^2$ corrects for attractive interactions: molecules pulling on each other reduce the force they exert on the walls, lowering the effective pressure.
The constants $a$ and $b$ are specific to each gas and must be measured; when $a = b = 0$ the ideal gas law is recovered.


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


## Real solutions: activity

In solution, the role of pressure is played by concentration, and the same logic applies.
The ideal chemical potential $\mu_i = \mu_i^\circ + RT\ln(c_i/c^\circ)$ assumes no interactions between solute molecules.
The **activity** $a_i$ replaces the dimensionless concentration $c_i/c^\circ$ with a corrected quantity that accounts for interactions:

$$
\mu_i = \mu_i^\circ(T) + RT\ln a_i
$$

The activity is written as:

$$
a_i = \gamma_i \frac{c_i}{c^\circ}
$$

where $\gamma_i$ is the **activity coefficient**.
In the dilute limit, molecules are far apart and interactions are negligible, so $\gamma_i \to 1$ and $a_i \to c_i/c^\circ$.
At higher concentrations, $\gamma_i$ deviates from 1: attractive interactions give $\gamma_i < 1$; repulsive or excluded-volume interactions give $\gamma_i > 1$.


## Phase equilibrium and activity: Raoult's and Henry's laws

Activity coefficients connect the chemical potential formalism to liquid and vapour phase equilibrium for mixtures.
At equilibrium between liquid and vapour phases, $\mu_i(\text{liquid}) = \mu_i(\text{vapour})$ (see @chemical-potential).
For an ideal vapour phase (valid at low pressure) this gives $P_i \propto x_i \gamma_i$.
Two limiting laws follow, corresponding to different choices of the activity coefficient model.

**Raoult's law** applies when components are chemically similar (e.g., water + ethanol).
It assumes ideal mixing in the liquid: $\gamma_i = 1$, giving:

$$P_i = x_i P_i^*$$

The vapour pressure of each component is proportional to its mole fraction, with proportionality constant $P_i^*$ (the vapour pressure of the pure component).
Raoult's law is exact for ideal solutions and provides a symmetric reference state: $\gamma_i \to 1$ as $x_i \to 1$.

**Henry's law** applies at infinite dilution when the solute is chemically unlike the solvent (e.g., O$_2$ or CO$_2$ in water).
It uses an asymmetric reference state: $\gamma_i \to 1$ as $x_i \to 0$.
At low solute concentration, the vapour pressure follows:

$$P_i = K_{H,i}\, x_i$$

where $K_{H,i}$ is the **Henry's law constant**.
Henry's law is the appropriate limit for trace components in a mixture where the solvent dominates.

Both Raoult's and Henry's laws are special cases of $\mu_i(\text{liq}) = \mu_i(\text{vap})$ with different activity coefficient models.
In practice, the choice between them depends on the system and the composition range of interest.


## The Gibbs-Duhem constraint on activity coefficients

The Gibbs-Duhem equation from @laws-of-thermodynamics states that the intensive variables $T$, $P$, and $\{\mu_i\}$ cannot all vary independently:

$$
S\,dT - V\,dP + \sum_i n_i\,d\mu_i = 0
$$

At constant $T$ and $P$ this reduces to $\sum_i n_i\,d\mu_i = 0$.
Substituting $\mu_i = \mu_i^\circ + RT\ln a_i$:

$$
\sum_i n_i\,d\ln a_i = 0
$$

This means the activity coefficients of different species in the same solution are linked: the activity coefficients cannot be chosen independently for each species.
If the activity coefficient of the solute increases with concentration, the activity coefficient of the solvent must respond accordingly.
In practice this constrains how activity coefficient models are constructed and tested.
The underlying thermodynamic potential is the **excess Gibbs energy** $G^E = \Delta G^\text{mix} - \Delta G^\text{mix,ideal}$: activity coefficients arise as its partial molar derivatives, $RT\ln\gamma_i = (\partial G^E/\partial n_i)_{T,P,n_{j\neq i}}$, in exact analogy with $\mu_i = (\partial G/\partial n_i)_{T,P,n_{j\neq i}}$.
Models such as NRTL or UNIQUAC (for non-electrolyte mixtures) and Pitzer (for concentrated electrolytes) are built by specifying $G^E$ and deriving $\gamma_i$ from it; the Debye-Hückel and Davies models used in this book bypass $G^E$ and specify $\gamma_i$ directly from the physics of the ionic atmosphere.


## Debye-Hückel theory for electrolyte solutions

For ionic solutions, the dominant source of non-ideality is long-range electrostatic interactions between ions.
Debye and Hückel (1923) computed these analytically by considering the **ionic atmosphere**: each ion is surrounded by a diffuse cloud of oppositely charged ions that partially screens its charge.

The starting point is the linearised Poisson-Boltzmann equation, which relates the electric potential $\psi$ around a central ion to the charge distribution of the surrounding ions.
Its solution shows that $\psi$ decays as $e^{-\kappa r}/r$, where $\kappa$ is the inverse **Debye length**:

$$
\kappa^2 = \frac{2N_A e^2 I}{\varepsilon_0 \varepsilon_r RT}
$$

Here $e$ is the elementary charge, $\varepsilon_r$ is the relative permittivity of the solvent, and $I$ is the **ionic strength**:

$$
I = \frac{1}{2}\sum_i z_i^2 \frac{c_i}{c^\circ}
$$

where $z_i$ is the charge number of ion $i$ and $c^\circ = 1\ \text{mol/L}$ is the standard concentration, making $I$ dimensionless but numerically equal to the ionic strength in mol/L.

The excess free energy from the ionic atmosphere gives the **Debye-Hückel limiting law** for a single ionic species $i$ with charge $z_i$:

$$
\log_{10} \gamma_i = -A z_i^2 \sqrt{I},
$$

where $A \approx 0.509\ \text{mol}^{-1/2}\text{L}^{1/2}$ at 25 °C.
The law is exact as $I \to 0$ and accurate to $I \approx 0.1$.

At higher ionic strengths the point-charge approximation breaks down.
The **Davies equation** adds an empirical linear correction that extends the useful range to $I \approx 0.5$:

$$
\log_{10} \gamma_i = -A z_i^2 \left(\frac{\sqrt{I}}{1 + \sqrt{I}} - 0.3\,I\right).
$$

The electrochemistry convention uses $c^\circ = 1\ \text{mol/L}$, so $I$ above is numerically equal to the molar ionic strength.
CADET and this library use SI units throughout: concentrations in mol/m³, so $1\ \text{mol/L} = 10^3\ \text{mol/m}^3$.
The activity coefficient models convert internally; the figure below uses mol/m³ on the horizontal axis, where $1\ \text{mol/m}^3 = 1\ \text{mM}$.

Both expressions scale as $z_i^2$: a divalent ion at the same ionic strength has a correction four times larger than a monovalent one (@fig-activity).

```{admonition} Intuition
:class: tip

In solution, each ion is surrounded by a diffuse cloud of oppositely charged ions.
This cloud partially cancels the ion's electric field at larger distances, so other ions interact with a weakened effective charge.
As salt concentration increases, these clouds become denser and screening becomes stronger, leading to $\gamma_i < 1$.
At higher ionic strengths the simple dilute-solution picture begins to break down, which is why more detailed activity models such as Davies introduce empirical corrections.
```

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-activity

import numpy as np
import matplotlib.pyplot as plt
from reactions.api import (
    ActivityCoefficientDebyeHuckel,
    ActivityCoefficientDavies,
    PhysicalState,
)

dh  = ActivityCoefficientDebyeHuckel()
dav = ActivityCoefficientDavies()

I_dh  = np.linspace(1, 100, 200)    # mol/m³, DH valid range (~0–100 mM)
I_dav = np.linspace(1, 500, 400)    # mol/m³, Davies valid range (~0–500 mM)

def gamma_series(model, I_values, z):
    return np.array([
        model.activity(
            PhysicalState(c=np.array([1.0]), T=298.15, I=float(I)),
            np.array([float(z)]),
        )[0]
        for I in I_values
    ])

fig, axes = plt.subplots(1, 2, figsize=(9, 3.8), sharey=True)

for ax, z, title in [
    (axes[0], 1, r"$|z| = 1$  (monovalent)"),
    (axes[1], 2, r"$|z| = 2$  (divalent)"),
]:
    ax.plot(I_dh,  gamma_series(dh,  I_dh,  z), color="C0", label="Debye-Hückel")
    ax.plot(I_dav, gamma_series(dav, I_dav, z), color="C1", label="Davies")
    ax.axvline(100, color="gray", lw=0.7, ls=":", label="DH valid limit")
    ax.set_xlabel(r"$I$  [mol/m³]")
    ax.set_xlim(0, 550)
    ax.set_title(title)
    ax.legend(fontsize=8)

axes[0].set_ylabel(r"$\gamma_i$")
fig.tight_layout()
```

```{figure} #cell-activity
:name: fig-activity

Activity coefficient $\gamma_i$ vs ionic strength at 25 °C for monovalent (left) and divalent (right) species.
Concentrations are in mol/m³ (SI); $1\ \text{mol/m}^3 = 1\ \text{mM}$.
Debye-Hückel is accurate below $\approx 100\ \text{mol/m}^3$ (dotted line); Davies extends the range to $\approx 500\ \text{mol/m}^3$.
Divalent species experience markedly stronger suppression because both models scale as $z_i^2$.
```

---

## From thermodynamics to reactions

The framework developed in Parts 1 and 2 is now complete.
Statistical mechanics provided the microscopic foundation: microstates, entropy, the Boltzmann factor, and the ideal gas law.
Thermodynamics built on that foundation: the laws, the fundamental relation, the potentials, and the chemical potential.
Non-ideality corrections via activity coefficients complete the toolkit for real solutions.

Part 3 puts this machinery to work.
The central quantity is the **reaction Gibbs energy** $\Delta_r G = \sum_i \nu_i \mu_i$: the slope of $G$ as the reaction progresses.
Setting $\Delta_r G = 0$ defines the equilibrium constant $K$, and $\Delta_r G^\circ = -RT \ln K$ connects it to tabulated thermodynamic data.
Activities appear throughout: the equilibrium constant and all subsequent results are written in terms of $a_i = \gamma_i c_i / c^\circ$, inheriting the corrections introduced above.
