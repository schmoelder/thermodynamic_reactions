---
kernelspec:
  name: python3
  display_name: Python 3 (ipykernel)
---

(mixing)=
# Mixing and Excess Gibbs Energy

@nonidealities established that pure-fluid non-ideality is captured by fugacity: the real chemical potential replaces $P$ with an effective pressure $f$ derived from the equation of state.
That correction applies to each component individually; it says nothing about what happens when unlike species are combined.
Mixing introduces a distinct second layer of non-ideality: the interactions between unlike species can differ from those in the pure components, so the Gibbs energy of mixing deviates from the ideal baseline even when each pure component is itself ideal.
Both corrections modify the chemical potential: fugacity addresses deviations in the $P$–$V$–$T$ sector; activity addresses deviations in the composition sector.
The two are formally separable projections of the same state function $G(T, P, \{n_i\})$.

The chemical potential of a species in an ideal solution is $\mu_i = \mu_i^* + RT\ln x_i$, where $\mu_i^*(T,P)$ is the pure-component value from @chemical-potential.
That expression rests on the assumption of no unlike-species interactions.
Real mixtures violate it.
Mixing equal volumes of water and ethanol yields a solution roughly $4\,\%$ smaller than the sum of their pure volumes; the process is exothermic.
Mixing acetone and chloroform also releases heat.
Mixing hexane and heptane is nearly ideal.
The source of the difference is not geometry or accident: it is whether the unlike intermolecular interactions are stronger or weaker than the like ones.

The strategy here is not to replace the ideal framework but to correct it with one object: the **excess Gibbs energy** $G^E$.
Every departure from ideal mixing, in volume, enthalpy, entropy, and activity, is a derivative of $G^E$ with respect to the appropriate variable.


## Solutions, solvents, and solutes

A **solution** is a liquid mixture in which one component — the **solvent** — is present in large excess and sets the physical environment: the dielectric constant, the density, the standard state.
All other components are **solutes**, present at lower concentration.
In aqueous chemistry the solvent is water; in chromatography it is the mobile phase.
The solvent/solute distinction is not merely a convention: it determines which limiting law applies.
The solvent obeys Raoult's law ($\gamma \to 1$ as $x \to 1$, symmetric reference), while dilute solutes follow Henry's law ($\gamma \to \gamma^\infty$ as $x \to 0$, asymmetric reference); both emerge naturally from the $G^E$ framework developed below.


## Ideal mixing

Consider forming a mixture from $n_1$ moles of component 1 and $n_2$ moles of component 2 at constant $T$ and $P$, with mole fractions $x_i = n_i / n$ and total amount $n = \sum_i n_i$.
For an ideal solution, $\mu_i = \mu_i^*(T,P) + RT\ln x_i$, so the Gibbs energy of the mixture is:

$$
G = \sum_i n_i \mu_i^*(T,P) + RT\sum_i n_i\ln x_i
$$

The Gibbs energy of the unmixed pure components is $G^\text{pure} = \sum_i n_i \mu_i^*$.
Their difference, the **Gibbs energy of mixing**, is:

$$
\Delta_\text{mix} G^\text{id} = nRT\sum_i x_i\ln x_i
$$

Since $x_i \in (0,1)$, every $\ln x_i < 0$, so $\Delta_\text{mix} G^\text{id} < 0$ for any mixture.
Ideal mixing is always spontaneous; it is driven entirely by entropy.
This follows directly: $\Delta_\text{mix} S^\text{id} = -(\partial \Delta_\text{mix} G^\text{id}/\partial T)_P = -nR\sum_i x_i\ln x_i > 0$.

The ideal mixing enthalpy and volume are both zero:

$$
\Delta_\text{mix} H^\text{id} = 0, \qquad \Delta_\text{mix} V^\text{id} = 0
$$

These identities are the defining feature of ideal mixing: replacing like interactions with unlike ones introduces no additional energetic or volumetric preference, so mixing produces neither heat exchange nor volume change.

The entropy of mixing is nevertheless nonzero: mixing increases the number of accessible molecular arrangements and therefore increases entropy.
What vanishes in the ideal limit is the *excess* contribution beyond this combinatorial baseline:

$$
G^E = H^E = V^E = S^E = 0
$$


## Excess Gibbs energy

The **excess Gibbs energy** $G^E$ is the difference between the actual Gibbs energy of mixing and the ideal prediction:

$$
G^E = \Delta_\text{mix} G - \Delta_\text{mix} G^\text{id}
$$

When $G^E = 0$ the mixture is ideal.
When $G^E \neq 0$, interactions between unlike species cause measurable deviations.
Every excess property follows by differentiating $G^E$ with respect to the appropriate variable.
From the fundamental relation $dG = -S\,dT + V\,dP + \sum_i \mu_i\,dn_i$ (@thermodynamic-potentials):

$$
V^E = \left(\frac{\partial G^E}{\partial P}\right)_{T,n}, \qquad
S^E = -\left(\frac{\partial G^E}{\partial T}\right)_{P,n}, \qquad
H^E = G^E - T\left(\frac{\partial G^E}{\partial T}\right)_{P,n}
$$

The excess heat capacity is:

$$
C_p^E = -T\left(\frac{\partial^2 G^E}{\partial T^2}\right)_{P,n}
$$

One function, $G^E$, determines all excess properties simultaneously and self-consistently.
For water/ethanol, $V^E < 0$ and $H^E < 0$ from the same $G^E$: the mixture is denser than ideal and releases heat on mixing, both consequences of the water--ethanol hydrogen bond being stronger than the average of water--water and ethanol--ethanol interactions.

In practice, $G^E$ models are written in intensive form.
Define the molar excess Gibbs energy $g^E = G^E/n$; most models specify $g^E/RT$ as a dimensionless function of composition and temperature.


## Activity coefficient as partial molar excess Gibbs energy

The real chemical potential of species $i$ can be written as:

$$
\mu_i = \mu_i^*(T,P) + RT\ln x_i + RT\ln\gamma_i
$$

where $\gamma_i$ is the **activity coefficient**.
The first two terms are the ideal expression from @chemical-potential; the third is the excess contribution.
Since $G^E = \sum_i n_i(\mu_i - \mu_i^\text{id})$ and $RT\ln\gamma_i = \mu_i - \mu_i^\text{id}$, the activity coefficient is the partial molar excess Gibbs energy:

$$
RT\ln\gamma_i = \left(\frac{\partial G^E}{\partial n_i}\right)_{T,P,n_{j\neq i}}
$$

This is the central result: $\gamma_i$ is not an empirical correction appended to the ideal expression; it is the derivative of $G^E$ with respect to the moles of species $i$.
Any model for $G^E(x, T, P)$ yields all $\gamma_i$ by differentiation.

The **activity** in mole-fraction convention is then:

$$
a_i = \gamma_i\, x_i
$$

and the chemical potential takes the compact form $\mu_i = \mu_i^* + RT\ln a_i$, exact for any mixture.
The ideal limit $\gamma_i \to 1$ as $x_i \to 1$ is a consequence of the definition: when species $i$ is the dominant component it interacts mainly with itself, recovering the pure-component reference.

```{admonition} Intuition
:class: tip

$G^E$ is the extra free energy of the mixture relative to the ideal prediction: positive when unlike interactions are unfavourable, negative when they are favourable.
Differentiating with respect to the amount of species $i$ gives the per-mole version of that correction, the activity coefficient.
A negative $G^E$ (exothermic, contracting) gives $\gamma_i < 1$: the species prefers the mixture to the pure state, so its effective concentration is lower than its mole fraction.
A positive $G^E$ (endothermic) gives $\gamma_i > 1$: the species would rather be with its own kind, behaving as if present at a higher effective concentration.
```


## Gibbs-Duhem consistency

Activity coefficients of different species in the same solution are not independent.
The Gibbs-Duhem equation at constant $T$ and $P$ states $\sum_i n_i\,d\mu_i = 0$ (@thermodynamic-potentials).
Substituting $\mu_i = \mu_i^* + RT\ln(\gamma_i x_i)$ and using $\sum_i n_i\,dx_i = 0$:

$$
\sum_i x_i\,d\ln\gamma_i = 0
$$

This is the **Gibbs-Duhem consistency condition** for activity coefficients: if $\gamma_1$ increases with composition, $\gamma_2$ must adjust in a compensating way.
Any model for $G^E$ satisfies the condition automatically, because the $\gamma_i$ are partial derivatives of the same function.
Models that specify $\gamma_i$ directly, rather than through $G^E$, must be tested against this constraint.


## Raoult and Henry as limiting cases

Two classical limiting laws emerge from the $G^E$ framework.

**Raoult's law** corresponds to $G^E \to 0$ as $x_i \to 1$.
All $\gamma_i \to 1$, so the vapour pressure of each component is $P_i = x_i P_i^*$, with $P_i^*$ the vapour pressure of the pure liquid.
This **symmetric reference state** ($\gamma_i \to 1$ as $x_i \to 1$) holds exactly for ideal mixtures and approximately when the species are chemically similar.

**Henry's law** is the dilute limit $x_i \to 0$.
In this regime $\gamma_i \to \gamma_i^\infty$, the **activity coefficient at infinite dilution**.
The vapour pressure then follows $P_i = K_{H,i}\,x_i$ where $K_{H,i} = \gamma_i^\infty P_i^*$.
This **asymmetric reference state** defines the standard state so that $\gamma_i \to 1$ as $x_i \to 0$; it is the natural convention for dilute solutes and the basis for the concentration-based activity used in Part 3.

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-mixing

import numpy as np
import matplotlib.pyplot as plt
from reactions.plots import setup_figure, COLORS

x1 = np.linspace(0.001, 0.999, 400)
x2 = 1 - x1

fig, axes = setup_figure(1, 2)

for A, color, label in [
    (0.0, "C0", r"ideal  ($A = 0$)"),
    (1.8, "C1", r"positive  ($A = 1.8$)"),
    (-1.8, "C2", r"negative  ($A = -1.8$)"),
]:
    dG = x1 * np.log(x1) + x2 * np.log(x2) + A * x1 * x2
    g1 = np.exp(A * x2**2)
    axes[0].plot(x1, dG, color=color, label=label)
    axes[1].plot(x1, g1, color=color, label=label)

axes[0].axhline(0, color="gray", lw=0.5, ls="--")
axes[0].set_xlabel(r"$x_1$")
axes[0].set_ylabel(r"$\Delta_\mathrm{mix}G \;/\; nRT$")
axes[0].set_xlim(0, 1)
axes[0].set_xticks([0, 0.5, 1])
axes[0].legend(fontsize=8)

axes[1].axhline(1, color="gray", lw=0.5, ls="--")
axes[1].set_xlabel(r"$x_1$")
axes[1].set_ylabel(r"$\gamma_1$")
axes[1].set_xlim(0, 1)
axes[1].set_xticks([0, 0.5, 1])
axes[1].legend(fontsize=8)

fig.tight_layout()
```

```{figure} #cell-mixing
:name: fig-mixing

Left: Gibbs energy of mixing $\Delta_\text{mix}G/(nRT)$ vs mole fraction $x_1$ for the Margules one-suffix model $g^E/RT = Ax_1x_2$.
Right: corresponding activity coefficient $\gamma_1 = \exp(Ax_2^2)$.
The ideal case ($A = 0$) corresponds to enthalpy-neutral mixing, so the Gibbs energy of mixing is purely entropic and always negative.
Positive $A$ represents unfavourable unlike interactions, producing $\gamma_1 > 1$ and positive deviations from ideality.
Negative $A$ represents favourable unlike interactions, giving $\gamma_1 < 1$ and stabilised mixing.
In all cases the Raoult limit is recovered as $x_1 \to 1$, where $\gamma_1 \to 1$ (dashed line).
```

## Activity in solution: concentration convention

The mole-fraction activity $a_i = \gamma_i x_i$ from the previous sections uses the pure-component standard state.
For dilute solutions the natural reference shifts: the solute is rarely close to pure, and $c^\circ = 1\ \mathrm{mol/L}$ is a more convenient standard.
In the concentration convention the chemical potential is:

$$
\mu_i = \mu_i^\circ(T) + RT\ln a_i, \qquad a_i = \gamma_i \frac{c_i}{c^\circ}
$$

where $\gamma_i \to 1$ as $c_i \to 0$ (Henry-law reference state).
This is the convention used throughout Part 3 and Part 4.
In the dilute limit $a_i \to c_i/c^\circ$; at higher concentrations $\gamma_i$ deviates from unity according to the $G^E$ model appropriate for the system.

Models such as NRTL and UNIQUAC (non-electrolyte mixtures) and Pitzer (concentrated electrolytes) are built by specifying a $G^E$ model and deriving $\gamma_i$ as partial molar derivatives; they are mentioned here but not developed further in this book.


## Debye-Hückel theory for electrolyte solutions

For ionic solutions, the dominant source of $G^E$ is long-range electrostatic interactions between ions.
{cite:t}`debye1923` derived the electrostatic contribution to $G^E$ by modelling the **ionic atmosphere**: each ion is surrounded by a diffuse cloud of oppositely charged ions that partially screens its charge.
Solving the linearised Poisson-Boltzmann equation gives an electrostatic potential that decays as $e^{-\kappa r}/r$, where $\kappa^{-1}$ is the **Debye length**:

$$
\kappa^2 = \frac{2N_A e^2 I}{\varepsilon_0 \varepsilon_r RT}
$$

Here $e$ is the elementary charge, $\varepsilon_r$ is the relative permittivity of the solvent, and $I$ is the **ionic strength**:

$$
I = \frac{1}{2}\sum_i z_i^2 \frac{c_i}{c^\circ}
$$

where $z_i$ is the charge number of ion $i$.
The electrostatic excess free energy of the ionic atmosphere gives the **Debye-Hückel limiting law**:

$$
\log_{10} \gamma_i = -A z_i^2 \sqrt{I}
$$

where $A \approx 0.509\ \mathrm{mol^{-1/2}\,L^{1/2}}$ at 25 °C.
The law is exact as $I \to 0$ and accurate to $I \approx 0.1\ \mathrm{mol/L}$.

At higher ionic strengths the point-charge approximation breaks down.
The **Davies equation** {cite:p}`davies1938` adds an empirical linear correction that extends the useful range to $I \approx 0.5\ \mathrm{mol/L}$:

$$
\log_{10} \gamma_i = -A z_i^2 \left(\frac{\sqrt{I}}{1 + \sqrt{I}} - 0.3\,I\right)
$$

Both expressions scale as $z_i^2$: a divalent ion has a correction four times larger than a monovalent one at the same ionic strength.

```{admonition} Implementation note
:class: note

The activity coefficient models in this library use mol/m³ as the internal unit for ionic strength; $1\ \mathrm{mol/m}^3 = 1\ \mathrm{mM}$.
The horizontal axis of @fig-activity uses this unit convention.
```

```{admonition} Intuition
:class: tip

In solution, each ion is surrounded by a diffuse cloud of oppositely charged ions.
This cloud partially cancels the ion's electric field at larger distances, so other ions interact with a weakened effective charge.
The Debye-Hückel result is the $G^E$ of the ionic atmosphere: differentiating with respect to $n_i$ gives $\gamma_i < 1$ — ions are stabilised by their surroundings relative to the ideal dilute limit.
As salt concentration increases the clouds become denser, screening becomes stronger, and $\gamma_i$ decreases further.
Above $I \approx 0.1\ \mathrm{mol/L}$ the simple point-charge picture breaks down, which is why Davies adds an empirical correction.
```

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-activity

import numpy as np
import matplotlib.pyplot as plt
from reactions.plots import setup_figure, COLORS
from reactions.api import (
    ActivityCoefficientDebyeHuckel,
    ActivityCoefficientDavies,
    PhysicalState,
)

dh = ActivityCoefficientDebyeHuckel()
dav = ActivityCoefficientDavies()

I_dh = np.linspace(1, 100, 200)
I_dav = np.linspace(1, 500, 400)


def gamma_series(model, I_values, z):
    return np.array(
        [
            model.activity(
                PhysicalState(c=np.array([1.0]), T=298.15, I=float(I)),
                np.array([float(z)]),
            )[0]
            for I in I_values
        ]
    )


fig, axes = setup_figure(1, 2, sharey=True)

for ax, z, title in [
    (axes[0], 1, r"$|z| = 1$  (monovalent)"),
    (axes[1], 2, r"$|z| = 2$  (divalent)"),
]:
    ax.plot(I_dh, gamma_series(dh, I_dh, z), color="C0", label="Debye-Hückel")
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
Concentrations are in mol/m³; $1\ \mathrm{mol/m}^3 = 1\ \mathrm{mM}$.
Debye-Hückel is accurate below $\approx 100\ \mathrm{mol/m}^3$ (dotted line); Davies extends the range to $\approx 500\ \mathrm{mol/m}^3$.
```

---

The framework developed in Parts 1 and 2 is now complete.
Statistical mechanics provided the microscopic foundation: microstates, entropy, the Boltzmann factor, and the ideal gas law.
Thermodynamics built on that foundation: the laws, the fundamental relation, the potentials, the chemical potential, and two layers of non-ideality — fugacity correcting $\mu$ in the $P$–$V$–$T$ sector, activity correcting $\mu$ in the composition sector — both projections of the same Gibbs energy $G(T, P, \{n_i\})$.

Part 3 puts this machinery to work.
The central quantity is the **reaction Gibbs energy** $\Delta_r G = \sum_i \nu_i \mu_i$: the slope of $G$ as the reaction progresses.
Setting $\Delta_r G = 0$ defines the equilibrium constant $K$, and $\Delta_r G^\circ = -RT \ln K$ connects it to tabulated thermodynamic data.
Activities appear throughout: the equilibrium constant and all subsequent results are written in terms of $a_i = \gamma_i c_i / c^\circ$, inheriting the $G^E$-derived corrections introduced here.
