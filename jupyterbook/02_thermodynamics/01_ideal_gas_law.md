---
kernelspec:
  name: python3
  display_name: Python 3 (ipykernel)
---

(ideal-gas)=
# The Ideal Gas Law

If you worked through Part 1, you already know that $PV = nRT$ follows from molecular motion and equipartition (@pressure).
This chapter approaches the same law from the opposite direction: the historical experiments that established it empirically, long before anyone understood why it was true.
If you skipped Part 1, $PV = nRT$ is the one result to take as a starting point before the thermodynamic framework begins.

## Four experiments, one law

Between 1662 and 1811, four independent lines of experiment each revealed a simple proportionality between the state variables of a gas.
None of the experimentalists had a molecular picture; they were simply measuring.

**Boyle's law (1662):**
Robert Boyle measured the pressure and volume of a fixed amount of gas at constant temperature and found, doubling the pressure halves the volume:

$$PV = \text{const} \quad (T,\, n \text{ fixed})$$

The product $PV$ is independent of how it is distributed between $P$ and $V$.


**Charles's law (1787):**
Jacques Charles found that the volume of a gas at constant pressure is proportional to its temperature:

$$\frac{V}{T} = \text{const} \quad (P,\, n \text{ fixed})$$

Cooling a gas shrinks it; heating expands it, linearly in $T$.

**Gay-Lussac's law (1808):**
Joseph Louis Gay-Lussac measured gases at constant volume and found pressure proportional to temperature:

$$\frac{P}{T} = \text{const} \quad (V,\, n \text{ fixed})$$

A sealed container whose contents are heated sees its pressure rise in proportion to the absolute temperature.

**Avogadro's law (1811):**
Amedeo Avogadro proposed that equal volumes of any gas at the same temperature and pressure contain equal numbers of molecules, so volume is proportional to the amount of substance:

$$\frac{V}{n} = \text{const} \quad (T,\, P \text{ fixed})$$

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-gas-laws

import numpy as np
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 4, figsize=(14, 3.5))

V = np.linspace(0.5, 3, 200)
for C, label in [(1, "$T_1$"), (2, "$T_2$"), (3, "$T_3$")]:
    axes[0].plot(V, C / V, label=label)
axes[0].set_xlabel("Volume $V$")
axes[0].set_ylabel("Pressure $P$")
axes[0].set_title("Boyle's law\n$PV = $ const")
axes[0].legend(fontsize=8)

T = np.linspace(100, 400, 200)
for C, label in [(0.5, "$P_1$"), (1.0, "$P_2$"), (1.5, "$P_3$")]:
    axes[1].plot(T, C * T / 200, label=label)
axes[1].set_xlabel("Temperature $T$  [K]")
axes[1].set_ylabel("Volume $V$")
axes[1].set_title("Charles's law\n$V/T = $ const")
axes[1].legend(fontsize=8)

for C, label in [(0.5, "$V_1$"), (1.0, "$V_2$"), (1.5, "$V_3$")]:
    axes[2].plot(T, C * T / 200, label=label)
axes[2].set_xlabel("Temperature $T$  [K]")
axes[2].set_ylabel("Pressure $P$")
axes[2].set_title("Gay-Lussac's law\n$P/T = $ const")
axes[2].legend(fontsize=8)

n = np.linspace(0, 3, 200)
for C, label in [(1, "$T_1, P_1$"), (2, "$T_2, P_2$"), (3, "$T_3, P_3$")]:
    axes[3].plot(n, C * n / 3, label=label)
axes[3].set_xlabel("Amount $n$  [mol]")
axes[3].set_ylabel("Volume $V$")
axes[3].set_title("Avogadro's law\n$V/n = $ const")
axes[3].legend(fontsize=8)

for ax in axes[1:]:
    ax.set_ylabel("")

fig.tight_layout()
```

```{figure} #cell-gas-laws
:name: fig-gas-laws

The four empirical gas laws. Each panel holds two variables fixed and varies the third, revealing a simple proportionality in every case.
```

## Unification: $PV = nRT$

The four laws are not independent.
Combining them into a single equation:

$$PV = nRT$$

where $R = 8.314\ \mathrm{J/(mol\,K)}$ is the **gas constant**.
Every one of the four laws is a special case: fixing $T$ and $n$ recovers Boyle's law; fixing $P$ and $n$ recovers Charles's law; and so on.

$R$ is related to Boltzmann's constant $k_B = 1.380 \times 10^{-23}$ J/K by $R = N_A k_B$, where $N_A = 6.022 \times 10^{23}\ \text{mol}^{-1}$ is Avogadro's number.
$k_B$ is the per-particle version; $R$ is the per-mole version.

```{admonition} Moles: the chemist's counting unit
:class: note

At the molecular scale, statistical mechanics counts $N$ discrete particles; $N_A$ converts particle counts to moles (the chemist's counting unit); $M_i$ converts moles to mass (what a balance measures).
The mole is the interface between the quantum/statistical world and bulk engineering quantities; this is why molar concentrations are natural throughout this book.
```

## Composition variables for mixtures

The ideal gas law for a single component becomes a set of relations for a mixture.
Several variables describe the composition of gas mixtures (and, by extension, solutions).

**Mole fraction.** For a mixture of $N$ components, the **mole fraction** of species $i$ is:

$$x_i = \frac{n_i}{\sum\limits_{j=1}^N n_j}$$

It is dimensionless and satisfies $0 \leq x_i \leq 1$ with $\sum_i x_i = 1$.

**Dalton's law (partial pressure).** For an ideal gas mixture, each component behaves as if it alone occupied the container at pressure $P_i$:

$$P_i = x_i P$$

The **partial pressure** $P_i$ is the pressure component $i$ would exert if it occupied the same volume at the same temperature.
The total pressure is the sum of all partial pressures: $P = \sum_i P_i$.

**Molar volume.** The molar volume is the volume per mole of mixture:

$$V_m = \frac{V}{n} = \frac{V}{\sum_i n_i} = \frac{RT}{P}$$ 

For an ideal gas at fixed $T$ and $P$, $V_m$ is the same for any mixture and equals the molar volume of the pure components at the same conditions.

**Molar concentration.** The molar concentration of species $i$ is the amount per unit volume:

$$c_i = \frac{n_i}{V} \quad \mathrm{[mol/m^3\ or\ mol/L]}$$

Molar concentrations are the primary composition variables used throughout this book.
Where "concentration" appears without qualification, it means molar concentration.

**Mass concentration.** The mass concentration of species $i$ is:

$$\rho_i = \frac{m_i}{V} = c_i M_i \quad \mathrm{[kg/m^3\ or\ g/L]}$$

where $M_i$ is the molar mass of species $i$ (kg/mol).
Mass concentration is the engineering variant and appears where mass balance or density calculations are needed.

**Total mass density.** The total mass density of the mixture is:

$$\rho = \sum_i \rho_i = \sum_i c_i M_i$$ 

For liquids, $\rho$ is nearly constant; for gases, it depends on $T$ and $P$ via the ideal gas law.

---

The ideal gas law describes how $P$, $V$, $T$, and $n$ relate at equilibrium, but not how the energy of a system changes or what drives a system toward equilibrium.
It also rests on the assumption that molecules do not interact; where that assumption fails, corrections are needed; a later chapter on non-idealities covers this in detail.
The next chapter introduces the thermodynamic framework that answers both questions, with $PV = nRT$ as one of its key building blocks.
