---
kernelspec:
  name: python3
  display_name: Python 3 (ipykernel)
---

(kinetics-temperature)=
# Temperature Dependence of Rate Constants

@kinetics established how mechanism determines the form of a rate law.
@mass-action-law showed that thermodynamic consistency requires $k_f/k_r = K$.
This chapter addresses how the rate constant $k$ itself depends on temperature, and what that implies for the barrier heights.

## Activation energy

The connection between thermodynamics and kinetics runs through the Maxwell-Boltzmann distribution.
@maxwell-boltzmann showed that at temperature $T$, the fraction of particles with energy exceeding a threshold $E_a$ is proportional to $e^{-E_a/RT}$, where $R = N_A k_B$ is the gas constant.
For a reaction to occur, colliding molecules must collectively carry enough energy to overcome the **activation energy** $E_a$: the energy barrier separating reactants from products.

@fig-arrhenius-barrier illustrates this barrier on a potential energy diagram.
Reactants approach from the left and must surmount the barrier to reach products.
At the top sits the **transition state** $[\mathrm{AB}]^\ddagger$, an unstable arrangement through which all reactive trajectories must pass.
The reverse barrier is $E_a' = E_a - \Delta_r H^\circ$: smaller than $E_a$ for an exothermic reaction and larger for an endothermic one.

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-arrhenius-barrier

import numpy as np
import matplotlib.pyplot as plt

Ea_fwd = 60.0
drH = -25.0

n = 500
x_left  = np.linspace(0, 0.5, n // 2)
x_right = np.linspace(0.5, 1.0, n // 2)

t_left  = np.linspace(0, np.pi, n // 2)
t_right = np.linspace(0, np.pi, n // 2)

y_left  = Ea_fwd * (1 - np.cos(t_left)) / 2
y_right = Ea_fwd + (drH - Ea_fwd) * (1 - np.cos(t_right)) / 2

x = np.concatenate([x_left, x_right[1:]])
y = np.concatenate([y_left, y_right[1:]])

fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(x, y, color="#1c4f8a", linewidth=2.5)

ax.hlines(0.0, -0.1, 0.32, colors="gray", linestyles="--", linewidth=1.0)
ax.hlines(drH, 0.68, 1.1, colors="gray", linestyles="--", linewidth=1.0)

ax.annotate(
    "",
    xy=(0.30, Ea_fwd),
    xytext=(0.30, 0.0),
    arrowprops=dict(arrowstyle="<->", color="C1", lw=1.5),
)
ax.text(0.28, Ea_fwd / 2, r"$E_a$",
        ha="right", va="center", fontsize=12, color="C1")

ax.annotate(
    "",
    xy=(0.70, drH),
    xytext=(0.70, 0.0),
    arrowprops=dict(arrowstyle="<->", color="C2", lw=1.5),
)
ax.text(0.72, drH / 2, r"$\Delta_r H^\circ$",
        ha="left", va="center", fontsize=11, color="C2")

ax.text(0.5, Ea_fwd + 3, r"$[\mathrm{AB}]^\ddagger$",
        ha="center", va="bottom", fontsize=12)
ax.text(-0.05, 2, "Reactants", ha="right", va="bottom", fontsize=10, color="gray")
ax.text(1.05, drH + 2, "Products", ha="left", va="bottom", fontsize=10, color="gray")

ax.set_xlim(-0.18, 1.18)
ax.set_ylim(drH - 15, Ea_fwd + 18)
ax.set_xlabel("Reaction coordinate", fontsize=11)
ax.set_ylabel("Potential energy", fontsize=11)
ax.set_xticks([])
ax.set_yticks([])
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
fig.tight_layout()
```

```{figure} #cell-arrhenius-barrier
:name: fig-arrhenius-barrier

Potential energy along the reaction coordinate for an exothermic reaction.
Reactants must supply energy $E_a$ to reach the transition state $[\mathrm{AB}]^\ddagger$.
Products sit below reactants by $|\Delta_r H^\circ|$, so the reverse barrier is $E_a' = E_a - \Delta_r H^\circ > E_a$.
```

Raising the temperature inflates the high-energy tail of the Maxwell-Boltzmann distribution, increasing the fraction of collisions that are reactive.
This is the microscopic origin of the observation that reaction rates increase with temperature.

## Arrhenius equation

Svante Arrhenius (1889) recognised that rate constants follow:

$$
k = A\,e^{-E_a/RT}
$$

where $A$ is the **pre-exponential factor** (frequency factor).
The exponential term is the Boltzmann-weighted fraction of collisions that are sufficiently energetic; the Arrhenius equation follows directly from the Maxwell-Boltzmann distribution and the concept of an activation barrier.

Taking the logarithm:

$$
\ln k = \ln A - \frac{E_a}{R} \cdot \frac{1}{T}
$$

A plot of $\ln k$ against $1/T$ is an **Arrhenius plot** (@fig-arrhenius): linear with slope $-E_a/R$ and intercept $\ln A$.
Rate constants measured at several temperatures determine $E_a$ from the slope.

```{admonition} Intuition
:class: tip
The Arrhenius equation links activation energy to the temperature sensitivity of a rate constant.
For a typical reaction with $E_a \approx 50\,\mathrm{kJ/mol}$, a 10 °C increase near room temperature roughly doubles the rate.
For $E_a \approx 100\,\mathrm{kJ/mol}$, the same increase leads to an approximate fourfold rise.
This exponential dependence reflects the Boltzmann-weighted fraction of molecules that have enough energy to cross the barrier, as derived in @maxwell-boltzmann.
Thermodynamic consistency further constrains forward and reverse reactions, linking their activation energies to the reaction enthalpy, as discussed in the section below.
```

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-arrhenius

import numpy as np
import matplotlib.pyplot as plt

R = 8.314
T = np.linspace(300, 1200, 400)
T_inv = 1000 / T

Ea_vals = [20_000, 40_000, 60_000, 80_000]
colors  = ["C0", "C1", "C2", "C3"]

fig, ax = plt.subplots(figsize=(7, 4.5))
for Ea, color in zip(Ea_vals, colors):
    ln_k  = -Ea / R / T
    label = rf"$E_a = {Ea // 1000:.0f}\ \mathrm{{kJ/mol}}$"
    ax.plot(T_inv, ln_k, color=color, linewidth=2.0, label=label)

ax.set_xlabel(r"$1000/T\ [\mathrm{K}^{-1}]$", fontsize=11)
ax.set_ylabel(r"$\ln(k/A)$", fontsize=11)
ax.legend(fontsize=9)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
fig.tight_layout()
```

```{figure} #cell-arrhenius
:name: fig-arrhenius

Arrhenius plot: $\ln(k/A)$ versus $1000/T$ for four activation energies.
Each curve is linear with slope $-E_a/R$; steeper slopes indicate stronger temperature sensitivity.
```

## Transition State Theory

The pre-exponential factor $A$ receives a mechanistic interpretation in Transition State Theory (TST, Eyring-Evans-Polanyi, 1935).
TST postulates a quasi-equilibrium between reactants and an activated complex at the transition state:

$$
k = \frac{k_BT}{h}\,e^{-\Delta G^\ddagger/RT}
$$

where $h$ is Planck's constant and $\Delta G^\ddagger = \Delta H^\ddagger - T\Delta S^\ddagger$ is the **activation Gibbs energy**.
The factor $k_BT/h \approx 6\times10^{12}\ \mathrm{s}^{-1}$ at room temperature sets the attempt frequency.
The entropic term $e^{\Delta S^\ddagger/R}$ accounts for how ordered the transition state is relative to the reactants and contributes to the empirical pre-exponential $A$.
TST explains why $E_a$ is an enthalpy, not a Gibbs energy: entropy enters separately through $A$.

## Barrier heights and reaction enthalpy

@mass-action-law established that thermodynamic consistency requires $k_f/k_r = K$.
With both rate constants following Arrhenius, $k_f = A_f e^{-E_{a,f}/RT}$ and $k_r = A_r e^{-E_{a,r}/RT}$, their ratio is $\propto e^{-(E_{a,f} - E_{a,r})/RT}$.
The van't Hoff equation gives $K(T) \propto e^{-\Delta H^\circ_r/RT}$ when $\Delta H^\circ_r$ is constant.
Matching the two expressions requires:

$$
E_{a,f} - E_{a,r} = \Delta H^\circ_r
$$

The barrier heights must differ by exactly the reaction enthalpy.
This is visible in @fig-arrhenius-barrier: the forward and reverse barriers are separated by $|\Delta_r H^\circ|$, and the exothermic direction has the lower barrier.
How to maintain $k_f/k_r = K(T)$ automatically at every temperature is taken up in @implementation-kinetics.

---

The kinetics arc ends with one more phenomenon that power-law rate laws cannot describe: saturation, where adding more reactant no longer increases the rate because the catalytic capacity is fully occupied.
