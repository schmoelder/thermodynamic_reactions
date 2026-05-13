---
kernelspec:
  name: python3
  display_name: Python 3 (ipykernel)
---

(equilibrium-temperature)=
# Temperature Dependence of Equilibrium

The equilibrium constant $K$ is tabulated at 298 K, but equilibrium systems rarely operate at exactly 25 °C.
Chromatographic adsorption columns are commonly run at 4--37 °C, and pKa values of buffers shift measurably between lab and physiological temperature.
This chapter derives how $K(T)$ varies via the **Kirchhoff relations** and the **van't Hoff equation**.
The result is general: it applies uniformly to any equilibrium constant derived from $\Delta_r G^\circ = -RT \ln K$, including the adsorption constant $K_\text{ads}$ from @adsorption and the acid dissociation constant $K_a$ introduced in @acid-base.

## Standard thermodynamic data

The standard reaction Gibbs energy is computed from tabulated standard quantities:

$$\Delta_r G^\circ = \Delta_r H^\circ - T\Delta_r S^\circ$$

$$\Delta_r H^\circ = \sum_i \nu_i \Delta_f H^\circ_i, \qquad \Delta_r S^\circ = \sum_i \nu_i S^\circ_i$$

Standard enthalpies of formation $\Delta_f H^\circ$ and standard molar entropies $S^\circ$ are measured at $T^\circ = 298.15\,\text{K}$ and $P^\circ = 1\,\mathrm{bar}$.
The third law (@laws-of-thermodynamics) provides absolute values of $S^\circ$ via integration of $C_p / T$ from 0 K; unlike enthalpies, entropies have a natural zero.

```{admonition} Intuition
:class: tip
$\Delta_r H^\circ$ links a reaction to its heat release or absorption.
Combustion of one mole of methane releases 890 kJ, enough to heat roughly 3 L of water from room temperature to boiling.
Reactions span orders of magnitude: weak interactions (e.g., van der Waals) involve ~1–10 kJ/mol, while strong acid-base neutralizations release ~50–60 kJ/mol.
The van't Hoff equation quantifies this: a small $|\Delta_r H^\circ|$ yields a shallow slope in a $\ln K$ vs $1/T$ plot, making $K$ nearly temperature-independent.
```

## Kirchhoff relations

Tabulated values are given at 298.15 K.
To find $\Delta_r H^\circ$ at another temperature, integrate the heat capacity difference:

$$\Delta_r H^\circ(T_2) = \Delta_r H^\circ(T_1) + \int_{T_1}^{T_2} \Delta_r C_p\, dT$$

where $\Delta_r C_p = \sum_i \nu_i C_{p,i}$.
Similarly, $\Delta_r S^\circ(T_2) = \Delta_r S^\circ(T_1) + \int_{T_1}^{T_2} (\Delta_r C_p / T)\, dT$.
When $\Delta_r C_p$ is small or the temperature range is narrow, $\Delta_r H^\circ$ and $\Delta_r S^\circ$ are approximately constant and both integrals can be dropped.
The Kirchhoff relations correct the thermodynamic data; the van't Hoff equation (below) then uses that corrected $\Delta_r H^\circ$ to predict how $K$ shifts with temperature.

## Van't Hoff equation

Since $\Delta_r G^\circ = -RT \ln K$, the temperature dependence of $K$ follows directly from the Gibbs-Helmholtz equation derived in @chemical-potential.
Differentiating $\ln K = -\Delta_r G^\circ / RT$ with respect to $T$:

$$\frac{d \ln K}{dT} = \frac{\Delta_r H^\circ}{RT^2}$$

This is the **van't Hoff equation** {cite:p}`vanthoff1884`.
In integrated form, assuming $\Delta_r H^\circ$ is approximately constant over the temperature range:

$$\ln \frac{K(T_2)}{K(T_1)} = -\frac{\Delta_r H^\circ}{R}\left(\frac{1}{T_2} - \frac{1}{T_1}\right)$$

An exothermic reaction ($\Delta_r H^\circ < 0$) has $K$ decreasing with temperature: Le Chatelier's principle in quantitative form.
A van't Hoff plot of $\ln K$ versus $1/T$ has slope $-\Delta_r H^\circ / R$, which is how $\Delta_r H^\circ$ is measured experimentally when direct calorimetry is impractical.

```{admonition} Intuition
:class: tip
The van't Hoff equation links temperature to the position of equilibrium: heating favours the endothermic direction, cooling the exothermic one.
The Haber process ($\ce{N2 + 3 H2 <=> 2 NH3}$, $\Delta_r H^\circ = -92\,\mathrm{kJ/mol}$) illustrates the scale: between 25 °C and 500 °C, $K$ drops by several orders of magnitude, drastically reducing equilibrium yield.
Industrial plants operate at 400--500 °C as a compromise between reaction rate and equilibrium yield.
The same equation governs $K_\text{ads}$ for adsorption (@adsorption) and $K_a$ for acid-base equilibria (@acid-base), with each reaction's $\Delta_r H^\circ$ setting the slope.
```

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-vant-hoff

import numpy as np
import matplotlib.pyplot as plt
from reactions.plots import setup_figure, COLORS

R = 8.314  # J/(mol·K)
T = np.linspace(280, 400, 300)

cases = [
    (r"exothermic  ($\Delta_r H^\circ = -40\,\mathrm{kJ/mol}$)", -40e3, -100.0, "C0"),
    (r"endothermic  ($\Delta_r H^\circ = +40\,\mathrm{kJ/mol}$)", +40e3, +100.0, "C3"),
]

fig, axes = setup_figure(1, 2)

for ax, (label, dH, dS, color) in zip(axes, cases):
    dG = dH - T * dS
    K = np.exp(-dG / (R * T))
    ax.plot(1e3 / T, np.log(K), color=color, lw=2)
    ax.set_xlabel(r"$1/T$  [$10^{-3}$ K$^{-1}$]")
    ax.set_ylabel(r"$\ln K$")
    ax.set_title(label)
    slope = -dH / R
    ax.annotate(
        f"slope $= -\\Delta_r H^\\circ / R$\n$= {slope / 1e3:.1f}$ kK",
        xy=(1e3 / T[150], np.log(K[150])),
        xytext=(20, -30 if dH < 0 else 30),
        textcoords="offset points",
        fontsize=9,
        arrowprops=dict(arrowstyle="->", color="gray"),
    )

fig.tight_layout()
```

```{figure} #cell-vant-hoff
:name: fig-vant-hoff

Van't Hoff plots ($\ln K$ versus $1/T$) for an exothermic (left) and endothermic (right) reaction.
The slope equals $-\Delta_r H^\circ / R$: exothermic reactions have $K$ decreasing with temperature (negative slope), while endothermic reactions have $K$ increasing (positive slope).
```

---

The next chapter applies the same equilibrium framework to acid-base reactions, where $K_a$ is the equilibrium constant for proton transfer and the van't Hoff equation governs its shift with temperature (@acid-base).
