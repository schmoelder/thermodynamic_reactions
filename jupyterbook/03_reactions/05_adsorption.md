---
kernelspec:
  name: python3
  display_name: Python 3 (ipykernel)
---

(adsorption)=
# Adsorption Equilibrium

This chapter applies the equilibrium framework to **adsorption**: the distribution of solute between a mobile phase and a stationary phase.
In chromatography, adsorption is treated as a quasi-equilibrium process: solute molecules rapidly partition between the flowing mobile phase and the porous stationary phase.

## Chromatographic adsorption equilibrium

Chromatographic processes rely on the equilibrium between solute in the mobile phase and solute adsorbed onto the stationary phase:

$$\text{A}_\text{mobile} \rightleftharpoons \text{A}_\text{adsorbed}$$

**Heat of adsorption.**
The retention factor $k'$ (measurable from chromatogram peak positions) is proportional to the adsorption equilibrium constant $K_\text{ads}$.
Measuring $k'$ at several temperatures and plotting $\ln k'$ versus $1/T$ gives a van't Hoff plot with slope $-\Delta H_\text{ads} / R$ and intercept $\Delta S_\text{ads} / R$.
This provides the heat of adsorption $\Delta H_\text{ads}$ and the adsorption entropy without any calorimetric measurement.

**Langmuir isotherm.**
At a single temperature, $K_\text{ads}$ determines how the adsorbed amount $q$ depends on solution concentration $c$.
The stationary phase has a **finite number of adsorption sites** with maximum loading $q_\text{max}$.
Each site is either empty or occupied by one solute molecule; the two states are in equilibrium with the solution.
Treating the site occupancy as a chemical equilibrium (exactly as in @equilibrium) gives the **Langmuir isotherm**:

$$q = \frac{q_\text{max}\, K_\text{ads}\, c}{1 + K_\text{ads}\, c}$$

At low loading ($K_\text{ads} c \ll 1$) this is linear (**Henry's law region**, $q \approx H c$ with $H = q_\text{max} K_\text{ads}$); at high loading it saturates at $q_\text{max}$ because all sites are occupied.
The saturation is an equilibrium consequence of finite site capacity, not a kinetic phenomenon.

```{admonition} Intuition
:class: tip

The Langmuir isotherm behaves like seats on a train.
In the Henry region the train is nearly empty: every additional passenger (solute molecule) finds a seat immediately, and the number seated is proportional to the number arriving ($q \propto c$).
As loading increases, fewer seats remain and each additional molecule is less likely to adsorb.
At saturation the train is full: $q = q_\text{max}$ regardless of how many more molecules are in solution.
A large $K_\text{ads}$ means the train fills at lower platform concentration.
```

The same finite-site logic reappears in kinetics: when an enzyme or cell population acts as the catalytic pool, the rate saturates in exactly the same way (@saturation).

**Bridge between phase and reaction equilibrium.**
Adsorption equilibrium is a special case that sits between pure phase equilibrium (@phase-equilibrium) and chemical reaction equilibrium (@equilibrium):
- Like phase equilibrium: it describes partitioning between two regions (mobile and stationary phases)
- Like reaction equilibrium: it can be written as a quasi-reaction with its own equilibrium constant $K_\text{ads}$


---

The next chapter shifts to acid-base equilibria (@acid-base), the most practically important chemical equilibria for buffer design and pH control in chromatography.
