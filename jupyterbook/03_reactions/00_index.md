---
kernelspec:
  name: python3
  display_name: Python 3 (ipykernel)
---

(reactions)=
# Part 3: Chemical Reactions

Parts 1 and 2 built the thermodynamic toolkit: entropy, Gibbs energy, chemical potential, activity.
Part 3 applies it to chemical reactions and to any phenomenon that reduces to Gibbs energy minimization under a structural constraint.

Two constraint types organize all phenomena here.
**Stoichiometric constraints** tie composition to a single scalar $\xi$; equilibrium is $\Delta_r G = 0$.
**Finite-site constraints** partition a conserved pool of sites between free and occupied states; the same occupancy logic produces saturation in both equilibrium (Langmuir isotherm) and kinetics (Michaelis-Menten).

Stage I covers equilibrium under both constraint types, then their applications (temperature dependence, acid-base, speciation).
Stage II covers dynamics: how composition relaxes toward equilibrium, and how the same finite-site constraint that gives the Langmuir isotherm yields Michaelis-Menten kinetics.

```{admonition} Key results
:class: tip

**Reaction driving force.**
$$\Delta_r G = \sum_i \nu_i \mu_i = \Delta_r G^\circ + RT \ln Q$$
Equilibrium is reached when $\Delta_r G = 0$, giving $Q = K$.

**Equilibrium constant.**
$$\Delta_r G^\circ = -RT \ln K$$

**Thermodynamic consistency.**
Forward and reverse rate constants are not independent:
$$\frac{k_f}{k_r} = K$$

**Temperature dependence.**
Both equilibrium constants and rate constants shift with temperature:
$$\frac{d \ln K}{dT} = \frac{\Delta_r H^\circ}{RT^2} \qquad k = A\,e^{-E_a/RT}$$
Consistency requires $E_{a,f} - E_{a,r} = \Delta_r H^\circ_r$.

**Saturation (finite-site constraints).**
The same occupancy constraint governs both equilibrium and kinetic saturation:
$$q = \frac{q_\text{max}\,K_\text{ads}\,c}{1 + K_\text{ads}\,c} \qquad r = \frac{V_\text{max}\,[S]}{K_m + [S]}$$
```


**Chapters in this part:**

*Stage I — Equilibrium*

- **@reaction-coordinates**: extent of reaction $\xi$, stoichiometric coefficients, physical range
- **@reaction-gibbs-energy**: $G(\xi)$, reaction Gibbs energy $\Delta_r G$, reaction quotient $Q$, Le Chatelier's principle
- **@equilibrium**: equilibrium constant $K$, mass-action form, tabulated data
- **@adsorption**: finite-site equilibrium, Langmuir isotherm
- **@equilibrium-temperature**: Kirchhoff relations, van't Hoff equation
- **@acid-base**: proton transfer, $K_a$, pH, Henderson-Hasselbalch
- **@speciation-buffers**: Bjerrum diagrams, buffer capacity, ionic strength corrections

*Stage II — Dynamics*

- **@kinetics**: reaction rates, elementary steps, mechanisms, rate-determining step
- **@mass-action-law**: mass-action rate law, thermodynamic consistency $k_f/k_r = K$
- **@kinetics-temperature**: Arrhenius equation, transition state theory, barrier-enthalpy relation
- **@saturation**: Michaelis-Menten kinetics, Monod equation, connection to Langmuir
