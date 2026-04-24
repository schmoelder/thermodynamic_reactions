---
kernelspec:
  name: python3
  display_name: Python 3 (ipykernel)
---

(reactions)=
# Part 3: Chemical Reactions

Parts 1 and 2 built the thermodynamic toolkit from first principles.
Part 3 applies it to chemical reactions, in two arcs.
Part 4 implements the resulting framework in CADET.

**Equilibrium arc.**
The central new quantity is the **extent of reaction** $\xi$, a single coordinate that tracks how far a reaction has proceeded.
The Gibbs energy is a function of $\xi$: its slope $\Delta_r G = dG/d\xi$ is the reaction Gibbs energy, and its minimum defines equilibrium.
Setting $\Delta_r G = 0$ yields the equilibrium constant $K$; the relation $\Delta_r G^\circ = -RT \ln K$ connects $K$ to tabulated standard enthalpies, entropies, and heat capacities.
The van't Hoff equation gives how $K$ shifts with temperature, and two concrete applications follow: adsorption equilibria (Langmuir isotherm) and acid-base equilibria (pH, buffers, speciation).

**Kinetics arc.**
Kinetics adds the time dimension.
A mechanism decomposes a reaction into elementary steps and determines which step limits the overall rate.
Thermodynamic consistency ties forward and reverse rate constants to $K$: the same equilibrium that the equilibrium arc predicts is the attractor of the kinetic dynamics.
The Arrhenius equation, rooted in the Boltzmann factor from @maxwell-boltzmann, describes how each rate constant varies with temperature.
When a finite catalytic pool is involved, the rate saturates at high substrate concentrations via the same finite-site logic as the Langmuir isotherm.

```{admonition} Key results
:class: tip

**Reaction Gibbs energy.**
The driving force for a reaction is the slope of $G$ along the extent of reaction $\xi$:

$$\Delta_r G = \sum_i \nu_i \mu_i = \Delta_r G^\circ + RT \ln Q$$

Equilibrium is reached when $\Delta_r G = 0$, giving $Q = K$.

**Equilibrium constant.**
The standard Gibbs energy of reaction is related to the equilibrium constant by:

$$\Delta_r G^\circ = -RT \ln K$$

**Van't Hoff equation.**
The temperature dependence of $K$ is governed by the reaction enthalpy:

$$\frac{d \ln K}{dT} = \frac{\Delta_r H^\circ}{RT^2}$$

**Kinetics.**
Reaction rates are governed by reversible rate laws of the form $r = k_f \prod_i c_i^{\nu_i^f} - k_r \prod_i c_i^{\nu_i^b}$.
Thermodynamic consistency requires that the forward and reverse rate constants are not independent:

$$\frac{k_f}{k_r} = K$$

**Arrhenius equation.**
Each rate constant follows an activated temperature dependence:

$$k = A\,e^{-E_a/RT}$$

Consistency then requires $E_{a,f} - E_{a,r} = \Delta_r H^\circ$.
```

**Chapters in this part:**

*Equilibrium arc*

- @reaction-coordinates: extent of reaction $\xi$, stoichiometric coefficients, physical range
- @reaction-gibbs-energy: $G(\xi)$, reaction Gibbs energy $\Delta_r G$, standard state, reaction quotient $Q$
- @equilibrium: $K$, mass-action form, tabulated data
- @equilibrium-temperature: Kirchhoff relations, van't Hoff equation
- @adsorption: chromatographic adsorption, Langmuir isotherm as finite-site equilibrium
- @acid-base: proton transfer, $K_a$, pH, Henderson-Hasselbalch
- @speciation-buffers: Bjerrum diagrams, buffer capacity, ionic strength corrections

*Kinetics arc*

- @kinetics: reaction rates, elementary steps, mechanisms, rate-determining step
- @mass-action-law: mass-action rate law, thermodynamic consistency condition $k_f/k_r = K$
- @kinetics-temperature: Arrhenius equation, transition state theory, barrier-enthalpy relation
- @saturation: Michaelis-Menten kinetics, Monod equation, connection to Langmuir
