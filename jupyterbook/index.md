---
kernelspec:
  name: python3
  display_name: Python 3 (ipykernel)
---

# Thermodynamic Reactions in CADET

CADET simulates chromatographic separations by solving coupled mass balances in packed-bed columns.
Reactions between components (protonation, complex formation, aggregation) are currently described using mass-action kinetics:

$$
r = k^f \prod_i c_i^{\nu_{f,i}} - k^r \prod_i c_i^{\nu_{b,i}},
$$

where $c_i$ are species concentrations and $\nu_{f,i}$, $\nu_{b,i}$ are the kinetic orders in the forward and reverse directions.
In the standard mass action law these equal the stoichiometric coefficients; CADET also allows them to be set independently, which matters when the apparent kinetic order differs from stoichiometry.
The ratio $k^f/k^r$ is an empirical parameter, and non-ideal effects are absorbed into the rate constants.
This works well for dilute, isothermal systems. Real chromatographic conditions routinely violate both assumptions.

**Ionic strength.**
Buffer concentrations of 50 to 1000 mM cause strong electrostatic screening between solutes, stationary phase, and solvent.
The effective driving force of a reaction is not the same as the measured bulk concentration.

**pH.**
As buffer fronts propagate through the column, local acid-base equilibria shift and change the protonation state of solutes.
Reaction rates couple to transport and local composition.
The effect can be represented with kinetic formulations based on conserved moieties, but these get stiff for high pKa values and obscure the underlying equilibrium structure.

**Temperature.**
Temperature influences both thermodynamic equilibrium (adsorption capacity, reaction extent) and kinetics (rate constants, activation energy barriers).
Equilibrium constants shift with $T$ through the van't Hoff relation; rate constants follow Arrhenius exponentially.
Constant parameters capture neither.

The current CADET implementation cannot express any of this.
This book builds the thermodynamic foundation needed to fix it.
Not as a patch, but as a principled extension that makes the reaction framework self-consistent.

The key insight is that thermodynamics is not a collection of special-case equations.
It is a coherent mathematical structure: a state function, its natural variables, and rules for differentiating under constraints.
Named results — van't Hoff, Henderson-Hasselbalch, Debye-Hückel — are not independent empirical patches but consequences of differentiating the same underlying function under different conditions.
That structure is what makes the extension principled: thermodynamic consistency follows automatically from building the framework correctly, rather than being enforced case by case.

In thermodynamics, equilibrium is defined through chemical potentials, not rate constants.
At equilibrium, reactions satisfy

$$
\sum_i \nu_i \mu_i = 0,
$$

which links stoichiometry to the underlying driving forces.
For ideal and non-ideal solutions, this defines an equilibrium constant that depends on temperature and composition through activities rather than concentrations.
The ratio of forward and reverse rate coefficients is not arbitrary but must satisfy

$$
\frac{k^f(T)}{k^r(T)} = K(T),
$$

where $K(T)$ is determined by thermodynamics.
This does not remove kinetic parameter estimation: it restricts kinetic models to a thermodynamically consistent class, while preserving flexibility in rate law form and prefactor estimation.
The target rate law is

$$
r = k^f(T) \prod_i a_i^{e_i^f} - k^r(T) \prod_i a_i^{e_i^r}, \qquad \frac{k^f(T)}{k^r(T)} = K(T),
$$

where $a_i = \gamma_i c_i / c^\circ$ are activities, $e_i^f$ and $e_i^r$ are kinetic exponents (equal to $|\nu_i|$ for elementary reactions), $K(T)$ is the thermodynamic equilibrium constant, and the ratio $k^f / k^r$ is never a free parameter.
pH enters not through this equation directly, but through speciation: acid-base reactions are modelled as fast-equilibrium constraints, and the activities that feed into the rate law already reflect the local protonation state.

The book is structured in four parts.
Part 1 arrives at thermodynamics from statistical mechanics, deriving entropy, temperature, and the Boltzmann factor from first principles, so the design decisions in the later parts are motivated rather than assumed.
Part 2 builds the thermodynamic framework: Gibbs energy, chemical potential, and activity for non-ideal solutions.
Part 3 develops the theory of chemical reactions, from equilibrium constants and acid-base speciation through kinetics and the mass action law.
Part 4 implements the complete framework as a Python library and connects it to the CADET solver interface.

This is not a gentle introduction to chemistry.
It is a principle-first reconstruction of reaction thermodynamics for computational modeling, aimed at graduate students and researchers who already have some exposure to calculus, equilibrium chemistry, and the concept of Gibbs energy.
The goal is a structurally correct and self-consistent foundation, not a minimal explanation.
Readers who want a first exposure to thermodynamics should start at Part 1; those who already know the Gibbs energy can enter at Part 2 or later (see reading paths below).

```{admonition} Reading paths
:class: tip

The book supports multiple entry points depending on background.

| Background                             | Recommended entry point                             |
| -------------------------------------- | --------------------------------------------------- |
| New to thermodynamics                  | {ref}`statistical-mechanics`                        |
| Know ideal gas law                     | {ref}`laws-of-thermodynamics`                       |
| Know Gibbs energy                      | {ref}`nonidealities`                                |
| Know MAL, want thermodynamic grounding | {ref}`equilibrium` → {ref}`multicomponent` → Part 4 |
| CADET developer, know MAL              | {ref}`multicomponent` → Part 4                      |

Note for the Part 4 entry reader: {ref}`nonidealities` (activity, Debye-Hückel, Davies) is a prerequisite for the activity and acid-base chapters in Part 4.
```

