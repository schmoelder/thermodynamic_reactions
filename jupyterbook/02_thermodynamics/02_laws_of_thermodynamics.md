---
kernelspec:
  name: python3
  display_name: Python 3 (ipykernel)
---

(laws-of-thermodynamics)=
# Laws of Thermodynamics

@ideal-gas introduced the equation of state $PV = nRT$, which relates the macroscopic variables of a gas at equilibrium.
This chapter goes deeper: what does temperature mean, how does energy change, and what determines the direction of spontaneous change?
The four laws of thermodynamics answer these questions and form the foundation of chemical thermodynamics.

## The zeroth law

Before energy and entropy can be discussed, temperature must be well defined.
The **zeroth law of thermodynamics** provides this foundation: if system A is in thermal equilibrium with system B, and system B is in thermal equilibrium with system C, then A and C are also in thermal equilibrium with each other.

This transitivity is what makes a thermometer meaningful and what allows temperature to be assigned to a system rather than only to a pair of systems.
In @entropy, temperature emerged from the condition $\partial S_1/\partial U_1 = \partial S_2/\partial U_2$ for two systems sharing energy; the zeroth law guarantees that this condition is consistent across any number of systems in contact.
Without it, the concept of a uniform temperature throughout a system would have no logical foundation.
The law was recognized after the first and second laws were already named, hence the retroactive numbering.

## The first law

The internal energy $U$ of a closed system can change in three ways:

- **Heat** $\delta Q$: energy flows across the boundary because of a temperature difference.
- **Work** $\delta W$: the system changes volume against an external pressure.
- **Matter**: particles of species $i$ enter or leave, each carrying chemical energy $\mu_i$.

This gives the **first law of thermodynamics**:

$$
dU = \delta Q - \delta W + \sum_i \mu_i\,dn_i
$$

Energy is conserved: any change in $U$ is accounted for by heat, work, or matter transfer.
The first law sets the bookkeeping but says nothing about which direction a process goes.

## The second law

Two processes can both conserve energy yet only one is observed to happen spontaneously: heat flows from hot to cold, gases expand into vacuum, reactions proceed toward equilibrium rather than away from it.
The second law captures this directionality.

For an **isolated system** (no exchange of energy or matter with surroundings):

$$
dS \geq 0
$$

Entropy never decreases.
Equality holds for a reversible process; strict inequality holds for a spontaneous, irreversible one.

```{admonition} Connection to statistical mechanics
:class: note

For readers of Part 1: the second law $dS \geq 0$ and the statistical result $S = k_B \ln \Omega$ are the same principle.
In @entropy, the isolated system evolves toward the macrostate with the most microstates (maximising $\mathcal{H} = \ln \Omega$).
Here, $S = k_B \mathcal{H}$ converts the dimensionless statistical count to thermodynamic units (J/K), and $dS \geq 0$ expresses the same idea: nature moves toward the most probable configuration.
```

The statistical picture from @entropy makes this concrete: $S = k_B \ln \Omega$ counts the number of microstates consistent with a given macrostate.
All microstates are equally probable, but macrostates with more microstates are overwhelmingly more likely, exponentially so.
A system evolves spontaneously toward the macrostate with the highest entropy not because entropy is a mysterious driving force, but because that outcome is the most probable one by an enormous margin.
Irreversibility is statistics, not magic.

For a system in thermal contact with its surroundings, the **Clausius inequality** relates entropy change to heat:

$$
dS \geq \frac{\delta Q}{T}
$$

Equality holds for reversible processes: $\delta Q_\text{rev} = T\,dS$.
For irreversible processes, entropy is also generated internally, so $dS > \delta Q / T$.

## The third law

The first and second laws establish how energy and entropy change; neither fixes the absolute value of entropy.
The **third law of thermodynamics** provides that reference point:

$$
S \to 0 \quad \text{as} \quad T \to 0
$$

For a perfect crystal at absolute zero, there is only one microstate ($W = 1$), so $S = k_B \ln 1 = 0$.
This grounds entropy on an absolute scale: the standard molar entropies $S^\circ$ in thermodynamic tables are absolute values, measured by integrating $C_p / T$ from 0 K to the temperature of interest.
Without the third law, only entropy differences could be tabulated; with it, absolute entropies can be used directly to compute $\Delta S$ for any process.

---

The next chapter builds the mathematical framework that follows from these laws: the fundamental relation, the Gibbs-Duhem equation, and the thermodynamic potentials.
