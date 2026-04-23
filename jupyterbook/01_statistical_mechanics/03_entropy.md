---
kernelspec:
  name: python3
  display_name: Python 3 (ipykernel)
---

(entropy)=
# Entropy

The previous chapter showed that collisions drive the energy distribution toward a fixed shape.
Understanding why that shape and not another requires a way to measure how many microscopic arrangements are compatible with a given distribution.
That measure is entropy.

## Macrostates and microstates

A **microstate** is the complete specification of the system: the position and velocity of every particle.
A **macrostate** is a coarser description; for example, the number of particles $n_i$ with energy in each bin $[\varepsilon_i,\, \varepsilon_i + \Delta\varepsilon)$.

Many different microstates are compatible with the same macrostate.
The number of microstates corresponding to a given macrostate is denoted $\Omega$:

$$
\Omega = \frac{N!}{n_1!\, n_2!\, n_3! \cdots}
$$

This is the multinomial coefficient: the number of ways to assign $N$ distinguishable particles to bins with those occupations.

## Statistical entropy

$\Omega$ counts the number of microstates compatible with a given macrostate.
For two independent systems, such as two boxes of gas brought into contact, each microstate of the combined system is obtained by pairing a microstate of system 1 with a microstate of system 2.
As a result, the total number of microstates multiplies: $\Omega = \Omega_1 \Omega_2$.

A thermodynamic state function should be additive for independent systems (extensive).
$\Omega$ itself is not additive because it multiplies when systems are combined, and it grows extremely rapidly with system size (typically exponentially in $N$).
Taking the logarithm turns this multiplicative structure into an additive one:
$\ln(\Omega_1 \Omega_2) = \ln \Omega_1 + \ln \Omega_2$.
The resulting quantity scales linearly with system size and behaves like other extensive variables such as energy or volume.

The **statistical entropy** is therefore defined as

$$
\mathcal{H} = \ln \Omega
$$

$\mathcal{H}$ is a dimensionless number, a property of the distribution like its mean or variance.
Where variance measures how spread values are around their mean, $\mathcal{H}$ measures how spread a macrostate is across microstates.
A macrostate with only one compatible microstate has $\Omega = 1$, so $\mathcal{H} = 0$.
Every additional accessible microstate increases $\mathcal{H}$.
This is also Shannon's entropy from information theory: the same formula arises when quantifying how much is unknown about which microstate the system currently occupies.

```{code-cell} ipython3
:tags: [remove-cell]

import numpy as np
import matplotlib.pyplot as plt
from math import lgamma

def log_W(occupations):
    N = sum(occupations)
    return lgamma(N + 1) - sum(lgamma(n + 1) for n in occupations)
```

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-macrostates

macrostates = [
    [ 0, 20,  0,  0,  0],  # spike: Omega=1
    [15,  0,  0,  0,  5],  # bimodal extremes
    [10,  0, 10,  0,  0],  # bimodal centre
    [12,  2,  2,  2,  2],  # concentrated with equal tails
    [10,  4,  3,  2,  1],  # Boltzmann-like (exponential decrease)
]
# All have N=20 particles and total energy U = sum(n_i * eps_i) = 20

eps = np.arange(len(macrostates[0]))  # energy levels: 0, 1, 2, 3, 4

fig, axes = plt.subplots(1, len(macrostates), figsize=(13, 3), sharey=True)

for ax, occ in zip(axes, macrostates):
    lw = log_W(occ)
    ax.bar(eps, occ, color="C0", edgecolor="white")
    ax.set_xlabel(r"particle energy $\varepsilon_i$")
    ax.set_title(f"$\\mathcal{{H}} = {lw:.1f}$")
    ax.set_xticks(eps)

axes[0].set_ylabel("$n_i$")
plt.tight_layout()
plt.show()
```

```{figure} #cell-macrostates
:name: fig-macrostates

Five ways to distribute $N = 20$ particles across five energy levels, all with the same total energy $U = 20$.
The title of each panel shows the statistical entropy $\mathcal{H} = \ln \Omega$.
The Boltzmann-like panel (rightmost) has the highest $\mathcal{H}$ and by far the largest number of compatible microstates.
```

As seen in @fig-macrostates, $\mathcal{H}$ ranges from $0$ for the spike to over $20$ for the Boltzmann-like distribution, a factor of $e^{20} \approx 5 \times 10^8$ in $\Omega$.
For $N \sim 10^{23}$ particles the gap is so extreme that the highest-$\mathcal{H}$ macrostate is, for all practical purposes, the only one ever observed.

The occupation fractions $p_i = n_i / N$ visible in the bar charts are probabilities: $p_i$ is the fraction of particles found in energy bin $i$.
Applying Stirling's approximation ($\ln n! \approx n \ln n - n$) to $\mathcal{H} = \ln N! - \sum_i \ln n_i!$:

$$
\mathcal{H}
\approx N \ln N - \sum_i n_i \ln n_i
= -N \sum_i \frac{n_i}{N} \ln \frac{n_i}{N}
= -N \sum_i p_i \ln p_i
$$

Since $N$ is fixed, maximising $\mathcal{H}$ is equivalent to maximising $-\sum_i p_i \ln p_i$.
This is the **Gibbs entropy**, the same formula Shannon derived independently in information theory as a measure of uncertainty over a probability distribution.
Both forms are the same quantity; the probability form is more convenient when searching for the optimal distribution in the next chapter.

## The maximum-entropy principle

An isolated system evolves by moving between microstates.
If all microstates are equally likely (which follows from the symmetry of the laws of mechanics), the system spends most of its time in the macrostates with the most microstates, i.e., those with the highest $\mathcal{H}$.
This is the statistical foundation of the second law of thermodynamics:

> An isolated system at equilibrium is in the macrostate that maximises $\mathcal{H}$.

The initial state from the previous chapter, where every particle has the same energy, has $\mathcal{H} = 0$.
Any redistribution of energy among particles increases $\mathcal{H}$.
The system evolves not because anything pushes it, but because there are vastly more microstates in the spread-out macrostates than in the spike.

## Temperature and thermodynamic entropy

$\mathcal{H}$ is dimensionless and unit-free.
We can use it to define temperature directly from the statistics, before introducing any units.

Consider two systems that can exchange energy, with fixed total energy $U_1 + U_2$.
At equilibrium, total entropy is maximised.
Transferring a small amount $\delta U$ from system 2 to system 1 (so $\delta U_1 = -\delta U_2 = \delta U$) must leave total entropy unchanged to first order:

$$
\frac{\partial \mathcal{H}_1}{\partial U_1}\,\delta U - \frac{\partial \mathcal{H}_2}{\partial U_2}\,\delta U = 0
\quad\Longrightarrow\quad
\frac{\partial \mathcal{H}_1}{\partial U_1} = \frac{\partial \mathcal{H}_2}{\partial U_2}
$$

This shared quantity is temperature in disguise: a system with high $\partial \mathcal{H}/\partial U$ gains a lot of entropy per unit energy; it is cold and readily absorbs energy.
Low $\partial \mathcal{H}/\partial U$ means hot.
We define temperature as the inverse of this quantity:

$$
\frac{\partial \mathcal{H}}{\partial U}\bigg|_{N,V} \equiv \frac{1}{k_B T}
$$

where $k_B = 1.380 \times 10^{-23}$ J/K is a unit conversion factor, a consequence of measuring temperature in Kelvin rather than in Joules.
Multiplying $\mathcal{H}$ by $k_B$ gives thermodynamic entropy with conventional units:

$$
S = k_B \mathcal{H} = k_B \ln \Omega, \qquad \frac{\partial S}{\partial U} = \frac{1}{T}
$$

---

The next chapter uses the maximum entropy principle to identify the equilibrium distribution.
