---
kernelspec:
  name: python3
  display_name: Python 3 (ipykernel)
---

(kinetics)=
# Reaction Mechanisms and Rate Laws

@equilibrium established *when* a system stops changing: the composition where $\Delta_r G = 0$.
It gives no information about *how quickly* that state is reached.
Kinetics addresses the complementary question: at what rate does a reaction proceed, and how does that rate depend on concentration and mechanism?

## Reaction rate

The rate of a reaction is the time derivative of the extent of reaction per unit volume:

$$
r = \frac{1}{V}\frac{d\xi}{dt}
$$

giving $r$ units of mol/(L$\cdot$s).
From $dn_i = \nu_i\,d\xi$, the rate connects to the concentration change of any species:

$$
r = \frac{1}{\nu_i}\frac{dc_i}{dt}
$$

**Rate laws.** The rate law relates $r$ to the current composition.
For many reactions it takes the power-law form:

$$
r = k \prod_{\text{reactants}} c_i^{n_i}
$$

where $n_i$ are the **reaction orders** and $k$ is the **rate constant**.
The rate constant $k$ is independent of composition but depends strongly on temperature (see @kinetics-temperature).

## Elementary reactions

A reaction is **elementary** if it occurs in a single molecular event with no intermediate species.
For an elementary step, the reaction orders equal the stoichiometric coefficients by definition: the rate is proportional to the probability of the required molecules meeting simultaneously.

**Example: SN2 substitution.**
The nucleophilic substitution $\text{CH}_3\text{Br} + \text{OH}^- \to \text{CH}_3\text{OH} + \text{Br}^-$ is a single bimolecular collision.
The rate law is exactly

$$
r = k\,[\text{CH}_3\text{Br}][\ce{OH-}]
$$

Order equals stoichiometry because there is one elementary step and two molecules must meet for it to occur.
Termolecular elementary steps (three simultaneous collisions) are extremely rare; practically all elementary steps are uni- or bimolecular.

## Complex mechanisms: order from experiment

Most reactions written as balanced equations are **not** elementary; they proceed through a sequence of elementary steps called the **mechanism**.
The observed rate law reflects the mechanism, not the stoichiometry, and must be determined experimentally.

```{important}
For a non-elementary reaction, the exponents $n_i$ in $r = k\prod c_i^{n_i}$ bear
no required relationship to the stoichiometric coefficients $\nu_i$.
Stoichiometric coefficients are fixed by conservation of atoms; reaction orders are
properties of the mechanism and must be measured.
Only for elementary steps does $n_i = |\nu_i|$ hold by definition.
```

**Example: SN1 substitution.**
The hydrolysis of tert-butyl bromide appears bimolecular from the balanced equation:

$$
(\text{CH}_3)_3\text{CBr} + \text{OH}^- \to (\text{CH}_3)_3\text{COH} + \text{Br}^-
$$

The mechanism, however, has two steps:

$$
(\text{CH}_3)_3\text{CBr} \xrightarrow{k_1,\;\text{slow}} (\text{CH}_3)_3\text{C}^+ + \text{Br}^-
$$
$$
(\text{CH}_3)_3\text{C}^+ + \text{OH}^- \xrightarrow{k_2,\;\text{fast}} (\text{CH}_3)_3\text{COH}
$$

The first step (ionization) is rate-determining; the carbocation is consumed almost instantly.
The observed rate law is:

$$
r = k_1\,[\,(\text{CH}_3)_3\text{CBr}\,]
$$

First order, despite the bimolecular stoichiometry.
OH$^-$ does not appear in the rate law because it acts only after the slow step.

**Example: fractional order from a chain mechanism.**
The overall reaction $\text{H}_2 + \text{Br}_2 \to 2\,\text{HBr}$ suggests a bimolecular rate law.
The experimental rate law is instead:

$$
r = \frac{k[\ce{H2}][\ce{Br2}]^{1/2}}{1 + k'[\ce{HBr}]/[\ce{Br2}]}
$$

The half-order in Br$_2$ and the inhibition by product HBr both arise from the underlying radical chain mechanism (initiation, propagation, termination steps).
No balanced equation can reproduce this structure; only the mechanism does.

## Rate-determining step

When one elementary step in a mechanism is much slower than all others, it controls the overall rate: the **rate-determining step** (RDS).
The rate law of the overall reaction then equals the rate law of the RDS, which may involve reactive intermediates that must be eliminated.

**Elimination of intermediates.** Steps before the RDS are fast and can be treated as equilibrated.
If an intermediate X appears in the RDS rate expression, its concentration is expressed in terms of observable species using the equilibrium condition of the fast preceding steps.
This approach is used to derive the Michaelis-Menten rate law in @saturation.

## Pseudo-first-order

For $\text{A} + \text{B} \to \text{products}$ with $r = k[\text{A}][\text{B}]$, if species B is in large excess, its concentration barely changes during the reaction.
In that regime:

$$
[\text{B}] \approx [\text{B}]_0 \quad\Rightarrow\quad
r \approx k_\text{obs}[\text{A}], \qquad k_\text{obs} = k[\text{B}]_0
$$

The reaction appears first-order in A even though the true rate law is second-order overall.
Pseudo-first-order conditions are widely exploited experimentally to isolate the concentration dependence of a single species; $k$ is then recovered from $k_\text{obs}$ and the known $[\text{B}]_0$.

---

The next chapter derives the mass action law, which connects rate constants directly to the equilibrium constant and imposes thermodynamic consistency on any kinetic model.
