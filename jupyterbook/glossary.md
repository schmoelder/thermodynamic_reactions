---
kernelspec:
  name: python3
  display_name: Python 3 (ipykernel)
---

(glossary-page)=
# Glossary

```{glossary}

Activation energy
: The minimum energy $E_a$ that colliding particles must possess for a reaction to occur.
  The fraction of particles exceeding $E_a$ is proportional to $e^{-E_a/k_BT}$, giving the Arrhenius dependence of reaction rates on temperature.
  Introduced in @maxwell-boltzmann.

Acid dissociation constant
: $K_a = a_{\text{H}^+} a_{\text{A}^-} / a_\text{HA}$: the equilibrium constant for the proton-transfer reaction $\text{HA} \rightleftharpoons \text{H}^+ + \text{A}^-$.
  It is a special case of the general equilibrium constant $K = \prod_i a_i^{\nu_i}$ (see @equilibrium).
  Introduced in @acid-base.

Adsorption
: The distribution of solute between a mobile phase and a stationary phase, treated as a quasi-equilibrium process in chromatography.
  Characterised by an adsorption equilibrium constant $K_\text{ads}$ and described by isotherms such as Langmuir.
  Introduced in @adsorption.

Adsorption equilibrium constant
: $K_\text{ads} = \exp(-\Delta G^\circ_\text{ads}/RT)$: the equilibrium constant for the adsorption reaction $\text{A}_\text{mobile} \rightleftharpoons \text{A}_\text{adsorbed}$.
  Determines the partitioning of a species between mobile and stationary phases in chromatography.
  Introduced in @adsorption.

Activity
: The thermodynamic quantity $a_i = \gamma_i c_i / c^\circ$ that replaces the dimensionless concentration in the chemical potential of a real solution: $\mu_i = \mu_i^\circ + RT\ln a_i$.
  In the dilute limit $\gamma_i \to 1$ and $a_i \to c_i/c^\circ$.
  Introduced in @nonidealities.

Activity coefficient
: The factor $\gamma_i$ that corrects the ideal concentration term for intermolecular interactions: $a_i = \gamma_i c_i/c^\circ$.
  Activity coefficients quantify deviations from ideal solution behaviour and are constrained collectively by the Gibbs-Duhem equation.
  Introduced in @nonidealities.

Apparent pKa
: $\text{p}K_a^\text{app}$: the negative logarithm of the effective acid dissociation constant that includes activity coefficient corrections.
  For a species with charge $z$, the Davies-based shift is $\text{p}K_a^\text{app} = \text{p}K_a - 2A(1-z)(\sqrt{I}/(1+\sqrt{I}) - 0.3I)$.
  Introduced in @speciation-buffers.

Arrhenius equation
: The empirical rate law $k = A\,e^{-E_a/k_BT}$ relating a reaction rate constant to temperature through the activation energy $E_a$.
  The exponential factor is the Boltzmann-weighted fraction of particles energetic enough to react.
  Introduced in @kinetics-temperature.

Avogadro's number
: $N_A = 6.022 \times 10^{23}\ \text{mol}^{-1}$.
  The number of particles per mole; converts between particle-level quantities (using $k_B$) and molar quantities (using $R = N_A k_B$).

Bjerrum diagram
: A plot of species speciation fractions $f_i$ as a function of pH for an acid-base system.
  Each species traces a sigmoid curve whose inflection point coincides with the corresponding pKa.
  For a polyprotic acid the diagram shows one sigmoid per dissociation step, spanning the full pH axis.
  Introduced in @acid-base.

Boltzmann constant
: $k_B = 1.380 \times 10^{-23}\ \text{J/K}$.
  The conversion factor between the statistical temperature $\partial H/\partial U$ and the Kelvin scale.
  Introduced in @entropy.

Buffer capacity
: $\beta = dc_b/d(\text{pH})$: the concentration of strong base required to raise the pH by one unit.
  For a monoprotic buffer at concentration $c_\text{tot}$, the maximum capacity $\beta_\text{max} = (\ln 10 / 4)\,c_\text{tot}$ occurs at $\text{pH} = \text{p}K_a$.
  Introduced in @acid-base.

Boltzmann factor
: The factor $e^{-\varepsilon/k_BT}$ giving the relative probability of a state with energy $\varepsilon$ in a system at thermal equilibrium.
  States with higher energy are exponentially less probable, with the energy scale set by $k_BT$.
  Derived in @maxwell-boltzmann.

Chemical equilibrium
: The state in which each chemical reaction satisfies $\sum_i \nu_i \mu_i = 0$.
  This ensures no net driving force for any reaction in the system.
  Introduced in @thermodynamic-potentials.

Chemical potential
: $\mu_i = \partial G/\partial n_i\big|_{T,P,n_{j\neq i}}$: the free energy cost of adding one mole of species $i$ at constant temperature, pressure, and all other amounts.
  Species flow from high to low $\mu_i$; equilibrium requires uniformity of chemical potentials across space and reactions.
  Introduced in @chemical-potential.

Clausius inequality
: $dS \geq \delta Q / T$: for any process, the entropy change of a system is at least as large as the heat absorbed divided by temperature.
  Equality holds for reversible processes; strict inequality for irreversible ones.
  Introduced in @laws-of-thermodynamics.

Clausius-Clapeyron equation
: $dP/dT = \Delta H_\text{trans}/(T \cdot \Delta V)$: the slope of the coexistence curve between two phases.
  Derived from the phase equilibrium condition $\mu_i(\alpha) = \mu_i(\beta)$; predicts how boiling point or melting point changes with pressure.
  Introduced in @chemical-potential.

Characteristic speeds
: The typical velocities that characterise a Maxwell-Boltzmann distribution: most probable speed $v_p$, mean speed $\langle v \rangle$, and root-mean-square speed $v_\text{rms}$.
  For an ideal gas: $v_p = \sqrt{2k_BT/m}$, $\langle v \rangle = \sqrt{8k_BT/\pi m}$, $v_\text{rms} = \sqrt{3k_BT/m}$. 
  Introduced in @maxwell-boltzmann.

Conjugate pairs
: The intensive-extensive variable pairs $(T, S)$, $(P, V)$, and $(\mu_i, n_i)$ that appear in the fundamental relation.
  Each product has units of energy.
  Introduced in @thermodynamic-potentials.

Davies equation
: An empirical extension of Debye-Hückel for ionic strengths up to $I \approx 0.5\ \text{mol/L}$:
  $\log_{10}\gamma_i = -Az_i^2\!\left(\sqrt{I}/(1+\sqrt{I}) - 0.3\,I\right)$, where $I$ is in mol/L.
  Introduced in @nonidealities.

Dalton's law
: In a mixture of ideal gases, the partial pressure of each component equals its mole fraction times the total pressure: $P_i = x_i P$.
  A direct consequence of the ideal gas law applied to each species independently.
  Introduced in @ideal-gas.

Debye-Hückel theory
: A theory for electrolyte solutions that accounts for long-range electrostatic interactions via the ionic atmosphere concept.
  In the limiting law, $\log_{10}\gamma_i = -Az_i^2\sqrt{I}$, exact as $I \to 0$.
  Introduced in @nonidealities.

Debye length
: $\kappa^{-1}$: the characteristic decay length of the electrostatic potential around an ion in solution, set by the ionic strength and the permittivity of the solvent.
  Shorter Debye lengths indicate stronger screening of ionic interactions.
  Introduced in @nonidealities.

Entropy (statistical)
: $H = \ln W$: the logarithm of the number of microstates compatible with a given macrostate.
  An additive, dimensionless measure of how spread a macrostate is across microstates.
  Introduced in @entropy.

Entropy (thermodynamic)
: $S = k_B \ln W$: statistical entropy scaled by the Boltzmann constant to give conventional units of J/K.
  Satisfies $\partial S/\partial U\big|_{N,V} = 1/T$ and $dS \geq \delta Q/T$.
  Introduced in @entropy.

Equation of state
: A relation between the macroscopic state variables $P$, $V$, $T$, and $n$ of a substance.
  The ideal gas law $PV = nRT$ is the simplest; more realistic equations (van der Waals, etc.) correct for intermolecular interactions.
  Introduced in @nonidealities.

Elastic collision
: A collision between particles in which kinetic energy is conserved (no energy loss to internal degrees of freedom or deformation).
  In a gas, elastic collisions between particles drive the redistribution of energy toward the Maxwell-Boltzmann distribution.
  Introduced in @collisions.

Elementary reaction
: A single molecular event in which reactants collide and form products in one step.
  For elementary reactions, the reaction order equals the stoichiometric coefficient exactly.
  Introduced in @kinetics.

Complex reaction
: A reaction whose balanced chemical equation hides a multi-step mechanism.
  The reaction order is determined experimentally and need not match the stoichiometric coefficients.
  Introduced in @kinetics.

Equipartition theorem
: Each quadratic degree of freedom contributes $\frac{1}{2}k_BT$ to the mean energy.
  For a monatomic ideal gas, three translational degrees give $\langle\varepsilon\rangle = \frac{3}{2}k_BT$.
  Derived in @maxwell-boltzmann.

Equilibrium constant
: $K = \prod_i a_i^{\nu_i}$ evaluated at equilibrium: the value of the reaction quotient when $\Delta_r G = 0$.
  Related to the standard reaction Gibbs energy by $\Delta_r G^\circ = -RT\ln K$.
  Introduced in @equilibrium.

Equilibrium mode
: A reaction modeling mode for fast reactions treated as instantaneous equilibria.
  Replaces an ODE with an algebraic constraint $\ln Q_j - \ln K_j = 0$, turning the transport PDE into a PDAE.
  Introduced in @implementation-source-term.

Extent of reaction
: The single coordinate $\xi$ (SI unit: mol) that tracks how far a reaction has proceeded.
  All composition changes follow $dn_i = \nu_i\,d\xi$; specifying $\xi$ and the initial amounts determines the full composition.
  Introduced in @equilibrium.

First law of thermodynamics
: Energy is conserved: $dU = \delta Q - \delta W + \sum_i \mu_i\,dn_i$.
  Any change in internal energy is accounted for by heat, mechanical work, and chemical work.
  Introduced in @laws-of-thermodynamics.

Fugacity
: The effective pressure $f$ of a real gas, defined so that $\mu = \mu^\circ + RT\ln(f/P^\circ)$ holds exactly.
  The fugacity coefficient $\varphi = f/P$ measures deviation from ideal behaviour; $\varphi \to 1$ as $P \to 0$.
  Introduced in @nonidealities.

Fugacity coefficient
: $\varphi = f/P$: the ratio of fugacity to pressure, measuring how much a real gas deviates from ideal behaviour.
  Introduced in @nonidealities.

Fundamental relation
: The combined first and second law of thermodynamics: $dU = T\,dS - P\,dV + \sum_i \mu_i\,dn_i$.
  The starting point for all thermodynamic potentials via Legendre transforms.
  Introduced in @thermodynamic-potentials.

Gas constant
: $R = N_A k_B = 8.314\ \text{J/(mol\,K)}$.
  The per-mole version of the Boltzmann constant; appears in the ideal gas law $PV = nRT$ and the ideal chemical potential.

Gibbs free energy
: $G = U - TS + PV$: the thermodynamic potential minimised at constant temperature and pressure.
  Satisfies $G = \sum_i \mu_i n_i$ as an identity; $dG = -S\,dT + V\,dP + \sum_i \mu_i\,dn_i$.
  Introduced in @thermodynamic-potentials.

Gibbs-Duhem equation
: $S\,dT - V\,dP + \sum_i n_i\,d\mu_i = 0$: the intensive variables $T$, $P$, and $\{\mu_i\}$ are not all independent.
  Derived from the fundamental relation and its integrated form; constrains activity coefficients.
  Introduced in @thermodynamic-potentials.

Henderson-Hasselbalch equation
: $\text{pH} = \text{p}K_a + \log_{10}([\text{A}^-]/[\text{HA}])$: the pH of a weak acid buffer at equilibrium in dilute ideal solution.
  Derived directly from the equilibrium condition $Q = K_a$; valid within roughly two pH units of $\text{p}K_a$ before activity corrections become significant.
  Introduced in @acid-base.

Henry's law
: For a dilute solute in equilibrium with its vapor, the partial pressure is proportional to mole fraction: $P_i = K_{H,i} x_i$.
  A asymmetric activity model for the dilute limit, contrasted with Raoult's law for the pure-component reference.
  Introduced in @nonidealities.

Heat capacity
: $C_p = (\partial H/\partial T)_{P,n}$: the heat required to raise the temperature of a system by one kelvin at constant pressure.
  Tabulated $C_p$ values allow extrapolation of $\Delta H$ and $\Delta S$ to other temperatures via the Kirchhoff relations.
  Introduced in @thermodynamic-potentials.

Helmholtz energy
: $A = U - TS$: the thermodynamic potential minimised at constant temperature and volume.
  Natural setting for molecular simulation; not used further in this book, where pressure (not volume) is the controlled variable.
  Introduced in @thermodynamic-potentials.

Ideal gas
: A gas whose molecules do not interact: no intermolecular forces, no excluded volume.
  Obeys $PV = nRT$ exactly.
  Introduced in @ideal-gas.

Internal energy
: $U$: the total microscopic energy of a system.
  For an ideal gas, purely kinetic; in general, includes potential energy from intermolecular interactions.
  The natural potential minimised at constant entropy and volume.
  Introduced in @particles.

Ionic strength
: $I = \frac{1}{2}\sum_i z_i^2\,c_i/c^\circ$: a measure of the total ion concentration weighted by the square of the charge number.
  Controls the strength of electrostatic interactions and the Debye length.
  Introduced in @nonidealities.

Kirchhoff relations
: Equations giving the temperature dependence of $\Delta_r H^\circ$ and $\Delta_r S^\circ$:
  $\Delta_r H^\circ(T_2) = \Delta_r H^\circ(T_1) + \int_{T_1}^{T_2} \Delta_r C_p\,dT$.
  Allow tabulated 298 K data to be extrapolated to other temperatures.
  Introduced in @equilibrium.

Kinetic energy
: The energy of a particle due to its motion, $\varepsilon = \frac{1}{2}mv^2$ for a single particle.
  The sum of kinetic energies of all particles comprises the internal energy $U$ of an ideal gas.
  Introduced in @particles.

Kinetic mode
: A reaction modeling mode in which reaction rates evolve on timescales comparable to transport.
  Each reaction contributes a flux $\varphi_j$; the source term is $f_{\text{react},i} = \sum_j \nu_{ij}\varphi_j$, integrated as an ODE.
  Introduced in @implementation-source-term.

Legendre transform
: A mathematical operation that replaces an extensive natural variable with its conjugate intensive one, producing a new thermodynamic potential.
  Subtracting $TS$ from $U$ swaps $(S, T)$; adding $PV$ swaps $(V, P)$.
  Introduced in @thermodynamic-potentials.

Le Chatelier's principle
: When a system at equilibrium is perturbed, it responds in a way that counteracts the perturbation.
  A consequence of the minimisation of Gibbs energy: changing $T$ shifts entropy, changing $P$ shifts volume, changing composition shifts species amounts, all to restore $dG = 0$. 
  Introduced in @thermodynamic-potentials.

Langmuir isotherm
: $q = q_\text{max}\, K_\text{ads}\, c / (1 + K_\text{ads}\, c)$: adsorption isotherm for a finite number of binding sites.
  Describes $q$ (adsorbed amount) as function of $c$ (solution concentration); linear at low loading, saturates at $q_\text{max}$.
  Introduced in @adsorption.

Macrostate
: A coarse-grained description of a system, such as the number of particles in each energy bin.
  Many microstates are compatible with the same macrostate; the number of compatible microstates is $W$.
  Introduced in @entropy.

Mass action law
: The kinetic rate expression $r = k_f\prod_\text{reactants} c_i^{|\nu_i|} - k_r\prod_\text{products} c_j^{\nu_j}$ that assigns reaction orders equal to stoichiometric coefficients.
  Exact for elementary reactions; consistent with thermodynamics when $k_f/k_r = K$ is enforced.
  Introduced in @mass-action-law.

Maxwell-Boltzmann distribution
: The equilibrium energy distribution of a classical ideal gas: $p(\varepsilon) \propto \varepsilon^{1/2}\,e^{-\varepsilon/k_BT}$.
  The unique distribution maximising entropy at fixed $N$ and $U$.
  Derived in @maxwell-boltzmann.

Mass concentration
: $\rho_i = c_i M_i$: the mass of species $i$ per unit volume.
  Engineering variant of concentration; related to molar concentration $c_i$ by the molar mass $M_i$.
  Introduced in @ideal-gas.

Mass density
: $\rho = \sum_i \rho_i$: the total mass per unit volume of a mixture.
  For an ideal gas, $\rho = \sum_i c_i M_i$.
  Introduced in @ideal-gas.

Maximum entropy principle
: In an isolated system with fixed $U$, $V$, and $N$, the equilibrium macrostate is the one that maximises the entropy $\mathcal{H} = \ln \Omega$.
  Equivalently, it maximises the Gibbs/Shannon form $\mathcal{H} = -\sum_i p_i \ln p_i$ subject to constraints.
  Introduced in @entropy.

Mechanical equilibrium
: The state in which pressure is uniform throughout the system, so no net volume changes occur.
  One of the three conditions for full thermodynamic equilibrium.
  Introduced in @thermodynamic-potentials.

Michaelis constant
: $K_m = (k_{-1} + k_2)/k_1$: the substrate concentration at which an enzyme-catalysed reaction proceeds at half its maximum rate.
  Sets the transition between first-order (substrate-limited) and zeroth-order (enzyme-saturated) kinetic regimes.
  Introduced in @saturation.

Michaelis-Menten kinetics
: The rate law $r = V_\text{max}[\text{S}]/(K_m + [\text{S}])$ for enzyme-catalysed reactions, derived from a quasi-steady-state approximation on the enzyme-substrate complex.
  Exhibits first-order behaviour at low substrate and zeroth-order saturation at high substrate, independent of the reaction stoichiometry.
  Introduced in @saturation.

Microstate
: The complete microscopic specification of a system: positions and velocities (or quantum numbers) of every particle.
  Introduced in @entropy.

Monod equation
: The empirical growth-rate expression $\mu_g = \mu_{g,\text{max}}[\text{S}]/(K_s + [\text{S}])$ for microbial cultures limited by a single nutrient.
  Mathematically identical to Michaelis-Menten kinetics; reflects the same saturation of catalytic capacity.
  Introduced in @saturation.

Mole fraction
: $x_i = n_i / \sum_j n_j$: the dimensionless fraction of species $i$ in a mixture.
  For ideal gas mixtures, the partial pressure is $P_i = x_i P$ (Dalton's law).
  Introduced in @ideal-gas.

Molar concentration
: $c_i = n_i / V$: the amount of species $i$ per unit volume.
  Primary concentration measure throughout this book; SI unit mol/m$^3$, often expressed as mol/L.
  Introduced in @ideal-gas.

Molar volume
: $V_m = V/n$: the volume per mole of substance.
  For an ideal gas, $V_m = RT/P$.
  Introduced in @ideal-gas.

Pseudo-first-order
: A kinetic regime in which a reaction that is formally higher-order appears first-order because one or more reactants are in large excess.
  For $\ce{A + B -> P}$ with $[B] \gg [A]$, the rate $r \approx k'[A]$ where $k' = k_f[B]$ is the pseudo-first-order rate constant.
  Introduced in @kinetics.

Phase
: A homogeneous region of matter that is uniform in chemical composition and physical state, separated from other regions by a phase boundary.
  Examples: liquid water, water vapour, ice, oil in an oil-water mixture.
  Introduced in @chemical-potential.

Phase equilibrium
: The state where two or more phases coexist at constant temperature and pressure.
  At phase equilibrium, the chemical potential of each species is equal in all phases: $\mu_i(\alpha) = \mu_i(\beta)$.
  Introduced in @chemical-potential.

pH
: $\text{pH} = -\log_{10} a_{\text{H}^+}$: a logarithmic measure of proton activity in solution.
  In dilute ideal solution, $\text{pH} \approx -\log_{10}([\text{H}^+]/c^\circ)$.
  Neutral pH at 25 °C is 7.0; it decreases to 6.5 near 60 °C because $K_w$ is temperature-dependent.
  Introduced in @acid-base.

pKa
: $\text{p}K_a = -\log_{10} K_a$: the negative logarithm of the acid dissociation constant.
  Provides a convenient scale for comparing acid strengths; the pH at which equal amounts of the acid and conjugate base are present.
  Introduced in @acid-base.

Rate-determining step
: The slowest step in a multi-step reaction mechanism that controls the overall rate.
  All other steps are assumed to be at quasi-equilibrium; only the RDS contributes to the rate law.
  Introduced in @kinetics.

Reaction rate
: $r = (1/V)(d\xi/dt)$: the rate of change of extent of reaction per unit volume.
  For a power-law rate expression, $r = k\prod_i c_i^{n_i}$ where $n_i$ is the reaction order with respect to species $i$.
  Introduced in @kinetics.

Rate constant
: The factor $k$ in a rate law $r = k\prod_i c_i^{n_i}$: independent of composition but strongly temperature-dependent.
  Its temperature dependence is captured by the Arrhenius equation $k = A\,e^{-E_a/RT}$.
  Introduced in @kinetics.

Reaction Gibbs energy
: $\Delta_r G = dG/d\xi = \sum_i \nu_i \mu_i$: the slope of the Gibbs energy with respect to the extent of reaction.
  At equilibrium $\Delta_r G = 0$; for $\Delta_r G < 0$ the reaction proceeds forward; for $\Delta_r G > 0$ it proceeds backward.
  Introduced in @reaction-gibbs-energy.

Reaction order
: The exponent $n_i$ in a power-law rate expression $r = k\prod_i c_i^{n_i}$.
  For elementary steps the order equals the stoichiometric coefficient; for complex mechanisms it is determined experimentally and need not be integer.
  Introduced in @kinetics.

Reaction quotient
: $Q = \prod_i a_i^{\nu_i}$: the current composition expressed as the product of activities raised to stoichiometric powers.
  At equilibrium $Q = K$; away from equilibrium $\Delta_r G = \Delta_r G^\circ + RT\ln Q$.
  Introduced in @reaction-gibbs-energy.

Reaction flux
: $\varphi_j$: the rate of reaction $j$ in the canonical form $\mathbf{S}\boldsymbol{\varphi}$.
  For kinetic mode, $\varphi_j = k_j\prod_i a_i^{\nu_{ij}}$; for equilibrium mode, fluxes are replaced by algebraic constraints.
  Introduced in @implementation-source-term.

Reaction source term
: $\mathbf{f}_\text{react}$: the contribution of chemical reactions to the species balance $\partial \mathbf{c}/\partial t = -\nabla\cdot\mathbf{J} + \mathbf{f}_\text{react}$ in CADET.
  In canonical form: $f_{\text{react},i} = \sum_j \nu_{ij}\varphi_j$.
  Introduced in @implementation-source-term.

Raoult's law
: For an ideal mixture, the partial pressure of each component equals its mole fraction times the vapor pressure of the pure component: $P_i = x_i P_i^*$.
  A symmetric activity model for the pure-component reference; $\gamma_i \to 1$ as $x_i \to 1$.
  Introduced in @nonidealities.

Second law of thermodynamics
: The entropy of an isolated system never decreases: $dS \geq 0$.
  Spontaneous processes increase entropy because macrostates with more microstates are overwhelmingly more probable.
  Introduced in @laws-of-thermodynamics.

Spontaneous process
: A process that proceeds without external driving force, characterised by $\Delta G < 0$ at constant $T$ and $P$.
  A negative $\Delta G$ establishes thermodynamic possibility; rate is a separate question governed by kinetics.
  Introduced in @chemical-potential.

Speciation fraction
: $f_i$: the fraction of a total element or species present in a particular form.
  In a monoprotic system: $f_\text{HA} = 1/(1 + 10^{\text{p}K_a - \text{pH}})$, $f_\text{A}^- = 1 - f_\text{HA}$. 
  Introduced in @speciation-buffers.

Standard chemical potential
: $\mu_i^\circ(T)$: the chemical potential of species $i$ at the standard state ($P^\circ$ or $c^\circ$).
  Absorbs all temperature-dependent contributions; tabulated in thermodynamic databases.
  Introduced in @chemical-potential.

Standard state
: The reference condition used to make the argument of the logarithm dimensionless in $\mu_i = \mu_i^\circ + RT\ln(\cdots)$.
  By convention: $P^\circ = 1\ \text{bar}$ for gases, $c^\circ = 1\ \text{mol/L}$ for solutions.
  Only differences $\Delta\mu$ are measurable, so the choice does not affect predictions.
  Introduced in @chemical-potential.

Stoichiometric coefficient
: The signed integer $\nu_i$ in the balanced reaction $\sum_i \nu_i \text{A}_i = 0$: positive for products, negative for reactants.
  Determines how the amount of each species changes with extent of reaction: $dn_i = \nu_i\,d\xi$.
  Introduced in @reaction-coordinates.

Stoichiometric matrix
: $\mathbf{S}$: the $n_\text{species} \times n_\text{reactions}$ matrix with entries $S_{ij} = \nu_{ij}$.
  In canonical form it maps reaction fluxes to source terms: $\mathbf{f}_\text{react} = \mathbf{S}\boldsymbol{\varphi}$.
  Introduced in @implementation-source-term.

Temperature
: The quantity that is equal in two systems when they are in thermal equilibrium.
  Defined statistically as $\partial S/\partial U\big|_{N,V} = 1/T$; the zeroth law guarantees consistency.
  Introduced in @entropy.

Thermal equilibrium
: The state in which temperature is uniform throughout and no net heat flows.
  One of the three conditions for full thermodynamic equilibrium.
  Introduced in @laws-of-thermodynamics.

Thermodynamic equilibrium
: The state in which all three equilibrium conditions hold simultaneously: thermal, mechanical, and chemical equilibrium.
  The relevant thermodynamic potential is minimised.
  Introduced in @thermodynamic-potentials.

Third law of thermodynamics
: The entropy of a perfect crystal approaches zero as temperature approaches absolute zero: $S \to 0$ as $T \to 0$.
  Grounds entropy on an absolute scale, enabling tabulation of standard entropies.
  Introduced in @laws-of-thermodynamics.

Transition state
: The unstable arrangement of atoms at the top of the activation energy barrier, denoted $[\mathrm{AB}]^\ddagger$.
  In Transition State Theory it is treated as a quasi-equilibrium species; the activation Gibbs energy $\Delta G^\ddagger$ determines the rate.
  Introduced in @kinetics-temperature.

Van der Waals equation
: $\left(P + an^2/V^2\right)\left(V - nb\right) = nRT$: corrects the ideal gas law for finite molecular size ($b$) and attractions ($a$).
  Introduced in @nonidealities.

Van't Hoff equation
: $d\ln K/dT = \Delta_r H^\circ/RT^2$: temperature dependence of the equilibrium constant.
  A van't Hoff plot of $\ln K$ vs $1/T$ has slope $-\Delta_r H^\circ/R$.
  Introduced in @equilibrium-temperature.

Water autoionization constant
: $K_w = [\text{H}^+][\text{OH}^-] = 1.01 \times 10^{-14}$ at 25 °C: equilibrium constant for $\text{H}_2\text{O} \rightleftharpoons \text{H}^+ + \text{OH}^-$.
  Determines the pH of pure water.
  Introduced in @acid-base.

Zeroth law of thermodynamics
: If system A is in thermal equilibrium with B, and B with C, then A and C are also in thermal equilibrium.
  Establishes temperature as a well-defined property.
  Introduced in @laws-of-thermodynamics.
```
