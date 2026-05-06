---
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

(multicomponent)=
# Multicomponent Equilibria

Single-component Langmuir (@adsorption) and single-reaction equilibrium (@equilibrium) treat one species in isolation.
Real chromatographic systems involve many species competing for the same finite-site pool, coupled through shared counter-ions and a common solvent.
This chapter extends both stoichiometric and finite-site equilibria to multicomponent systems, recovering the competitive Langmuir isotherm, the Steric Mass Action model, and the Donnan equilibrium as special cases of the same $\mu$-equality framework.


## Competitive adsorption

When $N$ species compete for the same adsorption sites, the finite-site constraint becomes $\sum_i \theta_i + \theta_\text{free} = 1$.
Applying $\mu$-equality jointly under the shared site-balance constraint yields the competitive Langmuir isotherm:

$$
q_i = \frac{q_{\max,i}\,K_i\,c_i}{1 + \sum_j K_j\,c_j}.
$$

The denominator is shared across all species: a high concentration of any species $j$ suppresses $q_i$ even when $c_i$ is unchanged, because all species draw from the same finite site pool.
This global coupling, not pairwise competition, is what produces nonlinear displacement cascades and sharpening concentration fronts in overloaded chromatography.
Species-dependent $q_{\max,i}$ values reflect differences in steric footprint or effective accessibility within a single-site model, not necessarily distinct site populations.

The **separation factor** $\alpha_{ij} = K_i/K_j$ determines which species displaces which; $\alpha_{ij} > 1$ means species $i$ has higher affinity and will displace species $j$ from the stationary phase.


## Steric Mass Action

Ion-exchange adsorption involves charge displacement: a protein ion of characteristic charge $\nu_i$ displaces $\nu_i$ counter-ions from the resin to bind, and sterically shields $\sigma_i$ additional sites from other proteins.
The SMA equilibrium condition is

$$
q_i = K_i\, c_i \left(\frac{\bar{q}_0}{c_s}\right)^{\!\nu_i},
$$

where $c_s$ is the counter-ion concentration in solution and $\bar{q}_0$ is the number of free binding sites,

$$
\bar{q}_0 = \Lambda - \sum_j \left(\nu_j + \sigma_j\right) q_j.
$$

Here $\Lambda$ is the ionic capacity of the resin (total chargeable sites per unit volume); each adsorbed protein of species $j$ removes $\nu_j + \sigma_j$ sites from the available pool.
Unlike the Langmuir case, where the site balance closes the system independently, SMA is a coupled adsorption-ion balance problem: the counter-ion concentration $c_s$ is not an independent parameter but a state variable tied to the protein loadings by electroneutrality of the resin phase,

$$
q_s = \Lambda - \sum_j \nu_j\, q_j,
$$

which replaces "free counter-ion availability" with a coupled constraint.
The ratio $\bar{q}_0 / c_s$ is dimensionless for the same reason $c/c^\circ$ is: $c_s$ plays the role of the standard concentration $c^\circ$, making $K_i$ a true dimensionless equilibrium constant.
The factor $(\bar{q}_0 / c_s)^{\nu_i}$ is the available ionic capacity raised to the stoichiometric power: the IEX analogue of the shared denominator in competitive Langmuir.
In the single-component limit with $\nu_i = 1$ and $\sigma_i = 0$, electroneutrality gives $q_s = \Lambda - q_i$, so $\bar{q}_0 = q_s = \Lambda - q_i$, and the SMA reduces to a Langmuir isotherm with $q_\text{max} = \Lambda$ and $K_\text{eff} = K_i/c_s$; the saturation denominator arises from counter-ion depletion rather than a steric term.


## Donnan equilibrium

SMA describes charge uptake by the solid phase: how much counter-ion is displaced and how much protein is bound.
Donnan equilibrium describes the complementary question: how do mobile ions partition between the bulk solution and the solution inside the charged pore space of the resin?
This is a phase equilibrium problem, not an adsorption problem.
An ion-exchange resin carries a fixed charge density that cannot leave the resin phase; this creates an electrostatic potential difference $\psi_D$ (the **Donnan potential**) at the resin-solution boundary, which shifts the equilibrium partition of all mobile ions between the two fluid phases.

The equilibrium condition is $\mu$-equality across the Donnan boundary, now including an electrostatic term:

$$
\mu_i^\text{bulk} = \mu_i^\text{resin}
\quad\Longrightarrow\quad
RT \ln a_i^\text{bulk} = RT \ln a_i^\text{resin} + z_i F \psi_D,
$$

where $a_i^\text{bulk}$ and $a_i^\text{resin}$ are the activities of species $i$ in the bulk solution and resin phases, $z_i$ is the charge number, and $F$ is the Faraday constant.
The Donnan potential is determined self-consistently by the electroneutrality condition in the resin phase.
Donnan equilibrium is the mechanistic basis of the colloidal isotherm model in CADET, which treats the resin as a charged Donnan phase in contact with the bulk solution.


## Multicomponent reaction networks

Multiple simultaneous reactions are described by a stoichiometric matrix $\mathbf{S}$ with one column per reaction $j$ and one row per species $i$, so $S_{ij} = \nu_{ij}$.
The equilibrium condition $\Delta_r G_j = 0$ holds independently for each reaction:

$$
\sum_i \nu_{ij}\,\mu_i = 0 \quad \forall\, j
\quad\Longleftrightarrow\quad
Q_j(\mathbf{a}) = K_j(T) \quad \forall\, j.
$$

Species shared across reactions, most commonly $\ce{H+}$ in buffer systems, couple the constraints into a nonlinear system.
Equilibrium is the composition that simultaneously satisfies all $\mu$-equality conditions and all conservation constraints: it is the minimum of $G$ on the feasible composition manifold, equivalently the intersection of the $Q_j = K_j$ hypersurfaces with the mass and charge conservation constraints.
`solve_equilibrium` in Part 4 finds this point by Newton iteration on the full coupled residual.

---

The next chapter adds temperature dependence to all equilibrium constants via the van't Hoff equation, which applies uniformly to $K$, $K_\text{ads}$, $K_a$, and the SMA binding constants $K_i$ (@equilibrium-temperature).
