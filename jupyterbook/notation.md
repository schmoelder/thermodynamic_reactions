---
kernelspec:
  name: python3
  display_name: Python 3 (ipykernel)
---

(notation-page)=
# Notation and Key Equations

## Symbols

### Latin

| Symbol                  | Meaning                                                                                                                                         | Unit                         | First use                      |
| ----------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------- | ------------------------------ |
| $A$                     | Helmholtz energy ($= U - TS$)                                                                                                                   | J                            | @thermodynamic-potentials      |
| $A$                     | Pre-exponential factor                                                                                                                          | same as $k$                  | @kinetics-temperature          |
| $A^f$, $A^r$            | Forward and reverse pre-exponential factors                                                                                                     | same as $k$                  | @kinetics-temperature          |
| $a$, $b$                | Van der Waals constants                                                                                                                         | various                      | @nonidealities                 |
| $a_i$                   | Activity of species $i$                                                                                                                         |                              | @mixing                        |
| $c^\circ$               | Standard concentration $= 1$ mol/L                                                                                                              | mol/L                        | @chemical-potential            |
| $c_i$                   | Molar concentration of species $i$ $= n_i/V$                                                                                                    | mol/L                        | @ideal-gas                     |
| $c_s$                   | Counter-ion concentration in solution (SMA)                                                                                                     | mol/L                        | @multicomponent                |
| $C_p$                   | Heat capacity at constant pressure                                                                                                              | J/(mol K)                    | @thermodynamic-potentials      |
| $C_{p,k}$               | Molar heat capacity of solvent species $k$                                                                                                      | J/(mol K)                    | @implementation-energy-balance |
| $e$                     | Elementary charge $= 1.602\times10^{-19}$ C                                                                                                     | C                            | @mixing                        |
| $e_{ij}^f$, $e_{ij}^r$  | Kinetic exponents of species $i$ in reaction $j$, forward and reverse; $e_{ij}^f = \lvert\nu_{ij}\rvert$ for reactants for elementary reactions | —                            | @implementation-source-term    |
| $E_a$                   | Activation energy                                                                                                                               | J/mol                        | @maxwell-boltzmann             |
| $E_a^f$, $E_a^r$        | Forward and reverse activation energies; $E_a^f - E_a^r = \Delta_r H^\circ$                                                                     | J/mol                        | @kinetics-temperature          |
| $E_\text{tot}$          | Total enzyme concentration                                                                                                                      | mol/L                        | @implementation-enzyme         |
| $F$                     | Faraday constant $= 96\,485$ C/mol                                                                                                              | C/mol                        | @multicomponent                |
| $f$                     | Fugacity                                                                                                                                        | Pa                           | @nonidealities                 |
| $f_i$                   | Speciation fraction of species $i$ at given pH                                                                                                  |                              | @speciation-buffers            |
| $f_{\text{react},i}$    | Net production rate of species $i$ $= \sum_j \nu_{ij}\,\varphi_j$                                                                               | mol/(L$\cdot$s)              | @implementation-source-term    |
| $G$                     | Gibbs free energy ($= U - TS + PV$)                                                                                                             | J                            | @thermodynamic-potentials      |
| $G^E$                   | Excess Gibbs energy $= \Delta_\text{mix}G - \Delta_\text{mix}G^\text{id}$                                                                       | J                            | @mixing                        |
| $g^E$                   | Molar excess Gibbs energy $= G^E/n$; models typically specify $g^E/RT$ as a function of composition                                             | J/mol                        | @mixing                        |
| $\Delta_\text{mix}G$    | Gibbs energy of mixing                                                                                                                          | J                            | @mixing                        |
| $H^E$                   | Excess enthalpy $= G^E - T(\partial G^E/\partial T)_{P,n}$                                                                                      | J                            | @mixing                        |
| $S^E$                   | Excess entropy $= -(\partial G^E/\partial T)_{P,n}$                                                                                             | J/K                          | @mixing                        |
| $C_p^E$                 | Excess heat capacity $= -T(\partial^2 G^E/\partial T^2)_{P,n}$                                                                                  | J/K                          | @mixing                        |
| $V^E$                   | Excess molar volume $= (\partial G^E/\partial P)_{T,n}$                                                                                         | m$^3$                        | @mixing                        |
| $\Delta G^\ddagger$     | Activation Gibbs energy                                                                                                                         | J/mol                        | @kinetics-temperature          |
| $\Delta_r G$            | Reaction Gibbs energy $= dG/d\xi = \sum_i \nu_i \mu_i$                                                                                          | J/mol                        | @reaction-gibbs-energy         |
| $\Delta_r G^\circ$      | Standard reaction Gibbs energy                                                                                                                  | J/mol                        | @equilibrium                   |
| $H$                     | Enthalpy ($= U + PV$)                                                                                                                           | J                            | @thermodynamic-potentials      |
| $\Delta H_\text{trans}$ | Latent heat: enthalpy absorbed at constant pressure during a phase transition                                                                   | J/mol                        | @thermodynamic-potentials      |
| $\Delta_r H^\circ$      | Standard reaction enthalpy                                                                                                                      | J/mol                        | @equilibrium-temperature       |
| $\mathcal{H}$           | Statistical entropy ($= \ln \Omega = -N\sum_i p_i \ln p_i$)                                                                                     |                              | @entropy                       |
| $I$                     | Ionic strength $= \tfrac{1}{2}\sum_i c_i z_i^2$                                                                                                 | mol/m$^3$                    | @mixing                        |
| $\mathbf{J}$            | Transport flux vector (convection and dispersion) in the CADET PDE                                                                              | mol/(m$^2$ s)                | @implementation-source-term    |
| $k$                     | Rate constant                                                                                                                                   | mol$^{1-n}$L$^{n-1}$s$^{-1}$ | @kinetics                      |
| $k_B$                   | Boltzmann constant $= 1.380\times10^{-23}$ J/K                                                                                                  | J/K                          | @entropy                       |
| $k_\text{cat}$          | Catalytic rate constant (enzyme turnover number)                                                                                                | s$^{-1}$                     | @implementation-enzyme         |
| $k^f$, $k^r$            | Forward and reverse rate constants; with reaction index: $k_j^f$, $k_j^r$                                                                       | same as $k$                  | @mass-action-law               |
| $K$                     | Equilibrium constant $= Q\big                                                                                                                   | _\text{eq}$                  |                                | @equilibrium |
| $K_a$                   | Acid dissociation constant                                                                                                                      |                              | @acid-base                     |
| $K_{H,i}$               | Henry's law constant for species $i$; $P_i = K_{H,i}\,x_i$ in the dilute limit                                                                 | Pa                           | @mixing                        |
| $K_\text{ads}$          | Adsorption equilibrium constant                                                                                                                 |                              | @adsorption                    |
| $K_m$                   | Michaelis constant                                                                                                                              | mol/L                        | @saturation                    |
| $K_w$                   | Water autoionization constant $= [\ce{H+}][\ce{OH-}]$                                                                                           |                              | @acid-base                     |
| $m$                     | Particle mass                                                                                                                                   | kg                           | @particles                     |
| $M_i$                   | Molar mass of species $i$                                                                                                                       | kg/mol                       | @ideal-gas                     |
| $n$                     | Amount of substance                                                                                                                             | mol                          | @pressure                      |
| $N$                     | Number of particles                                                                                                                             |                              | @particles                     |
| $N_A$                   | Avogadro's number $= 6.022\times10^{23}$ mol$^{-1}$                                                                                             | mol$^{-1}$                   | @pressure                      |
| $p_i$                   | Occupation fraction of energy bin $i$ ($= n_i/N$)                                                                                               |                              | @entropy                       |
| $P$                     | Pressure                                                                                                                                        | Pa                           | @pressure                      |
| $P^\circ$               | Standard pressure $= 1$ bar                                                                                                                     | Pa                           | @chemical-potential            |
| $P_c$                   | Critical pressure                                                                                                                               | Pa                           | @nonidealities                 |
| $P_i$                   | Partial pressure of species $i$ $= x_i P$                                                                                                       | Pa                           | @ideal-gas                     |
| pH                      | Negative log proton activity $= -\log_{10} a_{\ce{H+}}$                                                                                         |                              | @acid-base                     |
| $\text{p}K_a$           | Negative log acid dissociation constant $= -\log_{10} K_a$                                                                                      |                              | @acid-base                     |
| $q$                     | Adsorbed amount per volume                                                                                                                      | mol/m$^3$                    | @adsorption                    |
| $\bar{q}_0$             | Free binding sites in SMA resin $= \Lambda - \sum_j(\nu_j+\sigma_j)q_j$                                                                         | mol/m$^3$                    | @multicomponent                |
| $q_\text{max}$          | Maximum adsorption capacity                                                                                                                     | mol/m$^3$                    | @adsorption                    |
| $q_s$                   | Counter-ion loading in resin phase (SMA)                                                                                                        | mol/m$^3$                    | @multicomponent                |
| $\delta Q$              | Infinitesimal heat                                                                                                                              | J                            | @laws-of-thermodynamics        |
| $\dot{Q}_\text{ext}$    | External heat input per unit volume                                                                                                             | W/m$^3$                      | @implementation-energy-balance |
| $P^f_j$, $P^r_j$        | Forward and reverse activity products for reaction $j$: $\prod_i a_i^{e_{ij}^f}$, $\prod_i a_i^{e_{ij}^r}$                                     | —                            | @implementation-energy-balance |
| $Q$                     | Reaction quotient $= \prod_i a_i^{\nu_i}$                                                                                                       |                              | @reaction-gibbs-energy         |
| $r$                     | Reaction rate                                                                                                                                   | mol/(L$\cdot$s)              | @kinetics                      |
| $R$                     | Gas constant $= 8.314$ J/(mol K)                                                                                                                | J/(mol K)                    | @pressure                      |
| $S$                     | Thermodynamic entropy ($= k_B \ln W$)                                                                                                           | J/K                          | @entropy                       |
| $\mathbf{S}$            | Stoichiometric matrix $(S_{ij} = \nu_{ij})$                                                                                                     |                              | @multicomponent                |
| $T$                     | Temperature                                                                                                                                     | K                            | @entropy                       |
| $T_c$                   | Critical temperature                                                                                                                            | K                            | @nonidealities                 |
| $U$                     | Internal energy                                                                                                                                 | J                            | @particles                     |
| $v$, $\mathbf{v}$       | Particle velocity / speed                                                                                                                       | m/s                          | @particles                     |
| $v_p$                   | Most probable speed                                                                                                                             | m/s                          | @maxwell-boltzmann             |
| $\langle v \rangle$     | Mean speed                                                                                                                                      | m/s                          | @maxwell-boltzmann             |
| $v_\text{rms}$          | Root-mean-square speed                                                                                                                          | m/s                          | @maxwell-boltzmann             |
| $V$                     | Volume                                                                                                                                          | m$^3$                        | @pressure                      |
| $\bar{V}_i$             | Partial molar volume $= (\partial V/\partial n_i)_{T,P,n_{j\neq i}}$                                                                            | m$^3$/mol                    | @chemical-potential            |
| $V_c$                   | Critical molar volume $= 3b$ (van der Waals)                                                                                                    | m$^3$/mol                    | @nonidealities                 |
| $V_m$                   | Molar volume $= V/n$                                                                                                                            | m$^3$/mol                    | @ideal-gas                     |
| $V_\text{max}$          | Maximum reaction rate (Michaelis-Menten)                                                                                                        | mol/(L$\cdot$s)              | @saturation                    |
| $\delta W$              | Infinitesimal work                                                                                                                              | J                            | @laws-of-thermodynamics        |
| $x_i$                   | Mole fraction of species $i$ $= n_i/\sum_j n_j$                                                                                                 |                              | @ideal-gas                     |
| $z_i$                   | Charge number of ion $i$                                                                                                                        |                              | @mixing                        |


### Greek


| Symbol               | Meaning                                                                                          | Unit            | First use                      |
| -------------------- | ------------------------------------------------------------------------------------------------ | --------------- | ------------------------------ |
| $\alpha_{ij}$        | Separation factor $= K_i/K_j$                                                                    |                 | @multicomponent                |
| $\beta$              | Inverse temperature $= 1/k_BT$                                                                   | J$^{-1}$        | @maxwell-boltzmann             |
| $\beta$              | Buffer capacity $= dc_b/d(\text{pH})$                                                            | mol/L           | @speciation-buffers            |
| $\gamma_i$           | Activity coefficient of species $i$                                                              |                 | @mixing                        |
| $\gamma_i^\infty$    | Activity coefficient at infinite dilution                                                        |                 | @mixing                        |
| $\varepsilon$        | Particle energy                                                                                  | J               | @particles                     |
| $\varepsilon_0$      | Permittivity of free space $= 8.854\times10^{-12}$ F/m                                          | F/m             | @mixing                        |
| $\varepsilon_r$      | Relative permittivity (dielectric constant) of the solvent                                       |                 | @mixing                        |
| $\kappa^{-1}$        | Debye length                                                                                     | m               | @mixing                        |
| $\Lambda$            | Ionic capacity of ion-exchange resin (total chargeable sites)                                    | mol/m$^3$       | @multicomponent                |
| $\mu_i$              | Chemical potential of species $i$                                                                | J/mol           | @laws-of-thermodynamics        |
| $\mu_i^*$            | Pure-component chemical potential (mole-fraction reference state)                                | J/mol           | @mixing                        |
| $\mu_i^\circ$        | Standard chemical potential                                                                      | J/mol           | @chemical-potential            |
| $\mu_g$              | Specific growth rate (Monod)                                                                     | s$^{-1}$        | @saturation                    |
| $\mu_{g,\text{max}}$ | Maximum specific growth rate                                                                     | s$^{-1}$        | @saturation                    |
| $\nu_i$              | Stoichiometric coefficient (positive: product, negative: reactant); characteristic charge in SMA |                 | @reaction-coordinates          |
| $\sigma_i$           | Steric shielding factor (SMA): sites blocked by adsorbed species $i$                             |                 | @multicomponent                |
| $\theta$             | Surface occupancy fraction $= q/q_\text{max}$                                                    |                 | @adsorption                    |
| $\xi$                | Extent of reaction                                                                               | mol             | @reaction-coordinates          |
| $\rho$               | Total mass density $= \sum_i \rho_i = \sum_i c_i M_i$                                            | kg/m$^3$        | @ideal-gas                     |
| $\rho_i$             | Mass concentration of species $i$ $= c_i M_i$                                                    | kg/m$^3$        | @ideal-gas                     |
| $\rho C_p$           | Volumetric heat capacity of the fluid $= \sum_k x_k (\rho_k/M_k) C_{p,k}$                        | J/(m$^3$ K)     | @implementation-energy-balance |
| $\dot{Q}_\text{ext}$ | External heat source rate per unit volume (wall transfer, jacket)                                | W/m$^3$         | @implementation-energy-balance |
| $\varphi$            | Fugacity coefficient $= f/P$                                                                     |                 | @nonidealities                 |
| $\varphi_j$          | Reaction flux of reaction $j$ (CADET)                                                            | mol/(L$\cdot$s) | @implementation-source-term    |
| $\psi_D$             | Donnan potential                                                                                 | V               | @multicomponent                |
| $\Omega$             | Number of microstates                                                                            |                 | @entropy                       |

**Note.** S
ome symbols are overloaded.
The context makes the meaning clear in each case:
- The symbol $\beta$ is used for the inverse temperature $\beta = 1/k_BT$ (@maxwell-boltzmann) and the buffer capacity $\beta = dc_b/d(\text{pH})$ (@speciation-buffers).
- The symbol $\varphi$ is used for the fugacity coefficient $\varphi = f/P$ (@nonidealities) and, with a reaction index, for the CADET reaction flux $\varphi_j$ (@implementation-source-term).
- The symbol $A$ is used for the Helmholtz energy $A = U - TS$ (@thermodynamic-potentials) and the Arrhenius pre-exponential factor (@kinetics-temperature).

---

## Key equations

### Statistical mechanics

**Kinetic energy of a particle**

$$\varepsilon = \tfrac{1}{2}mv^2$$

**Internal energy (ideal gas)**

$$U = \sum_{i=1}^N \varepsilon_i = \sum_{i=1}^N \tfrac{1}{2}mv_i^2$$

**Number of microstates**

$$\Omega = \frac{N!}{n_1!\,n_2!\,n_3!\cdots}$$

**Statistical entropy**

$$\mathcal{H} = \ln \Omega = -N\sum_i p_i \ln p_i, \qquad p_i = n_i/N$$

**Thermodynamic entropy and temperature**

$$S = k_B \ln W, \qquad \left.\frac{\partial S}{\partial U}\right|_{N,V} = \frac{1}{T}$$

**Maxwell-Boltzmann distribution**

$$p(\varepsilon) \propto \varepsilon^{1/2}\,e^{-\varepsilon/k_BT}$$

**Boltzmann factor**

$$n_j \propto e^{-\varepsilon_j/k_BT}$$

**Equipartition theorem** (monatomic ideal gas)

$$\langle\varepsilon\rangle = \tfrac{3}{2}k_BT$$

**Characteristic speeds**

$$v_p = \sqrt{\frac{2k_BT}{m}}, \qquad \langle v\rangle = \sqrt{\frac{8k_BT}{\pi m}}, \qquad v_\text{rms} = \sqrt{\frac{3k_BT}{m}}$$

**Arrhenius equation**

$$k = A\,e^{-E_a/k_BT}$$

---

### Ideal gas law

$$PV = nRT = Nk_BT, \qquad R = N_Ak_B$$

**Historical special cases**

| Law | Fixed | Result |
| --- | ----- | ------ |
| Boyle (1662) | $T$, $n$ | $PV = \text{const}$ |
| Charles (1787) | $P$, $n$ | $V/T = \text{const}$ |
| Gay-Lussac (1808) | $V$, $n$ | $P/T = \text{const}$ |
| Avogadro (1811) | $T$, $P$ | $V/n = \text{const}$ |

---

### Laws of thermodynamics

**First law**

$$dU = \delta Q - \delta W + \sum_i \mu_i\,dn_i$$

**Second law** (isolated system)

$$dS \geq 0$$

**Clausius inequality**

$$dS \geq \frac{\delta Q}{T}$$

**Third law**

$$S \to 0 \quad \text{as} \quad T \to 0$$

---

### Fundamental relation and potentials

**Fundamental relation**

$$dU = T\,dS - P\,dV + \sum_i \mu_i\,dn_i$$

**Integrated fundamental relation**

$$U = TS - PV + \sum_i \mu_i n_i$$

**Gibbs-Duhem equation**

$$S\,dT - V\,dP + \sum_i n_i\,d\mu_i = 0$$

**Thermodynamic potentials**

| Potential | Definition | Natural variables | Differential |
| --------- | ---------- | ----------------- | ------------ |
| $U$ | -- | $S$, $V$, $n_i$ | $T\,dS - P\,dV + \sum_i\mu_i\,dn_i$ |
| $H = U + PV$ | enthalpy | $S$, $P$, $n_i$ | $T\,dS + V\,dP + \sum_i\mu_i\,dn_i$ |
| $A = U - TS$ | Helmholtz | $T$, $V$, $n_i$ | $-S\,dT - P\,dV + \sum_i\mu_i\,dn_i$ |
| $G = U - TS + PV$ | Gibbs | $T$, $P$, $n_i$ | $-S\,dT + V\,dP + \sum_i\mu_i\,dn_i$ |

**Gibbs energy as sum of chemical potentials** (identity, holds generally)

$$G = \sum_i \mu_i n_i$$

**Heat capacity at constant pressure**

$$C_p = \left.\frac{\partial H}{\partial T}\right|_{P,n}$$

---

### Spontaneity and chemical potential

**Spontaneity criterion** (constant $T$, $P$)

$$\Delta G = \Delta H - T\Delta S < 0$$

**Chemical potential**

$$\mu_i = \left.\frac{\partial G}{\partial n_i}\right|_{T,P,n_{j\neq i}}$$

**Pressure dependence of $\mu$**

$$\left.\frac{\partial \mu}{\partial P}\right|_T = V_m$$

**Ideal chemical potential (gas)**

$$\mu_i = \mu_i^\circ(T) + RT\ln\frac{P}{P^\circ}$$

**Ideal chemical potential (solution)**

$$\mu_i = \mu_i^\circ(T) + RT\ln\frac{c_i}{c^\circ}$$

---

### Non-idealities

**Van der Waals equation of state**

$$\left(P + \frac{an^2}{V^2}\right)\left(V - nb\right) = nRT$$

**Chemical potential with fugacity**

$$\mu = \mu^\circ(T) + RT\ln\frac{f}{P^\circ}, \qquad \varphi = \frac{f}{P} \to 1 \text{ as } P \to 0$$

**Activity and activity coefficient**

$$\mu_i = \mu_i^\circ(T) + RT\ln a_i, \qquad a_i = \gamma_i\frac{c_i}{c^\circ}, \qquad \gamma_i \to 1 \text{ as } c_i \to 0$$

**Ionic strength**

$$I = \frac{1}{2}\sum_i z_i^2\,\frac{c_i}{c^\circ}$$

**Debye-Hückel limiting law**

$$\log_{10}\gamma_i = -A z_i^2\sqrt{I}$$

**Davies equation**

$$\log_{10}\gamma_i = -A z_i^2\!\left(\frac{\sqrt{I}}{1+\sqrt{I}} - 0.3\,I\right)$$

where $A \approx 0.509\ \mathrm{mol^{-1/2}\,L^{1/2}}$ in water at 25 °C and $I$ is ionic strength in mol/L.

---

### Chemical reactions and equilibrium

**Composition change with extent of reaction**

$$n_i(\xi) = n_{i,0} + \nu_i\,\xi$$

**Reaction Gibbs energy**

$$\Delta_r G = \frac{dG}{d\xi} = \sum_i \nu_i \mu_i = \Delta_r G^\circ + RT\ln Q$$

**Equilibrium condition**

$$\Delta_r G = 0 \quad\Longrightarrow\quad \Delta_r G^\circ = -RT\ln K, \qquad K = \prod_i a_i^{\nu_i}\bigg|_\text{eq}$$

**Standard reaction Gibbs energy from tabulated data**

$$\Delta_r G^\circ = \Delta_r H^\circ - T\Delta_r S^\circ, \qquad \Delta_r H^\circ = \sum_i \nu_i\,\Delta_f H^\circ_i, \qquad \Delta_r S^\circ = \sum_i \nu_i\,S^\circ_i$$

**Kirchhoff relations**

$$\Delta_r H^\circ(T_2) = \Delta_r H^\circ(T_1) + \int_{T_1}^{T_2}\!\Delta_r C_p\,dT$$

**Van't Hoff equation**

$$\frac{d\ln K}{dT} = \frac{\Delta_r H^\circ}{RT^2} \quad\Longrightarrow\quad \ln\frac{K(T_2)}{K(T_1)} = -\frac{\Delta_r H^\circ}{R}\!\left(\frac{1}{T_2} - \frac{1}{T_1}\right)$$

**Arrhenius equation**

$$k = A\,e^{-E_a/RT}$$

**Mass action law**

$$r = k^f \prod_\text{reactants} c_i^{|\nu_i|} - k^r \prod_\text{products} c_j^{\nu_j}, \qquad \frac{k^f}{k^r} = K$$

**Michaelis-Menten / Monod**

$$r = \frac{V_\text{max}[\ce{S}]}{K_m + [\ce{S}]}$$

---

### Acid-base equilibria

**Acid dissociation constant**

$$K_a = \frac{[\ce{H+}][\ce{A-}]}{[\ce{HA}]}, \qquad \text{p}K_a = -\log_{10} K_a$$

**Water autoionisation**

$$K_w = [\ce{H+}][\ce{OH-}] = 1.01\times10^{-14} \text{ at } 25^\circ\text{C}$$

**Henderson-Hasselbalch equation**

$$\text{pH} = \text{p}K_a + \log_{10}\!\frac{[\ce{A-}]}{[\ce{HA}]}$$

**Speciation fractions (monoprotic)**

$$f_{\ce{A-}} = \frac{1}{1 + 10^{\text{p}K_a - \text{pH}}}, \qquad f_{\ce{HA}} = 1 - f_{\ce{A-}}$$

**Buffer capacity**

$$\beta = \ln 10\left([\ce{H+}] + [\ce{OH-}] + c_\text{tot}\,f_{\ce{HA}}\,f_{\ce{A-}}\right), \qquad \beta_\text{max} = \frac{\ln 10}{4}\,c_\text{tot}$$

**Apparent pKa shift (Davies equation, acid charge $z$)**

$$\text{p}K_a^\text{app} = \text{p}K_a - 2A(1-z)\!\left(\frac{\sqrt{I}}{1+\sqrt{I}} - 0.3\,I\right)$$
