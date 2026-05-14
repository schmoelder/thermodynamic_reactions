---
kernelspec:
  name: python3
  display_name: Python 3 (ipykernel)
---

(saturation)=
# Enzyme Kinetics

Power-law rate laws assume that adding more reactant always increases the rate.
When a **catalyst** is involved, this breaks down: the catalyst can be fully occupied at high substrate concentrations, and the rate saturates.
The saturation curve has the same mathematical origin as the Langmuir isotherm (@adsorption): a finite number of sites is distributed between occupied and unoccupied states.
The difference is that Langmuir describes the equilibrium occupation fraction, while Michaelis-Menten describes the steady-state flux through those sites.

## Michaelis-Menten kinetics

**Mechanism** {cite:p}`michaelis1913,michaelis2011`**.** The substrate S binds reversibly to enzyme E, forming complex ES, which releases product P:

$$
\ce{E + S <=>[$k_1$][$k_{-1}$] ES ->[$k_2$] E + P}
$$

**Quasi-steady-state approximation.** When the enzyme concentration is much lower than the substrate concentration, the complex ES reaches a quasi-steady state rapidly and $[\ce{ES}]$ is approximately constant on the timescale of the overall reaction.
This is an application of the rate-determining step concept from @kinetics: the slow step is $k_2$ release, and the fast binding equilibrium gives $[\ce{ES}]$ in terms of observable quantities.

Setting $d[\ce{ES}]/dt \approx 0$ and using the conservation $[\ce{E}]_0 = [\ce{E}] + [\ce{ES}]$:

$$
[\ce{ES}] = \frac{[\ce{E}]_0 [\ce{S}]}{K_m + [\ce{S}]},
\qquad K_m = \frac{k_{-1} + k_2}{k_1}
$$

where $K_m$ is the **Michaelis constant**.
The rate of product formation is:

$$
r = k_2[\ce{ES}] = \frac{V_\text{max}[\ce{S}]}{K_m + [\ce{S}]},
\qquad V_\text{max} = k_2[\ce{E}]_0
$$

**Limiting cases** (@fig-michaelis-menten):

- $[\ce{S}] \ll K_m$: $r \approx (V_\text{max}/K_m)[\ce{S}]$, first order in S.
- $[\ce{S}] \gg K_m$: $r \approx V_\text{max}$, zeroth order; the enzyme is saturated.

The transition between regimes occurs near $[\ce{S}] = K_m$, where $r = V_\text{max}/2$.
These apparent orders arise from enzyme saturation, not stoichiometry; no balanced equation predicts them.

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-michaelis-menten

import numpy as np
import matplotlib.pyplot as plt
from reactions.plots import setup_figure, COLORS

Vmax = 1.0
Km = 1.5
S = np.linspace(0, 8, 400)

r = Vmax * S / (Km + S)
r_first = (Vmax / Km) * S
S_lo = S[S <= 2.0]

fig, ax = setup_figure()
ax.plot(S, r, color=COLORS["primary"], linewidth=2.5, label="Michaelis-Menten")
ax.plot(
    S_lo,
    r_first[S <= 2.0],
    color="C1",
    linewidth=1.5,
    linestyle="--",
    label=r"first-order limit: $r \approx (V_\mathrm{max}/K_m)[\mathrm{S}]$",
)
ax.axhline(Vmax, color="C2", linewidth=1.5, linestyle="--", label=r"$V_\mathrm{max}$")
ax.axvline(Km, color="gray", linewidth=1.0, linestyle=":")
ax.plot(Km, Vmax / 2, "o", color="gray", markersize=7, zorder=5)
ax.text(
    Km + 0.25,
    Vmax / 2,
    r"$[\mathrm{S}] = K_m,\ r = V_\mathrm{max}/2$",
    fontsize=9,
    va="center",
)
ax.text(Km, -0.08, r"$K_m$", ha="center", va="top", fontsize=10, color="gray")
ax.set_xlabel(r"$[\mathrm{S}]$ [a.u.]")
ax.set_ylabel(r"$r$ [a.u.]")
ax.set_ylim(-0.05, 1.35)
ax.set_xlim(0, 8)
ax.legend(fontsize=9, loc="lower right")
fig.tight_layout()
```

```{figure} #cell-michaelis-menten
:name: fig-michaelis-menten

Michaelis-Menten rate as a function of substrate concentration ($V_\text{max} = K_m = 1.5$).
At low $[\ce{S}]$ the rate is linear (first order); at high $[\ce{S}]$ it saturates at $V_\text{max}$ (zeroth order).
The half-maximum rate occurs at $[\ce{S}] = K_m$.
```

```{admonition} Intuition
:class: tip

Enzyme kinetics combine reversible binding with an irreversible conversion step.
At low substrate concentration most enzymes are free, and the reaction rate increases approximately linearly with $[\mathrm{S}]$.
As substrate increases, more enzymes form the intermediate complex $\mathrm{ES}$, and the availability of free enzyme becomes limiting.
At high substrate concentration nearly all enzymes are bound, and the rate approaches a maximum set by the turnover number $k_2$, which governs how fast the bound complex converts product.
The Michaelis constant $K_m$ characterizes the balance between binding and conversion; smaller values indicate that saturation is reached at lower substrate concentrations.
```

**Kinetic vs. equilibrium $K_m$.**
The Michaelis constant $K_m = (k_{-1} + k_2)/k_1$ is not in general an equilibrium dissociation constant.
Only when $k_2 \ll k_{-1}$ does $K_m \approx k_{-1}/k_1 = K_S$, the true binding equilibrium constant.
In that limit the mathematical form becomes identical to Langmuir: site occupancy at quasi-equilibrium, multiplied by the turnover rate $k_2$.


## Hill equation

Some enzymes show a sigmoidal rate-concentration curve steeper than Michaelis-Menten, arising when binding at one active site alters the affinity at another {cite:p}`hill1910`.
The Hill equation captures this phenomenologically:

$$
v_\text{Hill} = \frac{V_\text{max}\,[S]^n}{K_m^n + [S]^n},
$$

where $n$ is the Hill coefficient.
At $n = 1$ the equation reduces to Michaelis-Menten exactly.
For $n > 1$ the response is sigmoidal: the rate rises slowly at low $[S]$, steepens near $[S] = K_m$, then saturates at $V_\text{max}$; the transition is sharper and more switch-like than the $n = 1$ curve.
$K_m$ retains its meaning as the half-saturation concentration for all $n$.

The Hill coefficient is a phenomenological index, not a mechanistic parameter.
It quantifies the steepness of the sigmoidal transition and is fitted to rate-concentration data; its value reflects cooperativity without specifying its molecular origin.
Fractional values $n < 1$ describe negative cooperativity, where binding at one site reduces affinity at others.
The mechanistic basis, whether conformational coupling (Monod-Wyman-Changeux), sequential binding, or allosteric effects, determines the apparent $n$ but does not change the mathematical form.


## Monod equation

An empirically identical form describes the growth rate of microorganisms limited by a single nutrient S {cite:p}`monod1949`:

$$
\mu_g = \frac{\mu_{g,\text{max}}[\ce{S}]}{K_s + [\ce{S}]}
$$

where $\mu_g$ is the specific growth rate, $\mu_{g,\text{max}}$ the maximum growth rate, and $K_s$ the half-saturation constant.
The mathematical equivalence with Michaelis-Menten reflects the same physical idea: a limiting resource is processed by a finite pool of catalytic capacity (enzymes or cells), which saturates at high substrate concentrations.


(enzyme-ph-activity)=
## pH-dependent enzyme activity

Many enzymes require a specific protonation state at the active site to function {cite:p}`dixon1953`.
The standard model considers two ionizable residues: one must be deprotonated to bind substrate (controlling $\text{p}K_{a1}$) and one must remain protonated to stabilise the transition state (controlling $\text{p}K_{a2}$).
Three enzyme states result from sequential proton loss:

$$
\ce{EH2+ <=>[K_{a1}] EH^* <=>[K_{a2}] E^-}
$$

where $\ce{EH^*}$ is the sole active form.
Treating both proton-transfer equilibria as fast relative to catalysis (the same timescale-separation argument used for acid-base equilibria in chromatography, @acid-base), the fraction of enzyme in the active state is the quasi-equilibrium occupancy of $\ce{EH^*}$.
This is exactly the middle Bjerrum fraction of a diprotic system (@speciation-buffers):

$$
f_\text{active}(a_{\mathrm{H}^+})
= \frac{K_{a1}\,a_{\mathrm{H}^+}}{a_{\mathrm{H}^+}^2 + K_{a1}\,a_{\mathrm{H}^+} + K_{a1}\,K_{a2}},
$$

where $a_{\mathrm{H}^+} = \gamma_\pm\,c_{\mathrm{H}^+}/c^\circ$ is the proton activity.
The connection to speciation is not an analogy: $f_\text{active}$ is the speciation fraction of $\ce{EH^*}$, computed by the same denominator as any other Bjerrum diagram.

Setting $\partial f_\text{active}/\partial a_{\mathrm{H}^+} = 0$ gives the optimal proton activity $a_{\mathrm{H}^+}^* = \sqrt{K_{a1}K_{a2}}$, hence:

$$
\text{pH}_\text{opt} = \tfrac{1}{2}(\text{p}K_{a1} + \text{p}K_{a2}).
$$

The full rate law combines saturation kinetics with the pH modulation:

$$
r([\ce{S}],\, a_{\mathrm{H}^+})
= \frac{V_\text{max}\,f_\text{active}(a_{\mathrm{H}^+})\,[\ce{S}]}{K_m + [\ce{S}]}.
$$

$f_\text{active}$ enters as a multiplicative prefactor on $V_\text{max}$: at $\text{pH}_\text{opt}$ the enzyme operates at full capacity; away from the optimum the effective maximum rate is reduced by the fraction of enzyme not in the active state (@fig-ph-bell).

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-ph-bell

import numpy as np
import matplotlib.pyplot as plt
from reactions.plots import setup_figure

pKa1, pKa2 = 5.5, 8.5
Ka1 = 10.0 ** (-pKa1)
Ka2 = 10.0 ** (-pKa2)

pH = np.linspace(3, 12, 400)
a_H = 10.0 ** (-pH)
f = Ka1 * a_H / (a_H**2 + Ka1 * a_H + Ka1 * Ka2)

fig, ax = setup_figure()
ax.plot(pH, f, color="C2", linewidth=2.0)
ax.axvline(
    0.5 * (pKa1 + pKa2),
    color="gray",
    linewidth=0.9,
    linestyle="--",
    label=rf"$\mathrm{{pH}}_\mathrm{{opt}} = {0.5*(pKa1+pKa2):.1f}$",
)
ax.axvline(pKa1, color="gray", linewidth=0.6, linestyle=":")
ax.axvline(pKa2, color="gray", linewidth=0.6, linestyle=":")
ax.text(pKa1 - 0.15, 0.06, rf"p$K_{{a1}}={pKa1}$", fontsize=8, ha="right", color="gray")
ax.text(pKa2 + 0.15, 0.06, rf"p$K_{{a2}}={pKa2}$", fontsize=8, ha="left", color="gray")
ax.set_xlabel("pH")
ax.set_ylabel(r"$f_\mathrm{active}$")
ax.set_xlim(3, 12)
ax.set_ylim(-0.02, 1.08)
ax.legend(fontsize=9)
fig.tight_layout()
```

```{figure} #cell-ph-bell
:name: fig-ph-bell

Bell-shaped pH--activity profile ($\text{p}K_{a1} = 5.5$, $\text{p}K_{a2} = 8.5$) from the three-state model.
The optimum lies at $\text{pH}_\text{opt} = \tfrac{1}{2}(\text{p}K_{a1} + \text{p}K_{a2}) = 7.0$; the dotted lines mark the two $\text{p}K_a$ values where $f_\text{active}$ drops to half its maximum.
```

---

Part 4 turns from single reactions to the full reaction framework in CADET: how thermodynamic consistency is enforced numerically, how activity corrections are incorporated, and how the equilibrium constraint couples to transport in a chromatographic column (@implementation-source-term).
The combined model, coupling proton-transfer equilibria to pH-dependent saturating kinetics, appears in @implementation-enzyme.
