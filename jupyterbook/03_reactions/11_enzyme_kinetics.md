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

**Mechanism.** The substrate S binds reversibly to enzyme E, forming complex ES, which releases product P:

$$
\text{E} + \text{S}
\underset{k_{-1}}{\stackrel{k_1}{\rightleftharpoons}}
\text{ES}
\xrightarrow{k_2} \text{E} + \text{P}
$$

**Quasi-steady-state approximation.** When the enzyme concentration is much lower than the substrate concentration, the complex ES reaches a quasi-steady state rapidly and $[\text{ES}]$ is approximately constant on the timescale of the overall reaction.
This is an application of the rate-determining step concept from @kinetics: the slow step is $k_2$ release, and the fast binding equilibrium gives $[\text{ES}]$ in terms of observable quantities.

Setting $d[\text{ES}]/dt \approx 0$ and using the conservation $[\text{E}]_0 = [\text{E}] + [\text{ES}]$:

$$
[\text{ES}] = \frac{[\text{E}]_0 [\text{S}]}{K_m + [\text{S}]},
\qquad K_m = \frac{k_{-1} + k_2}{k_1}
$$

where $K_m$ is the **Michaelis constant**.
The rate of product formation is:

$$
r = k_2[\text{ES}] = \frac{V_\text{max}[\text{S}]}{K_m + [\text{S}]},
\qquad V_\text{max} = k_2[\text{E}]_0
$$

**Limiting cases** (@fig-michaelis-menten):

- $[\text{S}] \ll K_m$: $r \approx (V_\text{max}/K_m)[\text{S}]$, first order in S.
- $[\text{S}] \gg K_m$: $r \approx V_\text{max}$, zeroth order; the enzyme is saturated.

The transition between regimes occurs near $[\text{S}] = K_m$, where $r = V_\text{max}/2$.
These apparent orders arise from enzyme saturation, not stoichiometry; no balanced equation predicts them.

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-michaelis-menten

import numpy as np
import matplotlib.pyplot as plt

Vmax = 1.0
Km   = 1.0
S    = np.linspace(0, 8, 400)

r       = Vmax * S / (Km + S)
r_first = (Vmax / Km) * S
S_lo    = S[S <= 1.5]

fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(S, r, color="#1c4f8a", linewidth=2.5, label="Michaelis-Menten")
ax.plot(
    S_lo,
    r_first[S <= 1.5],
    color="C1",
    linewidth=1.5,
    linestyle="--",
    label=r"first-order limit: $r \approx (V_\mathrm{max}/K_m)[\mathrm{S}]$",
)
ax.axhline(Vmax, color="C2", linewidth=1.5, linestyle="--", label=r"$V_\mathrm{max}$")
ax.axvline(Km, color="gray", linewidth=1.0, linestyle=":")
ax.plot(Km, Vmax / 2, "o", color="gray", markersize=7, zorder=5)
ax.text(Km + 0.1, Vmax / 2,
        r"$[\mathrm{S}] = K_m,\ r = V_\mathrm{max}/2$", fontsize=9, va="center")
ax.text(Km, -0.08, r"$K_m$", ha="center", va="top", fontsize=10, color="gray")
ax.set_xlabel(r"$[\mathrm{S}]$ [a.u.]", fontsize=11)
ax.set_ylabel(r"$r$ [a.u.]", fontsize=11)
ax.set_ylim(-0.05, 1.35)
ax.set_xlim(0, 8)
ax.legend(fontsize=9, loc="lower right")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
fig.tight_layout()
```

```{figure} #cell-michaelis-menten
:name: fig-michaelis-menten

Michaelis-Menten rate as a function of substrate concentration ($V_\text{max} = K_m = 1$).
At low $[\text{S}]$ the rate is linear (first order); at high $[\text{S}]$ it saturates at $V_\text{max}$ (zeroth order).
The half-maximum rate occurs at $[\text{S}] = K_m$.
```

```{admonition} Intuition
:class: tip

Enzyme kinetics combine reversible binding with an irreversible conversion step.
At low substrate concentration most enzymes are free, and the reaction rate increases approximately linearly with $[\mathrm{S}]$.
As substrate increases, more enzymes form the intermediate complex $\mathrm{ES}$, and the availability of free enzyme becomes limiting.
At high substrate concentration nearly all enzymes are bound, and the rate approaches a maximum set by the turnover number $k_2$, which governs how fast the bound complex converts product.
The Michaelis constant $K_m$ characterises the balance between binding and conversion; smaller values indicate that saturation is reached at lower substrate concentrations.
```

**Kinetic vs. equilibrium $K_m$.**
The Michaelis constant $K_m = (k_{-1} + k_2)/k_1$ is not in general an equilibrium dissociation constant.
Only when $k_2 \ll k_{-1}$ does $K_m \approx k_{-1}/k_1 = K_S$, the true binding equilibrium constant.
In that limit the mathematical form becomes identical to Langmuir: site occupancy at quasi-equilibrium, multiplied by the turnover rate $k_2$.

## Monod equation

An empirically identical form describes the growth rate of microorganisms limited by a single nutrient S:

$$
\mu_g = \frac{\mu_{g,\text{max}}[\text{S}]}{K_s + [\text{S}]}
$$

where $\mu_g$ is the specific growth rate, $\mu_{g,\text{max}}$ the maximum growth rate, and $K_s$ the half-saturation constant.
The mathematical equivalence with Michaelis-Menten reflects the same physical idea: a limiting resource is processed by a finite pool of catalytic capacity (enzymes or cells), which saturates at high substrate concentrations.

---

Part 4 turns from single reactions to the full reaction framework in CADET: how thermodynamic consistency is enforced numerically, how activity corrections are incorporated, and how the equilibrium constraint couples to transport in a chromatographic column (@implementation-source-term).
