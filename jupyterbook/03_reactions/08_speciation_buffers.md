---
kernelspec:
  name: python3
  display_name: Python 3 (ipykernel)
---

(speciation-buffers)=
# Speciation, Buffers, and Ionic Strength

The Henderson-Hasselbalch equation (@acid-base) gives pH for a single acid at a known ratio of forms.
This chapter asks the inverse question: given a total concentration and a pH, what fraction of each species is present?
That answer enables buffer design and reveals how ionic strength shifts the apparent equilibrium.

## Speciation: Bjerrum diagrams

The total analytical concentration $c_\text{tot} = [\ce{HA}] + [\ce{A-}]$ is conserved.
Solving for the individual fractions at a given pH:

$$
f_{\ce{A-}} = \frac{1}{1 + 10^{\text{p}K_a - \text{pH}}},
\qquad f_{\ce{HA}} = 1 - f_{\ce{A-}}
$$

A plot of these fractions versus pH is called a **Bjerrum diagram**.
Each species traces a sigmoid; the inflection of $f_{\ce{A-}}$ lies at $\text{pH} = \text{p}K_a$.

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-bjerrum-acetic

import numpy as np
import matplotlib.pyplot as plt
from reactions.plots import setup_figure, COLORS

pKa = 4.756  # acetic acid, NIST
pH = np.linspace(1, 9, 400)

f_AcO = 1.0 / (1.0 + 10.0 ** (pKa - pH))
f_AcOH = 1.0 - f_AcO

fig, ax = setup_figure()
ax.plot(pH, f_AcO, color="C0", linewidth=2.0, label=r"$f_{\mathrm{AcO}^-}$")
ax.plot(pH, f_AcOH, color="C1", linewidth=2.0, label=r"$f_{\mathrm{AcOH}}$")
ax.axvline(pKa, color="gray", linestyle=":", linewidth=1.0)
ax.text(pKa + 0.08, 0.52, rf"$\mathrm{{p}}K_a = {pKa}$", fontsize=9, color="gray")
ax.set_xlabel("pH")
ax.set_ylabel("Speciation fraction")
ax.set_xlim(1, 9)
ax.set_ylim(-0.02, 1.05)
ax.legend(fontsize=10)
fig.tight_layout()
```

```{figure} #cell-bjerrum-acetic
:name: fig-bjerrum-acetic

Bjerrum diagram for acetic acid ($\text{p}K_a = 4.756$, NIST).
The two forms cross at $\text{pH} = \text{p}K_a$ where each fraction equals one half.
```

**Polyprotic acids.** For an acid with $n$ dissociation steps with constants $K_{a,1},\ldots,K_{a,n}$, define the cumulative products $\alpha_0 = 1$ and $\alpha_k = \prod_{j=1}^{k} K_{a,j}$.
The speciation fraction of the species that has lost $k$ protons is:

$$
f_k = \frac{\alpha_k\,[\ce{H+}]^{n-k}}{\displaystyle\sum_{j=0}^{n}
\alpha_j\,[\ce{H+}]^{n-j}}
$$

Each transition produces one inflection point in the diagram at the corresponding pKa.

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-bjerrum-phosphate

import numpy as np
import matplotlib.pyplot as plt
from reactions.plots import setup_figure, COLORS
from reactions.api import speciation_fractions

pKa1, pKa2, pKa3 = 2.148, 7.198, 12.350

pH = np.linspace(0, 14, 500)
fracs = speciation_fractions(pH, [pKa1, pKa2, pKa3])

fig, ax = setup_figure()
labels = [
    r"$\mathrm{H_3PO_4}$",
    r"$\mathrm{H_2PO_4^-}$",
    r"$\mathrm{HPO_4^{2-}}$",
    r"$\mathrm{PO_4^{3-}}$",
]
for f, lbl, c in zip(fracs, labels, ["C0", "C1", "C2", "C3"]):
    ax.plot(pH, f, color=c, linewidth=2.0, label=lbl)
for pka, lbl in zip(
    [pKa1, pKa2, pKa3],
    [r"p$K_{a1}$", r"p$K_{a2}$", r"p$K_{a3}$"],
):
    ax.axvline(pka, color="gray", linestyle=":", linewidth=0.8)
    ax.text(pka + 0.15, 1.01, lbl, fontsize=8, color="gray", ha="left", va="bottom")
ax.set_xlabel("pH")
ax.set_ylabel("Speciation fraction")
ax.set_xlim(0, 14)
ax.set_ylim(-0.02, 1.10)
ax.legend(fontsize=9, loc="upper right")
fig.tight_layout()
```

```{figure} #cell-bjerrum-phosphate
:name: fig-bjerrum-phosphate

Bjerrum diagram for phosphoric acid ($\text{p}K_{a1} = 2.15$, $\text{p}K_{a2} = 7.20$, $\text{p}K_{a3} = 12.35$, NIST).
The three transitions span the full biochemical pH range.
Near pH 7, $\ce{H2PO4-}$ and $\ce{HPO4^2-}$ coexist, making phosphate an effective physiological buffer.
```

## Buffer capacity

A buffer resists pH change by consuming added acid or base.
The **buffer capacity** $\beta$ quantifies this resistance:

$$
\beta = \frac{dc_b}{d(\text{pH})}
$$

where $c_b$ is the concentration of strong base added to the solution.

**Derivation (monoprotic).** For a monoprotic buffer at total concentration $c_\text{tot}$, the proton balance at any stage of titration is:

$$
[\ce{H+}] + c_b = [\ce{OH-}] + [\ce{A-}]
$$

Differentiating with respect to pH:

$$
\beta = \ln 10 \left(
  [\ce{H+}] + [\ce{OH-}] + c_\text{tot}\,f_{\ce{HA}}\,f_{\ce{A-}}
\right)
$$

The first two terms arise from water autoionization and dominate at extreme pH.
The buffer term $c_\text{tot}\,f_{\ce{HA}}\,f_{\ce{A-}}$ is maximized when $f_{\ce{HA}} = f_{\ce{A-}} = \tfrac{1}{2}$, i.e. at $\text{pH} = \text{p}K_a$:

$$
\beta_\text{max} = \frac{\ln 10}{4}\,c_\text{tot} \approx 0.576\,c_\text{tot}
$$

The buffer capacity falls to half its peak value one pH unit away from $\text{p}K_a$, defining a practical **buffer range** of $\text{p}K_a \pm 1$.

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-buffer-capacity

import numpy as np
import matplotlib.pyplot as plt
from reactions.plots import setup_figure, COLORS

pKa = 4.756
c_tot = 0.1  # mol/L
Kw = 1e-14

pH = np.linspace(1, 9, 500)
h = 10.0 ** (-pH)
oh = Kw / h

f_AcOH = 1.0 / (1.0 + 10.0 ** (pH - pKa))
f_AcO = 1.0 - f_AcOH

ln10 = np.log(10)
beta_water = ln10 * (h + oh)
beta_buffer = ln10 * c_tot * f_AcOH * f_AcO
beta_total = beta_water + beta_buffer

fig, ax = setup_figure()
ax.plot(
    pH, beta_total, color=COLORS["primary"], linewidth=2.5, label=r"$\beta$ (total)"
)
ax.plot(
    pH,
    beta_buffer,
    color="C1",
    linewidth=1.5,
    linestyle="--",
    label="buffer contribution",
)
ax.plot(
    pH,
    beta_water,
    color="C2",
    linewidth=1.5,
    linestyle="--",
    label="water contribution",
)
ax.axvline(pKa, color="gray", linestyle=":", linewidth=1.0)
ax.text(
    pKa + 0.1,
    beta_total.max() * 0.92,
    rf"$\mathrm{{p}}K_a = {pKa}$",
    fontsize=9,
    color="gray",
)
ax.set_xlabel("pH")
ax.set_ylabel(r"Buffer capacity $\beta$ [mol L$^{-1}$]")
ax.set_xlim(1, 9)
ax.legend(fontsize=9)
fig.tight_layout()
```

```{figure} #cell-buffer-capacity
:name: fig-buffer-capacity

Buffer capacity of a 0.1 mol/L acetic acid buffer ($\text{p}K_a = 4.756$).
The buffer contribution peaks at $\text{pH} = \text{p}K_a$; water contributions dominate below pH 3 and above pH 11.
```

**Polyprotic generalisation.** For a system with $n+1$ species carrying $0, 1, \ldots, n$ donated protons, the buffer contribution is the variance of the protonation state weighted by the speciation fractions (King, 1990):

$$
\beta_\text{buffer} = \ln 10 \cdot c_\text{tot}
\left[\sum_{k=0}^{n} k^2 f_k - \left(\sum_{k=0}^{n} k\,f_k\right)^2\right]
$$

For a monoprotic acid the sum reduces to $f_{\ce{HA}}\,f_{\ce{A-}}$, recovering the formula above.
For phosphate ($n = 3$), the dominant contribution near pH 7 comes from the $\ce{H2PO4-}$/$\ce{HPO4^2-}$ pair ($\text{p}K_{a2} = 7.20$), which is why phosphate is the standard physiological buffer.

**Common buffers.** The design rule is to choose a buffer whose $\text{p}K_a$ is within one unit of the target pH.
For coverage over a wider range, multiple buffers with staggered $\text{p}K_a$ values are combined.

| Buffer      | Relevant equilibrium       | p$K_a$ (25 °C)   | Useful range |
| ----------- | -------------------------- | ---------------- | ------------ |
| Citric acid | H$_2$Cit$^-$/HCit$^{2-}$   | 3.13, 4.76, 6.40 | 2.5--6.5     |
| Acetic acid | AcOH/AcO$^-$               | 4.76             | 3.8--5.8     |
| MES         |                            | 6.15             | 5.5--6.7     |
| MOPS        |                            | 7.20             | 6.5--7.9     |
| HEPES       |                            | 7.55             | 6.8--8.2     |
| Phosphate   | H$_2$PO$_4^-$/HPO$_4^{2-}$ | 7.20             | 6.2--8.2     |
| Tris        |                            | 8.06             | 7.0--9.0     |
| Glycine     | $^+$H$_3$N-CH$_2$-COOH     | 9.78             | 8.6--10.6    |

## Ionic strength and the apparent pKa

The equilibrium constants in thermodynamic tables are based on activities.
At nonzero ionic strength, activity coefficients shift the apparent constant that governs observable concentration ratios.

The activity coefficients $\gamma_i(I, z_i)$ follow from Debye-Hückel or Davies theory (@nonidealities, @fig-activity).
Substituting $a_i = \gamma_i c_i/c^\circ$ into $Q = \prod_i a_i^{\nu_i}$ factorises the reaction quotient as $Q = \bigl(\prod_i \gamma_i^{\nu_i}\bigr) K^\text{app}$, where $K^\text{app}$ collects the concentration-ratio terms.
Setting $Q = K$ and taking $-\log_{10}$ gives the general result:

$$
\text{p}K_a^\text{app} = \text{p}K_a + \sum_i \nu_i \log_{10} \gamma_i
$$

where $\nu_i$ are the signed stoichiometric coefficients (negative for reactants).
The Davies equation gives $\log_{10}\gamma_i = -Az_i^2\,f(I)$ with $f(I) = \frac{\sqrt{I}}{1+\sqrt{I}} - 0.3\,I$.
Denoting the charge of the undissociated acid form as $z$, the three species in $\ce{HA <=> A- + H+}$ carry charges $z$, $z-1$, and $+1$, so the sum becomes:

$$
\sum_i \nu_i \log_{10}\gamma_i
= A f(I)\bigl[z^2 - (z-1)^2 - 1\bigr]
= -2A(1-z)\,f(I)
$$

Substituting:

$$
\text{p}K_a^\text{app} = \text{p}K_a
- 2A(1 - z)\!\left(\frac{\sqrt{I}}{1+\sqrt{I}} - 0.3\,I\right)
$$

The factor $(1 - z)$ grows with the charge magnitude of the acid.
For a neutral acid ($z = 0$, meaning the undissociated acid form carries zero net charge, e.g. acetic acid AcOH): the shift is $-2A\,f(I)$, about $-0.24$ at physiological ionic strength ($I = 0.15\ \mathrm{mol/L}$).
For $\ce{H2PO4-}$ ($z = -1$): the factor is 4, giving a shift near $-0.5$.
For $\ce{HPO4^2-}$ ($z = -2$): the factor is 6, giving a shift near $-0.7$.

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-pka-ionic-strength

import numpy as np
import matplotlib.pyplot as plt
from reactions.plots import setup_figure, COLORS

A_davies = 0.509  # log10 Davies constant, water at 25 °C

I = np.linspace(0, 0.5, 400)
f = np.sqrt(I) / (1 + np.sqrt(I)) - 0.3 * I

cases = [
    (0, r"$z = 0$ (neutral acid, e.g. AcOH)"),
    (-1, r"$z = -1$ (e.g. $\mathrm{H_2PO_4^-}$)"),
    (-2, r"$z = -2$ (e.g. $\mathrm{HPO_4^{2-}}$)"),
]

fig, ax = setup_figure()
for z, label in cases:
    shift = -2 * A_davies * (1 - z) * f
    ax.plot(I * 1000, shift, linewidth=2.0, label=label)

ax.axhline(0, color="gray", linewidth=0.8, linestyle=":")
ax.axvline(150, color="gray", linewidth=0.8, linestyle="--")
ax.text(155, -0.05, r"$I = 150\ \mathrm{mmol/L}$", fontsize=8, color="gray")
ax.set_xlabel(r"Ionic strength $I$ [mmol/L]")
ax.set_ylabel(r"$\mathrm{p}K_a^\mathrm{app} - \mathrm{p}K_a$")
ax.legend(fontsize=9, loc="lower left")
fig.tight_layout()
```

```{figure} #cell-pka-ionic-strength
:name: fig-pka-ionic-strength

Apparent pKa shift as a function of ionic strength for three charge-type transitions (Davies equation, 25 °C).
At $I = 150\ \mathrm{mmol/L}$, a neutral acid shifts by $-0.24$; a dianion acid shifts by $-0.72$.
```

## Consequences for chromatography

In ion-exchange chromatography (IEX), buffer pH determines the net charge of every ionisable species, including the protein being separated.
Three effects couple the chemistry developed above to the process:

1. **pH shift on dilution or mixing.** A buffer prepared at low ionic strength delivers a different pH at the column inlet, where it mixes with the mobile phase at higher salt concentration.
2. **Salt gradient coupling.** A simultaneous NaCl gradient changes $I$ along the column, shifting all apparent pKa values and therefore the local pH even if no acid or base is added.
3. **Protein surface effects.** pKa values of surface residues on a resin-bound protein differ from those in free solution; the relevant equilibrium is between solution and surface, not between two solution species.

Designing a buffer system that delivers a predictable, linear pH gradient at the column outlet despite these couplings requires the full machinery of this chapter: Bjerrum diagrams to select components, buffer capacity to choose concentrations, and Davies corrections for the salt gradient.

---

The equilibrium arc is now complete: from the driving force $\Delta_r G$, through the equilibrium constant $K$, to adsorption, temperature dependence, and pH buffering.
The next chapter begins the kinetics arc, asking not where reactions end up but how fast they get there (@kinetics).
A later case study (@case-ph-gradient) shows how the acid-base framework developed here is used to model pH gradients in ion-exchange chromatography.
