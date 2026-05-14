---
kernelspec:
  name: python3
  display_name: Python 3 (ipykernel)
---

(acid-base)=
# Acid-Base Equilibria and pH

Proton transfer between a weak acid and water is fast: the bimolecular rate constant for H$^+$ transfer in aqueous solution is of order $10^{10}$ M$^{-1}$s$^{-1}$, many orders of magnitude faster than the transport timescale in any chromatographic column.
Proton-transfer reactions are therefore treated as instantaneous equilibria throughout this book, and pH is computed from $Q = K$ rather than from a kinetic rate law.

## Acids and bases

Three definitions are in common use, each generalising the previous.

**Arrhenius (1884).** An acid produces H$^+$ in water; a base produces OH$^-$.
Historically first, but limited to aqueous solution and excludes many important bases (e.g. ammonia, which does not contain OH$^-$).

**Brønsted-Lowry (1923).** An acid is a proton donor; a base is a proton acceptor.
Every acid has a conjugate base (the species left after donating H$^+$) and every base has a conjugate acid.
This is the standard definition for aqueous acid-base chemistry and the one used throughout this book.

**Lewis (1923).** An acid is an electron-pair acceptor; a base is an electron-pair donor.
The Brønsted-Lowry definition is a special case.
Lewis acid-base chemistry is important in organometallic and coordination chemistry but is not needed for the buffer and pH problems addressed here.

```{important}
**Strong vs weak acids.**
A strong acid (p$K_a \ll 0$, e.g. HCl) is essentially fully dissociated.
A weak acid (p$K_a > 0$, e.g. acetic acid) exists in partial equilibrium.

"Strong" does not mean "concentrated."
At $c = 10^{-8}$ mol/L, HCl contributes only $10^{-8}$ mol/L of H$^+$, while water autoionization contributes $\approx 10^{-7}$ mol/L.
The resulting pH is $\approx 6.98$, not 8.
Neglecting the water equilibrium at very low concentrations of strong acids is a common source of error.
```

## Proton transfer as chemical equilibrium

A Brønsted acid HA donates a proton to water:

$$
\ce{HA <=> H+ + A-}
$$

Its equilibrium constant, the **acid dissociation constant**, follows directly from the equilibrium condition $\Delta_r G = 0$ (see @equilibrium):

$$
K_a = \frac{a_{\ce{H+}}\,a_{\ce{A-}}}{a_{\ce{HA}}}
$$

In dilute ideal solution $a_i = c_i/c^\circ$ and water is treated as a pure solvent with $a_{\ce{H2O}} = 1$, giving:

$$
K_a = \frac{[\ce{H+}][\ce{A-}]}{[\ce{HA}]}
$$

where brackets denote concentration relative to the standard state $c^\circ = 1\ \mathrm{mol/L}$.

Water itself autoionises with equilibrium constant:

$$
K_w = [\ce{H+}][\ce{OH-}]
$$

At 25 °C, $K_w = 1.01 \times 10^{-14}$.
The van't Hoff equation (@equilibrium-temperature) with $\Delta H^\circ \approx +55.8\ \mathrm{kJ/mol}$ predicts that $K_w$ rises by roughly one order of magnitude between 0 and 60 °C, shifting the neutral pH from 7.47 to 6.51.

The **pKa** scale compresses the wide dynamic range of acid strengths:

$$
\text{p}K_a = -\log_{10} K_a
$$

Strong acids ($\text{p}K_a < 0$) are essentially fully dissociated.
Weak acids ($\text{p}K_a > 0$) exist in partially protonated form at analytical concentrations.
Biochemically relevant buffers typically have $\text{p}K_a$ in the range 2--12.

## pH and Henderson-Hasselbalch

The **pH** is the negative base-10 logarithm of the proton activity:

$$
\text{pH} = -\log_{10} a_{\ce{H+}}
$$

In dilute ideal solution, $a_{\ce{H+}} \approx [\ce{H+}]$.
At equilibrium, $Q = K_a$ gives $[\ce{H+}] = K_a\,[\ce{HA}]/[\ce{A-}]$.
Taking the negative logarithm:

$$
\text{pH} = \text{p}K_a + \log_{10}\!\frac{[\ce{A-}]}{[\ce{HA}]}
$$

This is the **Henderson-Hasselbalch equation** {cite:p}`henderson1908,hasselbalch1916`.
At $\text{pH} = \text{p}K_a$ the two forms are present in equal amounts.
The ratio $[\ce{A-}]/[\ce{HA}]$ changes by a factor of ten for each pH unit.

---

Henderson-Hasselbalch describes a single acid at a single pH.
The next chapter extends this to speciation across the full pH range (Bjerrum diagrams), buffer capacity, and the ionic-strength corrections that shift apparent pKa values in real solutions.
