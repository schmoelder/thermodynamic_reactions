---
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

(implementation-buffer)=
# Buffer Capacity

Buffer capacity $\beta = -dn_b / d\text{pH}$ measures how much strong base $n_b$ (mol/L) is needed to shift the pH by one unit (@speciation-buffers).
In equilibrium calculations, $\beta$ is computed from the sensitivity of the equilibrium composition to changes in proton activity; no explicit base addition is required; $n_b$ is implicitly represented by the proton mass balance.
Computing it numerically requires only the equilibrium model from @implementation-acid-base: sweep over target proton activities, re-solve equilibrium consistently at each point, then differentiate the resulting proton balance numerically.

```{code-cell} ipython3
import numpy as np

from reactions.api import (
    Species,
    Component,
    ActivityCoefficientDavies,
    IonicStrengthIdeal,
    IonicStrengthBackground,
    ThermodynamicReaction,
    ReactionModel,
    pKa,
    water,
    H_plus,
    OH_minus,
    buffer_capacity,
)
from reactions.solver import solve_equilibrium
```

## Phosphate buffer: pH curve and β

Three coupled algebraic constraints describe phosphate speciation.
Sweeping over target proton activities and re-solving equilibrium at each point gives the speciation curve;
numerical differentiation of the proton balance yields $\beta$:

```{code-cell} ipython3
:tags: [remove-cell]

pKa1, pKa2, pKa3 = 2.148, 7.198, 12.350
Ka1 = 10 ** (-pKa1)
Ka2 = 10 ** (-pKa2)
Ka3 = 10 ** (-pKa3)
c_phos = 100.0  # mol/m³ total phosphate

phosphate = Component(
    "phosphate",
    [
        Species("H3PO4", charge=0),
        Species("H2PO4-", charge=-1),
        Species("HPO4-2", charge=-2),
        Species("PO4-3", charge=-3),
    ],
)

model_phos = ReactionModel(
    components=[phosphate, H_plus, OH_minus, water],
    reactions=[
        ThermodynamicReaction(
            "H3PO4 <-> H2PO4- + H+", mode="equil", equilibrium_constant=pKa(pKa1)
        ),
        ThermodynamicReaction(
            "H2PO4- <-> HPO4-2 + H+", mode="equil", equilibrium_constant=pKa(pKa2)
        ),
        ThermodynamicReaction(
            "HPO4-2 <-> PO4-3 + H+", mode="equil", equilibrium_constant=pKa(pKa3)
        ),
        ThermodynamicReaction(
            "H2O <-> H+ + OH-", mode="equil", equilibrium_constant=pKa(14.00)
        ),
    ],
    T=298.15,
)


def phosphate_fractions(pH):
    h = 10.0 ** (-pH)
    D = h**3 + Ka1 * h**2 + Ka1 * Ka2 * h + Ka1 * Ka2 * Ka3
    return h**3 / D, Ka1 * h**2 / D, Ka1 * Ka2 * h / D, Ka1 * Ka2 * Ka3 / D


def solve_at_pH(model, pH, c_tot):
    a0, a1, a2, a3 = phosphate_fractions(pH)
    H = 10.0 ** (-pH) * H_plus.c_ref
    OH = 1e-14 * H_plus.c_ref**2 / H
    c0 = {
        "H3PO4": max(a0 * c_tot, 1e-10),
        "H2PO4-": max(a1 * c_tot, 1e-10),
        "HPO4-2": max(a2 * c_tot, 1e-10),
        "PO4-3": max(a3 * c_tot, 1e-10),
        "H+": H,
        "OH-": OH,
    }
    return solve_equilibrium(model, c0, T=298.15, prescribed={"H2O": water.c_ref})
```

```{code-cell} ipython3
pH_vals = np.linspace(3.0, 12.0, 200)
eq_phos = [solve_at_pH(model_phos, pH, c_phos) for pH in pH_vals]
pH_num = np.array([-np.log10(r["H+"] / H_plus.c_ref) for r in eq_phos])

Q = np.array([r["H+"] - r["OH-"] - r["H2PO4-"] - 2*r["HPO4-2"] - 3*r["PO4-3"] for r in eq_phos])
beta_num = np.abs(np.diff(Q) / np.diff(pH_num))
pH_mid = 0.5 * (pH_num[:-1] + pH_num[1:])
```

```{code-cell} ipython3
:tags: [remove-cell]

_beta_lib = buffer_capacity(model_phos, {"phosphate": c_phos}, pH_mid)
beta_ana = _beta_lib["phosphate"] + _beta_lib["water"]
```

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-beta-phosphate

from reactions.plots import setup_figure

fig, ax = setup_figure()
ax.plot(pH_mid, beta_num, color="C0", lw=2.0, label="numerical")
ax.plot(pH_mid, beta_ana, color="C1", lw=1.5, ls="--", label="analytical")
for pKa_val, lbl in [
    (pKa1, r"p$K_{a1}$"),
    (pKa2, r"p$K_{a2}$"),
    (pKa3, r"p$K_{a3}$"),
]:
    ax.axvline(pKa_val, color="gray", ls=":", lw=0.8)
    ax.text(
        pKa_val + 0.08,
        0.95,
        lbl,
        fontsize=8,
        color="gray",
        transform=ax.get_xaxis_transform(),
        va="top",
    )
ax.set_xlabel("pH")
ax.set_ylabel(r"$\beta$ / (mol m$^{-3}$ pH$^{-1}$)")
ax.set_xlim(3, 12)
ax.legend(fontsize=10)
fig.tight_layout()
```

```{figure} #cell-beta-phosphate
:name: fig-beta-phosphate

Buffer capacity $\beta$ of 100 mol/m³ phosphate: numerical finite-difference (solid) versus `buffer_capacity` from `reactions.analysis` (dashed).
All three peaks reach the same height because every adjacent phosphate transition is a unit charge step ($\Delta z = 1$), giving $\beta_\text{max} = \ln(10)\,c/4$ at each p$K_a$.
Small deviations at the pH extremes arise from the finite-difference step size.
```

The Van Slyke formula sums $(\Delta z)^2 \alpha_i \alpha_j$ over all species pairs; for phosphate every adjacent transition carries $\Delta z = 1$, so all three peaks are equal.
A buffer with non-unit charge steps between adjacent forms — for example a species jumping directly from $z = 0$ to $z = -2$ — would produce a taller peak at that transition because $(\Delta z)^2 = 4$.

## Mixed buffer: citrate + phosphate

No single buffer covers a wide pH range with flat capacity.
Citrate (three p$K_a$ values spanning pH 3--7) combined with phosphate gives a nearly flat $\beta$ from pH 3 to 8:

```{code-cell} ipython3
:tags: [remove-cell]

pKa_cit = [3.128, 4.761, 6.396]  # NIST
Ka_cit = [10 ** (-p) for p in pKa_cit]
c_cit = 50.0  # mol/m³
c_phos = 50.0  # mol/m³

citrate = Component(
    "citrate",
    [
        Species("H3Cit", charge=0),
        Species("H2Cit-", charge=-1),
        Species("HCit-2", charge=-2),
        Species("Cit-3", charge=-3),
    ],
)
```

```{code-cell} ipython3
model_mixed = ReactionModel(
    components=[citrate, phosphate, H_plus, OH_minus, water],
    reactions=[
        ThermodynamicReaction(
            "H3Cit <-> H2Cit- + H+", mode="equil", equilibrium_constant=pKa(3.128)
        ),
        ThermodynamicReaction(
            "H2Cit- <-> HCit-2 + H+", mode="equil", equilibrium_constant=pKa(4.761)
        ),
        ThermodynamicReaction(
            "HCit-2 <-> Cit-3 + H+", mode="equil", equilibrium_constant=pKa(6.396)
        ),
        ThermodynamicReaction(
            "H3PO4 <-> H2PO4- + H+", mode="equil", equilibrium_constant=pKa(pKa1)
        ),
        ThermodynamicReaction(
            "H2PO4- <-> HPO4-2 + H+", mode="equil", equilibrium_constant=pKa(pKa2)
        ),
        ThermodynamicReaction(
            "HPO4-2 <-> PO4-3 + H+", mode="equil", equilibrium_constant=pKa(pKa3)
        ),
        ThermodynamicReaction(
            "H2O <-> H+ + OH-", mode="equil", equilibrium_constant=pKa(14.00)
        ),
    ],
    T=298.15,
)
```

```{code-cell} ipython3
:tags: [remove-cell]

def citrate_fractions(pH):
    h = 10.0 ** (-pH)
    D = (
        h**3
        + Ka_cit[0] * h**2
        + Ka_cit[0] * Ka_cit[1] * h
        + Ka_cit[0] * Ka_cit[1] * Ka_cit[2]
    )
    return (
        h**3 / D,
        Ka_cit[0] * h**2 / D,
        Ka_cit[0] * Ka_cit[1] * h / D,
        Ka_cit[0] * Ka_cit[1] * Ka_cit[2] / D,
    )


def solve_mixed_at_pH(pH):
    a0c, a1c, a2c, a3c = citrate_fractions(pH)
    a0p, a1p, a2p, a3p = phosphate_fractions(pH)
    H = 10.0 ** (-pH) * H_plus.c_ref
    OH = 1e-14 * H_plus.c_ref**2 / H
    c0 = {
        "H3Cit": max(a0c * c_cit, 1e-10),
        "H2Cit-": max(a1c * c_cit, 1e-10),
        "HCit-2": max(a2c * c_cit, 1e-10),
        "Cit-3": max(a3c * c_cit, 1e-10),
        "H3PO4": max(a0p * c_phos, 1e-10),
        "H2PO4-": max(a1p * c_phos, 1e-10),
        "HPO4-2": max(a2p * c_phos, 1e-10),
        "PO4-3": max(a3p * c_phos, 1e-10),
        "H+": H,
        "OH-": OH,
    }
    return solve_equilibrium(model_mixed, c0, T=298.15, prescribed={"H2O": water.c_ref})


pH_mix = np.linspace(3.0, 10.0, 150)
eq_mix = [solve_mixed_at_pH(pH) for pH in pH_mix]
pH_m = np.array([-np.log10(r["H+"] / H_plus.c_ref) for r in eq_mix])
Q_mix = np.array([
    r["H+"] - r["OH-"]
    - r["H2Cit-"] - 2*r["HCit-2"] - 3*r["Cit-3"]
    - r["H2PO4-"] - 2*r["HPO4-2"] - 3*r["PO4-3"]
    for r in eq_mix
])
beta_m = np.abs(np.diff(Q_mix) / np.diff(pH_m))
pH_m_mid = 0.5 * (pH_m[:-1] + pH_m[1:])
```

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-beta-mixed

from reactions.plots import setup_figure

fig, ax = setup_figure()
ax.plot(pH_m_mid, beta_m, color="C0", lw=2.0)
for pKa_val, lbl in [
    (pKa_cit[0], r"p$K_{a1}^{\rm cit}$"),
    (pKa_cit[1], r"p$K_{a2}^{\rm cit}$"),
    (pKa_cit[2], r"p$K_{a3}^{\rm cit}$"),
    (pKa2, r"p$K_{a2}^{\rm phos}$"),
]:
    ax.axvline(pKa_val, color="gray", ls=":", lw=0.8)
    ax.text(
        pKa_val + 0.06,
        0.95,
        lbl,
        fontsize=7,
        color="gray",
        transform=ax.get_xaxis_transform(),
        va="top",
    )
ax.set_xlabel("pH")
ax.set_ylabel(r"$\beta$ / (mol m$^{-3}$ pH$^{-1}$)")
ax.set_xlim(3, 10)
fig.tight_layout()
```

```{figure} #cell-beta-mixed
:name: fig-beta-mixed

Buffer capacity of a mixed citrate (50 mol/m³) and phosphate (50 mol/m³) system.
The three citrate p$K_a$ values (pH 3.1, 4.8, 6.4) and the second phosphate p$K_a$ (pH 7.2) produce overlapping peaks that yield nearly flat $\beta$ from pH 3 to 8.
```

The apparent additivity of the two buffer systems holds because citrate and phosphate couple only weakly through the shared proton reservoir; the full equilibrium solve captures the residual coupling.

## Ionic strength effect on β

A salt background shifts apparent pKa values (@speciation-buffers, @implementation-activity),
which broadens and shifts the buffer capacity peak.
The ideal speciation is used as the initial guess; the non-ideal model then fully re-equilibrates under the activity-corrected residuals.
The same model with `IonicStrengthBackground` and `ActivityCoefficientDavies` shows the ionic strength correction.
The sweep is limited to pH 3--11: above pH 11.5, the Davies correction for $\ce{PO4^3-}$ ($z = -3$) becomes large enough to make the Newton Jacobian ill-conditioned, and the solver diverges.
This is outside the practical IEX range, so the truncation has no consequence here; a more robust solver or a different activity model (e.g. Pitzer) would extend the range:

```{code-cell} ipython3
model_phos_davies = ReactionModel(
    components=[phosphate, H_plus, OH_minus, water],
    reactions=[
        ThermodynamicReaction(
            "H3PO4 <-> H2PO4- + H+",
            mode="equil",
            equilibrium_constant=pKa(pKa1),
            activity_coefficient=ActivityCoefficientDavies(),
        ),
        ThermodynamicReaction(
            "H2PO4- <-> HPO4-2 + H+",
            mode="equil",
            equilibrium_constant=pKa(pKa2),
            activity_coefficient=ActivityCoefficientDavies(),
        ),
        ThermodynamicReaction(
            "HPO4-2 <-> PO4-3 + H+",
            mode="equil",
            equilibrium_constant=pKa(pKa3),
            activity_coefficient=ActivityCoefficientDavies(),
        ),
        ThermodynamicReaction(
            "H2O <-> H+ + OH-",
            mode="equil",
            equilibrium_constant=pKa(14.00),
            activity_coefficient=ActivityCoefficientDavies(),
        ),
    ],
    ionic_strength=IonicStrengthBackground(I_bg=150.0),
    T=298.15,
)
```

```{code-cell} ipython3
:tags: [remove-cell]

pH_dav = np.linspace(3.0, 11.0, 180)
eq_dav = [
    solve_equilibrium(
        model_phos_davies,
        solve_at_pH(model_phos, pH, c_phos),
        T=298.15,
        prescribed={"H2O": water.c_ref},
    )
    for pH in pH_dav
]
pH_dav = np.array([-np.log10(r["H+"] / H_plus.c_ref) for r in eq_dav])
Q_dav = np.array([r["H+"] - r["OH-"] - r["H2PO4-"] - 2*r["HPO4-2"] - 3*r["PO4-3"] for r in eq_dav])
beta_dav = np.abs(np.diff(Q_dav) / np.diff(pH_dav))
```

```{code-cell} ipython3
:tags: [remove-cell]
:label: cell-beta-davies

from reactions.plots import setup_figure

pH_dav_mid = 0.5 * (pH_dav[:-1] + pH_dav[1:])
_beta_lib_dav = buffer_capacity(model_phos, {"phosphate": c_phos}, pH_dav_mid)
beta_ideal_dav = _beta_lib_dav["phosphate"] + _beta_lib_dav["water"]

fig, ax = setup_figure()
ax.plot(pH_dav_mid, beta_ideal_dav, color="C0", lw=2.0, label=r"ideal ($I = 0$)")
ax.plot(pH_dav_mid, beta_dav, color="C1", lw=2.0, label=r"Davies ($I = 150\ \mathrm{mM}$)")
ax.axvline(pKa2, color="gray", ls=":", lw=0.8)
ax.text(
    pKa2 + 0.08,
    0.95,
    r"p$K_{a2}$",
    fontsize=8,
    color="gray",
    transform=ax.get_xaxis_transform(),
    va="top",
)
ax.set_xlabel("pH")
ax.set_ylabel(r"$\beta$ / (mol m$^{-3}$ pH$^{-1}$)")
ax.set_xlim(3, 11)
ax.legend(fontsize=10)
fig.tight_layout()
```

```{figure} #cell-beta-davies
:name: fig-beta-davies

Ionic strength effect on buffer capacity for 50 mol/m³ phosphate: ideal ($I = 0$) versus Davies correction at physiological ionic strength ($I = 150\ \mathrm{mM}$).
The activity correction shifts the p$K_{a2}$ peak to lower pH by roughly $0.1$--$0.3$ units and reduces the peak height slightly.
```

In IEX gradient design this shift is not negligible: a buffer prepared at low ionic strength will deliver a different pH on-column once the salt gradient is applied (@case-ph-gradient).

---

The next chapter extends the rate framework to saturation kinetics, where finite enzyme active sites produce a concentration-dependent rate that mass action cannot capture (@implementation-enzyme).
