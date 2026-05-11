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
    Species, Component,
    ActivityCoefficientDavies,
    IonicStrengthIdeal, IonicStrengthBackground,
    ThermodynamicReaction, ReactionModel,
    pKa,
    water, H_plus, OH_minus,
)
from reactions.solver import solve_equilibrium
```

## Phosphate buffer: pH curve and β

Three coupled algebraic constraints describe phosphate speciation.
Sweeping over target proton activities and re-solving equilibrium at each point gives the speciation curve;
numerical differentiation of the proton balance yields $\beta$:

```{code-cell} ipython3
pKa1, pKa2, pKa3 = 2.148, 7.198, 12.350
Ka1 = 10**(-pKa1); Ka2 = 10**(-pKa2); Ka3 = 10**(-pKa3)
c_phos = 100.0    # mol/m³ total phosphate

phosphate = Component("phosphate", [
    Species("H3PO4",  charge=0),
    Species("H2PO4-", charge=-1),
    Species("HPO4-2", charge=-2),
    Species("PO4-3",  charge=-3),
])

model_phos = ReactionModel(
    components=[phosphate, H_plus, OH_minus, water],
    reactions=[
        ThermodynamicReaction("H3PO4 <-> H2PO4- + H+",  mode="equil",
                              equilibrium_constant=pKa(pKa1)),
        ThermodynamicReaction("H2PO4- <-> HPO4-2 + H+", mode="equil",
                              equilibrium_constant=pKa(pKa2)),
        ThermodynamicReaction("HPO4-2 <-> PO4-3 + H+",  mode="equil",
                              equilibrium_constant=pKa(pKa3)),
        ThermodynamicReaction("H2O <-> H+ + OH-",        mode="equil",
                              equilibrium_constant=pKa(14.00)),
    ],
    T=298.15,
)


def phosphate_fractions(pH):
    h = 10.0**(-pH)
    D = h**3 + Ka1*h**2 + Ka1*Ka2*h + Ka1*Ka2*Ka3
    return h**3/D, Ka1*h**2/D, Ka1*Ka2*h/D, Ka1*Ka2*Ka3/D


def solve_at_pH(model, pH, c_tot):
    a0, a1, a2, a3 = phosphate_fractions(pH)
    H  = 10.0**(-pH) * H_plus.c_ref
    OH = 1e-14 * H_plus.c_ref**2 / H
    c0 = {
        "H3PO4":  max(a0*c_tot, 1e-10),
        "H2PO4-": max(a1*c_tot, 1e-10),
        "HPO4-2": max(a2*c_tot, 1e-10),
        "PO4-3":  max(a3*c_tot, 1e-10),
        "H+": H, "OH-": OH,
    }
    return solve_equilibrium(model, c0, T=298.15, prescribed={"H2O": water.c_ref})


pH_vals = np.linspace(3.0, 12.0, 200)
H_conc  = np.array([solve_at_pH(model_phos, pH, c_phos)["H+"] for pH in pH_vals])
pH_num  = -np.log10(H_conc / H_plus.c_ref)

dpH = np.diff(pH_num)
dH  = np.diff(H_conc)
beta_num = np.abs(dH / dpH) / 1000.0    # mol/L per pH unit

pH_mid  = 0.5 * (pH_num[:-1] + pH_num[1:])

def beta_analytical(pH):
    h  = 10.0**(-pH)
    kw = 1e-14
    a0, a1, a2, a3 = phosphate_fractions(pH)
    beta_water = np.log(10) * (h + kw/h)
    beta_buf   = np.log(10) * c_phos/1000 * (
        a0*a1 + a1*a2*4 + a2*a3*9
    )
    return beta_water + beta_buf

beta_ana = np.array([beta_analytical(pH) for pH in pH_mid])

print(f"{'pH':>6}  {'β num (mol/L/pH)':>18}  {'β ana':>10}  {'rel. diff':>10}")
for i, pH in enumerate(pH_mid[::20]):
    idx = np.searchsorted(pH_mid, pH)
    if idx < len(beta_num):
        diff = abs(beta_num[idx] - beta_ana[idx]) / beta_ana[idx]
        print(f"{pH:>6.2f}  {beta_num[idx]:>18.4f}  {beta_ana[idx]:>10.4f}  {diff:>10.3f}")
```

The numerical curve matches the analytical formula from @speciation-buffers; small deviations arise from finite-difference step size at the extremes of the pH range.
In the analytical expression, the coefficients 1, 4, 9 are the squared charge changes $(\Delta z)^2$ between successive protonation states; higher charges amplify the sensitivity of charge redistribution to proton activity.

## Mixed buffer: citrate + phosphate

No single buffer covers a wide pH range with flat capacity.
Citrate (three p$K_a$ values spanning pH 3--7) combined with phosphate gives a nearly flat $\beta$ from pH 3 to 8.
The apparent additivity holds because the two acid systems couple only weakly through the shared proton reservoir; the full equilibrium solve captures the residual coupling:

```{code-cell} ipython3
pKa_cit = [3.128, 4.761, 6.396]    # NIST
Ka_cit  = [10**(-p) for p in pKa_cit]

citrate = Component("citrate", [
    Species("H3Cit",  charge=0),
    Species("H2Cit-", charge=-1),
    Species("HCit-2", charge=-2),
    Species("Cit-3",  charge=-3),
])

c_cit  = 50.0    # mol/m³
c_phos = 50.0    # mol/m³

model_mixed = ReactionModel(
    components=[citrate, phosphate, H_plus, OH_minus, water],
    reactions=[
        ThermodynamicReaction("H3Cit <-> H2Cit- + H+",  mode="equil",
                              equilibrium_constant=pKa(pKa_cit[0])),
        ThermodynamicReaction("H2Cit- <-> HCit-2 + H+", mode="equil",
                              equilibrium_constant=pKa(pKa_cit[1])),
        ThermodynamicReaction("HCit-2 <-> Cit-3 + H+",  mode="equil",
                              equilibrium_constant=pKa(pKa_cit[2])),
        ThermodynamicReaction("H3PO4 <-> H2PO4- + H+",  mode="equil",
                              equilibrium_constant=pKa(pKa1)),
        ThermodynamicReaction("H2PO4- <-> HPO4-2 + H+", mode="equil",
                              equilibrium_constant=pKa(pKa2)),
        ThermodynamicReaction("HPO4-2 <-> PO4-3 + H+",  mode="equil",
                              equilibrium_constant=pKa(pKa3)),
        ThermodynamicReaction("H2O <-> H+ + OH-",        mode="equil",
                              equilibrium_constant=pKa(14.00)),
    ],
    T=298.15,
)


def citrate_fractions(pH):
    h = 10.0**(-pH)
    D = (h**3 + Ka_cit[0]*h**2 + Ka_cit[0]*Ka_cit[1]*h
         + Ka_cit[0]*Ka_cit[1]*Ka_cit[2])
    return h**3/D, Ka_cit[0]*h**2/D, Ka_cit[0]*Ka_cit[1]*h/D, Ka_cit[0]*Ka_cit[1]*Ka_cit[2]/D


def solve_mixed_at_pH(pH):
    a0c, a1c, a2c, a3c = citrate_fractions(pH)
    a0p, a1p, a2p, a3p = phosphate_fractions(pH)
    H  = 10.0**(-pH) * H_plus.c_ref
    OH = 1e-14 * H_plus.c_ref**2 / H
    c0 = {
        "H3Cit":  max(a0c*c_cit,  1e-10),
        "H2Cit-": max(a1c*c_cit,  1e-10),
        "HCit-2": max(a2c*c_cit,  1e-10),
        "Cit-3":  max(a3c*c_cit,  1e-10),
        "H3PO4":  max(a0p*c_phos, 1e-10),
        "H2PO4-": max(a1p*c_phos, 1e-10),
        "HPO4-2": max(a2p*c_phos, 1e-10),
        "PO4-3":  max(a3p*c_phos, 1e-10),
        "H+": H, "OH-": OH,
    }
    return solve_equilibrium(model_mixed, c0, T=298.15, prescribed={"H2O": water.c_ref})


pH_mix  = np.linspace(3.0, 10.0, 150)
H_mix   = np.array([solve_mixed_at_pH(pH)["H+"] for pH in pH_mix])
pH_m    = -np.log10(H_mix / H_plus.c_ref)
beta_m  = np.abs(np.diff(H_mix) / np.diff(pH_m)) / 1000.0
pH_m_mid = 0.5 * (pH_m[:-1] + pH_m[1:])

print("Mixed citrate + phosphate buffer capacity:")
print(f"{'pH':>6}  {'β (mol/L/pH)':>14}")
for pH, b in zip(pH_m_mid[::15], beta_m[::15]):
    print(f"{pH:>6.2f}  {b:>14.4f}")
```

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
        ThermodynamicReaction("H3PO4 <-> H2PO4- + H+",  mode="equil",
                              equilibrium_constant=pKa(pKa1),
                              activity_coefficient=ActivityCoefficientDavies()),
        ThermodynamicReaction("H2PO4- <-> HPO4-2 + H+", mode="equil",
                              equilibrium_constant=pKa(pKa2),
                              activity_coefficient=ActivityCoefficientDavies()),
        ThermodynamicReaction("HPO4-2 <-> PO4-3 + H+",  mode="equil",
                              equilibrium_constant=pKa(pKa3),
                              activity_coefficient=ActivityCoefficientDavies()),
        ThermodynamicReaction("H2O <-> H+ + OH-",        mode="equil",
                              equilibrium_constant=pKa(14.00),
                              activity_coefficient=ActivityCoefficientDavies()),
    ],
    ionic_strength=IonicStrengthBackground(I_bg=150.0),
    T=298.15,
)

pH_dav = np.linspace(3.0, 11.0, 180)
H_dav = np.array([
    solve_equilibrium(model_phos_davies, solve_at_pH(model_phos, pH, c_phos), T=298.15, prescribed={"H2O": water.c_ref})["H+"]
    for pH in pH_dav
])
pH_dav  = -np.log10(H_dav / H_plus.c_ref)
beta_dav = np.abs(np.diff(H_dav) / np.diff(pH_dav)) / 1000.0

idx_ideal  = np.argmax(beta_num)
idx_davies = np.argmax(beta_dav)
print(f"Peak β (ideal):  {beta_num[idx_ideal]:.4f} mol/L/pH  at pH {pH_mid[idx_ideal]:.2f}")
print(f"Peak β (Davies): {beta_dav[idx_davies]:.4f} mol/L/pH  at pH {0.5*(pH_dav[idx_davies]+pH_dav[idx_davies+1]):.2f}")
```

The Davies correction shifts the peak by roughly $0.1$--$0.3$ pH units at
physiological ionic strength.
In IEX gradient design this shift is not negligible: a buffer prepared at low ionic
strength will deliver a different pH on-column once the salt gradient is applied
(@case-ph-gradient).

---

The next chapter extends the rate framework to saturation kinetics, where finite enzyme active sites produce a concentration-dependent rate that mass action cannot capture (@implementation-enzyme).
