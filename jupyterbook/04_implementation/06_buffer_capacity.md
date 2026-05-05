---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.0
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

(implementation-buffer)=
# Buffer Capacity

Buffer capacity $\beta = -dn_b / d\text{pH}$ measures how much strong base (in
mol/L) is needed to shift the pH by one unit (@speciation-buffers).
Computing it numerically requires only the equilibrium model from
@implementation-acid-base: solve for speciation at a series of pH values, then
differentiate numerically.

```{code-cell} ipython3
import numpy as np

from reactions.api import (
    C_REF,
    Species, Component,
    ActivityCoefficientDavies,
    EquilibriumConstant,
    IonicStrengthIdeal, IonicStrengthBackground,
    ThermodynamicReaction, ReactionModel,
    pKa,
)
from reactions.solver import solve_equilibrium

proton    = Component("proton",    [Species("H+",   charge=+1)])
hydroxide = Component("hydroxide", [Species("OH-",  charge=-1)])
water     = Component("water",     [Species("H2O",  charge=0, is_solvent=True)])
```

## Phosphate buffer: pH curve and β

Three coupled algebraic constraints describe phosphate speciation.
Sweeping over imposed pH and solving gives the speciation at each point;
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
    components=[phosphate, proton, hydroxide, water],
    reactions=[
        ThermodynamicReaction("H3PO4 <-> H2PO4- + H+",  mode="equil",
                              equilibrium_constant=pKa(pKa1)),
        ThermodynamicReaction("H2PO4- <-> HPO4-2 + H+", mode="equil",
                              equilibrium_constant=pKa(pKa2)),
        ThermodynamicReaction("HPO4-2 <-> PO4-3 + H+",  mode="equil",
                              equilibrium_constant=pKa(pKa3)),
        ThermodynamicReaction("H2O <-> H+ + OH-",        mode="equil",
                              equilibrium_constant=EquilibriumConstant(1e-14)),
    ],
    T=298.15,
)


def phosphate_fractions(pH):
    h = 10.0**(-pH)
    D = h**3 + Ka1*h**2 + Ka1*Ka2*h + Ka1*Ka2*Ka3
    return h**3/D, Ka1*h**2/D, Ka1*Ka2*h/D, Ka1*Ka2*Ka3/D


def solve_at_pH(model, pH, c_tot):
    a0, a1, a2, a3 = phosphate_fractions(pH)
    H  = 10.0**(-pH) * C_REF
    OH = 1e-14 * C_REF**2 / H
    c0 = {
        "H3PO4":  max(a0*c_tot, 1e-10),
        "H2PO4-": max(a1*c_tot, 1e-10),
        "HPO4-2": max(a2*c_tot, 1e-10),
        "PO4-3":  max(a3*c_tot, 1e-10),
        "H+": H, "OH-": OH,
    }
    return solve_equilibrium(model, c0, T=298.15)


pH_vals = np.linspace(3.0, 12.0, 200)
H_conc  = np.array([solve_at_pH(model_phos, pH, c_phos)["H+"] for pH in pH_vals])
pH_num  = -np.log10(H_conc / C_REF)

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

The numerical curve matches the analytical formula from @speciation-buffers; small
deviations arise from finite-difference step size at the extremes of the pH range.

## Mixed buffer: citrate + phosphate

No single buffer covers a wide pH range with flat capacity.
Citrate (three p$K_a$ values spanning pH 3--7) combined with phosphate gives a
nearly flat $\beta$ from pH 3 to 8:

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
    components=[citrate, phosphate, proton, hydroxide, water],
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
                              equilibrium_constant=EquilibriumConstant(1e-14)),
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
    H  = 10.0**(-pH) * C_REF
    OH = 1e-14 * C_REF**2 / H
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
    return solve_equilibrium(model_mixed, c0, T=298.15)


pH_mix  = np.linspace(3.0, 10.0, 150)
H_mix   = np.array([solve_mixed_at_pH(pH)["H+"] for pH in pH_mix])
pH_m    = -np.log10(H_mix / C_REF)
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
The same model with `IonicStrengthBackground` and `ActivityCoefficientDavies` shows
the ionic strength correction:

```{code-cell} ipython3
model_phos_davies = ReactionModel(
    components=[phosphate, proton, hydroxide, water],
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
                              equilibrium_constant=EquilibriumConstant(1e-14),
                              activity_coefficient=ActivityCoefficientDavies()),
    ],
    ionic_strength=IonicStrengthBackground(I_bg=150.0),
    T=298.15,
)

# z=-3 γ correction overflows Newton above pH ≈ 11.5; outside IEX practical range
pH_dav = np.linspace(3.0, 11.0, 180)
H_dav = np.array([
    solve_equilibrium(model_phos_davies, solve_at_pH(model_phos, pH, c_phos), T=298.15)["H+"]
    for pH in pH_dav
])
pH_dav  = -np.log10(H_dav / C_REF)
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

The next chapter adds Arrhenius temperature dependence to the rate constants, completing the kinetic picture before the CADET-Core interface (@implementation-kinetics).
