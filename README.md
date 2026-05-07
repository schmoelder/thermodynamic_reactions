# Thermodynamic Reactions for CADET

**A thermodynamically consistent reaction library, built from first principles.**

## Why?

[CADET](https://cadet.github.io) is a process modeling framework for chromatographic separations.
Its current reaction models can handle pH-dependent speciation, but lack thermodynamic consistency: rate constant ratios are empirical, and non-ideal effects like ionic strength are not accounted for.

This library fixes that: equilibrium constants derive from ΔrG°, rate constants respect `k_f/k_r = K(T)` by construction, and activity corrections handle non-ideal solutions.
It works today as a standalone Python library and educational resource, and serves as a foundation for future CADET integration.

## What it provides

This project provides a reaction framework where:
- Equilibrium constants derive from ΔrG° (not fitted numbers)
- Rate constants respect `k_f/k_r = K` at all temperatures
- Activity corrections (Debye-Hückel, Davies) apply throughout
- A unified `Species/Component/ReactionModel` interface covers equilibrium, kinetic, and enzymatic reactions

Developed alongside an educational [JupyterBook](jupyterbook/) that reconstructs reaction thermodynamics from statistical mechanics to CADET implementation.

---

## 📚 Project Structure

```
reactions/          # pip-installable Python library
├── api.py          # Public API: Species, Component, ReactionModel, ...
└── ...             # ThermodynamicReaction, MassActionReaction, EnzymaticReaction

jupyterbook/        # MyST JupyterBook (6 parts, 20+ chapters)
├── 01_statistical_mechanics/  # Part 1: From particles to thermodynamics
├── 02_thermodynamics/           # Part 2: Fundamentals (G, μ, activity)
├── 03_reactions/                 # Part 3: Equilibrium & kinetics
├── 04_implementation/            # Part 4: CADET integration
└── 05_case_studies/              # Part 5: pH gradients, IEX, buffers

tests/              # pytest suite
pyproject.toml      # Package metadata
```

**Book Status:** Parts 1–3 complete, Part 4 in progress, Part 5 outlined.

---

## ⚙️ Installation

```bash
# Core library
pip install -e .

# With optional dependencies
pip install -e . --group test    # pytest
pip install -e . --group docs    # jupyter-book, griffe
```

---

## 🔬 Quick Start

```python
from reactions.api import Component, ThermodynamicReaction, ReactionModel
from reactions.api import EquilibriumConstantVantHoff, RateConstantArrhenius
from reactions.solver import simulate, solve_equilibrium

# Equilibrium mode: K(T) from van't Hoff; no rate constant needed
model = ReactionModel(
    components=[Component("A"), Component("B")],
    reactions=[
        ThermodynamicReaction(
            "A <-> B",
            mode="equil",
            equilibrium_constant=EquilibriumConstantVantHoff(dH=-20e3, dS=-50.0),
        ),
    ],
)

# Concentrations in mol/m³ (SI); 1 mol/L = 1000 mol/m³
result = solve_equilibrium(model, c0={"A": 1000.0}, T=298.15)
print(result)  # {"A": 200.0, "B": 800.0}

# Kinetic mode: add Arrhenius rate constant; integrate forward in time
model_kin = ReactionModel(
    components=[Component("A"), Component("B")],
    reactions=[
        ThermodynamicReaction(
            "A <-> B",
            mode="kinetic",
            equilibrium_constant=EquilibriumConstantVantHoff(dH=-20e3, dS=-50.0),
            rate_constant=RateConstantArrhenius(A=1e10, Ea=40e3),
        ),
    ],
)

result = simulate(model_kin, c0={"A": 1000.0}, t_span=(0, 10.0), T=298.15)
# T(t) callable for temperature programmes:
result = simulate(model_kin, c0={"A": 1000.0}, t_span=(0, 10.0),
                  T=lambda t: 298.15 + 2.0 * t)
```

---
## 📖 Key Concepts

| Concept                       | Implementation                             | Purpose                                                |
| ----------------------------- | ------------------------------------------ | ------------------------------------------------------ |
| **Thermodynamic Consistency** | `kr = kf/K(T)` enforced                    | Rate constants respect equilibrium at all T            |
| **Activity Corrections**      | `ActivityCoefficientDavies()`              | Non-ideal solutions (I ≤ 0.5 M)                        |
| **Equilibrium Mode**          | `mode="equil"`                             | Fast reactions as DAE constraints                      |
| **Kinetic Mode**              | `mode="kinetic"` + `RateConstantArrhenius` | Finite-rate relaxation to equilibrium                  |
| **Saturation Kinetics**       | `EnzymaticReaction`                        | Langmuir, Michaelis-Menten from finite-site constraint |
| **Temperature Programmes**    | `T=lambda t: ...` in `simulate()`          | Prescribed T(t); stored in `result.T_profile`          |

---
## 🧪 Testing

```bash
pytest tests/           # Run full test suite
myst build --html      # Build the book (from jupyterbook/)
```

Book output: `jupyterbook/_build/html/`

---
## 🎯 Design Principles

1. **Each concept introduced once**, in order needed (no forward references)
2. **One abstraction + one failure mode + one repair** per chapter
3. **Thermodynamics first**: All kinetics derive from Gibbs energy minimization under constraints

---
## 📞 Get Involved

- **Contributing:** Open issues for chapter requests or library features
