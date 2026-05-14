"""
Build API reference pages for the reactions library.

Parses reactions.api and reactions.solver with griffe and emits
MyST-compatible markdown into jupyterbook/api/.

Usage:
    python scripts/build_api_docs.py
"""

from __future__ import annotations

from pathlib import Path

import griffe
from griffe import DocstringSectionKind, Object

OUT_DIR = Path(__file__).parent.parent / "jupyterbook" / "api"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Groups: (page_stem, page_title, label, [(module, name), ...])
# Abstract bases and internal helpers are omitted intentionally.
# ---------------------------------------------------------------------------

GROUPS = [
    (
        "01_core",
        "Core Objects",
        "api-core",
        [
            ("reactions.api", "Species"),
            ("reactions.api", "Component"),
            ("reactions.api", "PhysicalState"),
            ("reactions.api", "ConservationReport"),
            ("reactions.api", "ReactionModel"),
        ],
    ),
    (
        "02_ionic_strength",
        "Ionic Strength Models",
        "api-ionic-strength",
        [
            ("reactions.api", "IonicStrengthIdeal"),
            ("reactions.api", "IonicStrengthBackground"),
            ("reactions.api", "IonicStrengthFixed"),
        ],
    ),
    (
        "03_activity",
        "Activity Coefficient Models",
        "api-activity",
        [
            ("reactions.api", "ActivityCoefficientIdeal"),
            ("reactions.api", "ActivityCoefficientDebyeHuckel"),
            ("reactions.api", "ActivityCoefficientDavies"),
            ("reactions.api", "ActivityCoefficientCustom"),
        ],
    ),
    (
        "04_equilibrium_constants",
        "Equilibrium Constants",
        "api-equilibrium-constants",
        [
            ("reactions.api", "EquilibriumConstant"),
            ("reactions.api", "EquilibriumConstantVantHoff"),
            ("reactions.api", "EquilibriumConstantVantHoffCp"),
            ("reactions.api", "EquilibriumConstantCustom"),
            ("reactions.api", "EquilibriumConstantTabulated"),
            ("reactions.api", "EquilibriumConstantPolynomial"),
            ("reactions.api", "pKa"),
        ],
    ),
    (
        "05_rate_constants",
        "Rate Constants",
        "api-rate-constants",
        [
            ("reactions.api", "RateConstantFixed"),
            ("reactions.api", "RateConstantArrhenius"),
            ("reactions.api", "RateConstantPolynomial"),
            ("reactions.api", "RateConstantTabulated"),
        ],
    ),
    (
        "06_rate_laws",
        "Saturation Rate Laws",
        "api-rate-laws",
        [
            ("reactions.api", "MichaelisMenten"),
            ("reactions.api", "HillRate"),
            ("reactions.api", "CustomRate"),
        ],
    ),
    (
        "07_reactions",
        "Reaction Types",
        "api-reactions",
        [
            ("reactions.api", "MassActionReaction"),
            ("reactions.api", "ThermodynamicReaction"),
            ("reactions.api", "EnzymaticReaction"),
            ("reactions.api", "CustomReaction"),
        ],
    ),
    (
        "08_solver",
        "Solver",
        "api-solver",
        [
            ("reactions.solver", "simulate"),
            ("reactions.solver", "solve_equilibrium"),
            ("reactions.solver", "SimulationResult"),
        ],
    ),
    (
        "09_formulation",
        "Formulation",
        "api-formulation",
        [
            ("reactions.api", "Solution"),
        ],
    ),
    (
        "10_common",
        "Common Buffers",
        "api-common",
        [
            ("reactions.common", "autoionization"),
            ("reactions.common", "acetic_acid_equilibria"),
            ("reactions.common", "phosphate_equilibria"),
            ("reactions.common", "citric_acid_equilibria"),
            ("reactions.common", "tris_equilibria"),
            ("reactions.common", "hepes_equilibria"),
            ("reactions.common", "mops_equilibria"),
        ],
    ),
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def annotation_str(annotation) -> str:
    if annotation is None:
        return ""
    return str(annotation)


def render_object(obj: Object) -> str:
    """Render a single class or function as MyST markdown."""
    lines = []

    is_func = obj.kind.name == "FUNCTION"

    # Heading
    lines.append(f"### `{obj.name}`")
    lines.append("")

    if not obj.docstring:
        lines.append("*No docstring.*")
        lines.append("")
        return "\n".join(lines)

    parsed = obj.docstring.parsed

    # For classes, parameters come from __init__ docstring if present,
    # otherwise fall back to class docstring.
    param_sections = []
    text_sections = []
    example_sections = []

    def collect_sections(sections) -> None:
        for section in sections:
            if section.kind == DocstringSectionKind.text:
                text_sections.append(section.value)
            elif section.kind == DocstringSectionKind.parameters:
                param_sections.extend(section.value)
            elif section.kind == DocstringSectionKind.examples:
                example_sections.append(section.value)

    collect_sections(parsed)

    # Also pull __init__ params for classes
    if not is_func and "__init__" in obj.members:
        init = obj.members["__init__"]
        if init.docstring:
            collect_sections(init.docstring.parsed)

    # Prose description
    for text in text_sections:
        lines.append(text.strip())
        lines.append("")

    # Parameters table
    if param_sections:
        lines.append("**Parameters**")
        lines.append("")
        lines.append("| Name | Type | Description |")
        lines.append("|------|------|-------------|")
        for p in param_sections:
            ann = annotation_str(p.annotation)
            desc = p.description.replace("\n", " ").strip()
            lines.append(f"| `{p.name}` | `{ann}` | {desc} |")
        lines.append("")

    # Examples
    for ex in example_sections:
        lines.append("**Example**")
        lines.append("")
        for kind, content in ex if isinstance(ex, list) else [ex]:
            if kind == DocstringSectionKind.text:
                # code block
                lines.append("```python")
                lines.append(content.strip())
                lines.append("```")
                lines.append("")

    return "\n".join(lines)


def render_page(title: str, label: str, objects: list[Object]) -> str:
    lines = [
        f"({label})=",
        f"# {title}",
        "",
    ]
    for obj in objects:
        lines.append(render_object(obj))
        lines.append("---")
        lines.append("")
    return "\n".join(lines)


def render_index(groups) -> str:
    lines = [
        "(api)=",
        "# API Reference",
        "",
        "Complete reference for all public classes and functions in the `reactions` library.",
        "",
        "**Modules in this reference:**",
        "",
    ]
    for _, _, label, _ in groups:
        lines.append(f"- @{label}")
    lines.append("")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    # Load modules once
    modules: dict[str, griffe.Module] = {}
    for _, _, _, members in GROUPS:
        for mod_name, _ in members:
            if mod_name not in modules:
                modules[mod_name] = griffe.load(mod_name, docstring_parser="numpy")

    # Write index
    (OUT_DIR / "00_index.md").write_text(render_index(GROUPS))
    print("wrote 00_index.md")

    # Write one page per group
    for stem, title, label, members in GROUPS:
        objects = []
        for mod_name, obj_name in members:
            mod = modules[mod_name]
            if obj_name in mod.members:
                objects.append(mod.members[obj_name])
            else:
                print(f"  WARNING: {obj_name} not found in {mod_name}")

        page = render_page(title, label, objects)
        out_path = OUT_DIR / f"{stem}.md"
        out_path.write_text(page)
        print(f"wrote {out_path.name}")

    print(f"\nDone. Output in {OUT_DIR}")


if __name__ == "__main__":
    main()
