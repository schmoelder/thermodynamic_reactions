#!/usr/bin/env python3
"""Format Python code cells in MyST Markdown files using ruff format.

Usage:
    python format_cells.py jupyterbook/          # format all .md files
    python format_cells.py path/to/file.md       # format one file
    python format_cells.py --check jupyterbook/  # check without modifying
"""

import re
import subprocess
import sys
from pathlib import Path


def ruff_format(code: str) -> str | None:
    """Run ruff format on a code string. Returns formatted code, or None on error."""
    result = subprocess.run(
        ["ruff", "format", "--stdin-filename", "cell.py", "-"],
        input=code,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        print(f"ruff error:\n{result.stderr}", file=sys.stderr)
        return None
    return result.stdout


_CELL_RE = re.compile(
    r"(```\{code-cell\}\s+ipython3[^\n]*\n)"  # group 1: opening fence line
    r"(.*?)"  # group 2: cell content
    r"(^```[ \t]*$)",  # group 3: closing fence
    re.MULTILINE | re.DOTALL,
)


def split_preamble(content: str) -> tuple[str, str]:
    """Split cell content into (preamble, code).

    Preamble: leading MyST option lines (:key: value) plus one blank separator.
    Code: the Python source.
    """
    lines = content.splitlines(keepends=True)
    i = 0
    while i < len(lines) and lines[i].lstrip().startswith(":"):
        i += 1
    if i < len(lines) and lines[i].strip() == "":
        i += 1
    return "".join(lines[:i]), "".join(lines[i:])


def process_file(path: Path, check: bool = False) -> bool:
    """Return True if the file was changed (or would change in check mode)."""
    original = path.read_text(encoding="utf-8")
    output_parts: list[str] = []
    pos = 0
    changed = False

    for m in _CELL_RE.finditer(original):
        output_parts.append(original[pos : m.start()])
        opening, content, closing = m.group(1), m.group(2), m.group(3)
        preamble, code = split_preamble(content)

        if code.strip():
            formatted = ruff_format(code)
            if formatted is not None and formatted != code:
                changed = True
                if not check:
                    content = preamble + formatted

        output_parts.append(opening + content + closing)
        pos = m.end()

    output_parts.append(original[pos:])
    new_text = "".join(output_parts)

    if changed and not check:
        path.write_text(new_text, encoding="utf-8")

    return changed


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(
        description="Format Python code cells in MyST Markdown files using ruff."
    )
    parser.add_argument("paths", nargs="+", type=Path, metavar="path")
    parser.add_argument(
        "--check",
        action="store_true",
        help="Check only; exit 1 if any file would change.",
    )
    args = parser.parse_args()

    any_changed = False
    for p in args.paths:
        files = sorted(p.rglob("*.md")) if p.is_dir() else [p]
        for f in files:
            if "_build" in f.parts:
                continue
            if process_file(f, check=args.check):
                verb = "Would reformat" if args.check else "Reformatted"
                print(f"{verb}: {f}")
                any_changed = True

    if args.check and any_changed:
        sys.exit(1)


if __name__ == "__main__":
    main()
