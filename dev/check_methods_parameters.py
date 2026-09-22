#!/usr/bin/env python3
"""Check the methods chapters of the implementation book against the parameter manifest.

Rules enforced (see website/implementation/implementation_conventions.md,
"Methods and parameter tables"):

1. Every ``parameters.html#<name>`` anchor in ``website/implementation/methods_*.md``
   names a ``param_name`` present in the public-defaults manifest snapshot
   ``website/data/public_defaults/cfg_pipeline_parameters.tsv``.
2. In every table row whose Status column is ``Configurable``, the Value column is empty,
   so defaults are never repeated outside the manifest.
3. In every table row whose Status column starts with ``Hardcoded``, the Value column is
   not empty.
4. Every manifest parameter is linked from at least one methods chapter (warning only).
"""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

ANCHOR_RE = re.compile(r"parameters\.html#([A-Za-z0-9_]+)")


def read_manifest(path: Path) -> set[str]:
    with path.open(newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        names = set()
        for row in rows:
            cleaned = {key.strip(): (value or "").strip() for key, value in row.items() if key}
            if cleaned.get("param_name"):
                names.add(cleaned["param_name"])
    return names


def table_rows(text: str):
    """Yield (line_number, cells) for Markdown table body rows."""
    header: list[str] | None = None
    for number, line in enumerate(text.splitlines(), start=1):
        stripped = line.strip()
        if not stripped.startswith("|"):
            header = None
            continue
        cells = [cell.strip() for cell in stripped.strip("|").split("|")]
        if header is None:
            header = cells
            continue
        if all(re.fullmatch(r":?-+:?", cell) for cell in cells):
            continue
        yield number, dict(zip(header, cells))


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    manifest = read_manifest(root / "website" / "data" / "public_defaults" / "cfg_pipeline_parameters.tsv")
    chapters = sorted((root / "website" / "implementation").glob("methods_*.md"))
    if not chapters:
        print("No methods chapters found.", file=sys.stderr)
        return 1

    errors: list[str] = []
    linked: set[str] = set()
    for chapter in chapters:
        rel = chapter.relative_to(root)
        text = chapter.read_text()
        for match in ANCHOR_RE.finditer(text):
            name = match.group(1)
            linked.add(name)
            if name not in manifest:
                line = text.count("\n", 0, match.start()) + 1
                errors.append(f"{rel}:{line}: unknown manifest parameter '{name}'")
        for line, row in table_rows(text):
            status = row.get("Status", "")
            value = row.get("Value", "")
            if status.startswith("Configurable") and value:
                errors.append(f"{rel}:{line}: configurable row must not list a value ({value!r})")
            if status.startswith("Hardcoded") and not value:
                errors.append(f"{rel}:{line}: hardcoded row must list its value")
            if status and not (status.startswith("Configurable") or status.startswith("Hardcoded")):
                errors.append(f"{rel}:{line}: unknown status {status!r}")

    for name in sorted(manifest - linked):
        print(f"warning: manifest parameter '{name}' is not linked from any methods chapter")

    for error in errors:
        print(error, file=sys.stderr)
    print(f"Checked {len(chapters)} chapters, {len(linked)} linked parameters, {len(errors)} errors.")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
