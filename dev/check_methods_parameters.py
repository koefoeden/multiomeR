#!/usr/bin/env python3
"""Check the methods chapters of the implementation book against the repository.

Rules enforced (see website/implementation/implementation_conventions.md,
"Methods and parameter tables"):

1. Every ``parameters.html#<name>`` anchor in ``website/implementation/methods_*.md``
   names a ``param_name`` in the public-defaults manifest snapshot.
2. Status is ``Configurable`` or ``Fixed``. Configurable rows leave Value empty and cite
   manifest anchors or ``cfg_GEM_wells.tsv`` columns; Fixed rows list a value.
3. The Source cell of a Fixed row contains only backticked symbols, each of which resolves:
   ``f()`` is defined in project code, ``pkg::f()`` is called in project code, a path exists,
   and any other name occurs in a target file.
4. Every manifest parameter is linked from at least one methods chapter (warning only).
"""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

ANCHOR_RE = re.compile(r"parameters\.html#([A-Za-z0-9_]+)")
SYMBOL_RE = re.compile(r"`([^`]+)`")


def read_manifest(path: Path) -> set[str]:
    with path.open(newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        return {(value or "").strip() for row in rows for key, value in row.items() if key and key.strip() == "param_name"} - {""}


def read_tsv_header(path: Path) -> set[str]:
    with path.open() as handle:
        return {name.strip() for name in handle.readline().split("\t")}


def read_code(root: Path, patterns: list[str]) -> str:
    return "\n".join(path.read_text() for pattern in patterns for path in sorted(root.glob(pattern)))


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


def unresolved_symbols(source: str, root: Path, project_code: str, target_code: str) -> list[str]:
    if not SYMBOL_RE.search(source) or SYMBOL_RE.sub("", source).strip(" ,"):
        return [f"source must contain only backticked symbols: {source!r}"]
    problems = []
    for symbol in SYMBOL_RE.findall(source):
        if symbol.endswith("()") and "::" in symbol:
            found = re.search(re.escape(symbol[:-2]) + r"\(", project_code)
        elif symbol.endswith("()"):
            found = re.search(r"\b" + re.escape(symbol[:-2]) + r"\s*(<-|=)\s*function\b", project_code)
        elif "/" in symbol:
            found = (root / symbol).exists()
        else:
            found = re.search(r"(?<![\w.])" + re.escape(symbol) + r"(?![\w.])", target_code)
        if not found:
            problems.append(f"unresolved source symbol `{symbol}`")
    return problems


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    manifest = read_manifest(root / "website" / "data" / "public_defaults" / "cfg_pipeline_parameters.tsv")
    GEM_well_columns = read_tsv_header(root / "configuration" / "cfg_GEM_wells.tsv")
    project_code = read_code(root, ["R/*.R", "packages/multiomeRCore/R/*.R", "extra_targets/*.R", "module_*/*.R", "_targets.R"])
    target_code = read_code(root, ["extra_targets/*.R", "module_*/*.R", "_targets.R"])
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
            linked.add(match.group(1))
            if match.group(1) not in manifest:
                line = text.count("\n", 0, match.start()) + 1
                errors.append(f"{rel}:{line}: unknown manifest parameter '{match.group(1)}'")
        for line, row in table_rows(text):
            status, value, source = row.get("Status", ""), row.get("Value", ""), row.get("Source", "")
            if status == "Configurable":
                if value:
                    errors.append(f"{rel}:{line}: configurable row must not list a value ({value!r})")
                if not ANCHOR_RE.search(source) and not set(SYMBOL_RE.findall(source)) & GEM_well_columns:
                    errors.append(f"{rel}:{line}: configurable row must cite a manifest anchor or GEM-well column")
            elif status == "Fixed":
                if not value:
                    errors.append(f"{rel}:{line}: fixed row must list its value")
                errors += [f"{rel}:{line}: {problem}" for problem in unresolved_symbols(source, root, project_code, target_code)]
            elif status:
                errors.append(f"{rel}:{line}: unknown status {status!r}")

    for name in sorted(manifest - linked):
        print(f"warning: manifest parameter '{name}' is not linked from any methods chapter")
    for error in errors:
        print(error, file=sys.stderr)
    print(f"Checked {len(chapters)} chapters, {len(linked)} linked parameters, {len(errors)} errors.")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
