#!/usr/bin/env python3
"""List unique Mermaid node-label regexes for manual bypassing.

The output is sorted alphabetically by node label. Each line is an exact-match
regex pattern for the displayed Mermaid node label. Patterns are commented out
by default with "#". Uncomment or edit lines to choose which nodes
mermaid_bypass_patterns.py should bypass.

By default, an input file named example.mmd writes:

    example_bypass_patterns.txt
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import Counter
from pathlib import Path

from mermaid_bypass_patterns import node_label, parse_edges


def default_output_path(input_path: Path) -> Path:
    return input_path.with_name(f"{input_path.stem}_bypass_patterns.txt")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mermaid_file", type=Path, help="Input Mermaid .mmd file.")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Output file. Defaults to <input_stem>_bypass_patterns.txt.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_path = args.output or default_output_path(args.mermaid_file)
    if output_path.exists():
        raise SystemExit(f"refusing to overwrite existing file: {output_path}")

    lines = args.mermaid_file.read_text().splitlines(keepends=True)
    _, node_exprs = parse_edges(lines)
    if not node_exprs:
        raise SystemExit(f"no Mermaid graph nodes found in {args.mermaid_file}")

    labels = {node_id: node_label(expr) for node_id, expr in node_exprs.items()}
    duplicate_labels = sorted(
        label for label, count in Counter(labels.values()).items() if count > 1
    )
    if duplicate_labels:
        preview = ", ".join(duplicate_labels[:5])
        raise SystemExit(
            "node labels are not unique, so exact label regexes cannot identify "
            f"individual nodes: {preview}"
        )

    output_path.write_text(
        "\n".join(
            f"# ^{re.escape(label)}$" for label in sorted(labels.values(), key=str.lower)
        )
        + "\n"
    )

    print(
        f"listed Mermaid bypass patterns: output={output_path}, nodes={len(node_exprs)}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
