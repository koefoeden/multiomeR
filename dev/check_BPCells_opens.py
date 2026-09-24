#!/usr/bin/env python3
"""Check that targets open BPCells directories through open_BPCells_dir().

A target whose command is a bare BPCells::open_matrix_dir() or
open_fragments_dir() returns an object that hashes identically after its
directory changes, so its dependents would silently keep stale results.
"""

import re
import sys
from pathlib import Path

BARE_OPEN_RE = re.compile(r"command\s*=\s*BPCells::open_(?:matrix|fragments)_dir\(")


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    files = [root / "_targets.R", *root.glob("extra_targets/*.R"), *root.glob("module_*/**/*.R")]
    errors = [
        f"{path.relative_to(root)}:{number}: use open_BPCells_dir() instead of a bare BPCells open"
        for path in files if path.exists()
        for number, line in enumerate(path.read_text().splitlines(), start=1)
        if BARE_OPEN_RE.search(line)
    ]
    print("\n".join(errors) if errors else f"Checked {len(files)} target files, 0 bare BPCells opens.")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
