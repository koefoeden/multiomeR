#!/usr/bin/env python3
"""Check repository source links in the documentation sources.

Every https://github.com/koefoeden/multiomeR/{blob,tree}/main/<path> link in
website/**/*.md must name an existing path, and the methods fragments in
website/implementation/_shared_methods/, which the manuscript supplement also
includes, must contain no Markdown links.
"""

import re
import sys
from pathlib import Path

SOURCE_LINK_RE = re.compile(r"https://github\.com/koefoeden/multiomeR/(?:blob|tree)/main/([^)\s#]+)")
MARKDOWN_LINK_RE = re.compile(r"\]\(")


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    errors = []
    pages = [path for path in sorted((root / "website").rglob("*.md")) if path.name != "multiomeR-manual-llm.md"]
    for page in pages:
        text = page.read_text()
        for match in SOURCE_LINK_RE.finditer(text):
            if not (root / match.group(1)).exists():
                line = text.count("\n", 0, match.start()) + 1
                errors.append(f"{page.relative_to(root)}:{line}: missing source path {match.group(1)}")
    for fragment in sorted((root / "website" / "implementation" / "_shared_methods").glob("*.md")):
        if MARKDOWN_LINK_RE.search(fragment.read_text()):
            errors.append(f"{fragment.relative_to(root)}: shared methods fragments must not contain links")
    for error in errors:
        print(error, file=sys.stderr)
    print(f"Checked {len(pages)} pages, {len(errors)} errors.")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
