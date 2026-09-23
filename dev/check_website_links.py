#!/usr/bin/env python3
"""Report internal links in the rendered documentation books whose page or anchor is missing."""

import argparse
from pathlib import Path
import re
import sys

LINK = re.compile(r'href="(?![a-z][a-z0-9+.-]*:)([^"#]*\.html)?(?:#([^"]*))?"')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("site", type=Path, nargs="?", default=Path("docs"), help="rendered site directory")
    site = parser.parse_args().site.resolve()
    ids = {}

    def page_ids(page):
        if page not in ids:
            ids[page] = set(re.findall(r'id="([^"]+)"', page.read_text())) if page.is_file() else None
        return ids[page]

    broken = set()
    for page in sorted(site.rglob("*.html")):
        for target, anchor in LINK.findall(page.read_text()):
            destination = (page.parent / target).resolve() if target else page
            found = page_ids(destination)
            # The parameter browser resolves its anchors in JavaScript.
            if found is None:
                broken.add((page.relative_to(site), f"{target}#{anchor}" if anchor else target, "missing page"))
            elif anchor and destination.name != "parameters.html" and anchor not in found:
                broken.add((page.relative_to(site), f"{target}#{anchor}", "missing anchor"))
    for page, link, problem in sorted(broken):
        print(f"{page}: {link} ({problem})")
    print(f"{len(broken)} broken internal links")
    return 1 if broken else 0


if __name__ == "__main__":
    sys.exit(main())
