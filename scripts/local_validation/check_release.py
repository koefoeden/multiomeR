#!/usr/bin/env python3
"""Fail unless both full-demo statuses passed for an exact release commit."""
import argparse
import json
import re
import subprocess

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("sha")
args = parser.parse_args()
if not re.fullmatch(r"[0-9a-f]{40}", args.sha):
    parser.error("Pass the exact 40-character commit SHA")
pages = json.loads(subprocess.check_output([
    "gh", "api", "--paginate", "--slurp",
    f"repos/koefoeden/multiomeR/commits/{args.sha}/statuses?per_page=100"
], text=True))
latest = {}
for page in pages:
    for status in page:
        latest.setdefault(status["context"], status["state"])
for context in ("full-demo/human", "full-demo/mouse", "local-validation/differential-analysis"):
    if latest.get(context) != "success":
        raise SystemExit(f"{context}: {latest.get(context, 'missing')}")
print(f"Full-demo validation passed for {args.sha}")
