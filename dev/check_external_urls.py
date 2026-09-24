#!/usr/bin/env python3
"""Check that the example-data download URLs still respond.

The download script reads example_data/*manifest.tsv directly, so these are
exactly the URLs the public examples download. Each request asks for one byte.
"""

import csv
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


def probe(url: str) -> str | None:
    request = urllib.request.Request(url, headers={"Range": "bytes=0-0", "User-Agent": "multiomeR-url-check"})
    error = None
    for attempt in range(3):
        try:
            with urllib.request.urlopen(request, timeout=30) as response:
                if response.status in (200, 206):
                    return None
                error = f"HTTP {response.status}"
        except urllib.error.HTTPError as exception:
            error = f"HTTP {exception.code}"
        except (urllib.error.URLError, OSError) as exception:
            error = str(getattr(exception, "reason", exception))
        time.sleep(2**attempt)
    return error


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    urls = sorted({
        row["url"]
        for manifest in sorted((root / "example_data").glob("*manifest.tsv"))
        for row in csv.DictReader(manifest.open(), delimiter="\t")
    })
    with ThreadPoolExecutor(max_workers=8) as pool:
        failures = [(url, error) for url, error in zip(urls, pool.map(probe, urls)) if error]
    for url, error in failures:
        print(f"{url}: {error}")
    print(f"Checked {len(urls)} URLs, {len(failures)} failed.")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
