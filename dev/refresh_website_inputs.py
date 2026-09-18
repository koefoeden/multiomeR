#!/usr/bin/env python3
"""Refresh shared documentation inputs from an explicitly selected public checkout."""

import argparse
from pathlib import Path
import shutil


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--public-repo", type=Path, required=True)
    args = parser.parse_args()
    source = args.public_repo.resolve()
    destination = Path(__file__).resolve().parents[1] / "website/data/public_defaults"
    destination.mkdir(parents=True, exist_ok=True)
    files = [source / "cfg_pipeline_parameters.tsv"]
    files += [source / "configuration" / name for name in (
        "cfg_aggregations.yaml", "cfg_module_differential_analyses.yaml",
        "cfg_module_genetic_enrichment.yaml", "cfg_module_peak_gene_correlation.yaml",
    )]
    for path in files:
        shutil.copyfile(path, destination / path.name)


if __name__ == "__main__":
    main()
