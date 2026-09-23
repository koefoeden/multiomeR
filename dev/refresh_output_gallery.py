#!/usr/bin/env python3
"""Refresh the output gallery from the saved plots of one built aggregation.

Every plot target of the aggregation, plus the GEM-well plot targets of its
first GEM well, contributes one WebP preview to website/gallery_assets/ and one
entry to website/output_gallery.yaml. An entry keeps its source_file while that
file still exists, so a hand-picked example survives a refresh.
"""

from pathlib import Path
import argparse
import datetime
import json
import re
import shutil
import subprocess
import tempfile

import yaml
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "website" / "output_gallery.yaml"
ASSETS = ROOT / "website" / "gallery_assets"

# Checkpoint tag -> (page group, section heading within the group).
SECTIONS = {
    "1_pre-aggregation-QC": ("Primary module", "1. Pre-aggregation QC"),
    "2_GEX-PCA-QC": ("Primary module", "2. GEX dimension reduction"),
    "3_GEX-QC": ("Primary module", "3. GEX clusters and cell types"),
    "4_peak-QC": ("Primary module", "4. Peak QC"),
    "5_pre-LSI-QC": ("Primary module", "5. ATAC filtering"),
    "6_ATAC-LSI-QC": ("Primary module", "6. ATAC dimension reduction"),
    "7_ATAC-QC": ("Primary module", "7. ATAC clusters and motifs"),
    "8_multimodal-QC": ("Primary module", "8. WNN integration"),
    "differential_analyses": ("Differential analyses", None),
    "genetic_enrichment": ("Genetic enrichment", None),
    "peak_gene_correlation": ("Peak–gene correlation", None),
}
# File-name patterns, in order of preference, for targets that save several plots.
PREFERRED_FILES = [r"cluster_cell_type", r"cell_type", r"_named", r"nCount_RNA", r"^top_", r"^01_"]
MAX_PIXELS = 2_000_000
MAX_SIDE = 2600

R_QUERY = r"""
args <- commandArgs(trailingOnly = TRUE)
aggregation <- args[[1]]
config <- read_aggregation_config_tibble()
first_well <- config$aggregation_GEM_well_IDs[[match(aggregation, config$aggregation)]][[1]]
scopes <- c(aggregation, first_well)
manifest <- targets::tar_manifest(fields = c(name, description), callr_function = NULL)
meta <- targets::tar_meta(fields = c(name, path, type, children), complete_only = FALSE)
plots_root <- normalizePath(file.path(targets::tar_config_get("store"), "plots"))
records <- list()
for (i in seq_len(nrow(manifest))) {
  name <- manifest$name[[i]]
  scope <- scopes[endsWith(name, paste0(".", scopes))][1]
  row <- meta[meta$name == name, ]
  if (is.na(scope) || !nrow(row)) next
  paths <- if (identical(row$type, "pattern")) {
    unlist(meta$path[meta$name %in% unlist(row$children)])
  } else unlist(row$path)
  paths <- normalizePath(paths[grepl("\\.png$", paths)], mustWork = FALSE)
  paths <- sort(paths[startsWith(paths, file.path(plots_root, scope, "")) & file.exists(paths)])
  if (!length(paths)) next
  records[[length(records) + 1L]] <- list(
    target = sub(paste0("\\.", scope, "$"), "", name),
    checkpoint = sub(".*\\[checkpoint:([^]]+)\\].*", "\\1", manifest$description[[i]]),
    description = trimws(gsub("\\s*\\[[^]]*\\]", "", manifest$description[[i]])),
    plots_dir = file.path(plots_root, scope),
    files = I(substring(paths, nchar(file.path(plots_root, scope, "")) + 1L))
  )
}
jsonlite::write_json(records, args[[2]], auto_unbox = TRUE)
"""


def query_plot_targets(aggregation):
    with tempfile.TemporaryDirectory() as directory:
        script, output = Path(directory) / "query.R", Path(directory) / "targets.json"
        script.write_text(R_QUERY)
        subprocess.run(["Rscript", str(script), aggregation, str(output)], cwd=ROOT, check=True)
        return json.loads(output.read_text())


def choose_file(files, previous):
    if previous in files:
        return previous
    for pattern in PREFERRED_FILES:
        matches = [file for file in files if re.search(pattern, Path(file).name)]
        if matches:
            return matches[0]
    return files[0]


def write_preview(source, destination):
    image = Image.open(source)
    scale = min(1.0, (MAX_PIXELS / (image.width * image.height)) ** 0.5, MAX_SIDE / max(image.size))
    size = (max(1, round(image.width * scale)), max(1, round(image.height * scale)))
    destination.parent.mkdir(parents=True, exist_ok=True)
    image.convert("RGB").resize(size, Image.LANCZOS).save(destination, "WEBP", quality=80, method=6)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("aggregation", help="built aggregation whose plots populate the gallery")
    aggregation = parser.parse_args().aggregation
    Image.MAX_IMAGE_PIXELS = None

    previous = {}
    if MANIFEST.exists():
        previous = {item["target"]: item.get("source_file") for item in yaml.safe_load(MANIFEST.read_text())["items"]}
    records = query_plot_targets(aggregation)
    unknown = sorted({record["checkpoint"] for record in records} - SECTIONS.keys())
    if unknown:
        raise SystemExit(f"Add gallery sections for checkpoints: {', '.join(unknown)}")

    order = list(SECTIONS)
    records.sort(key=lambda record: order.index(record["checkpoint"]))  # stable: keeps pipeline order
    shutil.rmtree(ASSETS, ignore_errors=True)
    items = []
    for record in records:
        target, checkpoint = record["target"], record["checkpoint"]
        source_file = choose_file(record["files"], previous.get(target))
        asset = Path("gallery_assets") / checkpoint / f"{target}.webp"
        write_preview(Path(record["plots_dir"]) / source_file, ROOT / "website" / asset)
        group, section = SECTIONS[checkpoint]
        items.append({
            "checkpoint": checkpoint,
            "group": group,
            "section": section,
            "target": target,
            "description": record["description"],
            "source_file": source_file,
            "asset": asset.as_posix(),
        })

    manifest = {
        "source": {"aggregation": aggregation, "refreshed": datetime.date.today().isoformat()},
        "items": items,
    }
    MANIFEST.write_text(yaml.safe_dump(manifest, sort_keys=False, allow_unicode=True, width=1000))
    print(f"Wrote {len(items)} previews from {aggregation}.")


if __name__ == "__main__":
    main()
