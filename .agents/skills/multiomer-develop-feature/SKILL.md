---
name: multiomer-develop-feature
description: Add or change targets, R helpers, configuration parameters, or analysis steps in the active multiomeR root workflow. Use for new pipeline functionality and behavior changes to existing target families or optional modules.
---

# multiomeR Develop Feature

Use this when adding or changing targets, helpers, or analysis steps.

## Core Workflow

1. Inspect the owning target fragment, helper, and configuration contract.
2. Put reusable or otherwise unreadable target logic in `R/<domain>_helpers.R`.
3. Add targets to the nearest `extra_targets/*.R` or `module_*/*.R` fragment;
   reserve `_targets.R` for graph composition and mapping.
4. Add or revise YAML parameters in `cfg_pipeline_parameters.tsv` before using
   them in config readers or targets.
5. Run R through Pixi with `multiomer-run-r-code` and validate with
   `multiomer-validation-workflow`.

## Current Design Preferences

- Prefer BPCells, matrices, `GRanges`, tibbles, and `SummarizedExperiment` over Seurat/Signac objects for new ATAC/GEX processing.
- Keep Seurat compatibility only at explicit legacy boundaries. Do not create intermediate Seurat objects just to pass data between targets.
- Use direct package APIs rather than compatibility wrappers when practical.
- Prefer direct target commands and existing machinery. Add helpers only for
  reused or otherwise unreadable logic, and targets only for meaningful
  computation, output, or cache/invalidation boundaries—not to rename, forward,
  or organize other targets.
- Defer expensive compatibility artifacts until their consumer needs them.
- Keep target names specific and modality-suffixed where relevant, e.g.
  `consensus_peak_BPCells_matrix.ATAC` and
  `metadata_w_cell_types_tibble.ATAC`.
- Add aggregation-level downstream analyses as optional modules when they consume aggregation outputs. Use `module_<module>/`, opt in with the `modules` field in `cfg_aggregations.yaml`, and put module config in `module_<module>/cfg.yaml`.

## Target Patterns

```r
tar_target(
  name = result_name.ATAC,
  description = "Short user-facing description [part_of_graph:<graph_id>]",
  command = helper_function(upstream_target, cfg_value)
)
```

- Use `tarchetypes::tar_file()` for paths written to disk.
- Use `tarchetypes::tar_map()` for per-dataset or per-parameter expansion.
- Always set `deployment = "main"` on `tarchetypes::tar_files()` so file
  discovery stays in the pipeline driver environment rather than a Crew worker.
- Use `packages = w_def("PackageName")` only when a target needs to extend the default package set.
- Add `get_tar_resources()` from measured needs or the closest comparable
  target; do not copy a large CPU/RAM request as a generic default.
- Avoid new cross-store reads for work that can live as a root aggregation module. Extracted workflows such as GWAS processing now live in sister repos.

## Validation

- Use `multiomer-validation-workflow` after changing helpers, target files, `_targets.R`, target names, or mapping tibbles.
- If a target errors, use `multiomer-fix-errors`.

## Dependency Changes

- Add conda-forge or Bioconda dependencies to `pixi.toml` using the existing
  constraint style and refresh `pixi.lock`.
- Add GitHub-only R packages to `scripts/install_r_github_packages.R` with a
  pinned commit and add any required system or indirect dependencies to
  `pixi.toml`.
- Do not recreate `DESCRIPTION`; this repository is not currently maintained as
  an R package.
