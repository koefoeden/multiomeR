---
name: multiomer-develop-feature
description: Guide for adding new targets or features to the R targets workflows in multiomeR. Use when the user asks to add a new analysis step, new target, new helper function, or extend existing multiomeR functionality.
---

# multiomeR Develop Feature

Use this when adding or changing targets, helpers, or analysis steps.

## Core Workflow

1. Read the existing target group and helper code first. Preserve local naming, data shapes, and return objects unless there is a clear reason to change them.
2. Only write helpers in `R/<domain>_helpers.R` when the logic is reusable, or if the target command otherwise becomes too complex or difficult to understand
3. Add targets in the nearest existing `extra_targets/*.R` or `module_*/*.R` group - otherwise use `_targets.R`.
4. Run R only through the pixi environment; see `multiomer-run-r-code`.
5. Validate narrowly first with a helper-level smoke test or a single target, then broaden to relevant test datasets if the change affects target dependencies.

## Current Design Preferences

- Prefer BPCells, matrices, `GRanges`, tibbles, and `SummarizedExperiment` over Seurat/Signac objects for new ATAC/GEX processing.
- Keep Seurat compatibility only at explicit legacy boundaries. Do not create intermediate Seurat objects just to pass data between targets.
- Use direct package APIs where possible, e.g. `motifmatchr::matchMotifs()` on `GRanges` rather than Signac wrappers.
- Defer expensive compatibility artifacts until the consumer needs them. Example: compute `chromVAR::getBackgroundPeaks()` in gchromVAR targets, not in the upstream betterChromVAR target.
- Keep target names specific and modality-suffixed where relevant, e.g. `peak_dgCmatrix.ATAC`, `metadata_w_cell_types_tibble.ATAC`.
- Add aggregation-level downstream analyses as optional modules when they consume aggregation outputs. Use `module_<module>/`, opt in with the `modules` field in `cfg_aggregations.yaml`, and put module config in `module_<module>/cfg.yaml`.

## Target Patterns

```r
tar_target(
  name = result_name.ATAC,
  description = "Short user-facing description",
  command = helper_function(upstream_target, cfg_value),
  resources = get_tar_resources(cores_req = 15, RAM_GB_req = 60)
)
```

- Use `tarchetypes::tar_file()` for paths written to disk.
- Use `tarchetypes::tar_map()` for per-dataset or per-parameter expansion.
- Use `packages = w_def("PackageName")` only when a target needs to extend the default package set.
- Avoid new cross-store reads for work that can live as a root aggregation module. Extracted workflows such as GWAS processing now live in sister repos.

## Validation

- Use `multiomer-validation-workflow` after changing helpers, target files, `_targets.R`, target names, or mapping tibbles.
- If a target errors, use `multiomer-fix-errors`.

## Dependency Changes

- Add new R dependencies to `DESCRIPTION`; add `pixi.toml` dependencies only when they are available through conda/bioconda and needed for the runtime environment.
- If a package must come from GitHub/Bioconductor source, document it with `Remotes:` in `DESCRIPTION`.
