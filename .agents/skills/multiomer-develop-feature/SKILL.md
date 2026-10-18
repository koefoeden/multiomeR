---
name: multiomer-develop-feature
description: Add or change targets, R helpers, configuration parameters, or analysis steps in the active multiomeR root workflow. Use for new pipeline functionality and behavior changes to existing target families or optional modules.
---

# multiomeR Develop Feature

Use this when adding or changing targets, helpers, or analysis steps.

## Core Workflow

1. Inspect the owning target fragment, helper, and configuration contract.
   [Change the pipeline](../../../website/change_the_pipeline.md) documents the
   contracts to preserve: the parameter manifest, mapping tibbles, target names
   and description tags.
2. Place helpers and target fragments as the code-change list in
   `.github/CONTRIBUTING.md` describes; reserve `_targets.R` for graph
   composition and mapping.
3. Add or revise YAML parameters in `manifests/cfg_pipeline_parameters.tsv` before using
   them in config readers or targets. A changed public parameter also needs the
   documentation snapshot that `multiomer-edit-documentation` describes.
4. Run R through Pixi with `multiomer-run-r-code` and validate with
   `multiomer-validation-workflow`. If a target errors, use `multiomer-fix-errors`.

## Design Rules

- For new or revised plots, follow [multiomer-annotate-plots](../multiomer-annotate-plots/SKILL.md): title identifies the subject, subtitle guides interpretation, and caption explains relevant methods.
- Compute per-cell ATAC results from fragments by mapping over the per-GEM-well
  fragment directories, not `combined_BPCells_fragment_obj.ATAC`: BPCells
  filters cells only after merging every well, so even a one-cell selection on
  the combined object takes minutes.
- Add helpers only for reused or otherwise unreadable logic, and targets only
  for meaningful computation, output, or cache/invalidation boundaries—not to
  rename, forward, or organize other targets.
- Add aggregation-level downstream analyses as optional modules when they
  consume aggregation outputs: an `extra_targets/module_<module>/` fragment
  maps over `build_module_tibble()`, which reads the module's
  `cfg_module_<module>.yaml` only for aggregations that list the module under
  `modules` in `cfg_aggregations.yaml`.
- Keep a download in its own target that stores only the rows or fields the
  analysis needs, apart from their parsing, so that a parsing change does not
  download again. Run downloads one at a time with `deployment = "main"` when
  a server refuses many parallel connections.
- Tag every user-facing output of a checkpoint or module, file exports as well
  as plots, with its `[checkpoint:<name>]` so documented selections build it.
- Use `packages = w_def("PackageName")` only when a target needs to extend the default package set.
- Add `get_tar_resources()` from measured needs or the closest comparable
  target; do not copy a large CPU/RAM request as a generic default. Request
  about 1.5 times the measured peak RAM; the request selects the smallest
  controller in the selected `crew_controllers.R` that fits. Keep the default
  path within 32 GB, and record the measurement as a
  `[resource_observation: ...]` description tag.
- Avoid new cross-store reads for work that can live as a root aggregation module.

## Dependency Changes

The checkout's Pixi environment is shared by other sessions and by the workers
of a running pipeline, also on other nodes. Replacing a package that another R
process has loaded breaks its lazy-load database there, so install or update
packages only while no driver (see `multiomer-run-pipeline`) and none of its
scheduler workers are running.

- Add conda-forge or Bioconda dependencies to `pixi.toml` using the existing
  constraint style and refresh `pixi.lock`.
- Add GitHub-only R packages to `scripts/github_packages.R` with a
  pinned commit and add any required system or indirect dependencies to
  `pixi.toml`.
- Do not create a root `DESCRIPTION`; the workflow itself is not maintained as
  an R package.
