---
name: multiomer-validation-workflow
description: Efficient validation workflow for multiomeR R helper and targets edits. Use after changing _targets.R, extra_targets/*.R, module target files, R helpers, target names, tar_map mapping tibbles, or target graph construction.
---

# multiomeR Validation Workflow

Run R through Pixi. Use this default path after code edits and add target
execution only when static graph construction cannot prove the behavior.

## Fast Default

Validate each edited file once through the strongest applicable check. Focused
source or behavior checks already parse helpers and scripts; manifest
construction parses and sources target graph files. Use standalone `parse()`
only when an edited file is not otherwise loaded or when isolating a syntax
failure.

For target changes, build the manifest without callr. This catches missing helpers, bad target commands, malformed `tar_map()` values, and most graph-construction failures without running targets:

```bash
pixi run --use-environment-activation-cache Rscript - <<'EOF'
targets::tar_manifest(callr_function = NULL)
cat("manifest ok\n")
EOF
```

The full manifest takes about a minute, so build it once per batch of edits and
save it with `saveRDS()` when several queries need it.

For configuration-directory changes, compare resolved settings and manifest
names before/after, check the default and selected directories independently,
and exercise invalid selections and disabled modules without changing the live
selector or store: give a worktree its own `configuration.local`, or pass
`configuration_path()` the `project_root` of a scratch copy. Shared
documentation examples use `configuration/`.

In a worktree that runs targets on a store snapshot, the GEM-well table
resolves to the worktree's own path, so `GEM_well_config_tsv` looks outdated
and runs that reach aggregation targets are refused. Run
`GEM_well_metadata_tibble.<aggregation>` first; it rebuilds to identical
hashes.

Finish code edits with:

```bash
git diff --check
```

## Target Execution

For checkpoint-tag changes, validate the selected targets' transitive
dependencies: a description tag does not prevent a downstream dependency from
crossing the intended review boundary. Keep `manifests/QC_checkpoint_manifest.tsv` and the
checkpoint sections of `website/main_running.md` consistent.

Run targets only when the changed behavior needs runtime proof, using the
preview-and-run patterns in `multiomer-run-pipeline`. Run the smallest target
that exercises the changed code, and avoid setup/download targets unless the
edit directly changed download/setup behavior.

To show that a performance change reproduces its outputs, call the changed
helper on a stored branch's inputs and compare the result with that branch's
stored value using `all.equal()`. `tar_read(branches = i)` of an
`iteration = "group"` input returns every group, so filter `tar_group == i`.
For a bit-identical comparison, set the BLAS threads the worker had
(`RhpcBLASctl::blas_set_num_threads()`): OpenBLAS sums matrix products in
another order with another thread count, or when a matrix gains columns.

WNN results cannot be reproduced by rebuilding: `BPCells::knn_hnsw()` with more
than one thread is nondeterministic. Validate code after the neighbour search on
a shared kNN input, and expect an adopted WNN change to cascade.

For a compiled kernel in `src/`, run its reference-parity test listed in
`website/algorithm_validation.md`, compare it with the committed
kernel on a stored branch's inputs, and try degenerate inputs such as collinear
or constant columns: a kernel error fails the whole branch.

## Invalidation Impact

Use `multiomer-estimate-rebuild-impact` to find which existing targets an edit
reruns and to confirm a commit's impact line.

## Known Waste

Do not use `tar_validate()` while the project sets `error = "trim"`: the
installed `targets` version rejects that option before useful validation.
Prefer manifest construction and targeted runtime checks.
