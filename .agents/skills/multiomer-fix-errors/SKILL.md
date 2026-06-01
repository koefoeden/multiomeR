---
name: multiomer-fix-errors
description: Inspect errors and interactively debug R targets workflows in multiomeR using list_distinct_errored_targets() and load_workspace(). Use when a multiomeR run reports errored targets, when investigating a failing target, or when the user asks to debug multiomeR.
---

# multiomeR Fix Errors

Always determine the target and dataset from context first. Current multiomeR work uses the root `_targets.R` project. Run R through pixi; see `multiomer-run-r-code`.

## Triage

Start with distinct errors:

```r
list_distinct_errored_targets()
```

Filter errored targets with `se()` when you only need matching target names and
stored error text:

```r
se(target_name_pattern = "<target_or_dataset>")
```

Add commands and tracebacks when the root cause is not obvious. This helper, not
`list_distinct_errored_targets()`, accepts `target_name_pattern`:

```r
list_distinct_errored_targets_w_tracebacks(target_name_pattern = "<target_or_dataset>")
```

Use the full target name from this output, including dataset or dynamic branch suffixes.
If the stored error is a truncated external-command failure, inspect the target
workspace and rerun the helper call directly to recover complete `stderr` /
`stdout`.

## Inspect The Workspace

For target-local debugging, load the failing target workspace:

```r
load_workspace("<full_target_name>")
str(object_or_dependency)
```

Then re-run the failing helper expression directly with smaller objects or `head()` subsets when possible. For upstream objects, use `targets::tar_read()` / `targets::tar_load()` by full target name.

## Fixing

- Prefer fixing the invariant in helper or target code over adding broad guards.
- Use explicit checks with informative errors for missing columns, mismatched names, empty sets, or unsupported config.
- Avoid `tryCatch()` unless the failure is expected and recoverable; document why recovery is valid.
- Preserve current data contracts unless the user asked for a redesign. For new ATAC work, prefer BPCells, matrices, `GRanges`, tibbles, and `SummarizedExperiment` over Seurat/Signac objects.

## Verify

Re-run as narrowly as possible:

```r
tar_make(names = matches("<target_name>.*<dataset>"))
```

If the fix changes shared helpers or target dependencies, also smoke-test a second small dataset. Parse edited files, check lints, and run `git diff --check`. Only commit when the user explicitly asks; use the multiomeR commit-format skill.
