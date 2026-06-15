---
name: multiomer-fix-errors
description: Inspect errors and interactively debug R targets workflows in multiomeR using list_distinct_errored_targets(), inspect_target_workspace(), and targeted tar_read_raw() probes. Use when a multiomeR run reports errored targets, when investigating a failing target, or when the user asks to debug multiomeR.
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

For target-local debugging, first inspect the failing target workspace without
assigning large objects into the session:

```r
inspect_target_workspace("<full_target_name>")
```

Use the dependency summaries to decide whether the issue is already present in
the target inputs. Zero-dimensional matrices, empty lists, absent design
columns, empty data frames, or obviously wrong configured labels usually point
upstream. If the inputs look plausible, re-run the failing helper expression
directly with smaller objects or `head()` subsets when possible.

For deeper probes, read upstream objects by the full target or dynamic branch
name:

```r
targets::tar_read_raw("<target_or_branch_name>")
```

Avoid wrapping `targets::tar_workspace(name)` directly when `name` is a
character variable; `tar_workspace()` uses non-standard evaluation and will look
for a target literally called `name`.

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
