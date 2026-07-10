---
name: multiomer-fix-errors
description: Inspect errors and interactively debug R targets workflows in multiomeR using list_distinct_errored_targets(), inspect_target_workspace(), and targeted tar_read_raw() probes. Use when a multiomeR run reports errored targets, when investigating a failing target, or when the user asks to debug multiomeR.
---

# multiomeR Fix Errors

Determine the failing target and configured scope first. The project uses the
root `_targets.R` graph and the configured `outputs` store. Run R through Pixi;
see `multiomer-run-r-code`.

## Triage

Start with the narrowest stored error query that fits the request:

```r
list_errored_targets(target_name_pattern = "<target_or_dataset>")
```

Use the distinct summary only when several targets repeat the same failure:

```r
list_distinct_errored_targets()
```

Add commands and tracebacks when the root cause is not obvious. This helper, not
`list_distinct_errored_targets()`, accepts `target_name_pattern`:

```r
list_distinct_errored_targets_w_tracebacks(target_name_pattern = "<target_or_dataset>")
```

Use the full target name, including mapped or dynamic suffixes. If an external
command error is truncated, inspect the saved workspace and rerun only the
failing helper expression to recover complete output.

## Inspect The Workspace

Inspect the saved workspace when the stored error and traceback do not identify
the cause:

```r
inspect_target_workspace("<full_target_name>")
```

Use the compact dependency summaries before reading large inputs. Probe upstream
objects only when their dimensions, names, or configuration look suspicious.

For deeper probes, read upstream objects by the full target or dynamic branch
name:

```r
targets::tar_read_raw("<target_or_branch_name>")
```

Do not wrap `targets::tar_workspace(name)` when `name` is a character variable;
its non-standard evaluation looks for a target literally called `name`.

## Fixing

- Prefer fixing the invariant in helper or target code over adding broad guards.
- Use explicit checks with informative errors for missing columns, mismatched names, empty sets, or unsupported config.
- Avoid `tryCatch()` unless the failure is expected and recoverable; document why recovery is valid.
- Preserve current data contracts unless the task changes them deliberately.

## Verify

Re-run as narrowly as possible:

```bash
pixi run Rscript - <<'EOF'
selection <- targets::tar_manifest(
  names = tidyselect::matches("<target-and-scope-pattern>"),
  fields = c(name, description),
  callr_function = NULL
)
stopifnot(nrow(selection) > 0L)
print(selection)
targets::tar_make(
  names = tidyselect::matches("<target-and-scope-pattern>")
)
EOF
```

Use `multiomer-validation-workflow` for parse and graph checks. Exercise a second
configured scope only when the changed contract crosses species, modalities, or
configuration shapes and the first run cannot establish that behavior.
