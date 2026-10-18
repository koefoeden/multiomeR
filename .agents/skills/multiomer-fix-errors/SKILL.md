---
name: multiomer-fix-errors
description: Inspect errors and interactively debug R targets workflows in multiomeR using stored errors, saved workspaces, and targeted dependency probes. Use when a multiomeR run reports errored targets, when investigating a failing target, or when the user asks to debug multiomeR.
---

# multiomeR Fix Errors

Determine the failing target and configured scope first. The project uses the
root `_targets.R` graph and configured store. Run R through Pixi; see
`multiomer-run-r-code`.

## Triage

Start with the stored errors, which `tar_meta()` reads without building the
graph:

```r
list_errored_targets(target_name_pattern = "<target_or_dataset>")
```

Use the full target name, including mapped or dynamic suffixes. If an external
command error is truncated, inspect the saved workspace and rerun only the
failing helper expression to recover complete output.

The helper keeps only targets that errored in the store's latest run: it filters
by `tar_errored()`, which reads the progress file that every run start empties,
including another session's. If a run has started since, read
`tar_meta(fields = error)` and filter `!is.na(error)` instead; a stored error
stays until the target succeeds. A failed target with no stored error had its
worker crash; read the controller's worker log.

## Inspect The Workspace

Inspect the saved workspace when the stored error does not identify the cause:

```r
inspect_target_workspace("<full_target_name>")
```

It prints the stored error, traceback and command, and builds the full manifest
to find the command, which takes about a minute. Use the compact dependency
summaries before reading large inputs. Probe upstream objects only when their
dimensions, names, or configuration look suspicious. The workspace loads the
current upstream values, so an old workspace may no longer reproduce the error.
Use these helpers and exported `targets` functions rather than guessing at
`targets:::` internals.

For deeper probes, read upstream objects by the full target or dynamic branch
name:

```r
targets::tar_read_raw("<target_or_branch_name>")
```

## Fixing

- Prefer fixing the invariant in helper or target code over adding broad guards.
- Use explicit checks with informative errors for missing columns, mismatched names, empty sets, or unsupported config.
- Avoid `tryCatch()` unless the failure is expected and recoverable; document why recovery is valid.
- Preserve current data contracts unless the task changes them deliberately.

## Verify

Use `multiomer-validation-workflow`, then re-run the narrowest failing target
through `multiomer-run-pipeline`. Exercise a second configured scope only when
the changed contract crosses species, modalities, or configuration shapes and
the first run cannot establish that behavior.
