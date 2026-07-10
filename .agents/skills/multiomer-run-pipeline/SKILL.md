---
name: multiomer-run-pipeline
description: Inspect or run the active multiomeR root targets workflow through Pixi. Use for explicit requests to run or rerun targets and for read-only target status, progress, outdatedness, or process-lock checks. Do not start target execution for a request that only asks for inspection.
---

# multiomeR Run Pipeline

The active project is `_targets.R`; `_targets.yaml` currently resolves the store
to `outputs`. Differential analyses and genetic enrichment are aggregation-level
modules in this same graph. Run R through Pixi; see `multiomer-run-r-code`.

## Read-only inspection

Use the live configuration instead of a hardcoded store. `tar_pid()` can retain
the PID of a completed run, so verify it with `ps` before reporting that a run is
active.

```bash
pixi run Rscript - <<'EOF'
cat("store:", targets::tar_config_get("store"), "\n")
cat("recorded pid:", targets::tar_pid(), "\n")
print(targets::tar_progress_summary())
EOF
```

For a target family, keep outdatedness checks narrow:

```bash
pixi run Rscript - <<'EOF'
targets::tar_outdated(
  names = tidyselect::matches("<target-or-dataset-pattern>"),
  callr_function = NULL
)
EOF
```

## Run patterns

Preview a new or regex-based selection before running it. Fail if it matches no
defined targets.

```bash
pixi run Rscript - <<'EOF'
selection <- targets::tar_manifest(
  names = tidyselect::matches("<target-or-dataset-pattern>"),
  fields = c(name, description),
  callr_function = NULL
)
stopifnot(nrow(selection) > 0L)
print(selection)

targets::tar_make(
  names = tidyselect::matches("<target-or-dataset-pattern>")
)
EOF
```

Run all targets only when explicitly requested:

```bash
pixi run Rscript - <<'EOF'
targets::tar_make()
EOF
```

For mapped targets, select the defined target stem and configured suffix rather
than a dynamic branch hash. Use names from `tar_manifest()` or `tar_meta()` and
prefer checkpoint description tags when they express the requested endpoint.

## Practical Notes

- Use `load_CFG("<dataset>")` only for interactive configuration probes outside
  target commands.
- If targets error, switch to `multiomer-fix-errors`.
