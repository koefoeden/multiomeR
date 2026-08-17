---
name: multiomer-run-pipeline
description: Inspect or run the active multiomeR root targets workflow through Pixi. Use for explicit requests to run or rerun targets and for read-only target status, progress, outdatedness, or process-lock checks. Do not start target execution for a request that only asks for inspection.
---

# multiomeR Run Pipeline

The active project is `_targets.R`; use the store resolved by `_targets.yaml`.
Differential analyses and genetic enrichment are aggregation-level modules in
this same graph. Run R through Pixi; see `multiomer-run-r-code`.

## Read-only inspection

Use the live configuration instead of a hardcoded store. `tar_pid()` can retain
the PID of a completed run, so verify it with `ps` before reporting that a run is
active.

```bash
pixi run --use-environment-activation-cache Rscript - <<'EOF'
cat("store:", targets::tar_config_get("store"), "\n")
pid <- targets::tar_pid()
cat("recorded pid:", pid, "\n")
process <- if (length(pid) == 1L && !is.na(pid)) {
  system2("ps", c("-p", pid, "-o", "pid=,stat=,cmd="), stdout = TRUE, stderr = TRUE)
} else {
  character()
}
cat("live process:", if (length(process)) process else "<none>", "\n")
print(targets::tar_progress_summary())
EOF
```

For a target family, keep outdatedness checks narrow:

```bash
pixi run --use-environment-activation-cache Rscript - <<'EOF'
targets::tar_outdated(
  names = tidyselect::matches("<target-or-dataset-pattern>"),
  callr_function = NULL
)
EOF
```

## Run patterns

Immediately before `tar_make()`, confirm that no live driver already owns the
configured store. Do not infer ownership from a recorded PID without the process
check above.

Preview a new or regex-based selection before running it. Fail if it matches no
defined targets. Run a known exact target directly instead of constructing the
graph once for preview and again for execution.

```bash
pixi run --use-environment-activation-cache Rscript - <<'EOF'
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

After `tar_make()` returns, inspect the selected endpoints. Report completion
only when each endpoint has data, no error, and terminal progress; partial
progress or a returned driver is not success.

Run all targets only when explicitly requested:

```bash
pixi run --use-environment-activation-cache Rscript - <<'EOF'
targets::tar_make()
EOF
```

For mapped targets, select the defined target stem and configured suffix rather
than a dynamic branch hash. Use names from `tar_manifest()` or `tar_meta()` and
prefer checkpoint description tags when they express the requested endpoint.

For a narrow selection whose upstream targets are known to be current, consider
`shortcut = TRUE` in `tar_outdated()` or `tar_make()`. This can avoid expensive
upstream file-metadata checks, especially on network filesystems, by trusting
stored upstream metadata. Do not use it when upstream code, configuration, data,
or file freshness may have changed: the shortcut can incorrectly skip work
because it assumes those dependencies are already up to date. Keep the default
`shortcut = FALSE` for production validation and whenever freshness is uncertain.

## Practical Notes

- Use `load_CFG("<dataset>")` only for interactive configuration probes outside
  target commands.
- If targets error, switch to `multiomer-fix-errors`.
