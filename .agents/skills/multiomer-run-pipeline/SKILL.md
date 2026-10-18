---
name: multiomer-run-pipeline
description: Inspect or run the active multiomeR root targets workflow through Pixi. Use for explicit requests to run or rerun targets and for read-only target status, progress, outdatedness, or process-lock checks. Do not start target execution for a request that only asks for inspection.
---

# multiomeR Run Pipeline

The active project is `_targets.R`; use the store resolved by `_targets.yaml`.
Optional analyses are aggregation-level modules in this same graph.
Configuration follows `configuration/README.md`; never change its selection
during a run. Run R through Pixi; see `multiomer-run-r-code`.

With the default callr child process, keep `names` expressions self-contained:
use literal names or patterns, or inject a vector with `!!` as in
`tidyselect::all_of(!!selection)`, rather than a selection variable defined
only in the calling R session. With `callr_function = NULL`, the session's globals
become imports instead: a variable such as `names` enters the hash of every
function that calls `names()`, so freshly built targets later look outdated.
Define no variables before an in-process call.

## Check for an active run

Run this before starting, stopping or relabelling a run, and before reporting
one as active. `tar_pid()` keeps the PID of a finished run, and a process is
visible only on its own node, so the check also reads the progress file. It
needs only `targets`, so it skips the project bootstrap:

```bash
pixi run --use-environment-activation-cache Rscript --no-init-file - <<'EOF'
store <- targets::tar_config_get("store")
pid <- targets::tar_pid(store = store)
progress <- file.path(store, "meta", "progress")
cat(
  "store:", store,
  "\npid", pid, "working directory on this node:", Sys.readlink(file.path("/proc", pid, "cwd")),
  "\nprogress written", format(difftime(Sys.time(), file.mtime(progress), units = "mins"), digits = 2), "ago;",
  sum(targets::tar_progress(store = store)$progress == "dispatched"), "targets dispatched\n"
)
EOF
```

A run owns the store when the PID's working directory is this checkout, or when
progress was written within the last few minutes with targets dispatched: the
driver then runs on another node. `NA` with older progress means no active run.
`tar_progress_summary()` reports the run's counts.

## Outdatedness

Keep outdatedness checks narrow for a target family. A full-graph check takes
minutes and suits an impact baseline:

```bash
pixi run --use-environment-activation-cache Rscript - <<'EOF'
targets::tar_outdated(
  names = tidyselect::matches("<target-or-dataset-pattern>"),
  callr_function = NULL
)
EOF
```

A check fails with "unequal lengths of vars in map(...)" or "index is out of
bounds in pattern = slice()" when a pattern's number of branches changed since
it was built, as after a change to an aggregation's GWAS list: the check
expands downstream patterns from the stale stored branches. `tar_make()`
rebuilds the upstream pattern first and is unaffected. Leave the downstream
patterns out of `names` until then; they are outdated anyway.

An outdated count is not the cost of a run. Before running after code,
configuration or input changes, estimate the cost with
`multiomer-estimate-rebuild-impact`.

## Run patterns

Immediately before `tar_make()`, run the active-run check: `targets` itself
refuses a second driver only on the same node. Use one driver per store.
Workers start in the checkout and rerun `.Rprofile`, so while a run is active,
leave no broken `R/` file, changed package pin or changed Pixi environment in
the checkout: new workers would fail to start.

Agent tool calls have time limits, including for background commands, and
stopping the driver also stops its running targets. Start a driver that may
run for hours detached (`setsid nohup <script> > /dev/null 2>&1 < /dev/null &`),
have the script log `tar_make()` output and its exit status to a file, and
follow that log.

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

For many trivial targets, such as file tracking, `tar_make(use_crew = FALSE)`
avoids worker startup.

After `tar_make()` returns, inspect the selected endpoints. Report completion
only when each endpoint has data, no error, and terminal progress; partial
progress or a returned driver is not success.

To verify dynamic branches, read metadata with `complete_only = FALSE` and
follow `children` only for selected rows with `type == "pattern"`. Other
targets can have slice identifiers in `children` that are not target names.
Use these current branch lists rather than a name-prefix search that can
include obsolete branches. A pattern's row is written only once all its
branches complete, so after a branch errors it is missing or stale: select the
branches by `parent` instead, and read an errored branch's inputs from its
workspace (`tar_workspace()`).

Run all targets only when explicitly requested:

```bash
pixi run --use-environment-activation-cache Rscript - <<'EOF'
targets::tar_make()
EOF
```

For mapped targets, select the defined target stem and configured suffix rather
than a dynamic branch hash. Use names from `tar_manifest()` or `tar_meta()`, or
a checkpoint description tag when it expresses the requested endpoint:

```r
names = targets::tar_described_as(tidyselect::contains("checkpoint:<tag>]")) &
  tidyselect::ends_with(".<aggregation>")
```

`shortcut = TRUE` in `tar_outdated()` or `tar_make()` trusts stored upstream
metadata instead of checking it. Use it only for a narrow selection whose
upstream is known to be current, such as iterating on a figure, never after
upstream code, configuration, data or files may have changed.

If targets error, switch to `multiomer-fix-errors`.
