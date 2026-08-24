---
name: multiomer-run-r-code
description: Execute R snippets and scripts in multiomeR through Pixi, using a one-off Rscript, a retained Pixi shell for repeated processes, or an optional stateful R REPL. Use whenever a task requires evaluating R code. This skill does not itself authorize running targets; use multiomer-run-pipeline when target execution is requested.
---

# multiomeR Run R Code

Run R from the repository root through Pixi. Always pass
`--use-environment-activation-cache` to `pixi run`. This skill does not
authorize target execution.

## Choose a Mode

- **One invocation:** use
  `pixi run --use-environment-activation-cache Rscript`.
- **Repeated work:** retain one verified Pixi shell and start a fresh `Rscript`
  process for each check.
- **Expensive reusable state:** use an R REPL inside the verified Pixi shell
  only when repeatedly reconstructing or loading objects is materially costly.

For substantial multi-step exploration, stage a scratch `.R` script and run it
with `Rscript`. Source it into a retained REPL only when its objects need to
persist. Keep temporary scripts outside tracked paths unless they are intended
deliverables.

For parallel work, either use separate cached `pixi run` commands or launch
independent R processes from one verified Pixi shell. Track each PID and exit
status, separate their output, and avoid concurrent writes to the same files or
targets store.

Use ordinary `Rscript` startup so the repository's `.Rprofile` initializes the
project runtime. Do not suppress repository startup: target and helper checks
may depend on that bootstrap even when they look structural. Never use bare `R`
or `Rscript` outside a verified Pixi shell. Confirm important REPL-derived
conclusions in a fresh `Rscript` process.

## One-Off Check

```bash
pixi run --use-environment-activation-cache Rscript - <<'EOF'
# R code
EOF
```

For a parse check:

```bash
pixi run --use-environment-activation-cache Rscript -e 'invisible(parse(file = "path/to/file.R"))'
```

## Retain Pixi Activation

Start a PTY-backed session and retain its session ID:

```bash
pixi shell --change-ps1 false --no-completions
```

Verify the environment once inside that session:

```bash
test "${PIXI_IN_SHELL:-}" = "1"
test "$(realpath -- "$PIXI_PROJECT_ROOT")" = "$(realpath -- "$PWD")"
command -v R Rscript
```

Submit sequential `Rscript` commands and require each process to exit
successfully. Fresh processes retain Pixi activation while isolating R state
and failures. The same shell may launch independent background processes when
their CPU, memory, output, and write targets do not conflict.

## Stateful Exploration

```bash
R --quiet --no-save --no-restore
```

Retain only deliberately reusable objects. Wrap transient work in `local()`,
send one complete expression at a time, and print a unique completion marker
inside the successful expression:

```r
local({
  result <- summary(expensive_object)
  print(result)
  assign("result", result, envir = .GlobalEnv)
  cat("CODEX_R_DONE_<unique>\n")
})
```

Absence of the marker means failure or incomplete work; a prompt alone is not a
success signal. Avoid large chunks of independent top-level expressions because
R can continue after an error and leave partial global state.

Treat REPL results as exploratory. Restart after environment, package,
configuration, `.Rprofile`, or bootstrap changes. Reload edited functions and
data explicitly. End with `q("no")`, close the Pixi shell, and confirm important
conclusions in a fresh `Rscript` process.

## Project Rules

The repository is one root `targets` project. Do not set `TAR_PROJECT` or call
legacy `activate_*()` helpers. Call `load_CFG("<dataset>")` only when an
interactive probe needs resolved dataset-level configuration. Use
`multiomer-run-pipeline` for target execution.

Do not use `targets::tar_config_set()` to redirect a scratch probe: it rewrites
`_targets.yaml`. Create a disposable scratch project instead.
