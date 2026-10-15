---
name: multiomer-run-r-code
description: Execute R snippets and scripts in multiomeR through Pixi, using one-off Rscript calls or scratch scripts. Use whenever a task requires evaluating R code. This skill does not itself authorize running targets; use multiomer-run-pipeline when target execution is requested.
---

# multiomeR Run R Code

Run R from the repository root through Pixi; never use bare `R` or `Rscript`.
This skill does not authorize target execution.

## Choose a Mode

- **One invocation:** use
  `pixi run --use-environment-activation-cache Rscript`.
- **Substantial multi-step work:** write a scratch `.R` script outside tracked
  paths and run it with `Rscript`.
- **Expensive reusable state:** use an interactive R session only when
  reloading objects is materially costly. Send one complete expression at a
  time, and confirm conclusions in a fresh `Rscript` process.

For parallel work, use separate cached `pixi run` commands. Track each PID and
exit status, separate their output, and avoid concurrent writes to the same
files or targets store.

Ordinary `Rscript` startup runs the repository's `.Rprofile`, which loads the
project runtime in about five seconds. Keep it for `_targets.R`, project
helpers and anything else that assumes the bootstrap. A probe that needs only
packages, such as a function signature or `targets` metadata, can skip it with
`Rscript --no-init-file`.

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


## Project Rules

The repository is one root `targets` project. Call `load_CFG("<dataset>")` only
when an interactive probe needs resolved dataset-level configuration. It
assigns into the global environment, so never follow it with an in-process
(`callr_function = NULL`) targets call in the same session. For aggregation
settings, filter `read_aggregation_config_tibble()` by `aggregation`. Use
`multiomer-run-pipeline` for target execution.

Do not use `targets::tar_config_set()` to redirect a scratch probe: it rewrites
`_targets.yaml`. Create a disposable scratch project instead.
