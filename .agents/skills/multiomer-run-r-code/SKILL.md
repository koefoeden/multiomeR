---
name: multiomer-run-r-code
description: Execute R code in multiomeR through the Pixi-managed environment. Use whenever a task requires evaluating an R snippet or script. This skill does not itself authorize running targets; use multiomer-run-pipeline when target execution is requested.
---

# multiomeR Run R Code

Run R through Pixi from the repository root. Do not use bare `R`, bare
`Rscript`, or `renv`.

## Default Pattern

Use `pixi run Rscript` with a heredoc for reproducible one-off checks:

```bash
pixi run Rscript - <<'EOF'
targets::tar_meta(fields = name) |>
  head()
EOF
```

For parse checks:

```bash
pixi run Rscript - <<'EOF'
files <- c("path/to/edited_file.R")
invisible(lapply(files, parse))
cat("parse ok\n")
EOF
```

The repository is one root `targets` project. Do not set `TAR_PROJECT` or call
legacy `activate_*()` helpers. Call `load_CFG("<dataset>")` only when an
interactive probe explicitly needs resolved dataset-level configuration.
