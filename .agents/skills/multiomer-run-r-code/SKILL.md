---
name: multiomer-run-r-code
description: Public-reusable guidance for executing R code in multiomeR through the pixi-managed R environment. Use whenever you need to run R code, test a snippet, or execute an R script for multiomeR in a public or local checkout.
---

# multiomeR Run R Code

Run R through pixi in the multiomeR repository. Do not use bare `R`, bare `Rscript`, or `renv`.

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

The active repo is a root `targets` project. Do not set `TAR_PROJECT` or call `activate_*()` for current multiomeR work. Call `load_CFG("<dataset>")` explicitly inside the script only when you need dataset-level config values in the global environment.
