---
name: multiomer-validation-workflow
description: Efficient validation workflow for multiomeR R helper and targets edits. Use after changing _targets.R, extra_targets/*.R, module target files, R helpers, target names, tar_map mapping tibbles, or target graph construction.
---

# multiomeR Validation Workflow

Run R through Pixi. Use this default path after code edits and add target
execution only when static graph construction cannot prove the behavior.

## Fast Default

Parse every edited R file in one R session:

```bash
pixi run Rscript - <<'EOF'
files <- c(
  "path/to/edited_file.R"
)
invisible(lapply(files, parse))
cat("parse ok\n")
EOF
```

For target changes, build the manifest without callr. This catches missing helpers, bad target commands, malformed `tar_map()` values, and most graph-construction failures without running targets:

```bash
pixi run Rscript - <<'EOF'
targets::tar_manifest(callr_function = NULL)
cat("manifest ok\n")
EOF
```

Finish code edits with:

```bash
git diff --check
```

## Target Execution

Run targets only when the changed behavior needs runtime proof. Preview the
selection and fail on an empty match:

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

If the target is not dataset-suffixed, run the smallest target that exercises the changed code. Avoid setup/download targets unless the edit directly changed download/setup behavior.

## Known Waste

Do not use `tar_validate()` while the project sets `error = "trim"`: the
installed `targets` version rejects that option before useful validation.
Prefer manifest construction and targeted runtime checks.
