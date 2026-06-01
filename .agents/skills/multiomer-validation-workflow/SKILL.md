---
name: multiomer-validation-workflow
description: Efficient validation workflow for multiomeR R helper and targets edits. Use after changing _targets.R, extra_targets/*.R, module target files, R helpers, target names, tar_map mapping tibbles, or target graph construction.
---

# multiomeR Validation Workflow

Run R through pixi only. Use this as the default validation path after code edits; add target execution only when the change needs runtime proof.

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

Always finish substantive edits with:

```bash
git diff --check
```

## Target Execution

Run targets only when parsing and graph construction cannot prove the behavior. Keep the run narrow:

```r
tar_make(names = matches("<target_name>.*muscle_test"))
```

If the target is not dataset-suffixed, run the smallest target that exercises the changed code. Avoid setup/download targets unless the edit directly changed download/setup behavior.

## Known Waste

Do not use `tar_validate()` as a routine check while this project has `settings$error = "trim"` and the installed `targets` version rejects it. It fails before useful validation. Prefer manifest construction plus direct in-memory target inspection.
