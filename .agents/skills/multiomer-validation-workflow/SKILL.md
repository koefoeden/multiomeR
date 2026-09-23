---
name: multiomer-validation-workflow
description: Efficient validation workflow for multiomeR R helper and targets edits. Use after changing _targets.R, extra_targets/*.R, module target files, R helpers, target names, tar_map mapping tibbles, or target graph construction.
---

# multiomeR Validation Workflow

Run R through Pixi. Use this default path after code edits and add target
execution only when static graph construction cannot prove the behavior.

## Fast Default

Validate each edited file once through the strongest applicable check. Focused
source or behavior checks already parse helpers and scripts; manifest
construction parses and sources target graph files. Use standalone `parse()`
only when an edited file is not otherwise loaded or when isolating a syntax
failure.

For target changes, build the manifest without callr. This catches missing helpers, bad target commands, malformed `tar_map()` values, and most graph-construction failures without running targets:

```bash
pixi run --use-environment-activation-cache Rscript - <<'EOF'
targets::tar_manifest(callr_function = NULL)
cat("manifest ok\n")
EOF
```

For configuration-directory changes, compare resolved settings and manifest
names before/after, check the default and selected directories independently,
and exercise invalid selections and disabled modules without changing the live
selector or store. Shared documentation examples use `configuration/`.

Finish code edits with:

```bash
git diff --check
```

## Target Execution

For checkpoint-tag changes, validate the selected targets' transitive
dependencies: a description tag does not prevent a downstream dependency from
crossing the intended review boundary. Keep `QC_checkpoint_manifest.tsv` and the
review guide consistent.

Run targets only when the changed behavior needs runtime proof, using the
preview-and-run patterns in `multiomer-run-pipeline`. Run the smallest target
that exercises the changed code, and avoid setup/download targets unless the
edit directly changed download/setup behavior.

## Invalidation Impact

`targets` hashes parsed code, so comment and formatting edits never invalidate,
and a refactor that reproduces its outputs stops the cascade at the targets that
call the changed code. For graph-construction refactors, save
`tar_manifest(fields = c(name, command, pattern))` before and after and compare
them by name; identical rows prove that no target reruns from that edit.

To confirm an impact line against an existing store, save `tar_outdated()`
before editing code (it loads the current code), apply the change, build only
the earliest changed targets, and compare a second `tar_outdated()` with the
baseline. File targets hash content only: a moved but byte-identical input keeps
its hash, while a reformatted copy does not.

## Known Waste

Do not use `tar_validate()` while the project sets `error = "trim"`: the
installed `targets` version rejects that option before useful validation.
Prefer manifest construction and targeted runtime checks.
