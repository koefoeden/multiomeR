---
name: multiomer-run-pipeline
description: Run the active multiomeR root targets workflow with the correct pixi R environment. Use when the user asks to run, execute, re-run, or check multiomeR targets, datasets, target families, or target locks.
---

# multiomeR Run Pipeline

Read `.agents/wiki/articles/root-targets-workflow.md` first if you need current layout details. Run R through pixi; see `multiomer-run-r-code`.

## Run Patterns

Run all targets only when explicitly requested:

```bash
pixi run Rscript - <<'EOF'
targets::tar_make()
EOF
```

Run specific datasets:

```bash
pixi run Rscript - <<'EOF'
targets::tar_make(names = matches("muscle_test|blood_human"))
EOF
```

Run specific target families or target/dataset combinations:

```bash
pixi run Rscript - <<'EOF'
targets::tar_make(names = matches("metadata_w_cell_types_tibble\\.ATAC.*muscle_test"))
EOF
```

Prefer narrow runs after edits. For dynamic or mapped targets, match the defined target stem plus dataset suffix rather than guessing branch hashes.

## Practical Notes

- Call `load_CFG("<dataset>")` inside the script only when you need dataset-level config values outside target commands.
- Use current target names from the code or `tar_meta()`, not old Seurat-era names.
- Run optional `differential_analyses` and `genetic_enrichment` work through the root workflow; they are aggregation-level modules selected by the `modules` field in `cfg_aggregations.yaml`.
- The active store remains `pipelines/processing_and_aggregation/outputs` during the root migration.
- Keep commands visible to the user for long or heavy runs.
- If targets error, switch to `multiomer-fix-errors` and inspect with `list_distinct_errored_targets()` / `inspect_target_workspace()`.
