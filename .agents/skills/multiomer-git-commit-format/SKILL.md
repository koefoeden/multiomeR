---
name: multiomer-git-commit-format
description: Create git commits for multiomeR with the correct message format, including impact keywords when relevant. Use when the user asks to commit changes, write a commit message, or stage and commit files in the multiomeR repository.
---

# multiomeR Git Commit Format

## Workflow

Use standard workflow to commit changes. However, be sure to follow the format below, by passing the message via a HEREDOC to `git commit -m` to preserve newlines:
```
<short imperative summary>

[optional body: more details on what changed and why]

[mandatory impact keyword line(s) — see `multiomer-impact-keyword-lines` skill ]
```

## Examples

**New feature — no existing target touched:**
```
add continuous SNN UMAP plot target

non_target_breaking
```

**Bug fix — changes output of one plot target:**
```
fix SCAVENGE UMAP column filter

str_starts("score_") returns logical; replace with str_subset("^score_")

contained_target_breaking: SCAVENGE_UMAP_plots
```

**Refactor — changes an intermediate target whose output flows downstream:**
```
replace SingleFeatureMatrix with BPCells peak matrix generation

Switches per-cluster ATAC peak matrix computation to the BPCells
backend; coerces to dgCMatrix at the ChromatinAssay boundary so all
downstream Signac code is unchanged.

cascading_target_breaking: ATAC_per_cluster_peak_BPCells_matrix_dir
```
