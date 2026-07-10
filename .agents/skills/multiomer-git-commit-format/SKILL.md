---
name: multiomer-git-commit-format
description: Create git commits for multiomeR with the correct message format, including impact keywords when relevant. Use when the user asks to commit changes, write a commit message, or stage and commit files in the multiomeR repository.
---

# multiomeR Git Commit Format

## Workflow

Inspect the staged diff and write the message through stdin so line breaks are
preserved exactly:

```bash
git commit -F - <<'EOF'
<short imperative summary>

[optional body: what changed and why]

<mandatory impact keyword line(s)>
EOF
```

Use `multiomer-impact-keyword-lines` to classify the commit. Every commit gets
at least one impact line; use `non_target_breaking` for changes that invalidate
no existing target.

Message shape:

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

**Bug fix — changes one terminal plot target:**
```
fix ATAC UMAP column filter

str_starts("score_") returns logical; replace with str_subset("^score_")

contained_target_breaking: categorical.UMAPs.ATAC
```

**Refactor — changes an intermediate target whose output flows downstream:**
```
revise accepted ATAC cell metadata

Change the accepted barcode set used by downstream ATAC processing.

cascading_target_breaking: metadata_w_cell_types_tibble.ATAC
```
