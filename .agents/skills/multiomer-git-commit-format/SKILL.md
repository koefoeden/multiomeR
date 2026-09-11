---
name: multiomer-git-commit-format
description: Create git commits for multiomeR with the correct message format, including target-impact lines and Codex task provenance. Use when the user asks to commit changes, write a commit message, or stage and commit files in the multiomeR repository.
---

# multiomeR Git Commit Format

## Workflow

Inspect the staged diff. Resolve the current Codex task ID before committing:

```bash
codex_task_id=${CODEX_THREAD_ID:-${CODEX_SESSION_ID:-}}
test -n "$codex_task_id" || {
  printf 'No Codex task ID is available; stop before committing.\n' >&2
  exit 1
}
printf 'Codex task ID: %s\n' "$codex_task_id"
```

Prefer `CODEX_THREAD_ID`; use `CODEX_SESSION_ID` when the thread variable is
unavailable. Do not guess an ID or copy one from another task. Write the
resolved literal value into a `Codex-Task-ID` Git trailer, then write the
message through stdin so line breaks are preserved exactly:

```bash
git commit -F - <<'EOF'
<short imperative summary>

[optional body: what changed and why]

<mandatory impact keyword line(s)>

Codex-Task-ID: <resolved task ID>
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

Codex-Task-ID: <resolved task ID>
```

Include exactly one `Codex-Task-ID` trailer on a new commit. When amending a
commit, preserve any distinct existing task-ID trailers and add the current ID
only when the current task materially contributed to the amended commit.

## Examples

**New feature — no existing target touched:**
```
add continuous SNN UMAP plot target

non_target_breaking

Codex-Task-ID: <resolved task ID>
```

**Bug fix — changes one terminal plot target:**
```
fix ATAC UMAP column filter

str_starts("score_") returns logical; replace with str_subset("^score_")

contained_target_breaking: categorical.UMAPs.7_ATAC_QC

Codex-Task-ID: <resolved task ID>
```

**Refactor — changes an intermediate target whose output flows downstream:**
```
revise accepted ATAC cell metadata

Change the accepted barcode set used by downstream ATAC processing.

cascading_target_breaking: metadata_w_cell_types_tibble.ATAC

Codex-Task-ID: <resolved task ID>
```
