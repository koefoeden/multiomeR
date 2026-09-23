---
name: multiomer-git-commit-format
description: Create git commits for multiomeR with the correct message format, including target-impact lines and agent task provenance. Use when the user asks to commit changes, write or amend a commit message, or classify target invalidation for commits, pull requests, or release notes.
---

# multiomeR Git Commit Format

## Workflow

Inspect the staged diff, then write the message through stdin so line breaks
are preserved exactly:

```bash
git commit -F - <<'EOF'
<short imperative summary>

[optional body: what changed and why]

<mandatory impact keyword line(s)>

<Agent>-Task-ID: <resolved task ID>
EOF
```

## Impact Keyword Lines

Every commit gets at least one impact line; summarize the applicable lines in
pull requests and release notes.

| Impact keyword line | When to use | Example(s) |
|---|---|---|
| `non_target_breaking` | No existing targets will be invalidated by this change | Code comments, white-space changes, or changes to the `resources` argument of a target |
| `contained_target_breaking: <target>` | The named target reruns, but its output hash is unchanged or it has no downstream consumers | Refactoring a target command, changing a terminal plot, or revising a standalone export |
| `cascading_target_breaking: <target> [<dataset-scope>]` | The named target's output changes and downstream targets consume it; append a configured scope when only some datasets or aggregations are affected | Adding a column to a consumed tibble or changing accepted cells in a matrix |

List multiple earliest affected targets when no single upstream target captures
the invalidation boundary. Do not replace target names with vague families.

Trace the boundary rather than assuming a cascade. A target reruns when its
command, a function or global object it calls, or an upstream data hash
changes. Code that only constructs the graph, such as `tar_map()` values, affects
only commands that reference it; descriptions and checkpoint tags invalidate
nothing. A rerun that reproduces its hash stops the cascade: when consumers of
a changed output all reproduce theirs, report them as contained.
`tar_outdated()` cannot see this; confirm with `multiomer-validation-workflow`.

For PR and release summaries, use commit impact lines as evidence and reconcile
them with the final diff against the destination base or previous release.
Drop effects from reverted or superseded changes, deduplicate surviving lines,
and retain the cascading classification when it applies to the same target.

## Agent Task Provenance

An agent that creates a commit records its own thread or session identifier in
a trailer named after the agent. Resolve the literal value from the current
environment; do not guess an ID or copy one from another task.

| Agent | Trailer | Identifier |
|---|---|---|
| Codex | `Codex-Task-ID` | `CODEX_THREAD_ID`, else `CODEX_SESSION_ID` |
| Claude Code | `Claude-Task-ID` | `CLAUDE_CODE_SESSION_ID` |

Other agents use their own thread or session identifier in a matching
`<Agent>-Task-ID` trailer. If none is available, ask the user before committing
without one. Include exactly one trailer for the current agent on a new commit.
When amending, preserve distinct existing task-ID trailers and add the current
one only when the current task materially contributed to the amended commit.

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

Claude-Task-ID: <resolved task ID>
```

**Refactor — changes an intermediate target whose output flows downstream:**
```
revise accepted ATAC cell metadata

Change the accepted barcode set used by downstream ATAC processing.

cascading_target_breaking: metadata_w_cell_types_tibble.ATAC

Codex-Task-ID: <resolved task ID>
```
