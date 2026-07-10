---
name: multiomer-impact-keyword-lines
description: Classify multiomeR target invalidation in commit messages, pull requests, and release notes. Use when writing or interpreting non_target_breaking, contained_target_breaking, or cascading_target_breaking lines. These are not GitHub issue-label requirements.
---

# multiomeR Impact Keyword Lines

Append at least one impact line to each multiomeR commit and summarize applicable
lines in pull requests or release notes.

| Impact keyword line | When to use | Example(s) |
|---|---|---|
| `non_target_breaking` | No existing targets will be invalidated by this change | Code comments, white-space changes, or changes to the `resources` argument of a target |
| `contained_target_breaking: <target>` | The named target reruns, but its output hash is unchanged or it has no downstream consumers | Refactoring a target command, changing a terminal plot, or revising a standalone export |
| `cascading_target_breaking: <target> [<dataset-scope>]` | The named target's output changes and downstream targets consume it; append a configured scope when only some datasets or aggregations are affected | Adding a column to a consumed tibble or changing accepted cells in a matrix |

List multiple earliest affected targets when no single upstream target captures
the invalidation boundary. Do not replace target names with vague families.
