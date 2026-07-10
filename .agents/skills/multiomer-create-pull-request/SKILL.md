---
name: multiomer-create-pull-request
description: Draft, create, or update a multiomeR pull request from the complete branch diff, including validation and deduplicated target-impact lines. Use for PR text drafts, explicitly requested draft PRs, new PRs, and refreshes of existing PRs. Do not publish a text-only draft request.
---

# multiomeR Create Or Update Pull Request

## Workflow

Inspect the merge base, full commit list, and complete diff against the intended
base branch. Summarize the current branch state, not only the newest commit.

When updating an existing pull request:

- Re-read the current PR metadata and complete included diff.
- Regenerate the title and body unless the user requests a narrow edit.
- Preserve useful discussion context, but remove claims no longer supported by
  the branch.

## Body structure

Use only populated change sections among `New Features`, `Changes &
Optimizations`, and `Bug Fixes`. Always include:

- `Validation`: commands or checks actually completed, plus anything not run.
- `Pipeline Impact`: target invalidation after merge.

Collect impact lines from every branch commit using
`multiomer-impact-keyword-lines`. Deduplicate identical lines. If the same
target is both contained and cascading, keep the cascading line. Do not collapse
different targets without verifying their dependency relationship. If a commit
lacks a valid impact line, inspect that commit before classifying it.

If no existing target is invalidated, write exactly:

> None — no existing targets are invalidated by this PR.

Return text only for a PR-body draft. Use a draft PR when the user explicitly
asks to open one as draft; otherwise create or update the requested PR and
return its URL.
