---
name: multiomer-accept-pull-request
description: Accepts a GitHub pull request for the multiomeR project. Use when the user asks to merge or accept a pull request.
---

# multiomeR Accept Pull Request

## Workflow

1. Read the PR, full diff, reviews, and checks. Confirm it targets `main`, is not
   a draft, has no unresolved blocking review, and is mergeable.
2. Add or update `NEWS.md` on the PR branch only when the change is user-facing
   and the entry is missing. Link the PR and commit the update with
   `multiomer-git-commit-format`. If the source branch is not writable, report
   the blocker rather than changing `main` directly.
3. Recheck the PR after any NEWS commit.
4. Merge with GitHub's merge-commit method so the reviewed branch commits are
   preserved. Do not squash or rebase.
5. Verify the PR is merged and report the merge commit. Do not delete the source
   branch unless the user requested it.
