---
name: multiomer-create-release
description: Assess release readiness and draft or create a multiomeR GitHub release only when the repository has an explicit current versioning contract. Use for release planning, release drafts, tags, or publication. Do not recreate DESCRIPTION or invent a version source.
---

# multiomeR Create Release

## Preflight

Before changing state:

1. Inspect the latest Git tags and GitHub releases, `NEWS.md`, and any tracked
   version source.
2. Confirm the release commit is on an up-to-date `main` branch with a clean
   worktree and completed validation.
3. Determine whether the user requested release notes, a draft GitHub release,
   or publication.

For a release-note draft, return text only and make no repository or GitHub
changes.

If no tracked version source or established tag history exists, stop and ask the
user to define the versioning contract. Do not restore the deleted
`DESCRIPTION` file merely to satisfy this skill.

## Release workflow

Once a version contract exists:

1. Compare the release commit with the previous release boundary.
2. Finalize the matching `NEWS.md` section from merged changes and PRs. Include
   target-impact lines where relevant.
3. Update the canonical version source if the contract requires it, then commit
   with `multiomer-git-commit-format`.
4. Create an annotated `v<x.y.z>` tag at the validated release commit and push
   the commit and tag.
5. Draft or publish the GitHub release exactly as requested, using the finalized
   NEWS section, and verify its tag and URL.

Do not perform an automatic post-release development-version bump unless the
defined versioning contract explicitly requires one.
