# Releases

multiomeR is a configurable analysis workflow, not a root R package. `NEWS.md`
is the version and migration record; no root `DESCRIPTION` or independent
version file is needed. The nested `multiomeRCore` package has its own version.

## Versions and tags

Use semantic versions and immutable annotated tags `vMAJOR.MINOR.PATCH` on
validated public `main` commits. While the workflow is below 1.0, increment the
minor version for changes to configuration, target names, output schemas or
scientific methods; increment the patch version for compatible fixes. After
1.0, incompatible changes require a major version. Target invalidation is
reported separately: even a compatible fix can require recomputation.

Use `MAJOR.MINOR.PATCH (unreleased)` at the top of `NEWS.md` while preparing a
release. Replace `unreleased` with the release date only when publishing. Do not
bump a development version automatically after release; start the next section
when the next change needs it. Optional release candidates use `vX.Y.Z-rc.N`.

Version 0.5.0 starts this convention after the untagged 0.4.1.9000 development
series. Archive tags are historical preservation points, not software releases.
Do not manufacture retrospective release tags for them.

## Preparation and publication

Prepare coherent topic commits on a branch based on current public `main`.
Review the final diff for private configuration, paths, outputs and credentials.
Document user-visible changes, migration steps, scientific interpretation and
precise target-invalidation boundaries in `NEWS.md`.

Validate the destination configuration, relevant numerical tests, generated
manifest and dependency diagrams, and both documentation books. Record whether
validation used synthetic fixtures, saved data or actual target execution. A
successful manifest alone is not evidence of numerical or biological validity.
Use isolated worktrees when another checkout is running a pipeline; never alter
that checkout's files, environment or targets store during release preparation.

After review, merge the prepared branch into current `main`, verify the resulting
commit and validation, create an annotated version tag and publish the GitHub
release from the matching notes. Do not tag a preparation branch as a completed
release. Publish only the intended tag explicitly, never all local tags.

Downstreams merge the public release forward, then record their own overlays
and configurations in separate release notes and a distinct tag namespace.
Keep their release metadata out of public history.
