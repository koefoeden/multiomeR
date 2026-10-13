# Releases

multiomeR is a configurable analysis workflow, not a root R package. Its GitHub
releases are the version and migration record; no root `DESCRIPTION`, version
file or `NEWS.md` is needed. The nested `multiomeRCore` package has its own
version.

## Versions and tags

Use semantic versions and immutable annotated tags `vMAJOR.MINOR.PATCH` on
validated public `main` commits. Increment the major version for incompatible
changes to configuration, target names or output schemas, the minor version for
compatible additions and scientific-method changes, and the patch version for
compatible fixes. Before 1.0, the minor version marked incompatible changes.
Target invalidation is reported separately: even a compatible fix can require
recomputation.

Optional release candidates use `vX.Y.Z-rc.N`.

Version 0.5.0 starts this convention after the untagged 0.4.1.9000 development
series; tags up to v1.2.1 also carry the notes in `NEWS.md`. Archive tags are
historical preservation points, not software releases. Do not manufacture
retrospective release tags for them.

## Preparation and publication

Prepare coherent topic commits on a branch based on current public `main`.
Review the final diff for private configuration, paths, outputs and credentials.
Commit messages record each change's user-visible effect and target impact. At
release time, draft the notes from them and the diff since the previous tag:
user-visible changes, migration steps, scientific interpretation and precise
target-invalidation boundaries.

Validate the destination configuration, relevant numerical tests, generated
manifest and dependency diagrams, and the manual. Record whether
validation used synthetic fixtures, saved data or actual target execution. A
successful manifest alone is not evidence of numerical or biological validity.
Use isolated worktrees when another checkout is running a pipeline; never alter
that checkout's files, environment or targets store during release preparation.

After review, merge the prepared branch into current `main`, verify the resulting
commit and validation, create an annotated version tag and publish the GitHub
release with the drafted notes. Do not tag a preparation branch as a completed
release. Publish only the intended tag explicitly, never all local tags.

Downstreams merge the public release forward, then record their own overlays
and configurations in separate release notes and a distinct tag namespace.
Keep their release metadata out of public history.

## Local candidate validation

Run `pixi run --use-environment-activation-cache validate-local` on the candidate
checkout with a clean working tree. It enables the three public validation
aggregations, checks external sources, runs the pipeline incrementally in the
existing store, and verifies outputs. Retain the report matching the exact
release commit as `dev/validation/reports/VALIDATION-<version>.md`. See [dev/validation/README.md](validation/README.md). No CI branch
or separate profile repository is needed.
