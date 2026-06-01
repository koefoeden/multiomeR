---
name: multiomer-create-release
description: Guidance on creating a new GitHub release for multiomeR. Use when the user asks to create, draft, or submit a new multiomeR release.
---

# multiomeR Create Release

## Purpose

This skill guides the creation of a new GitHub release for multiomeR, to ensure that the conventions for updating `DESCRIPTION`, `NEWS.md` and git tags are followed.

## Versioning Rules

Use only two kinds of versions in `DESCRIPTION`:

- Development version: `x.y.z.9000` (never tagged or published as a release on GitHub)
- Released version: `x.y.z` (tagged and published as a release on GitHub.)

Usually, the user will specify one of three types of releases that they want to create:
- Patch release (x.y.z+1): bug fix, hotfix, narrow correction, or small user-visible behavior fix.
- Minor release (x.y+1.z): new targets capability, new targets, substantial workflow expansion, or major user-visible output change.
- Major release (x+1.y.z): reserve for future compatibility-era breaks if needed.

## Creating the release

1. Finalize the top `NEWS.md` development section based on changes from last time we released, ideally using and referring to merged pull requests (and potentially main-specific-commits), while also changing the label header to the correct, new version.
2. Increment the version to the correct, new release version in `DESCRIPTION`.
3. Commit the changes, create an annotated git tag `v<x.y.z>` and push to origin
4. Create a GitHub Release for the annotated tag, and publish it, with the title `multiomeR x.y.z` and notes matching the newest, finalized NEWS.md section.
5. Immediately bump `DESCRIPTION` to the next development version and restore an empty top development section in `NEWS.md`.
6. Commit the development version change and push to origin.
