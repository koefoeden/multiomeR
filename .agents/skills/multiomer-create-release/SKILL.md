---
name: multiomer-create-release
description: Prepare or publish multiomeR releases using the repository release contract, destination-specific notes and explicit validation evidence.
---

# multiomeR Create Release

Read `RELEASES.md`, `NEWS.md`, any downstream release notes, current remote refs,
tags and GitHub releases before selecting a boundary. The workflow is not a root
R package; never add `DESCRIPTION` for versioning. When the user asks you to
establish or revise the convention, document that choice and proceed within the
authorized scope instead of asking them to supply a version source.

Distinguish preparation from publication. Preparation includes isolated topic
branches, migration notes and destination validation. Publishing requires the
user's authorization and a validated commit integrated into current `main`.
Keep proposed notes marked unreleased until then. Follow the documented tag
namespace; do not assume every downstream uses public `vX.Y.Z` tags.

Preserve running checkouts and their environments and stores. Review the final
public diff for private material; merge accepted public work forward into the
private destination and retain its own configuration and release record.

For code and graph changes use `multiomer-validation-workflow`; regenerate
manifest diagrams and render both books and the Markdown export when affected.
State exactly which tests and runtime checks passed and which were not run.

At publication, update the release date, verify the destination commit and
validation, create an annotated immutable tag, push that tag explicitly and
publish the matching GitHub notes. Verify the resulting tag and release URL.
Do not automatically bump a post-release development version.
