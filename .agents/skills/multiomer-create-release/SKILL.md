---
name: multiomer-create-release
description: Prepare or publish multiomeR releases using the repository release contract, destination-specific notes and explicit validation evidence.
---

# multiomeR Create Release

Read any downstream release notes, current remote refs, tags and GitHub releases
before selecting a boundary. Choose the version as the README's Updating section
defines; optional release candidates are `vX.Y.Z-rc.N`. Versioning starts at
0.5.0; never tag older or archive commits retrospectively. The workflow is not a
root R package; never add `DESCRIPTION` or a version file. When the user asks
you to establish or revise the convention, document that choice and proceed
within the authorized scope instead of asking them to supply a version source.

Distinguish preparation from publication. Preparation includes isolated topic
branches, migration notes and destination validation. Publishing requires the
user's authorization and a validated commit integrated into current `main`.
Draft the notes from the commit messages and diff since the previous tag.
Downstreams tag in their own namespace, not as public `vX.Y.Z`.

Preserve running checkouts and their environments and stores. Review the final
public diff for private material; merge accepted public work forward into the
private destination and retain its own configuration and release record.

For code and graph changes use `multiomer-validation-workflow`; regenerate
manifest diagrams and render the manual and the Markdown export when affected.
State exactly which tests and runtime checks passed and which were not run.

Before publication, run the documented `validate-local` command on the exact
candidate with a clean working tree, once with `configuration/` and once with
`configuration_dev/` selected. Retain the report as
`dev/validation/reports/VALIDATION-<version>.md` and require all preflight
and output checks to pass; cached target results may be reused.

At publication, verify the destination commit and validation, create an
annotated immutable tag, push that tag explicitly and publish the drafted
GitHub notes, ordered as `multiomer-git-commit-format` describes. Verify the
resulting tag and release URL.
