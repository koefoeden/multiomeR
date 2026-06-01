---
name: multiomer-create-github-issue
description: Instructions for creating GitHub issues in the multiomeR project. Use when the user asks to create, open, file, or draft a GitHub issue for multiomeR.
---

# multiomeR Creating GitHub Issues

## Labels

Every issue **must** receive at least one label. Use `gh label list` to see all available labels before creating an issue.

Three labels are especially important — apply exactly one of them to classify the scope of impact:
`non_breaking`, `contained_target_breaking`, `cascading_target_breaking`. See the `multiomer-impact-keyword-lines` skill for more information on these labels.

Apply additional labels (e.g. `bug`, `enhancement`, `documentation`) as appropriate alongside the impact labels above.
