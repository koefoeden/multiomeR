---
name: multiomer-create-github-issue
description: Draft or create GitHub issues for multiomeR using labels that exist in the repository. Use when the user asks to draft issue text or explicitly open, file, or create an issue. Drafting text does not authorize publishing it.
---

# multiomeR Creating GitHub Issues

## Workflow

1. Ground the title, problem statement, and acceptance criteria in the current
   repository. Remove private paths, sample identifiers, and data.
2. If the user asked for a draft, return the proposed title, body, and labels
   without creating an issue.
3. Before creating an issue, read the current labels with `gh label list` and
   apply at least one label that exists, such as `bug`, `enhancement`, or
   `documentation`.
4. Create the issue only when the user explicitly requested publication, then
   return its URL.

Target-impact keywords classify implemented changes in commits and pull
requests. Do not require them on speculative issues. Apply an impact label only
when that label exists and the affected target boundary is already known.
