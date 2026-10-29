---
name: multiomer-create-github-issue
description: Draft or create GitHub issues for multiomeR using labels that exist in the repository, under a matching parent issue where the repository groups issues that way. Use when the user asks to draft issue text or explicitly open, file, or create an issue, and when AGENTS.md says to file tabled work or ideas as issues without asking. Drafting text does not authorize publishing it.
---

# multiomeR Creating GitHub Issues

## Workflow

1. Ground the title, problem statement, and acceptance criteria in the current
   repository. Remove private paths, sample identifiers, and data.
2. If the user asked for a draft, return the proposed title, body, labels and
   parent without creating an issue.
3. Before creating an issue, read the current labels with `gh label list` and
   apply at least one label that exists, such as `bug`, `enhancement`, or
   `documentation`.
4. Find a parent. List the open top-level issues with
   `gh issue list --state open --search "no:parent-issue"` and read the likely
   matches, such as an umbrella for the same rebuild, module or decisions. File
   the issue under the best one with `gh issue create --parent <number>`, or
   leave it top-level if none fits.
5. Create the issue only when the user explicitly requested publication, or
   when the repository's `AGENTS.md` says to file such issues without asking,
   then return its URL.

Target-impact keywords classify implemented changes in commits and pull
requests. Do not require them on speculative issues. Apply an impact label only
when that label exists and the affected target boundary is already known.
