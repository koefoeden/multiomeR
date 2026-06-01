---
name: multiomer-impact-keyword-lines
description: Description of multiomeR impact keyword lines that are added to commit messages, PR bodies, or release notes. Use when information is needed on what multiomeR impact keyword lines are.
---

# multiomeR Impact Keyword Lines

If a commit, pull request, or release changes any existing multiomeR targets workflow code, always append at least **one** of the following impact keyword lines at the end of the commit message/PR body/release notes:

| Impact keyword line | When to use | Example(s) |
|---|---|---|
| `non_target_breaking` | No existing targets will be invalidated by this change | Code comments, white-space changes, or changes to the `resources` argument of a target |
| `contained_target_breaking: <target>` | The named target's will re-run, but its output is identical to before, or it does not have any downstream dependencies  | A plotting target, standalone file export, code refactoring, optimizations, etc |
| `cascading_target_breaking: <target> [<dataset-scope>]` | The named target's output changes **and** downstream targets consume it, so they must also re-run. In some cases, such as for dataset-specific bug-fixes, the output will only change in a limited number of datasets. If this is the case, please add the <dataset-scope> to the impact keyword line to indicate the specific datasets that are affected |  Addition of new column to tibble, matrix or graph to Seurat object, etc |
