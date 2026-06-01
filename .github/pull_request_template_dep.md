## New Features
- Describe new features added in this PR.

## Changes & Optimizations
- Describe changes and optimizations made in this PR.

## Bug Fixes
- Describe bug fixes included in this PR.

## Pipeline Impact

Targets that will be invalidated the next time the pipeline runs after merging:

| Impact keyword line | Reason |
|--------|-------|
| `<impact_keyword_line>` | `<brief reason, e.g. "added new column to tibble xyz">` |

Explanation:
| Impact keyword line | When to use | Example(s) |
|---|---|---|
| `non_target_breaking` | No existing targets will be invalidated by this change | Code comments, white-space changes, or changes to the `resources` argument of a target |
| `contained_target_breaking: <target>` | The named target's will re-run, but its output is identical to before, or it does not have any downstream dependencies  | A plotting target, standalone file export, code refactoring, optimizations, etc |
| `cascading_target_breaking: <target> [<dataset-scope>]` | The named target's output changes **and** downstream targets consume it, so they must also re-run. In some cases, such as for dataset-specific bug-fixes, the output will only change in a limited number of datasets. If this is the case, please add the <dataset-scope> to the impact keyword line to indicate the specific datasets that are affected |  Addition of new column to tibble, matrix or graph to Seurat object, etc |
