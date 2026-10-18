---
name: multiomer-estimate-rebuild-impact
description: Estimate which existing targets a code, configuration or input change, a merge, or a run will actually rebuild, and what that costs. Use before running targets after such changes, before telling anyone that a change will or will not rerun targets, when tar_outdated() reports many outdated targets, and when classifying target impact for commits, merges, pull requests or release notes.
---

# Estimate Rebuild Impact

A target reruns when its command, a function or global object it calls, or an
upstream data hash changes. `targets` hashes parsed code, so comment and
formatting edits never invalidate. File targets hash content only: a moved but
byte-identical input keeps its hash, while a reformatted copy does not. A rerun
that reproduces its value stops the cascade: `tar_make()` skips consumers whose
inputs keep their hashes.

## Outdated is an upper bound

`tar_outdated()` also lists every target downstream of an outdated one, since
it cannot predict whether a rerun reproduces its value. Never report its count
as the cost of a run. List the outdated roots instead, the outdated targets
with no outdated upstream target:

```bash
pixi run --use-environment-activation-cache Rscript - <<'EOF'
network <- targets::tar_network(
  names = tidyselect::matches("<target-or-dataset-pattern>"),
  targets_only = TRUE, callr_function = NULL, reporter = "silent"
)
outdated <- network$vertices$name[network$vertices$status == "outdated"]
print(setdiff(outdated, network$edges$to[network$edges$from %in% outdated]))
EOF
```

When the roots and their direct consumers are cheap, rebuild only them and
compare their `tar_meta()` data hashes before and after; consumers that
reproduce theirs stop the cascade, which a second root listing confirms. Until
then, a target whose own command changed stays hidden below an outdated root;
`tar_sitrep(names = ...)` shows it at once in its `command` column. Do
this before any run whose cascade could reach WNN: its neighbour search is
nondeterministic, so WNN and everything downstream would rebuild with new
values.

## Code and configuration changes

Configuration is read when the graph is built, so its values appear in target
commands. Save `tar_manifest(fields = c(name, command, pattern))` before and
after a configuration, merge or graph-construction change and compare the rows
by name: identical rows rule out reruns from changed commands. They cannot show
an edited function body, which reruns every target that calls it directly or
through other functions, whatever configuration branch runs; list the roots as
above for that. To change one consumer, edit its call site rather than a shared
helper. Descriptions and resources are not cues, so edits limited to them
rerun nothing. Also check the
content of changed file targets, such as manifests. If no manifest was saved
before the edit, build it with `tar_manifest(script = ...)` from a temporary
copy of `_targets.R` in the repository root whose `source()` calls read
`git show HEAD:` copies of the edited target files.

Adding or activating GEM well rows changes only the `GEM_well_config_tsv` file
target. Each aggregation's `GEM_well_metadata_tibble` reruns and reproduces its
value unless the aggregation lists a changed well; compare its new value with
`tar_read()` before calling the change contained. Rebuilding both keeps the
aggregation current.

Dynamic branches take their names from their upstream branches, and those of
an `iteration = "vector"` stem from each element's value, never from the
pattern: a shifted `slice()` index or a key added beside existing ones reruns
no existing branch. Mapping over a different or regrouped upstream renames
every downstream branch, which then reruns.

To confirm a commit's impact line against an existing store, save the outdated
targets of the affected families before editing code, which the check loads,
with the commit they ran at. A full-graph baseline takes minutes and fails on
any stale pattern elsewhere. Apply the change, build only the earliest changed
targets, and compare a second listing with the baseline. Never revert shared
files to recreate a missed baseline; use the `git show HEAD:` copies above.
Rebuilding a target overwrites its `command`, `depend` and `seconds` records,
so before a run whose reruns you may have to explain, or whose run times are
measured, save `tar_sitrep()` and `tar_meta()` for the selection.
