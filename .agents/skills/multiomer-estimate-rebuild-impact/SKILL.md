---
name: multiomer-estimate-rebuild-impact
description: Estimate which existing targets a code, configuration or input change, a merge, or a run will actually rebuild, and what that costs. Use before running targets after such changes, when tar_outdated() reports many outdated targets, and when classifying target impact for commits, merges, pull requests or release notes.
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
reproduce theirs stop the cascade, which a second root listing confirms. Do
this before any run whose cascade could reach WNN: its neighbour search is
nondeterministic, so WNN and everything downstream would rebuild with new
values.

## Code and configuration changes

Configuration is read when the graph is built, so its values appear in target
commands. Save `tar_manifest(fields = c(name, command, pattern))` before and
after a configuration, merge or graph-construction change and compare the rows
by name: identical rows rule out reruns from changed commands. They cannot show
an edited function body, which reruns every target that calls it directly or
through other functions; list the roots as above for that. Also check the
content of changed file targets, such as manifests.

To confirm a commit's impact line against an existing store, save
`tar_outdated()` before editing code (it loads the current code), apply the
change, build only the earliest changed targets, and compare a second
`tar_outdated()` with the baseline.
