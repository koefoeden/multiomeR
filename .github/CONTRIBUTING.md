# Contributing

multiomeR welcomes bug reports and feature requests that are broadly relevant to users analyzing 10x Genomics Multiome data.

The project is intentionally a framework, not a closed product. Users are expected to adapt config files, target selections, helper functions, and downstream modules to their own biological questions. Contributions are most useful when they improve the shared framework without adding unnecessary complexity for everyone else.

## What to Report

Please open a GitHub issue for:

- reproducible pipeline failures
- incorrect or surprising target outputs
- missing documentation for common 10x Multiome workflows
- feature requests that are likely to help multiple projects or datasets
- installation or configuration problems with the public example workflow

When reporting a bug, include:

- the target name or command that failed
- the relevant traceback or error message
- the config fields involved, with private paths or sample IDs removed
- the multiomeR commit SHA
- whether you are using local execution or a scheduler-backed `crew` controller

## Feature Requests

Good feature requests describe the biological or workflow need first. The best requests are reusable across datasets, tissues, or cohorts.

Examples of broadly useful requests:

- support for another common 10x Multiome input shape
- clearer QC summaries
- better BPCells-native output handling
- a configurable but simple downstream module
- documentation for a common analysis path

## Pull Requests

Pull requests are welcome, but leanness is a high priority. Prefer a small, direct change over a flexible abstraction unless the abstraction removes real repeated complexity.

For code changes:

- keep target graph changes explicit and easy to inspect
- add pipeline-specific helper code under `R/` when target code would otherwise become hard to read or duplicated
- add generally reusable cross-repository helpers under `packages/multiomeRCore/R`; multiomeR sources this code directly, while standalone pipelines can install the nested package from a pinned commit
- add core target fragments under `extra_targets/`
- add optional downstream analyses under a `module_*` directory
- avoid compatibility layers for retired workflow shapes unless there is a clear public need
- prefer BPCells-native matrices, `GRanges`, tibbles, and file targets over creating large in-memory objects as intermediate state

## Validation

Before opening a pull request, run the smallest validation that matches the change:

```bash
pixi shell
R
```

```r
targets::tar_manifest(callr_function = NULL)
targets::tar_make(names = matches("healthy_PBMC_human|blood_human|brain_mouse"))
```

```bash
git diff --check
```

For target behavior changes, also run a narrow `targets::tar_make()` selection that exercises the affected target family.

For BPCells-native scoring-helper changes, run the synthetic parity check:

```bash
pixi run test-scoring-parity
```
