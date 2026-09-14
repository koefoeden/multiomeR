# Local incremental validation

The `ci` branch enables the human and mouse demos plus six public ENCODE
left-ventricle samples. It uses the standard configs, normal `tar_make()` and
unchanged `_targets.yaml`. Switch branches in the same checkout and keep the
same `outputs` store: targets reuses valid results and reruns affected work.

## Inputs

Download the demo inputs using the normal walkthrough if they are not present:

```sh
bash example_data/download_10X_cellranger_count_data.sh
```

Create an ignored local link named `example_data/encode` to the directory
containing the six ENCODE count directories named by accession. These are
locally processed public data: ARC 2.1.0, GRCh38 2024-A, with each count's
`outs/cellbender-output_gex_bent_only_filtered.h5`. Raw-read processing is not
part of this test. Public donor metadata and source URLs are included here.

## Validate a candidate

Only switch branches when the working tree is clean and no pipeline is running.
From the candidate branch, save its exact commit, then merge it into `ci`:

```sh
candidate=$(git rev-parse HEAD)
git switch ci
git merge --no-edit "$candidate"
```

Resolve any config conflicts while retaining the three CI aggregations. Review
`git diff "$candidate" HEAD`: only CI configs, metadata, checks and this guide
should differ. Keep pipeline implementation changes on the development branch;
do not merge `ci` into `main`.

On a suitably allocated compute node, run:

```sh
pixi run --use-environment-activation-cache Rscript -e 'targets::tar_make()'
pixi run --use-environment-activation-cache Rscript validation/check.R "$candidate"
```

Run the second command only after the first exits successfully. It checks all
three aggregation outputs, all four ENCODE differential modalities, the sex
composition result, current target errors and freshness. A passing report in
`outputs/validation/<CI-commit>/summary.tsv` records both commits. Save the console
log alongside it if needed. No timer, GitHub status service or extra repository
is required. Switching back leaves the store and downloaded inputs in place.

This validates the current dependency graph incrementally; it is not a clean
rebuild or an installation test. External changes can only invalidate work if
they are tracked by targets. Do not use `shortcut = TRUE` or manually mark
outdated targets current to obtain a passing report.

## ENCODE comparison

The model tests male versus female in cardiomyocytes and sex-associated cell
composition. There are four female and two male donors, with imbalanced ages
and health status; this is an integration test and an unadjusted association,
not a causal sex-effect estimate. It retains all eligible features and the
source marker sets, using public multiomeR's annotation method.

Source configuration: `multiomeR-CBMR` commit
`624a72f50f56e2b1355fb48fc8d9edbced6e088d`; migrated through validation profile
commit `18864f976b87d990e67e678df7e579db589764b9`. The superseded profile repository
contains historical private configuration and must remain private.
