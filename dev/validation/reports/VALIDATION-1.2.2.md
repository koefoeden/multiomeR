# Version 1.2.2 release validation

Checked on 2026-10-14 on commit `1f45018fda644b346e492b81d602365507b725d8`
with a clean working tree. The release commit adds only this report.

1.2.2 changes no pipeline function or target. In the R code it adds only
comments, and every function of the 12 commented files deparses identically
before and after. It changes `configuration_dev/` (the islet endocrine WNN
signature and 26 instead of 39 GWAS for `mixed_human_31x`), the manuscript
figure code and results, the manual and the repository documentation. The
checkout's existing Pixi environment was used as is: `pixi.lock`, `pixi.toml`
and the pinned GitHub packages are unchanged since 1.2.0, and no packages were
installed or updated. This is not a fresh-clone install test.

## Checks on the release commit

- `validate-local` preflight: 10 of 10 checks passed in each run (pinned GitHub
  packages Seurat, BPCells, Signac and betterChromVAR; four Open Targets
  datasets; CollecTRI checksum and schema; AnnotationHub metadata and Ensembl
  resources). R printed two "stack imbalance" warnings while the preflight
  loaded the pinned packages; the checks themselves passed.
- Unqualified `tar_make()` of the validation selection in the existing stores
  on one node: with `configuration/`, `immune_human_2x` and `brain_mouse`, all
  1,048 targets current; with `configuration_dev/` and `TAR_PROJECT=dev`,
  `mixed_human_31x`, all 6,164 targets current; no errors.
- Earlier, on 2026-10-13 on commit `f43477b`, with the same pipeline and
  configuration files, `mixed_human_31x` rebuilt its WNN labels, genetic
  enrichment and everything downstream (1,390 targets in 10 min on node-local
  scratch). The manuscript targets of `pipeline_manuscript_figures/targets.R`
  were rebuilt in the dev store as the figure code changed, last on commit
  `f466ef0`. The `comparison_*` benchmark aggregations keep their measured
  runs; Supplementary Figure S2 and its frozen resource table were redrawn
  without rerunning them.
- Output checks passed: no errored or outdated targets; Seurat/Signac exports
  with RNA and ATAC assays and the configured GEM wells; and male-versus-female
  results for all four differential-analysis matrices and the cell-type
  composition model.
- Reference-parity tests (`pixi run test`): 70 passed, with no failures,
  warnings or skips.
- The manual rendered without warnings, `check-website-links` found no broken
  internal links, `check-source-links` found no errors in 47 pages, and the
  working tree stayed clean.

| Aggregation | Nuclei | GEM wells | WNN cell types |
|---|---|---|---|
| `immune_human_2x` | 7838 | 2 | unchanged |
| `brain_mouse` | 2751 | 1 | unchanged |
| `mixed_human_31x` | 213212 | 31 | 19 (islet endocrine added) |

Inputs as for 1.2.0: the downloaded 10x Genomics demo wells and locally
processed ARC 2.1.0 / GRCh38-2024-A counts for the 24 ENCODE and seven
reprocessed 10x wells of `mixed_human_31x`. Raw-read processing was not part
of this test.

Main commands (from the validated checkout):

```sh
pixi run --frozen --use-environment-activation-cache validate-local
TAR_PROJECT=dev pixi run --frozen --use-environment-activation-cache validate-local  # with configuration_dev selected
pixi run --frozen --use-environment-activation-cache test
pixi run --frozen --use-environment-activation-cache -e dev render-website
pixi run --frozen --use-environment-activation-cache -e dev check-website-links
pixi run --frozen --use-environment-activation-cache -e dev check-source-links
```

Reports: `outputs/validation/20261014T124945` and
`outputs_dev/validation/20261014T125121`.
