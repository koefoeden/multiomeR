# Version 1.1.0 release validation

Checked on 2026-09-29 on commit `7f1d70802e7cb406ce56d0f2d34909a6afbac48a`,
which completes the dated `NEWS.md`, with a clean working tree. The release
commit adds only this report.

The checkout's existing Pixi environment was used as is: `pixi.lock` is
unchanged since 1.0.0, and no packages were installed or updated. This is not a
fresh-clone install test.

## Checks on the release commit

- `validate-local` preflight: 9 of 9 checks passed (pinned GitHub packages, Open
  Targets datasets, CollecTRI checksum and schema, AnnotationHub metadata and
  Ensembl resources). Earlier attempts that day, one of them on this commit,
  failed only Open Targets checks, because `ftp.ebi.ac.uk` intermittently
  refused HTTPS connections from the cluster; the release candidate therefore
  retries each live check up to four times, also after refused connections.
- Unqualified `tar_make()` of `immune_human_2x`, `brain_mouse` and
  `mixed_human_31x` in the existing store on one node: all 7,633 targets
  current, no errors. An earlier run of the validation selection that day, on
  commit `f0f02440ec0ed51a92998408baea091c29f48493` with the same pipeline
  files, rebuilt 461 targets and reused 7,172 in 15 min with no errors: the new
  demo aggregation and its GEM wells, and the GEM-well configuration projections
  of the other aggregations, whose values were unchanged. Before it, the stored
  BPCells objects of `healthy_PBMC_human`, whose paths belonged to the checkout
  that built the store, were invalidated so they reopened from this checkout.
- Output checks passed: no errored or outdated targets; Seurat/Signac exports
  with RNA and ATAC assays and the configured GEM wells; and male-versus-female
  results for all four differential-analysis matrices and the cell-type
  composition model.
- Reference-parity tests (`pixi run test`): 71 passed, with no failures,
  warnings or skips.
- The manual rendered without warnings, `check-website-links` found no broken
  internal links, `check-source-links` found no errors in 47 pages, and the
  working tree stayed clean.

| Aggregation | Nuclei | GEM wells |
|---|---|---|
| `immune_human_2x` | 2248 | 2 |
| `brain_mouse` | 3534 | 1 |
| `mixed_human_31x` | 197330 | 31 |

Every GEX, ATAC and WNN cluster of `immune_human_2x` received a cell-type
label; its WNN nuclei are 1,404 T, 598 Mono, 136 B and 110 NK.

Inputs: the downloaded 10x Genomics demo wells, of which `unsorted_PBMC_human`
is the published `pbmc_unsorted_3k` output linked into `example_data/`, and
locally processed ARC 2.1.0 / GRCh38-2024-A counts for the 24 ENCODE and seven
reprocessed 10x wells of `mixed_human_31x`. Raw-read processing was not part of
this test.

Main commands (from the validated checkout):

```sh
pixi run --frozen --use-environment-activation-cache validate-local
pixi run --frozen --use-environment-activation-cache test
pixi run --frozen --use-environment-activation-cache -e dev render-website
pixi run --frozen --use-environment-activation-cache -e dev check-website-links
pixi run --frozen --use-environment-activation-cache -e dev check-source-links
```
