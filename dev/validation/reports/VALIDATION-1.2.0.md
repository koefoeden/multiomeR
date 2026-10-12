# Version 1.2.0 release validation

Checked on 2026-10-12 on commit `e255eee22e0ea73e39c51286a8e11b707ff01f1e`,
which completes the dated `NEWS.md`, with a clean working tree. The release
commit adds only this report.

The checkout's existing Pixi environment was used as is: `pixi.lock` is
unchanged since 1.1.0, `pixi.toml` changes only task paths, and the Seurat 5.6.0
pin added to `scripts/github_packages.R` was already installed. No packages were
installed or updated. This is not a fresh-clone install test.

## Checks on the release commit

- `validate-local` preflight: 10 of 10 checks passed in each run (pinned GitHub
  packages Seurat, BPCells, Signac and betterChromVAR; four Open Targets
  datasets; CollecTRI checksum and schema; AnnotationHub metadata and Ensembl
  resources).
- Unqualified `tar_make()` of the validation selection in the existing stores
  on one node: with `configuration/`, `immune_human_2x` and `brain_mouse`, all
  1,041 targets current; with `configuration_dev/` and `TAR_PROJECT=dev`,
  `mixed_human_31x`, all 6,205 targets current; no errors.
- Earlier runs that day rebuilt both demos with the revised methods (848
  targets in 30 min, on commits `6ac6637` and `e6910c8`; `mixed_human_31x` had been
  rebuilt during release preparation). They found three defects, fixed before
  this commit: cluster-marker targets that branched over an empty table in
  `brain_mouse` (`7237596`), a reference check that compared stored paths from
  before the move of the reference metadata (`14c4520`), and `donor_id` plotted
  twice in the demo, which left one categorical UMAP target outdated
  (`8a66925`).
- Output checks passed: no errored or outdated targets; Seurat/Signac exports
  with RNA and ATAC assays and the configured GEM wells; and male-versus-female
  results for all four differential-analysis matrices and the cell-type
  composition model.
- Reference-parity tests (`pixi run test`): 70 passed, with no failures,
  warnings or skips.
- The manual rendered without warnings, `check-website-links` found no broken
  internal links, `check-source-links` found no errors in 47 pages, and the
  working tree stayed clean.

| Aggregation | Nuclei | GEM wells |
|---|---|---|
| `immune_human_2x` | 7838 | 2 |
| `brain_mouse` | 2734 | 1 |
| `mixed_human_31x` | 213212 | 31 |

`immune_human_2x`'s WNN nuclei are 4,942 T, 1,417 Mono, 331 B, 319 NK and 104
MAIT; three small clusters (725 nuclei) match no marker set above background
and stay unassigned. `brain_mouse`, an embryonic brain annotated with adult
marker sets, leaves all seven WNN clusters unassigned; it is an inactive example
and its markers are revised after this release.

Inputs: the downloaded 10x Genomics demo wells, `healthy_PBMC_human` and
`controller_PBMC_human` (the 10k PBMC Chromium Controller output, linked into
`example_data/`), and locally processed ARC 2.1.0 / GRCh38-2024-A counts for the
24 ENCODE and seven reprocessed 10x wells of `mixed_human_31x`. Raw-read
processing was not part of this test.

Main commands (from the validated checkout):

```sh
pixi run --frozen --use-environment-activation-cache validate-local
TAR_PROJECT=dev pixi run --frozen --use-environment-activation-cache validate-local  # with configuration_dev selected
pixi run --frozen --use-environment-activation-cache test
pixi run --frozen --use-environment-activation-cache -e dev render-website
pixi run --frozen --use-environment-activation-cache -e dev check-website-links
pixi run --frozen --use-environment-activation-cache -e dev check-source-links
```

Reports: `outputs/validation/20261012T100741` and
`outputs_dev/validation/20261012T100855`.
