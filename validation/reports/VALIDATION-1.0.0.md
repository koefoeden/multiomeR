# Version 1.0.0 release validation

Checked on 2026-09-29 on commit `a0da7008b1768e410ce46128d40d6745a4656eb4`,
which dates `NEWS.md`, with a clean working tree. The release commit adds only
this report.

The checkout's own Pixi environment was installed from `pixi.lock` with
`pixi install --locked --run-post-link-scripts`, updating that checkout's
earlier environment, and `install-r-github-packages` rebuilt BPCells, Signac
and betterChromVAR at their pinned revisions. No other packages were installed.
This is not a fresh-clone install test.

## Checks on the release commit

- `validate-local` preflight: 9 of 9 checks passed (pinned GitHub packages,
  Open Targets datasets, CollecTRI checksum and schema, AnnotationHub metadata
  and Ensembl resources).
- Unqualified `tar_make()` of `immune_human_2x`, `brain_mouse` and
  `mixed_human_31x` in the existing store on one 16-core node: all 7,625 targets
  current, no errors. An earlier run of the same files that day rebuilt six
  targets, the Open Targets and CollecTRI downloads and the GEM-well
  configuration file, whose recorded paths belonged to the checkout that built
  the store; their contents were unchanged and no dependent target reran.
- Output checks passed: no errored or outdated targets; Seurat/Signac exports
  with RNA and ATAC assays and the configured GEM wells; and male-versus-female
  results for all four differential-analysis matrices and the cell-type
  composition model.
- Reference-parity tests (`pixi run test`): 71 passed, with no failures,
  warnings or skips.
- Both documentation books rendered, `check-source-links` found no errors in 38
  pages, and the working tree stayed clean.

| Aggregation | Nuclei | GEM wells |
|---|---|---|
| `immune_human_2x` | 9785 | 2 |
| `brain_mouse` | 3534 | 1 |
| `mixed_human_31x` | 197330 | 31 |

## Earlier runs that built the store

Both used an existing environment whose lockfile adds five packages to this one.

- Commit `42de65127b33aa00fcf7638ecfcf7eb2dc075173` (2026-09-24/25), from a fresh
  store with the public local controllers on one 16-core node: 9 h 32 min,
  6,639 targets completed, no errors, and no outdated targets afterwards.
- Commit `830ccd5e3c3ae45d78bf708dda526d1931653481` (2026-09-26), after the
  per-GEM-well ATAC counts, the compiled UCell kernel, image-only plot saving
  and the gallery refresh: 1 h 0 min, 1,296 targets rebuilt and 6,329 reused,
  no errors, with unchanged nuclei and GEM-well counts.

Between `830ccd5` and the release commit, only `scripts/github_packages.R`
(an unused package removed) and `NEWS.md` changed.

Inputs: the downloaded 10x Genomics demo wells, and locally processed ARC 2.1.0 /
GRCh38-2024-A counts for the 24 ENCODE and seven reprocessed 10x wells of
`mixed_human_31x`. Raw-read processing was not part of this test.

Main commands (from the validated checkout):

```sh
pixi install --locked --run-post-link-scripts
pixi run --frozen --use-environment-activation-cache install-r-github-packages
pixi run --frozen --use-environment-activation-cache validate-local
pixi run --frozen --use-environment-activation-cache test
pixi install -e dev --locked
pixi run --frozen --use-environment-activation-cache -e dev render-website
pixi run --frozen --use-environment-activation-cache -e dev check-source-links
```
