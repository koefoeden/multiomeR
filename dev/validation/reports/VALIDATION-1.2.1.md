# Version 1.2.1 release validation

Checked on 2026-10-12 on commit `8422c42f9ce902103d90fe18a817b06525e763aa`,
which completes the dated `NEWS.md`, with a clean working tree. The release
commit adds only this report.

The checkout's existing Pixi environment was used as is: `pixi.lock`,
`pixi.toml` and the pinned GitHub packages are unchanged since 1.2.0, and no
packages were installed or updated. This is not a fresh-clone install test.

## Checks on the release commit

- `validate-local` preflight: 10 of 10 checks passed in each run (pinned GitHub
  packages Seurat, BPCells, Signac and betterChromVAR; four Open Targets
  datasets; CollecTRI checksum and schema; AnnotationHub metadata and Ensembl
  resources).
- Unqualified `tar_make()` of the validation selection in the existing stores
  on one node: with `configuration/`, `immune_human_2x` and `brain_mouse`, all
  1,048 targets current; with `configuration_dev/` and `TAR_PROJECT=dev`,
  `mixed_human_31x`, all 6,205 targets current; no errors.
- Earlier runs that day on commit `1152385`, with the same pipeline files,
  rebuilt `brain_mouse` with its new marker sets and the demos' UMAP sweeps
  (443 targets in 14 min) and `mixed_human_31x`'s UMAP sweeps and TRS interval
  plot (22 targets in 52 min). The manuscript figure targets of
  `pipeline_manuscript_figures/targets.R` were rebuilt separately in the dev
  store during release preparation.
- Output checks passed: no errored or outdated targets; Seurat/Signac exports
  with RNA and ATAC assays and the configured GEM wells; and male-versus-female
  results for all four differential-analysis matrices and the cell-type
  composition model.
- Reference-parity tests (`pixi run test`): 70 passed, with no failures,
  warnings or skips.
- The manual rendered without warnings, `check-website-links` found no broken
  internal links, `check-source-links` found no errors in 47 pages, and the
  working tree stayed clean.
- Separately, the public demo built from an empty store in a 16-CPU Slurm job
  limited to 32 GiB; its processes used at most 17 GiB (details in `NEWS.md`).

| Aggregation | Nuclei | GEM wells |
|---|---|---|
| `immune_human_2x` | 7838 | 2 |
| `brain_mouse` | 2751 | 1 |
| `mixed_human_31x` | 213212 | 31 |

`immune_human_2x`'s WNN labels are unchanged from 1.2.0. `brain_mouse`'s WNN
nuclei are 876 intermediate progenitors, 674 upper-layer neurons, 444
inhibitory neurons, 329 radial glia, 196 deep-layer neurons and 126
hippocampal neurons; one small mixed cluster (106 nuclei) stays unassigned.
Its GEX clusters differ slightly between builds, because the multithreaded
nearest-neighbor search is not deterministic.

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

Reports: `outputs/validation/20261012T144409` and
`outputs_dev/validation/20261012T144544`.
