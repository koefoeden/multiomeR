# Version 1.4.0 release validation

Checked on 2026-10-30 on commit `ddb03315d3a026184f93bb4a10f18ead6e958ad6`
with a clean working tree. Later commits add only this report and the
fresh-installation fix below.

1.4.0 changes results for every aggregation after its per-well steps, so both
stores rebuilt most of their targets before the final checks. The final
`validate-local` runs used the checkout's existing Pixi environment as is,
installed from the release's lock (MACS3 3.1.0, BPCells 5f2a5d21,
betterChromVAR d64e3b5d). The `configuration_dev/` rebuilds ran in a new
environment installed from the same lock on node-local disk; the three
BSgenome data packages, whose post-link download did not run there, were
copied from the existing environment at the locked versions.

## Checks on the release commit

- `validate-local` preflight: 11 of 11 checks passed in each run.
- Unqualified `tar_make()` of the validation selection:
  - with `configuration/`, `immune_human_2x` and `brain_mouse`: at the
    frozen code (c12aeb3b), 805 targets built in 9 min on 16 cores, from a
    store holding only the demo's targets, and 192 were current. On the
    release commit, all 997 were current.
  - with `configuration_dev/` and `TAR_PROJECT=dev`, `mixed_human_31x`: two
    runs on node-local disk at 81b29d7a, which differs from the release
    commit only in documentation, rebuilt 3,913 targets for the manuscript
    results (with `comparison_5x`'s parity fixtures) in 59 min and the rest
    of the validation selection, 1,614 targets, in 17 min, with no errors.
    On the release commit, all 6,026 were current.
- Output checks passed: no errored or outdated targets; Seurat/Signac exports
  with RNA and ATAC assays and the configured GEM wells; and male-versus-female
  results for the differential-analysis matrices and the cell-type
  composition model.
- Reference-parity tests (`pixi run test`, with `configuration/` selected): 61
  passed, with no failures, warnings or skips. The validation-mode test
  expects `configuration/`; with `configuration_dev/` selected in
  `configuration.local`, three of its checks fail.
- The manual rendered, `check-website-links` found no broken internal links,
  `check-source-links` found no errors in 50 pages, and the working tree
  stayed clean.
- The demo command of `website/demo_running.md` built from an empty store in
  23.5 min on 16 cores, with nothing outdated and no errors.

| Aggregation | Nuclei | GEM wells |
|---|---|---|
| `immune_human_2x` | 7805 | 2 |
| `brain_mouse` | 3529 | 1 |
| `mixed_human_31x` | 214643 | 31 |

Inputs as for 1.2.0. Main commands as for 1.3.0, with
`MULTIOMER_CONFIGURATION` selecting each configuration directory.

Reports: `outputs/validation/20261030T121318` and
`outputs_dev/validation/20261030T122556`.

## Fresh installation

A fresh environment lacked the three BSgenome data packages (hg38, mm10,
mm39). Pixi ran their post-link scripts, but each script's `R CMD INSTALL`
read the project's `.Rprofile`, whose runtime check stopped because the
pinned BPCells was not yet installed. Pixi reports a failed post-link script
only as a warning that its default log level hides. Since b8f466c6, R CMD
sessions skip the project runtime.

A fresh clone of b8f466c6 on node-local disk, with Pixi 0.82.1 and no
global Pixi configuration, followed the README's two install commands.
`pixi install --locked --run-post-link-scripts` installed all three packages
in 3 min. They loaded in R before and after
`install-r-github-packages`, and with the project runtime; the working tree
stayed clean.
