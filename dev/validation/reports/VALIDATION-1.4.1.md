# Version 1.4.1 release validation

Checked on 2026-11-06 on commit `db2af427b7f87c01ad1c1e0bc519ba39f6ca2be4`
with a clean working tree. The later commit adds only this report.

1.4.1 reruns every UMAP and the outputs that hold one, and the genetic-enrichment
tables without the most-accessible-peaks strategy. Both stores were rebuilt on
the release code before the final checks. All runs used an environment with
the same 623 linux-64 packages as this release's lock and the same pinned
GitHub packages. It lacked susieR and its two dependencies, which no target
uses since 1.4.0, and held packages that the pipeline does not load. The lock
is unchanged since 1.4.0.

## Checks on the release commit

- `validate-local` preflight: 11 of 11 checks passed in each run.
- Unqualified `tar_make()` of the validation selection:
  - with `configuration/`, `immune_human_2x` and `brain_mouse`: on the
    release code, 374 targets rebuilt in 8.1 min on 4 cores with no errors,
    all of them UMAP embeddings, the plots and tables that hold UMAP
    coordinates, and the Seurat/Signac exports. Targets that read those tables
    but not the coordinates, such as cluster composition tables and violin
    plots, reran with identical results. On the release commit, all 997 were
    current.
  - with `configuration_dev/` and `TAR_PROJECT=dev`, `mixed_human_31x`: on the
    release code, 382 targets rebuilt in 14 min on its Slurm workers with no
    errors: the same UMAP closure, and the genetic-enrichment tables and plots
    without the removed strategy. The targets of the manuscript's Figure 1, S5
    and S6 then rebuilt after each figure edit, in 2 min or less. On the
    release commit, all 5,998 were current.
- Output checks passed: no errored or outdated targets; Seurat/Signac exports
  with RNA and ATAC assays and the configured GEM wells; and male-versus-female
  results for the differential-analysis matrices and the cell-type
  composition model.
- The rebuilt WNN UMAP of `mixed_human_31x` equals the 1,000-epoch layout
  evaluated for the fix (largest coordinate difference 0).
- Reference-parity tests (`pixi run test`, with `configuration/` selected): 61
  passed, with no failures, warnings or skips.
- On `ee0fe65c`, whose manual the later figure commits leave unchanged, the
  manual rendered, `check-website-links` found no broken internal links,
  `check-source-links` found no errors in 50 pages, and the working tree
  stayed clean.

| Aggregation | Nuclei | GEM wells |
|---|---|---|
| `immune_human_2x` | 7805 | 2 |
| `brain_mouse` | 3529 | 1 |
| `mixed_human_31x` | 214643 | 31 |

Inputs as for 1.2.0. Main commands as for 1.3.0, with
`MULTIOMER_CONFIGURATION` selecting each configuration directory.

Reports: `outputs/validation/20261106T140339` and
`outputs_dev/validation/20261106T140452`.
