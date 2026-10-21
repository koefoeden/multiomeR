# Version 1.3.1 release validation

Checked on 2026-10-21 on commit `484292fe6d1dc1f478566ec58a2a401631502555`
with a clean working tree. The release commit adds only this report.

1.3.1 changes results only for aggregations with several Harmony covariates or
a count model with a one-column design, which the validation aggregations do
not have. The environment is unchanged since 1.3.0. The checkout's existing
Pixi environment was used as is; this is not a fresh-clone install test.

## Checks on the release commit

- `validate-local` preflight: 11 of 11 checks passed in each run.
- Unqualified `tar_make()` of the validation selection in the existing stores:
  - with `configuration/`, `immune_human_2x` and `brain_mouse`: 8 targets
    rebuilt, the Harmony embeddings and Seurat/Signac exports of both demos,
    and 980 were current;
  - with `configuration_dev/` and `TAR_PROJECT=dev`, `mixed_human_31x`: 9
    targets rebuilt in 9 min, the Harmony embeddings, Seurat/Signac exports,
    differential-analysis models and their four feature fits, and 6,072 were
    current. No errors.

  These are the targets that the release's impact lines name. Every one
  reproduced its value, so nothing downstream reran.
- The rebuilt Seurat/Signac exports recorded Seurat's usual warnings: feature
  names with underscores get dashes, and the approximate neighbour search did
  not find 21 to 201 nuclei as their own neighbour.
- Output checks passed: no errored or outdated targets; Seurat/Signac exports
  with RNA and ATAC assays and the configured GEM wells; and male-versus-female
  results for all four differential-analysis matrices and the cell-type
  composition model.
- Reference-parity tests (`pixi run test`): 66 passed, with no failures,
  warnings or skips.
- The manual rendered without warnings, `check-website-links` found no broken
  internal links, `check-source-links` found no errors in 48 pages, and the
  working tree stayed clean.

| Aggregation | Nuclei | GEM wells |
|---|---|---|
| `immune_human_2x` | 7887 | 2 |
| `brain_mouse` | 3532 | 1 |
| `mixed_human_31x` | 213885 | 31 |

Inputs as for 1.2.0. Main commands as for 1.3.0.

Reports: `outputs/validation/20261021T140210` and
`outputs_dev/validation/20261021T140847`.
