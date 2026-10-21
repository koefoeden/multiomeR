# Version 1.3.0 release validation

Checked on 2026-10-21 on commit `e923d26dd568823d849e6af53248de3382ae4758`
(tag `v1.3.0`) with a clean working tree. The release was published while this
validation ran, so this report was added to `main` afterwards and is not part
of the tag.

1.3.0 changes results, as its release notes describe. Since 1.2.2, `pixi.lock`
changes only the CA certificate bundle, `pixi.toml` adds `requires-pixi`, and
the pinned GitHub packages are unchanged. The checkout's existing Pixi
environment was used as is; this is not a fresh-clone install test.

## Checks on the release commit

- `validate-local` preflight: 11 of 11 checks passed in each run (pinned GitHub
  packages Seurat, BPCells, Signac and betterChromVAR; four Open Targets
  datasets; the S-LDSC reference archive on Zenodo; CollecTRI checksum and
  schema; AnnotationHub metadata and Ensembl resources).
- Unqualified `tar_make()` of the validation selection in the existing stores:
  - with `configuration/`, `immune_human_2x` and `brain_mouse`: all 988
    targets current;
  - with `configuration_dev/` and `TAR_PROJECT=dev`, `mixed_human_31x`: 91
    targets rebuilt in 34 min and 5,990 were current, with no errors. The
    rebuilt targets are the cell-type pseudobulk counts of genetic enrichment
    and the S-LDSC chain after them. They follow the S-LDSC fixes and the
    100-nucleus cell-type minimum committed on 2026-10-19 and 2026-10-20,
    after the store's last build.
- Four rebuilt targets recorded expected warnings. Three gene-annotation
  branches of `SLDSC_matched_enrichment_records` combined ranges without some
  chromosomes. `SLDSC_matched_size_comparison_plot` dropped six points with
  zero values from its log scale.
- Output checks passed: no errored or outdated targets; Seurat/Signac exports
  with RNA and ATAC assays and the configured GEM wells; and male-versus-female
  results for all four differential-analysis matrices and the cell-type
  composition model.
- Reference-parity tests (`pixi run test`): 66 passed, with no failures,
  warnings or skips.
- The manual rendered without warnings, `check-website-links` found no broken
  internal links, `check-source-links` found no errors in 48 pages, and the
  working tree stayed clean.

| Aggregation | Nuclei | GEM wells | Nuclei in 1.2.2 |
|---|---|---|---|
| `immune_human_2x` | 7887 | 2 | 7838 |
| `brain_mouse` | 3532 | 1 | 2751 |
| `mixed_human_31x` | 213885 | 31 | 213212 |

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

Reports: `outputs/validation/20261021T123639` and
`outputs_dev/validation/20261021T123812`.
