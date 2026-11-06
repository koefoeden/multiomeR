# Version 1.4.2 release validation

Checked on 2026-11-06 on commit `13203a034d23598bf5ea3b6d025e73f85979e05f`
with a clean working tree. The later commit adds only this report.

1.4.2 changes only the GWAS studies of `mixed_human_31x` in
`configuration_dev/`, so only that store rebuilt, and only its genetic
enrichment. The pipeline code, `configuration/` and the lock are unchanged
since 1.4.1, and all runs used the environment of the 1.4.1 validation.

## Checks on the release commit

- `validate-local` preflight: 11 of 11 checks passed in each passing run.
  Two earlier runs with `configuration/` stopped at preflight, before the
  pipeline started, because ftp.ebi.ac.uk did not finish its TLS handshakes
  within the 15-second connect timeout. The third run passed.
- Unqualified `tar_make()` of the validation selection:
  - with `configuration/`, `immune_human_2x` and `brain_mouse`: 996 of 997
    targets were current. `GEM_well_config_tsv` reran because the store was
    last built from another checkout path.
  - with `configuration_dev/` and `TAR_PROJECT=dev`, `mixed_human_31x`: on
    the release code, the downloads and S-LDSC fits of the seven new studies
    (37 targets) took 7 min, and the 270 genetic-enrichment targets after
    them 15 min on its Slurm workers, with no errors. The S-LDSC fits of the
    17 unchanged studies stayed current. The targets of the manuscript's
    Figure 1, S5 and S6 then rebuilt. On the release commit, all 5,938 were
    current, as were the figure targets and everything upstream of them.
- Output checks passed: no errored or outdated targets; Seurat/Signac exports
  with RNA and ATAC assays and the configured GEM wells; and male-versus-female
  results for the differential-analysis matrices and the cell-type
  composition model.
- Reference-parity tests (`pixi run test`, with `configuration/` selected): 61
  passed, with no failures, warnings or skips.
- The manual rendered, `check-website-links` found no broken internal links,
  `check-source-links` found no errors in 50 pages, and the working tree
  stayed clean.

| Aggregation | Nuclei | GEM wells |
|---|---|---|
| `immune_human_2x` | 7805 | 2 |
| `brain_mouse` | 3529 | 1 |
| `mixed_human_31x` | 214643 | 31 |

Inputs as for 1.2.0. Main commands as for 1.3.0, with
`MULTIOMER_CONFIGURATION` selecting each configuration directory.

Reports: `outputs/validation/20261106T164223` and
`outputs_dev/validation/20261106T161508`.
