# Parked cell-type annotation perturbation diagnostics

The primary module's UCell cluster annotation used to compute three diagnostics
alongside each decision, without ever changing an assignment:

- marker stability: the fraction of leave-one-marker-out variants, each removing
  one candidate marker and its matched control, in which the candidate keeps a
  positive advantage of at least the configured minimum;
- cell stability: the fraction of leave-one-block-out replicates, with nuclei
  split into ten blocks stratified by cluster and GEM well, that assign the same
  candidate;
- GEM-well agreement: the fraction of per-well subgroups with at least 25 nuclei
  that assign it.

They were only written to review tables, so they now live here instead of in the
target graph. `helpers.R` keeps the original one-pass implementation on top of
the production helpers in `R/cluster_annotation_helpers.R`. For a built
aggregation, compute them from the stored counts, cluster metadata, and matched
controls:

```bash
pixi run --use-environment-activation-cache Rscript \
  dev/cluster_annotation_diagnostics/run_diagnostics.R my_aggregation GEX
```

The script writes `clusters.tsv`, with the stability and agreement columns, and
`GEM_well_agreement.tsv` below
`<store>/files/my_aggregation/cluster_UCell_perturbation_diagnostics/GEX/`. Use
`ATAC` or `WNN` for the other cluster partitions.

## Tested alternatives to pooled control matching

Each marker's controls are matched on abundance and detection pooled over up to
50 reference nuclei per GEM well, so wells of different depth score them
differently. In the version 1.1.0 demo, whose PBMC and lymph-node wells had
median 1,743 and 1,068 genes, the controls scored 0.51 in PBMC T nuclei and
0.27 in lymph-node T nuclei. No T-cell marker set then cleared the background
in both. Two alternatives were tested on 2026-09-29 and not adopted:

- Matching on each well's own abundance and detection fixed the demo's T cells
  and was neutral on a one-tissue aggregation of 52 wells. On
  `mixed_human_31x`, with many tissues, it halved the assigned WNN nuclei (194k
  to 100k) and lost clear labels such as cardiomyocytes and hepatocytes,
  because the controls inherit each marker's tissue pattern.
- Matching within depth bins of the reference nuclei (2, 3 or 5 bins) labelled
  the PBMC T cells but not the lymph-node ones.

A better scheme must separate technical from biological differences between
wells.
