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
