# Output files and metadata

The [running guide](main_running.md#steps) lists the files to open after each
checkpoint. Read each plot's subtitle and caption for interpretation and
method details.

## Output folders and metric selection {#output-folders}

Plots live under `<store>/plots/<scope>/<checkpoint>/`. The store comes from
`_targets.yaml`; scope is a GEM well or aggregation. Existing files are not a
record of which targets are current: use `targets::tar_outdated()` before
reviewing results after a configuration change.

`QC_metric_manifest.tsv` selects metrics, display labels and plotting quantiles.
`do_plot = FALSE` hides a metric; plotting quantiles change the displayed range,
not the cells retained by the pipeline. Filtering is configured separately in
`cfg_GEM_wells.tsv` and `cfg_aggregations.yaml`.

GEX, ATAC and WNN parameter sweeps belong to checkpoints 3, 7 and 8.
The GEX compatibility object belongs to checkpoint 3 and the multimodal object
to checkpoint 8. Peak–gene correlation is a separate optional module.

## Cell retention tables {#cell-retention}

Append `.my_aggregation` to these target names when reading them with
`targets::tar_read()`:

| Target | Checkpoint |
|---|---|
| `cell_retention_tibble.GEX_input` | 1 |
| `cell_retention_tibble.GEX` | 3 |
| `cell_retention_tibble.ATAC_input` | 5 |
| `cell_retention_tibble.ATAC` | 7 |
| `cell_retention_tibble.WNN` | 8 |

Tables accumulate stages and retain wells with zero surviving cells. Rows identify
an exclusion action and its order, input count, excluded count and retained count.

## Cluster numbering {#cluster-numbering}

Leiden clusters receive size-ranked IDs before minimum-size filtering. Filtering
does not renumber survivors, so IDs can have gaps.

## Cluster annotation {#cluster-annotation}

Cluster labelling uses adjusted-score advantage over matched random controls.
Marker lists accept unsigned genes, positive `+` suffixes and negative `-` suffixes.
Signed UCell scores are clipped at zero per cell before averaging, for observed
and matched-control signatures alike. The method uses the same
GEX control reference for GEX, ATAC and WNN clusters, and returns `Assigned`,
or `Unassigned` annotation status. The adjusted score is the mean UCell score
minus that label's matched-control 95th percentile. The highest adjusted score
nominates the candidate. Set `aggregation_cluster_annotation_min_advantage`
(default **0.05**) to the required lead over both zero background and the
next-best adjusted score. Exact ties and scores at or below background remain
unassigned even at a zero threshold. Raising the threshold only withdraws
assignments; it never switches the nominated label. An unassigned cluster
retains its candidate and a reason (weak background advantage or competition)
in the diagnostics. This status does not remove cells. Mixture detection is
not included. The former module-score labeller and its
`aggregation_allow_multiple_cell_types` setting have been removed. Remove that
setting and `aggregation_cluster_annotation_method` from older configuration
files; there is now one annotation method. Signed marker panels remain supported.

Annotation exports include `clusters.tsv`, `marker_evidence.tsv`,
`control_gene_matching.tsv` and, where eligible, `GEM_well_agreement.tsv`.
The `cluster_UCell_evidence` targets cache scoring separately from the
threshold-dependent annotation. `marker_set_UCell_summary.3_GEX_QC` exports
`marker_sets.tsv` and `method.txt` before GEX doublet filtering.

## Further methods

See [algorithm validation](implementation/algorithm_validation.html) for reference
comparisons and deviations, and the [peak–gene module](downstream_peak_gene_correlation.md)
for its model, support filters and output paths.
