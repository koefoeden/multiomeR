# Output files and metadata

The [running guide](main_running.md#steps) lists the files to open after each checkpoint. Read each plot's subtitle and caption for interpretation and method details.

## Output folders and metric selection {#output-folders}

Plots live under `<store>/plots/<scope>/<checkpoint>/`. The store comes from `_targets.yaml`; scope is a GEM well or aggregation. Existing files are not a record of which targets are current: use `targets::tar_outdated()` before reviewing results after a configuration change.

`QC_metric_manifest.tsv` selects metrics, display labels and plotting quantiles. `do_plot = FALSE` hides a metric; plotting quantiles change the displayed range, not the cells retained by the pipeline. Filtering is configured separately in `cfg_GEM_wells.tsv` and `cfg_aggregations.yaml`.

GEX, ATAC and WNN parameter sweeps belong to checkpoints 3, 7 and 8. The GEX compatibility object belongs to checkpoint 3 and the multimodal object to checkpoint 8. Peak–gene correlation is a separate optional module.

## Cell retention tables {#cell-retention}

Append `.my_aggregation` to these target names when reading them with `targets::tar_read()`:

| Target | Checkpoint |
|---|---|
| `cell_retention_tibble.GEX_input` | 1 |
| `cell_retention_tibble.GEX` | 3 |
| `cell_retention_tibble.ATAC_input` | 5 |
| `cell_retention_tibble.ATAC` | 7 |
| `cell_retention_tibble.WNN` | 8 |

Tables accumulate stages and retain wells with zero surviving cells. Rows identify an exclusion action and its order, input count, excluded count and retained count.

## Cluster numbering {#cluster-numbering}

Leiden clusters receive size-ranked IDs before minimum-size filtering. Filtering does not renumber survivors, so IDs can have gaps.

## Cluster annotation {#cluster-annotation}

Cluster labelling uses adjusted-score advantage over matched random controls. Marker lists accept unsigned genes, positive `+` suffixes and negative `-` suffixes. The same GEX control reference annotates GEX, ATAC and WNN clusters, and each cluster receives an `Assigned` or `Unassigned` status; an unassigned cluster retains its candidate label and the reason in the diagnostics, and no cells are removed. Set [`aggregation_cluster_annotation_min_advantage`](parameters.html#aggregation_cluster_annotation_min_advantage) to the required lead over both background and the next-best label; raising it only withdraws assignments. The scoring rule, control construction and all fixed constants are documented in [Cell-type annotation and motif accessibility](implementation/methods_annotation_and_motifs.html). The former module-score labeller and its `aggregation_allow_multiple_cell_types` and `aggregation_cluster_annotation_method` settings have been removed; delete them from older configuration files.

Annotation exports include `clusters.tsv`, `marker_evidence.tsv`, `control_gene_matching.tsv` and, where eligible, `GEM_well_agreement.tsv`. The `cluster_UCell_evidence` targets cache scoring separately from the threshold-dependent annotation. `marker_set_UCell_summary.3_GEX_QC` exports `marker_sets.tsv` and `method.txt` before GEX doublet filtering.

## Further methods

See [algorithm validation](implementation/algorithm_validation.html) for reference comparisons and deviations, and the [peak–gene module](downstream_peak_gene_correlation.md) for its model, support filters and output paths.
