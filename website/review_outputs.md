# Output reference

This page collects facts about the saved outputs. For a guided first look, see [the demo results](demo_outputs.md); for what to review at each checkpoint, see the [running guide](main_running.md#steps). Plot subtitles and captions explain how to read each plot, and the [output gallery](gallery.md) shows an example of each.

## Output folders and metric selection {#output-folders}

Plots and files are saved under `<store>/plots/` and `<store>/files/`, in folders given by the target name read from right to left. For example, `cross.UMAPs.3_GEX_QC.my_aggregation` saves to `<store>/plots/my_aggregation/3_GEX_QC/UMAPs/cross/`. The first folder is the GEM well or aggregation. Existing files do not show which results are current: after a configuration change, run `targets::tar_outdated()` before reviewing them.

`QC_metric_manifest.tsv` in the repository root selects the plotted QC metrics, their display labels and plotting quantiles. `do_plot = FALSE` hides a metric. Plotting quantiles change the displayed range, not the nuclei retained; filters are set in `cfg_GEM_wells.tsv` and `cfg_aggregations.yaml`.

`UMAPs/cross/` at checkpoints 3, 7 and 8 redraws the cell-type UMAP for several numbers of dimensions and neighbours (neighbours only at checkpoint 8). Clusters and labels stay fixed, so these plots show whether the layout depends on the UMAP settings.

## Seurat objects {#seurat-objects}

Two targets export Seurat objects for work outside the pipeline:

- `GEX_Seurat_object.3_GEX_QC.my_aggregation`: GEX data only, available after checkpoint 3, before peak calling.
- `multimodal_Seurat_object.8_multimodal_QC.my_aggregation`: the final Seurat/Signac object with GEX, ATAC and WNN results.

Both keep the raw GEX counts in the `RNA` assay with a log-normalized `data` layer (counts scaled to 10,000 per cell, then log1p), computed lazily by BPCells. The Pearson residuals that the GEX PCA used, with [`aggregation_SCT_regress_vars`](parameters.html#aggregation_SCT_regress_vars) regressed out, are stored according to [`aggregation_GEX_PCA_backend`](parameters.html#aggregation_GEX_PCA_backend):

| Backend | GEX assay | Residuals |
|---|---|---|
| `BPCells_native` | `RNA` | `scale.data` of `RNA`: the variable genes, clipped to [−10, 10], computed lazily |
| `Seurat_SCT` | `SCT` | the `SCTransform()` assay, which also holds corrected counts and their log1p `data` |

The GEX assay is the default assay and carries the GEX reductions and graphs; the `misc$normalization` entry of each assay records how its layers were computed. The cell metadata also holds the Seurat cell-cycle scores `S.Score`, `G2M.Score`, `Phase` and `CC.Difference`, computed from the log-normalized counts with the 2019 Seurat gene sets whether or not they are regressed.

## Cell retention tables {#cell-retention}

Append `.my_aggregation` to these target names when reading them with `targets::tar_read()`. The `cell_retention_flow_plot.png` of each checkpoint draws its table.

| Target | Checkpoint |
|---|---|
| `cell_retention_tibble.GEX_input` | 1 |
| `cell_retention_tibble.GEX` | 3 |
| `cell_retention_tibble.ATAC_input` | 5 |
| `cell_retention_tibble.ATAC` | 7 |
| `cell_retention_tibble.WNN` | 8 |

Each row is one GEM well and one exclusion action (`discard_action`, in `action_order` within its `stage`), with `input_cells`, `excluded_cells`, `retained_cells` and `retained_fraction`. A nucleus matching several actions counts under the first. Each table includes the stages of the tables before it, and GEM wells with no remaining nuclei keep their rows.

## Cluster numbering {#cluster-numbering}

Leiden clusters are numbered by size, starting with 1 for the largest, before clusters smaller than [`aggregation_cluster_min_barcodes`](parameters.html#aggregation_cluster_min_barcodes) are removed. The remaining clusters keep their numbers, so IDs can have gaps.

## Cluster annotation {#cluster-annotation}

GEX, ATAC and WNN clusters are all labelled from GEX data. Each cluster is scored against every marker set in [`aggregation_GEX_marker_genes`](parameters.html#aggregation_GEX_marker_genes), relative to random control genes matched for expression level and detection rate. A cluster is `Assigned` its best-scoring label when that label leads both the control background and the next-best label by at least [`aggregation_cluster_annotation_min_advantage`](parameters.html#aggregation_cluster_annotation_min_advantage); otherwise it is `Unassigned`. Raising the threshold can only withdraw assignments, and unassigned clusters keep their nuclei. The scores are cached separately, so changing the threshold does not recompute them. The scoring rule, control construction and fixed constants are described in [Cell-type annotation and motif accessibility](implementation/methods_annotation_and_motifs.html).

In the cell metadata, the cluster IDs are in `PCA_harmony_SNN_cluster` (GEX), `LSI_harmony_SNN_cluster` (ATAC) and `WNN_harmony_SNN_cluster` (WNN). Each has companion columns with the suffixes `_cell_type` (the label, or `Unassigned`), `_named` (cluster ID and label) and `_annotation_status` (`Assigned` or `Unassigned`).

The `cluster_UCell_diagnostics` targets write these files to `<store>/files/my_aggregation/3_GEX_QC/cluster_UCell_diagnostics/` and to the matching folders under `7_ATAC_QC/` and `8_multimodal_QC/`:

| File | Contents |
|---|---|
| `clusters.tsv` | One row per cluster: status, label, leading candidate, runner-up and the reason for an `Unassigned` status. |
| `marker_evidence.tsv` | Each label's score, control background and adjusted score in each cluster. |
| `control_gene_matching.tsv` | Expression level and detection rate of each marker gene and its controls. |
| `settings.rds`, `method.txt` | The settings used and a short description of the rule. |

`<store>/files/my_aggregation/3_GEX_QC/marker_set_UCell_summary/marker_sets.tsv` summarizes each marker set across the GEX clusters before doublet filtering: how many clusters it leads or is assigned, and its closest competing set.

## Further methods

See [algorithm validation](implementation/algorithm_validation.html) for reference comparisons and deviations, and the [peak–gene module](downstream_peak_gene_correlation.md) for its model, support filters and output paths.
