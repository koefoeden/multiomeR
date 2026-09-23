# Run your own analysis

Configure your data once, then work through the eight checkpoints of the main pipeline in order. At each checkpoint, run its targets, review its plots and revise its settings until you accept the result. Start with the defaults and revise as the plots suggest: each plot's subtitle and caption say what to look for, and the [output gallery](gallery.md) shows an example of every plot.

Run the R commands from the repository root, as in the demo. Replace `my_aggregation` and `my_GEM_well` with your own names. `<store>` is the targets store set in `_targets.yaml`, normally `outputs/` in the repository root. Edit the files in the [selected configuration directory](main_overview.md#configuration-directory); linked parameters belong in your aggregation's entry in `cfg_aggregations.yaml`.

Each command also builds the earlier results it depends on. After changing a setting, rerun the checkpoint: `targets` rebuilds only the results the change affects, and later checkpoints pick up the change when they run. The methods behind the checkpoints are described in [Preprocessing and nucleus QC](implementation/methods_preprocessing_and_QC.html), [GEX, ATAC, batch correction and WNN](implementation/methods_GEX_ATAC_and_WNN.html) and [Cell-type annotation and motif accessibility](implementation/methods_annotation_and_motifs.html).

## Configure your data {#steps}

1. **GEM wells:** add one row per GEM well to `cfg_GEM_wells.tsv`, with `GEM_well_QC_exclude_list` set to `NA` for the first run; see [GEM well table](reference_GEM_wells.md).
2. **Donors:** prepare the donor metadata TSV, with one row for each donor in these GEM wells; see [Donor metadata table](reference_donor_metadata.md).
3. **Aggregation:** in `cfg_aggregations.yaml`, copy the `template_aggregation` entry and rename the copy `my_aggregation`. Set `is_active: true`, list your GEM wells in `aggregation_GEM_well_IDs`, give the path to the donor table in `aggregation_donor_id_metadata_tsv`, and replace the placeholder `aggregation_GEX_marker_genes` with markers for the cell types you expect. Optionally, list the donor or GEM well columns to show in the plots in [`aggregation_categorical_vars`](parameters.html#aggregation_categorical_vars) and [`aggregation_continuous_vars`](parameters.html#aggregation_continuous_vars). See [Aggregation configuration](reference_aggregations.md).
4. **Demo:** set `is_active: false` for the `immune_human_2x` aggregation and `GEM_well_is_active` to `FALSE` for its GEM wells, `healthy_PBMC_human` and `lymphoma_lymph_human`. Change both together, because an active aggregation cannot use an inactive GEM well. A later `targets::tar_make()` without `names` then builds only your data.

## Checkpoint 1: Pre-aggregation QC {#checkpoint-1}

Computes QC metrics and donor assignments for each GEM well, applies the GEM well QC filters and compares the GEM wells of the aggregation.

**Run pipeline:**

```{.r filename="R"}
targets::tar_make(
  names = tidyselect::ends_with("1_pre_aggregation_QC.my_aggregation")
)
```

**Review plots** ([examples](gallery.md#1-pre-aggregation-qc)):

```text
<store>/plots/my_aggregation/1_pre_aggregation_QC/
├── per_aggregation_GEM_well_QC_comparisons/
├── demultiplexing_assignment_bars.png
├── nuclei_per_donor_id_bars.png
├── aggregation_excluded_barcodes_by_type_upset.png
├── aggregation_excluded_cellranger_only_barcodes_by_type_upset.png
└── cell_retention_flow_plot.png
```

For the exclusion overlaps within one GEM well, run its checkpoint 1 targets:

```{.r filename="R"}
targets::tar_make(
  names = tidyselect::ends_with("1_pre_aggregation_QC.my_GEM_well")
)
```

```text
<store>/plots/my_GEM_well/1_pre_aggregation_QC/
├── excluded_barcodes_by_type_upset.png
└── excluded_cellranger_only_barcodes_by_type_upset.png
```

**Revise:** choose cutoffs from the `per_aggregation_GEM_well_QC_comparisons/` distributions and enter them in `GEM_well_QC_exclude_list` for each GEM well; see [Set QC filters after the first run](reference_GEM_wells.md#set-qc-filters-after-the-first-run). After the rerun, these plots mark numeric cutoffs, and the UpSet and retention plots show the nuclei each filter removes. Remove a GEM well that fails QC from `aggregation_GEM_well_IDs`.

## Checkpoint 2: GEX dimension reduction {#checkpoint-2}

Normalizes the combined GEX data, selects variable genes and computes PCA, with Harmony correction when covariates are configured.

**Run pipeline:**

```{.r filename="R"}
targets::tar_make(
  names = tidyselect::ends_with("2_GEX_PCA_QC.my_aggregation")
)
```

**Review plots** ([examples](gallery.md#2-gex-pca-qc)):

```text
<store>/plots/my_aggregation/2_GEX_PCA_QC/
├── variable_feature_plot.png
├── VizDimLoadings_plots/
├── PCA_singular_values_elbow_plot.png
├── PCA_embedding_sdev_plot.png
└── PCA_metadata_association_barplots.png
```

**Revise:** choose the PCs used for clustering with [`aggregation_GEX_data_PCs`](parameters.html#aggregation_GEX_data_PCs) and for the UMAP with [`aggregation_UMAP_GEX_PCs`](parameters.html#aggregation_UMAP_GEX_PCs). When components track a batch variable such as `GEM_well_ID`, list it in [`aggregation_harmony_correction_metadata_col_names`](parameters.html#aggregation_harmony_correction_metadata_col_names); Harmony then corrects both modalities, and the embedding and association plots add Harmony panels.

## Checkpoint 3: GEX clusters and cell types {#checkpoint-3}

Clusters the GEX data, labels the clusters from your marker genes and removes GEX doublets.

**Run pipeline:**

```{.r filename="R"}
targets::tar_make(
  names = tidyselect::ends_with("3_GEX_QC.my_aggregation")
)
```

**Review plots** ([examples](gallery.md#3-gex-qc)):

```text
<store>/plots/my_aggregation/3_GEX_QC/
├── UMAPs/
│   ├── categorical/{harmony,non_harmony}/
│   ├── continuous/{harmony,non_harmony}/
│   └── cross/
├── markers_by_cluster_dot_plot.png
├── markers_by_cell_type_dot_plot.png
├── module_scores_by_cluster_dot_plot.png
├── module_scores_by_cell_type_dot_plot.png
├── cluster_UCell_advantage_plots/
├── cluster_marker_volcano_plots/
├── cell_type_marker_volcano_plots.png
├── categorical_by_cluster_bars_plots/
├── categorical_by_cell_type_bars_plots/
├── continuous_by_cluster_violin_plot/
├── continuous_by_cell_type_violin_plot/
└── cell_retention_flow_plot.png
```

**Revise:**

- [`aggregation_GEX_cluster_res`](parameters.html#aggregation_GEX_cluster_res) sets how finely the nuclei are clustered; judge it by marker consistency and cluster sizes. By default, the GEX clusters are also the peak-calling groups of checkpoint 4.
- [`aggregation_GEX_marker_genes`](parameters.html#aggregation_GEX_marker_genes) and [`aggregation_cluster_annotation_min_advantage`](parameters.html#aggregation_cluster_annotation_min_advantage) determine the cell-type labels; see [Cluster annotation](review_outputs.md#cluster-annotation).
- `UMAPs/cross/` compares layouts across PC and neighbor counts, for choosing [`aggregation_UMAP_GEX_PCs`](parameters.html#aggregation_UMAP_GEX_PCs) and [`aggregation_UMAP_nNNs`](parameters.html#aggregation_UMAP_nNNs).

## Checkpoint 4: Peak QC {#checkpoint-4}

Calls peaks within groups of nuclei, by default the GEX clusters, merges them into one consensus peak set and computes peak-based ATAC QC metrics.

**Run pipeline:**

```{.r filename="R"}
targets::tar_make(
  names = tidyselect::ends_with("4_peak_QC.my_aggregation")
)
```

**Review plots** ([examples](gallery.md#4-peak-qc)):

```text
<store>/plots/my_aggregation/4_peak_QC/
├── peaks_QC_violins_plot/
└── peaks_similarity_tiles_plot.png
```

**Revise:** choose ATAC filters from the `peaks_QC_violins_plot/` distributions and enter them in [`aggregation_QC_exclude_list_combined_object`](parameters.html#aggregation_QC_exclude_list_combined_object); see [QC filters after peak calling](reference_aggregations.md#qc-filters-after-peak-calling). Checkpoint 5 applies them.

## Checkpoint 5: ATAC filtering {#checkpoint-5}

Applies the ATAC filters from checkpoint 4 to the GEX-retained nuclei.

**Run pipeline:**

```{.r filename="R"}
targets::tar_make(
  names = tidyselect::ends_with("5_pre_LSI_QC.my_aggregation")
)
```

**Review plots** ([examples](gallery.md#5-pre-lsi-qc)):

```text
<store>/plots/my_aggregation/5_pre_LSI_QC/
├── QC_excluded_upset_plot.png
└── cell_retention_flow_plot.png
```

**Revise:** if the filters remove more nuclei than you can justify, or remove one GEM well disproportionately, revise [`aggregation_QC_exclude_list_combined_object`](parameters.html#aggregation_QC_exclude_list_combined_object).

## Checkpoint 6: ATAC dimension reduction {#checkpoint-6}

Computes LSI embeddings of the filtered peak matrix, with Harmony correction when covariates are configured.

**Run pipeline:**

```{.r filename="R"}
targets::tar_make(
  names = tidyselect::ends_with("6_ATAC_LSI_QC.my_aggregation")
)
```

**Review plots** ([examples](gallery.md#6-atac-lsi-qc)):

```text
<store>/plots/my_aggregation/6_ATAC_LSI_QC/
├── LSI_singular_values_elbow_plot.png
├── LSI_embedding_sdev_plot.png
├── VizDimLoadings_plots/
└── LSI_metadata_association_barplots.png
```

**Revise:** choose the LSI components used for clustering with [`aggregation_ATAC_data_PCs`](parameters.html#aggregation_ATAC_data_PCs) and for the UMAP with [`aggregation_UMAP_ATAC_PCs`](parameters.html#aggregation_UMAP_ATAC_PCs). Both start at component 2 by default; check the association plots for components that track sequencing depth or other technical metrics. Add ATAC-only Harmony covariates with [`aggregation_extra_harmony_covars_ATAC`](parameters.html#aggregation_extra_harmony_covars_ATAC).

## Checkpoint 7: ATAC clusters and motifs {#checkpoint-7}

Clusters the ATAC data, names the clusters from GEX marker scores, removes ATAC doublets, and computes gene activity and motif accessibility.

**Configure:** optionally, list expected transcription factors per cell type in [`aggregation_ATAC_marker_TFs`](parameters.html#aggregation_ATAC_marker_TFs) to organize the motif-accessibility plots.

**Run pipeline:**

```{.r filename="R"}
targets::tar_make(
  names = tidyselect::ends_with("7_ATAC_QC.my_aggregation")
)
```

**Review plots** ([examples](gallery.md#7-atac-qc)):

```text
<store>/plots/my_aggregation/7_ATAC_QC/
├── UMAPs/{categorical,continuous,cross}/
├── marker_gene_activity_dot_plot.png
├── motif_family_accessibility_by_ATAC_cluster_heatmap.png
├── motif_family_accessibility_by_GEX_cluster_heatmap.png
├── motif_family_accessibility_by_GEX_cell_type_heatmap.png
├── motif_family_accessibility_marker_volcano_plots.png
├── coverage_tracks_plots/
├── cluster_UCell_advantage_plots/
├── confusion_matrices_plots.png
├── categorical_by_cluster_bars_plots/
├── categorical_by_cell_type_bars_plots/
├── continuous_by_cluster_violin_plot/
├── continuous_by_cell_type_violin_plot/
└── cell_retention_flow_plot.png
```

**Revise:** choose [`aggregation_ATAC_cluster_res`](parameters.html#aggregation_ATAC_cluster_res) from the motif patterns and the agreement with GEX labels in `confusion_matrices_plots.png`.

## Checkpoint 8: WNN integration {#checkpoint-8}

Integrates the GEX and ATAC embeddings with weighted nearest neighbors (WNN), clusters and labels the integrated graph, and builds the final multimodal Seurat/Signac object.

**Run pipeline:**

```{.r filename="R"}
targets::tar_make(
  names = tidyselect::ends_with("8_multimodal_QC.my_aggregation")
)
```

**Review plots** ([examples](gallery.md#8-multimodal-qc)):

```text
<store>/plots/my_aggregation/8_multimodal_QC/
├── UMAPs/
│   ├── categorical/
│   ├── continuous/
│   ├── cross/
│   ├── cluster_named_dim_tri_plot.png
│   └── cluster_cell_type_dim_tri_plot.png
├── markers_by_cluster_dot_plot.png
├── module_scores_by_cluster_dot_plot.png
├── cluster_UCell_advantage_plots/
├── confusion_matrices_plots/
├── WNN_weight_metadata_associations_plot/
├── categorical_by_cluster_bars_plots/
├── categorical_by_cell_type_bars_plots/
├── continuous_by_cluster_violin_plot/
├── continuous_by_cell_type_violin_plot/
└── cell_retention_flow_plot.png
```

**Revise:** choose [`aggregation_WNN_cluster_res`](parameters.html#aggregation_WNN_cluster_res). The WNN clusters and cell types are the populations used in downstream summaries and comparisons.

The final object is `multimodal_Seurat_object.8_multimodal_QC.my_aggregation`; read it as in [Inspect the demo results](demo_outputs.md).

## Choose your next analysis

To build every remaining output, run `targets::tar_make()` without `names`. Then continue with an optional module:

- [Differential analyses](downstream_differential_analyses.md), if you have many donors and a condition of interest.
- [Genetic enrichment](downstream_genetic_enrichment.md), to relate fine-mapped GWAS variants for a human trait to accessibility in cell types.
- [Peak–gene correlation](downstream_peak_gene_correlation.md), to test associations between peak accessibility and gene expression.
