# Run your own analysis

After completing the demo, replace its configuration and work through **Configure → Run pipeline → Review plots**, repeating each step as needed. The plots explain what to inspect and which settings to revise.

Run commands from the repository-root R session. Replace `my_GEM_well` with your configured identifiers. `<store>` means the folder selected in `_targets.yaml`, which is normally `outputs` in the root of the repository.

## 1. Pre-process the cellranger-arc count dirs (GEM wells) {#steps}

**Initial configuration:**

- Add entries for each GEM-well (cellranger-arc count dir) you want to process inside `cfg_GEM_wells.tsv`. See [GEM well table](reference_GEM_wells.md) for more information

- Locate the `template_aggregation` inside `cfg_aggregations.yaml` and replace the GEM_well_ID-placeholders with the GEM-well IDs that you just added - these GEM-wells are now officially part of this aggregation. Finish by renaming the aggregation-entry using a short, descriptive name. Throughout the rest of this page, replace the \<my_aggregation\>-placeholders with your custom name.

**Run pipeline:**

``` {.r filename="R"}
targets::tar_make(
  names = tidyselect::ends_with("1_pre_aggregation_QC.<my_aggregation>")
)
```

**Review plots:**

``` text
<store>/plots/my_aggregation/1_pre_aggregation_QC/
├── per_aggregation_GEM_well_QC_comparisons/
├── aggregation_excluded_cellranger_only_barcodes_by_type_upset.png
├── aggregation_excluded_barcodes_by_type_upset.png
├── nuclei_per_donor_id_bars.png
└── cell_retention_flow_plot.png
```

For exclusion overlaps within an individual GEM well, you can also request its plots:

``` {.r filename="R"}
targets::tar_make(
  names = tidyselect::ends_with("1_pre_aggregation_QC.<my_GEM_well_ID>")
)
```

``` text
<store>/plots/my_GEM_well/1_pre_aggregation_QC/
├── excluded_barcodes_by_type_upset.png
└── excluded_cellranger_only_barcodes_by_type_upset.png
```

## 2. Normalize, reduce dimensions, and batch-correct the gene-expression data

**Initial configuration:**

- Start with default settings and revise as described in the plots if necessary.

**Run pipeline:**

``` {.r filename="R"}
targets::tar_make(
  names = tidyselect::ends_with("2_GEX_PCA_QC.my_aggregation")
)
```

**Review plots:**

``` text
<store>/plots/my_aggregation/2_GEX_PCA_QC/
├── variable_feature_plot.png
├── VizDimLoadings_plots/
├── PCA_singular_values_elbow_plot.png
├── PCA_embedding_sdev_plot.png
└── PCA_metadata_association_barplots/
```

## 3. Cluster and label cell types with GEX

**Initial configuration:**

Set [`aggregation_GEX_marker_genes`](parameters.html#aggregation_GEX_marker_genes), [`aggregation_categorical_vars`](parameters.html#aggregation_categorical_vars) & [`aggregation_continuous_vars`](parameters.html#aggregation_continuous_vars)

``` {.r filename="R"}
targets::tar_make(
  names = tidyselect::ends_with("3_GEX_QC.my_aggregation")
)
```

**Review plots:**

``` text
<store>/plots/my_aggregation/3_GEX_QC/
├── UMAPs/
│   ├── categorical/{harmony,non_harmony}/
│   ├── continuous/{harmony,non_harmony}/
│   └── cross/
├── categorical_by_cell_type_bars_plots/
├── categorical_by_cluster_bars_plots/
├── continuous_by_cluster_violin_plot/
├── continuous_by_cell_type_violin_plot/
├── markers_by_cluster_dot_plot.png
├── markers_by_cell_type_dot_plot.png
├── module_scores_by_cluster_dot_plot.png
├── module_scores_by_cell_type_dot_plot.png
├── cluster_UCell_advantage_plots/
├── cluster_marker_volcano_plots/
├── cell_type_marker_volcano_plots.png
└── cell_retention_flow_plot.png
```

## 4. Call peaks, inspect ATAC quality, and filter nuclei

**Initial configuration:**

- Start with default settings and revise as described in the plots if necessary.

**Run pipeline:**

``` {.r filename="R"}
targets::tar_make(
  names = tidyselect::ends_with(c(
    "4_peak_QC.my_aggregation",
    "5_pre_LSI_QC.my_aggregation"
  ))
)
```

**Review plots:**

``` text
<store>/plots/my_aggregation/
├── 4_peak_QC/
│   ├── peaks_QC_violins_plot/
│   └── peaks_similarity_tiles_plot.png
└── 5_pre_LSI_QC/
    ├── QC_excluded_upset_plot.png
    └── cell_retention_flow_plot.png
```

## 5. Normalize, reduce dimensions and batch-correct the ATAC-data

Computes LSI on the retained nuclei with optional Harmony correction.

**Initial configuration:**

- Start with default settings and revise as described in the plots if necessary.

**Run pipeline:**

``` {.r filename="R"}
targets::tar_make(
  names = tidyselect::ends_with("6_ATAC_LSI_QC.my_aggregation")
)
```

**Review plots:**

``` text
<store>/plots/my_aggregation/6_ATAC_LSI_QC/
├── LSI_singular_values_elbow_plot.png
├── LSI_embedding_sdev_plot.png
├── VizDimLoadings_plots/
└── LSI_metadata_association_barplots/
```

## 6. Generate ATAC-clusters and motif accessibility

**Initial configuration:**

- `Set aggregation_ATAC_marker_TFs`

**Run pipeline:**

``` {.r filename="R"}
targets::tar_make(
  names = tidyselect::ends_with("7_ATAC_QC.my_aggregation")
)
```

**Review plots:**

``` text
<store>/plots/my_aggregation/7_ATAC_QC/
├── UMAPs/{categorical,continuous,cross}/
├── categorical_by_cell_type_bars_plots/
├── categorical_by_cluster_bars_plots/
├── marker_gene_activity_dot_plot.png
├── motif_family_accessibility_by_ATAC_cluster_heatmap.png
├── motif_family_accessibility_by_GEX_cluster_heatmap.png
├── motif_family_accessibility_by_GEX_cell_type_heatmap.png
├── motif_family_accessibility_marker_volcano_plots.png
├── coverage_tracks_plots/
├── cluster_UCell_advantage_plots/
├── confusion_matrices_plots.png
└── cell_retention_flow_plot.png
```

## 7. Integrate GEX and ATAC

**Initial configuration:**

- Start with default settings and revise as described in the plots if necessary.

**Run pipeline:**

``` {.r filename="R"}
targets::tar_make(
  names = tidyselect::ends_with("8_multimodal_QC.my_aggregation")
)
```

**Review plots:**

``` text
<store>/plots/my_aggregation/8_multimodal_QC/
├── UMAPs/{categorical,continuous,cross}/
├── categorical_by_cell_type_bars_plots/
├── categorical_by_cluster_bars_plots/
├── UMAPs/cluster_named_dim_tri_plot.png
├── UMAPs/cluster_cell_type_dim_tri_plot.png
├── markers_by_cluster_dot_plot.png
├── module_scores_by_cluster_dot_plot.png
├── continuous_by_cell_type_violin_plot/
├── continuous_by_cluster_violin_plot/
├── WNN_weight_metadata_associations_plot.png
├── confusion_matrices_plots/
├── cluster_UCell_advantage_plots/
└── cell_retention_flow_plot.png
```

The final object is `multimodal_Seurat_object.8_multimodal_QC.my_aggregation`.

## 8. Choose your next analysis!

Continue to one of the three modules:

- [Differential analyses](downstream_differential_analyses.md), if you have many donors and a condition of interest

- [Genetic enrichment](downstream_genetic_enrichment.md), if you are interested in pinpointing cell-type-level genetic enrichment

- [Peak–gene correlation](downstream_peak_gene_correlation.md), to test associations between peak accessibility and gene expression.