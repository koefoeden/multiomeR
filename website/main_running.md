# Run your own analysis

This section assumes that you have followed the demo, such that things are correctly installed, and you are comfortable with where to locate outputs and how to run the pipeline. Here, we will replace the relevant configuration for the demo and run a small part of the pipeline one step at a time, through a repeated series of actions following the same pattern: **Configure** -\> **Run** **pipeline**-\> **Review** **plots** (and Repeat if needed).

## 1. Pre-process the cellranger-arc count dirs (GEM wells) {#steps}

**Initial configuration:**

- Add entries for each GEM-well (cellranger-arc count dir) you want to process inside `cfg_GEM_wells.tsv`. See [GEM well table](reference_GEM_wells.md) for more information

- Add a bare-bones entry in `cfg_aggregations.yaml`. See [aggregation configuration](reference_aggregations.md) for more information.

**Run pipeline:**

``` {.r filename="R"}
targets::tar_make(
  names = tidyselect::ends_with("1_pre_aggregation_QC.my_aggregation")
)
```

**Review plots:**

``` text
<store>/plots/my_aggregation/1_pre_aggregation_QC/
├── per_aggregation_GEM_well_QC_comparisons/
├── aggregation_excluded_cellranger_only_barcodes_by_type_upset.png
├── aggregation_excluded_barcodes_by_type_upset.png
├── nuclei_per_donor_id_bars.png
└── cell_retention_flow_plot.png <to implement the reaction-level upsets here>
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

- `aggregation_GEX_marker_genes,`categorical and continuous metadata variables for the review plots: \< format as link to parameter reference\>

**Run pipeline:**

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

## 4. Call peaks and inspect ATAC quality \<combine 4 and 5 into a single checkpoint\>

**Initial configuration:**

- Start with default settings and revise as described in the plots if necessary.

``` {.r filename="R"}
targets::tar_make(
  names = tidyselect::ends_with("4_peak_QC.my_aggregation")
)
```

**Review plots**

``` text
<store>/plots/my_aggregation/4_peak_QC/
├── peaks_QC_violins_plot/
└── peaks_similarity_tiles_plot.png
```

## 5. Filter nuclei on ATAC quality

Applies the expressions from step 5 and shows what they removed.

**Initial configuration:**

**Run pipeline:**

``` {.r filename="R"}
targets::tar_make(
  names = tidyselect::ends_with("5_pre_LSI_QC.my_aggregation")
)
```

**Review plots:**

``` text
<store>/plots/my_aggregation/5_pre_LSI_QC/
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

- `aggregation_ATAC_marker_TFs`: marker transcription factors per cell type;

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

-  [Genetic enrichment](downstream_genetic_enrichment.md), if you are interested in pinpointing cell-type-level genetic enrichment

- Peak gene correlation, if you have many cells \<format as link\>