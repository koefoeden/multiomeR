PROCESSING_QC_FEATURE_SETS <- list(
  continuous = list(
    non_peak_QC = c(
      "log10_nCount_RNA",
      "RNA_mito_percent",
      "nFeature_RNA",
      "novelty",
      "TSS.enrichment",
      "nucleosome_signal",
      "amulet_q.value",
      "scDblFinder.score_GEX"
    ),
    peak_QC = c(
      "log10_nCount_ATAC",
      "atac_peak_counts_blacklist_frac",
      "atac_peak_counts_frac",
      "atac_peak_count_enrichment",
      "scDblFinder.score_ATAC"
    ),
    WNN = "ATAC.weight"
  ),
  categorical = list(
    aggregation_proj_spec_categorical_vars = character(),
    aggregation_GEX_categorical_vars = c(
      "PCA_harmony_SNN_cluster_named",
      "PCA_harmony_SNN_cluster_cell_type",
      "scDblFinder.class_GEX"
    ),
    aggregation_ATAC_categorical_vars = c(
      "LSI_harmony_SNN_cluster_named",
      "LSI_harmony_SNN_cluster_cell_type",
      "scDblFinder.class_ATAC"
    ),
    aggregation_WNN_categorical_vars = c(
      "WNN_harmony_SNN_cluster_named",
      "WNN_harmony_SNN_cluster_cell_type"
    ),
    aggregation_GEX_single_categorical_vars = "PCA_harmony_SNN_cluster",
    aggregation_ATAC_single_categorical_vars = "LSI_harmony_SNN_cluster",
    aggregation_WNN_single_categorical_vars = "WNN_harmony_SNN_cluster"
  )
)

PROCESSING_CONTINUOUS_QC_FEATURES <- PROCESSING_QC_FEATURE_SETS$continuous
