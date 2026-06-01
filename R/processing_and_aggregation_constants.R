PROCESSING_QC_FEATURE_SETS <- list(
  # Per-reaction and per-dataset QC plots use the metadata available before
  # aggregation-level GEX/ATAC processing has added modality-specific metrics.
  per_reaction = c(
    "nCount_RNA",
    "log10_nCount_RNA",
    "RNA_mito_percent",
    "nFeature_RNA",
    "novelty",
    "TSS.enrichment",
    "nucleosome_signal",
    "vireo_max_prob_singlet",
    "vireo_max_prob_doublet",
    "amulet_q.value",
    "ATAC_TSS_fragments_frac",
    "ATAC_peak_region_frac",
    "ATAC_peak_region_cutsites_frac",
    "GEX_exonic_to_intronic_umis_frac"
  ),
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
    GEX_cell_cycle = c("G2M.Score", "S.Score"),
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
PROCESSING_PER_REACTION_QC_VARS <- PROCESSING_QC_FEATURE_SETS$per_reaction
