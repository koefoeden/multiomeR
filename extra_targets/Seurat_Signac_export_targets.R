rlang::list2(
  targets::tar_target(
    name = signac_annotation_GRanges,
    description = "Build gene annotation ranges for the Signac compatibility assay",
    command = build_signac_annotation_GRanges(marker_validated_Ensembl_annotations_GRanges_list),
    packages = w_def("Signac")
  ),
  targets::tar_target(
    name = signac_fragment_records_tibble,
    description = "Map prefixed pipeline barcodes to original CellRanger fragment files for Signac",
    command = make_signac_fragment_records(
      metadata_tibble = metadata_w_cell_types_tibble.WNN,
      reaction_ID_vec = aggregation_reaction_IDs,
      cellranger_summary_files = aggregation_cellranger_summary_file_syms
    )
  ),
  targets::tar_target(
    name = multimodal_Seurat_object,
    description = "Build a Seurat/Signac compatibility object backed by BPCells matrices where possible. [checkpoint:multimodal]",
    command = build_seurat_signac_convenience_object(
      GEX_counts_matrix = aggregated_counts_BPCells_matrix.GEX,
      ATAC_peak_matrix = peak_QC_filtered_BPCells_matrix.ATAC,
      metadata_tibble = metadata_w_cell_types_tibble.WNN,
      ATAC_peak_GRanges = consensus_peak_GRanges.ATAC,
      ATAC_annotated_peak_GRanges = consensus_peak_annotated_GRanges.ATAC,
      PCA_results = PCA_BPCells.GEX,
      GEX_harmony_embeddings = harmony_embeddings_matrix.GEX,
      ATAC_LSI_results = LSI_BPCells.ATAC,
      ATAC_harmony_embeddings = harmony_embeddings_matrix.ATAC,
      WNN_results = WNN_results,
      GEX_UMAP_embeddings_tibble = UMAP_embeddings_tibble.GEX,
      GEX_non_harmony_UMAP_embeddings_tibble = UMAP_embeddings_tibble.GEX_non_harmony,
      ATAC_UMAP_embeddings_tibble = UMAP_embeddings_tibble.ATAC,
      WNN_UMAP_embeddings_tibble = UMAP_embeddings_tibble.WNN,
      TF_activity_matrix = TF_activity_BPCells_matrix.ATAC,
      gene_features_df = gene_features_df,
      signac_annotation_GRanges = signac_annotation_GRanges,
      fragment_records_tibble = signac_fragment_records_tibble,
      GEX_dims = aggregation_GEX_data_PCs,
      ATAC_dims = aggregation_ATAC_data_PCs,
      data_nNNs = aggregation_data_nNNs,
      graph_threads = 15
    ),
    packages = w_def(c("Seurat", "SeuratObject", "Signac", "BPCells"))
  )
)
