rlang::list2(
  targets::tar_target(
    name = GWAS_peak_weight_matrix,
    description = "Combine per-GWAS peak posterior-probability weights into one peak-by-GWAS annotation matrix [part_of_graph:genetic_enrichment_cell_type_contributions]",
    command = get_GWAS_chromVAR_peak_weight_matrix(
      peak_weight_records = GWAS_peak_weight_records
    ),
    resources = get_tar_resources(RAM_GB_req = 40)
  ),
  tarchetypes::tar_file(
    name = cell_type_pseudobulk_counts_BPCells_matrix_dir.ATAC,
    description = "Write ATAC counts summed per WNN cell type to BPCells for descriptive GWAS chromVAR plots [part_of_graph:genetic_enrichment_cell_type_contributions]",
    command = {
      out_dir <- get_structured_file_path()
      get_BPCells_group_pseudobulk_matrix(
        feature_matrix = peak_QC_filtered_BPCells_matrix.ATAC,
        metadata_tibble = metadata_w_cell_types_tibble.WNN,
        group_col = "PCA_harmony_SNN_cluster_cell_type",
        threads = 6
      ) |>
        methods::as("dgCMatrix") |>
        BPCells::write_matrix_dir(out_dir, overwrite = TRUE)
      out_dir
    },
    resources = get_tar_resources(cores_req = 6, RAM_GB_req = 60)
  ),
  targets::tar_target(
    name = cell_type_pseudobulk_counts_matrix.ATAC,
    description = "Open BPCells-backed ATAC counts summed per WNN cell type",
    command = BPCells::open_matrix_dir(cell_type_pseudobulk_counts_BPCells_matrix_dir.ATAC)
  ),
  targets::tar_target(
    name = cell_type_pseudobulk_support_tibble.ATAC,
    description = "Compute cell counts and ATAC depth support per WNN cell type",
    command = get_group_pseudobulk_depth_tibble(cell_type_pseudobulk_counts_matrix.ATAC) |>
      dplyr::left_join(
        get_group_cell_count_tibble(
          metadata_tibble = metadata_w_cell_types_tibble.WNN,
          group_col = "PCA_harmony_SNN_cluster_cell_type"
        ),
        by = "cluster"
      ) |>
      dplyr::mutate(modality = "ATAC")
  ),
  targets::tar_target(
    name = chromVAR_background_record.cell_type_pseudobulk,
    description = "Fit the cell-type pseudobulk betterChromVAR background used by GWAS score contributions [part_of_graph:genetic_enrichment_cell_type_contributions]",
    command = get_pseudobulk_chromVAR_background_record(
      pseudobulk_ATAC_data_matrix = cell_type_pseudobulk_counts_matrix.ATAC,
      chromVAR_obj = chromVAR_obj.ATAC
    ),
    resources = get_tar_resources(RAM_GB_req = 60)
  ),
  tarchetypes::tar_file(
    name = raw_deviation_unscaled.chromVAR_deviation_heatmaps.cell_type_pseudobulk,
    description = "Save ordinary GWAS deviations without within-GWAS scaling. [checkpoint:genetic_enrichment]",
    command = plot_GWAS_chromVAR_deviation_heatmap(
      chromVAR_deviation_tibble.cell_type_pseudobulk,
      GWAS_metadata_tracks_plot = GWAS_metadata_tracks_plot,
      compartments_patterns = genetic_enrichment_compartment_patterns,
      standardize = FALSE,
      beta_weighted = FALSE
    ) |>
      save_plots_structured(
        filetype = "png", width = 17,
        height = max(5.5, 0.3 * dplyr::n_distinct(chromVAR_deviation_tibble.cell_type_pseudobulk$GWAS_ID) + 3.5)
      )
  ),
  tarchetypes::tar_file(
    name = raw_deviation_scaled.chromVAR_deviation_heatmaps.cell_type_pseudobulk,
    description = "Save ordinary GWAS deviations standardized within each GWAS. [checkpoint:genetic_enrichment]",
    command = plot_GWAS_chromVAR_deviation_heatmap(
      chromVAR_deviation_tibble.cell_type_pseudobulk,
      GWAS_metadata_tracks_plot = GWAS_metadata_tracks_plot,
      compartments_patterns = genetic_enrichment_compartment_patterns,
      standardize = TRUE,
      beta_weighted = FALSE
    ) |>
      save_plots_structured(
        filetype = "png", width = 17,
        height = max(5.5, 0.3 * dplyr::n_distinct(chromVAR_deviation_tibble.cell_type_pseudobulk$GWAS_ID) + 3.5)
      )
  )
)
