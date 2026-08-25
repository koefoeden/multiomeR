rlang::list2(
  tarchetypes::tar_file(
    name = aggregated_GEX_BPCells_matrix_dir.GEX,
    description = "Write the combined GEX BPCells count matrix for all GEM wells in the dataset [part_of_graph:GEX] [part_of_graph:parallel] [part_of_graph:seurat_export]",
    command = {
      combined_counts_matrix <- purrr::reduce(aggregation_GEX_counts_BPCells_matrix_syms, cbind)

      out_dir <- get_structured_file_path()
      BPCells::write_matrix_dir(combined_counts_matrix, out_dir, overwrite = TRUE)
      out_dir
    },
    resources = get_tar_resources(RAM_GB_req = 32)
  ),
  targets::tar_target(
    name = aggregated_counts_BPCells_matrix.GEX,
    description = "Open the combined GEX BPCells count matrix",
    command = BPCells::open_matrix_dir(aggregated_GEX_BPCells_matrix_dir.GEX),
    resources = get_tar_resources(RAM_GB_req = 8)
  ),
  targets::tar_target(
    name = PCA_BPCells.GEX,
    description = "Run Seurat SCTransform on BPCells-backed GEX counts and compute PCA [part_of_graph:GEX] [part_of_graph:seurat_export]",
    command = run_GEX_PCA_BPCells(
      GEX_counts_matrix = aggregated_counts_BPCells_matrix.GEX,
      metadata_tibble = GEX_cellranger_kept_metadata_tibble |>
        dplyr::left_join(donor_id_SCT_metadata_tibble, by = "donor_id") |>
        dplyr::left_join(GEM_well_SCT_metadata_tibble, by = "GEM_well_ID"),
      organism_chr = organism_chr,
      GEX_PCA_backend = aggregation_GEX_PCA_backend,
      SCT_regress_vars = aggregation_SCT_regress_vars,
      n_components = utils::tail(aggregation_GEX_data_PCs, 1),
      threads = 6
    ),
    resources = get_tar_resources(cores_req = 6, RAM_GB_req = 60) # temporary increase for large datasets using SCT_backend until we optimize this.
  ),
  targets::tar_target(
    name = metadata_tibble.GEX,
    description = "Prepare CellRanger-kept GEX metadata aligned to the native PCA cells [part_of_graph:GEX] [part_of_graph:seurat_export]",
    command = prepare_GEX_metadata_tibble(
      metadata_tibble = GEX_cellranger_kept_metadata_tibble,
      barcode_vec = rownames(PCA_BPCells.GEX$cell_embeddings)
    )
  ),
  targets::tar_target(
    name = metadata_harmony_tibble.GEX,
    description = "Join only configured Harmony covariates onto aligned GEX metadata",
    command = prepare_GEX_metadata_tibble(
      metadata_tibble = metadata_tibble.GEX,
      barcode_vec = metadata_tibble.GEX$barcode_w_prefix,
      donor_id_metadata_tibble = donor_id_harmony_metadata_tibble,
      GEM_well_metadata_tibble = GEM_well_harmony_metadata_tibble
    )
  ),
  targets::tar_target(
    name = metadata_analysis_tibble.GEX,
    description = "Join configured analysis variables onto aligned GEX metadata",
    command = prepare_GEX_metadata_tibble(
      metadata_tibble = metadata_tibble.GEX,
      barcode_vec = metadata_tibble.GEX$barcode_w_prefix,
      donor_id_metadata_tibble = donor_id_analysis_metadata_tibble,
      GEM_well_metadata_tibble = GEM_well_analysis_metadata_tibble
    )
  ),
  targets::tar_target(
    name = PCA_embeddings_tibble.GEX,
    description = "Cell embeddings from SCTransform GEX PCA",
    command = PCA_BPCells.GEX$cell_embeddings |>
      tibble::as_tibble(rownames = "barcode_w_prefix")
  ),
  targets::tar_target(
    name = PCA_loadings_tibble.GEX,
    description = "Gene loadings from SCTransform GEX PCA",
    command = PCA_BPCells.GEX$feature_loadings |>
      tibble::as_tibble(rownames = "gene")
  ),
  targets::tar_target(
    name = top_variable_genes.GEX,
    description = "Extract the top 10 most variable GEX genes from SCTransform residuals",
    command = utils::head(PCA_BPCells.GEX$variable_features, 10)
  ),
  targets::tar_target(
    name = harmony_embeddings_matrix.GEX,
    description = "Harmony-corrected SCTransform GEX PCA embeddings [part_of_graph:GEX] [part_of_graph:WNN] [part_of_graph:seurat_export]",
    command = run_harmony_on_embedding_matrix(
      embedding_matrix = PCA_BPCells.GEX$cell_embeddings,
      metadata_tibble = metadata_harmony_tibble.GEX,
      harmony_correction_metadata_col_names = aggregation_harmony_correction_metadata_col_names,
      dims = aggregation_GEX_data_PCs,
      cores = 6,
      dim_prefix = "PCA_"
    ),
    resources = get_tar_resources(cores_req = 6, RAM_GB_req = 60)
  ),
  targets::tar_target(
    name = UMAP_embeddings_tibble.GEX,
    description = "UMAP coordinates from BPCells-native Harmony-corrected GEX PCA embeddings [part_of_graph:GEX] [part_of_graph:seurat_export]",
    command = run_UMAP_from_embedding_matrix(
      embedding_matrix = harmony_embeddings_matrix.GEX,
      dims = aggregation_UMAP_GEX_PCs,
      n_neighbors = aggregation_UMAP_nNNs,
      min_dist = aggregation_UMAP_min_dist,
      dim_prefix = "PCA_",
      col_prefix = "GEX_UMAP",
      threads = 6
    ) |>
      tibble::as_tibble(rownames = "barcode_w_prefix"),
    resources = get_tar_resources(cores_req = 6, RAM_GB_req = 32)
  ),
  targets::tar_target(
    name = UMAP_embeddings_tibble.GEX_non_harmony,
    description = "UMAP coordinates from uncorrected BPCells-native GEX PCA embeddings [part_of_graph:GEX] [part_of_graph:seurat_export]",
    command = run_UMAP_from_embedding_matrix(
      embedding_matrix = PCA_BPCells.GEX$cell_embeddings,
      dims = aggregation_UMAP_GEX_PCs,
      n_neighbors = aggregation_UMAP_nNNs,
      min_dist = aggregation_UMAP_min_dist,
      dim_prefix = "PCA_",
      col_prefix = "GEX_non_harmony_UMAP",
      threads = 6
    ) |>
      tibble::as_tibble(rownames = "barcode_w_prefix"),
    resources = get_tar_resources(cores_req = 6, RAM_GB_req = 32)
  ),
  tarchetypes::tar_file(
    name = VizDimLoadings_plots.GEX,
    description = "Plot top gene loadings for each GEX PCA dimension. [checkpoint:GEX]",
    command = PCA_loadings_tibble.GEX |>
      plot_embedding_loadings_from_tibble(
        dims = aggregation_GEX_data_PCs,
        feature_col = "gene",
        dim_prefix = "PCA_",
        nfeatures = 50
      ) |>
      save_plots_structured(),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  tarchetypes::tar_file(
    name = variable_feature_plot.GEX,
    description = "Plot residual variance of top variable genes selected from SCTransform residuals. [checkpoint:GEX]",
    command = {
      plot <- PCA_BPCells.GEX$variable_feature_stats |>
        utils::head(50) |>
        dplyr::mutate(gene = forcats::fct_reorder(.data$gene, .data$residual_variance)) |>
        ggplot2::ggplot(ggplot2::aes(
          x = .data$residual_variance,
          y = .data$gene,
          color = .data$PCA_weighted_loading_strength
        )) +
        ggplot2::geom_point(size = 2) +
        ggplot2::labs(x = "SCT residual variance", y = NULL, color = "Weighted PCA loading strength")
      save_plots_structured(plot)
    },
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  tarchetypes::tar_file(
    name = PCA_singular_values_elbow_plot.GEX,
    description = "Elbow plot of native GEX PCA singular values. [checkpoint:GEX]",
    command = {
      plot <- plot_embedding_singular_values(
        singular_values = PCA_BPCells.GEX$singular_values,
        dims = aggregation_GEX_data_PCs
      )
      save_plots_structured(plot)
    },
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  tarchetypes::tar_file(
    name = PCA_embedding_sdev_plot.GEX,
    description = "Non-Harmony and Harmony GEX PCA embedding coordinate SD plot. [checkpoint:GEX]",
    command = {
      plot <- plot_embedding_sdev(
        embedding_matrix = PCA_BPCells.GEX$cell_embeddings,
        dims = aggregation_GEX_data_PCs,
        harmony_embedding_matrix = if (length(aggregation_harmony_correction_metadata_col_names %||% character()) > 0) harmony_embeddings_matrix.GEX else NULL,
        dim_prefix = "PCA_"
      )
      save_plots_structured(plot)
    },
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  tarchetypes::tar_file(
    name = PCA_metadata_association_barplots.GEX,
    description = "Non-Harmony and Harmony GEX PCA metadata association bar plots. [checkpoint:GEX]",
    command = {
      plots <- plot_embedding_metadata_association_barplots(
        embedding_matrix = PCA_BPCells.GEX$cell_embeddings,
        harmony_embedding_matrix = if (length(aggregation_harmony_correction_metadata_col_names %||% character()) > 0) harmony_embeddings_matrix.GEX else NULL,
        metadata_tibble = metadata_analysis_tibble.GEX,
        dims = aggregation_GEX_data_PCs,
        dim_prefix = "PCA_",
        continuous_technical_cols = c(
          "log10_nCount_RNA",
          "nFeature_RNA",
          "RNA_mito_percent",
          "GEX_exonic_to_intronic_umis_frac",
          "TSS.enrichment",
          "nucleosome_signal",
          "ATAC_peak_region_frac",
          "vireo_max_prob_singlet",
          "GEM_well_pre_amp_cycles"
        ),
        categorical_technical_cols = c(
          "GEM_well_ID",
          "GEM_well_multiplex_batch",
          "GEM_well_multiplex_pool",
          "GEM_well_run_harmony",
          "GEM_well_run_harmony_batched",
          "vireo_type",
          "discarded_by_cellranger",
          "discarded_by_amulet"
        ),
        continuous_biological_cols = aggregation_continuous_vars %||% character(),
        categorical_biological_cols = aggregation_categorical_vars %||% character()
      )
      save_plots_structured(
        plots,
        height = purrr::map_dbl(plots, \(plot) max(5, 0.75 * get_num_facet_rows(plot)))
      )
    },
    resources = get_tar_resources(RAM_GB_req = 16)
  )
)
