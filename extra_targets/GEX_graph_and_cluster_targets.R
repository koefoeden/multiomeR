rlang::list2(
  targets::tar_target(
    name = PCA_clusters.GEX,
    description = "Leiden clusters from a BPCells nearest-neighbour graph on Harmony-corrected GEX PCA embeddings",
    command = cluster_embedding_matrix_BPCells(
      embedding_matrix = harmony_embeddings_matrix.GEX,
      dims = aggregation_GEX_data_PCs,
      k = aggregation_data_nNNs,
      resolution = aggregation_GEX_cluster_res,
      dim_prefix = "PCA_",
      threads = 15,
      min_barcodes = aggregation_cluster_min_barcodes
    ),
    resources = get_tar_resources(cores_req = 15, RAM_GB_req = 32)
  ),
  targets::tar_target(
    name = metadata_w_clusters_tibble.GEX,
    description = "BPCells-native GEX metadata with PCA clusters and UMAP coordinates",
    command = metadata_tibble.GEX |>
      dplyr::filter(.data$barcode_w_prefix %in% names(PCA_clusters.GEX)) |>
      dplyr::left_join(UMAP_embeddings_tibble.GEX, by = "barcode_w_prefix") |>
      dplyr::left_join(UMAP_embeddings_tibble.GEX_non_harmony, by = "barcode_w_prefix") |>
      dplyr::mutate(PCA_harmony_SNN_cluster = PCA_clusters.GEX[.data$barcode_w_prefix]),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = metadata_w_cell_types_unfiltered_tibble.GEX,
    description = "Assign GEX cell types from native marker-module scores before doublet removal",
    command = metadata_w_clusters_tibble.GEX |>
      add_GEX_module_scores_to_metadata(
        GEX_counts_matrix = aggregated_counts_BPCells_matrix.GEX,
        named_marker_genes_list = UCell_GEX_marker_genes_list
      ) |>
      add_cell_types_to_metadata_from_module_scores(
        named_marker_genes_list = UCell_GEX_marker_genes_list,
        allow_multiple_cell_types = aggregation_allow_multiple_cell_types,
        cluster_column = "PCA_harmony_SNN_cluster"
      ),
    resources = get_tar_resources(RAM_GB_req = 32)
  ),
  targets::tar_target(
    name = scDblFinder_reaction_tibble.GEX,
    description = "Prepare per-reaction GEX barcodes and cluster labels for scDblFinder",
    command = prepare_scDblFinder_reaction_tibble(
      metadata_tibble = metadata_w_cell_types_unfiltered_tibble.GEX,
      cluster_collapse_list = aggregation_scDblFinder_GEX_cell_type_collapse_list,
      cluster_col = "PCA_harmony_SNN_cluster_cell_type"
    ),
    iteration = "vector"
  ),
  targets::tar_target(
    name = scDblFinder_results_by_reaction_tibble.GEX,
    description = "Run GEX scDblFinder independently per 10x reaction from BPCells count slices",
    command = run_scDblFinder_BPCells_reaction(
      feature_matrix = aggregated_counts_BPCells_matrix.GEX,
      scDblFinder_reaction_tibble = scDblFinder_reaction_tibble.GEX,
      output_suffix = "GEX",
      dbr.sd = 1.0
    ),
    pattern = map(scDblFinder_reaction_tibble.GEX),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = scDblFinder_results_df.GEX,
    description = "Combine per-reaction GEX scDblFinder classifications",
    command = scDblFinder_results_by_reaction_tibble.GEX |>
      dplyr::bind_rows() |>
      dplyr::select(-dplyr::any_of("TENX_reaction_ID")) |>
      tibble::column_to_rownames("barcode_w_prefix"),
    resources = get_tar_resources(RAM_GB_req = 8)
  ),
  tarchetypes::tar_file(
    name = scDblFinder_score_violins_plot.GEX,
    description = "Violin plots of scDblFinder doublet scores per cluster and cell type. [checkpoint:GEX]",
    command = {
      scDblFinder_metadata <- metadata_w_cell_types_unfiltered_tibble.GEX |>
        dplyr::left_join(scDblFinder_results_df.GEX |> tibble::rownames_to_column("barcode_w_prefix"), by = "barcode_w_prefix")
      plot <- scDblFinder_metadata |>
        dplyr::select(PCA_harmony_SNN_cluster, PCA_harmony_SNN_cluster_cell_type, scDblFinder.score_GEX, scDblFinder.class_GEX, TENX_reaction_ID) |>
        tidyr::pivot_longer(cols = dplyr::all_of(c("PCA_harmony_SNN_cluster", "PCA_harmony_SNN_cluster_cell_type")), names_to = "cluster_type", values_to = "cluster_id") %>%
        ggplot2::ggplot(ggplot2::aes(x = cluster_id, y = scDblFinder.score_GEX)) +
        ggplot2::geom_violin(scale = "width") +
        ggplot2::geom_jitter(
          data = \(d) dplyr::filter(d, scDblFinder.class_GEX == "doublet"),
          ggplot2::aes(color = TENX_reaction_ID),
          size = 0.3,
          alpha = 0.5,
          width = 0.2
        ) +
        ggplot2::facet_wrap(~cluster_type, scales = "free_x") +
        ggplot2::labs(subtitle = "Points are cells classified as doublet by scDblFinder, colored by reaction.") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
        ggplot2::theme(legend.position = "none")

      save_plots_structured(plot)
    },
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = metadata_w_cell_types_tibble.GEX,
    description = "Remove GEX doublets and high-doublet clusters, then annotate remaining cells with scDblFinder scores",
    command = {
      scDblFinder_tibble <- scDblFinder_results_df.GEX |>
        tibble::rownames_to_column("barcode_w_prefix")
      excluded_doublets_barcodes <- scDblFinder_tibble |>
        dplyr::filter(.data$scDblFinder.class_GEX == "doublet") |>
        dplyr::pull(.data$barcode_w_prefix)

      high_doublet_clusters_vec <- metadata_w_cell_types_unfiltered_tibble.GEX |>
        dplyr::mutate(is_excluded = .data$barcode_w_prefix %in% excluded_doublets_barcodes) |>
        dplyr::group_by(.data$PCA_harmony_SNN_cluster) |>
        dplyr::summarise(frac_excluded = mean(is_excluded), .groups = "drop") |>
        dplyr::filter(frac_excluded > 0.5) |>
        dplyr::pull(.data$PCA_harmony_SNN_cluster)

      excluded_high_doublet_cluster_barcodes_vec <- metadata_w_cell_types_unfiltered_tibble.GEX |>
        dplyr::filter(.data$PCA_harmony_SNN_cluster %in% high_doublet_clusters_vec) |>
        dplyr::pull(.data$barcode_w_prefix)

      metadata_w_cell_types_unfiltered_tibble.GEX |>
        dplyr::filter(!.data$barcode_w_prefix %in% base::union(excluded_doublets_barcodes, excluded_high_doublet_cluster_barcodes_vec)) |>
        dplyr::left_join(scDblFinder_tibble, by = "barcode_w_prefix")
    },
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = SNN_cluster_marker_tibbles.GEX,
    description = "Compute BPCells marker genes for each GEX SNN cluster against other clusters within the same cell type",
    command = get_BPCells_markers_within_parent_groups_from_matrix(
      feature_matrix = aggregated_counts_BPCells_matrix.GEX,
      metadata_tibble = metadata_w_cell_types_tibble.GEX,
      group_col = "PCA_harmony_SNN_cluster",
      parent_group_col = "PCA_harmony_SNN_cluster_cell_type"
    ),
    resources = get_tar_resources(RAM_GB_req = 32)
  ),
  targets::tar_target(
    name = cell_type_marker_tibbles.GEX,
    description = "Compute BPCells marker genes for each GEX cell type against all other cell types",
    command = get_BPCells_markers_from_matrix(
      feature_matrix = aggregated_counts_BPCells_matrix.GEX,
      metadata_tibble = metadata_w_cell_types_tibble.GEX,
      group_col = "PCA_harmony_SNN_cluster_cell_type"
    ),
    resources = get_tar_resources(RAM_GB_req = 32)
  ),
  targets::tar_target(
    name = UMAP_n_dims_seq.GEX,
    description = "Generate a sequence of PC counts for cross-parameter UMAP exploration",
    command = round(seq(5, length(aggregation_GEX_data_PCs), length.out = 3))
  ),
  targets::tar_target(
    name = UMAP_neighbors_seq,
    description = "Generate a sequence of neighbour counts for cross-parameter UMAP exploration",
    command = round(seq(10, aggregation_UMAP_nNNs, length.out = 3))
  ),
  tarchetypes::tar_file(
    name = cross.UMAPs.GEX,
    description = "Compute GEX UMAPs across a sweep of PC counts and neighbour counts. [checkpoint:GEX]",
    command = {
      sweep_umap <- run_UMAP_from_embedding_matrix(
        embedding_matrix = harmony_embeddings_matrix.GEX,
        dims = seq_len(UMAP_n_dims_seq.GEX),
        n_neighbors = UMAP_neighbors_seq,
        min_dist = aggregation_UMAP_min_dist,
        dim_prefix = "PCA_",
        col_prefix = "GEX_UMAP"
      ) |>
        tibble::as_tibble(rownames = "barcode_w_prefix")

      metadata_w_cell_types_tibble.GEX |>
        dplyr::select(-dplyr::any_of(c("GEX_UMAP_1", "GEX_UMAP_2"))) |>
        dplyr::left_join(sweep_umap, by = "barcode_w_prefix") |>
        plot_UMAP_from_metadata(metadata_cols = "PCA_harmony_SNN_cluster_cell_type", umap_cols = c("GEX_UMAP_1", "GEX_UMAP_2"))
    } |>
      save_plots_structured(dyn_suffix_in_subdir = TRUE, override_suffix = paste0(UMAP_n_dims_seq.GEX, "_", UMAP_neighbors_seq)),
    pattern = cross(UMAP_n_dims_seq.GEX, UMAP_neighbors_seq),
    resources = get_tar_resources(cores_req = 15, RAM_GB_req = 32)
  ),
  tarchetypes::tar_file(
    name = harmony.categorical.UMAPs.GEX,
    description = "UMAPs colored by categorical metadata variables. [checkpoint:GEX]",
    command = metadata_w_cell_types_tibble.GEX |>
      plot_UMAP_from_metadata(metadata_cols = aggregation_GEX_categorical_vars, umap_cols = c("GEX_UMAP_1", "GEX_UMAP_2")) |>
      save_plots_structured(),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  tarchetypes::tar_file(
    name = harmony.continuous.UMAPs.GEX,
    description = "UMAPs colored by continuous QC and gene expression features. [checkpoint:GEX]",
    command = plot_UMAP_from_metadata(
      metadata_tibble = metadata_w_cell_types_tibble.GEX,
      metadata_cols = aggregation_non_peak_based_continuous_vars,
      feature_matrix = aggregated_counts_BPCells_matrix.GEX,
      feature_rows = c(GEX_marker_genes_vec, interesting_genes, top_variable_genes.GEX),
      umap_cols = c("GEX_UMAP_1", "GEX_UMAP_2")
    ) |>
      save_plots_structured(),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  tarchetypes::tar_file(
    name = non_harmony.categorical.UMAPs.GEX,
    description = "UMAPs colored by categorical metadata variables on the uncorrected PCA embedding. [checkpoint:GEX]",
    command = metadata_w_cell_types_tibble.GEX |>
      plot_UMAP_from_metadata(metadata_cols = aggregation_GEX_categorical_vars, umap_cols = c("GEX_non_harmony_UMAP_1", "GEX_non_harmony_UMAP_2")) |>
      save_plots_structured(),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  tarchetypes::tar_file(
    name = non_harmony.continuous.UMAPs.GEX,
    description = "UMAPs colored by continuous features on the uncorrected PCA embedding. [checkpoint:GEX]",
    command = plot_UMAP_from_metadata(
      metadata_tibble = metadata_w_cell_types_tibble.GEX,
      metadata_cols = aggregation_non_peak_based_continuous_vars,
      feature_matrix = aggregated_counts_BPCells_matrix.GEX,
      feature_rows = c(GEX_marker_genes_vec, interesting_genes),
      umap_cols = c("GEX_non_harmony_UMAP_1", "GEX_non_harmony_UMAP_2")
    ) |>
      save_plots_structured(),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  tarchetypes::tar_file(
    name = markers_dot_plot.GEX,
    description = "Dot plot of marker gene expression per GEX cell type. [checkpoint:GEX]",
    command = {
      plot <- plot_marker_expression_dot_BPCells(
        feature_matrix = aggregated_counts_BPCells_matrix.GEX,
        metadata_tibble = metadata_w_cell_types_tibble.GEX,
        features = GEX_marker_genes_vec,
        group_col = "PCA_harmony_SNN_cluster_cell_type"
      )
      save_plots_structured(plot)
    },
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  tarchetypes::tar_file(
    name = continuous_by_cell_type_violin_plot.GEX,
    description = "Violin plots of continuous QC and cell-cycle features per cell type. [checkpoint:GEX]",
    command = {
      plot <- metadata_w_cell_types_tibble.GEX |>
        dplyr::select(PCA_harmony_SNN_cluster_cell_type, dplyr::any_of(aggregation_continuous_features_vec.GEX)) |>
        tidyr::pivot_longer(cols = -PCA_harmony_SNN_cluster_cell_type, names_to = "feature", values_to = "value") |>
        dplyr::filter(
          value >= stats::quantile(value, probs = 0.02, na.rm = TRUE),
          value <= stats::quantile(value, probs = 0.98, na.rm = TRUE),
          .by = c(PCA_harmony_SNN_cluster_cell_type, feature),
        ) %>%
        ggplot2::ggplot(ggplot2::aes(x = PCA_harmony_SNN_cluster_cell_type, y = value, fill = PCA_harmony_SNN_cluster_cell_type)) +
        ggplot2::geom_violin(scale = "width") +
        ggplot2::facet_wrap(~feature, scales = "free") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), axis.title.x = ggplot2::element_blank(), legend.position = "none")
      save_plots_structured(plot)
    }
  ),
  tarchetypes::tar_file(
    name = continuous_by_cluster_violin_plot.GEX,
    description = "Violin plots of continuous QC and cell-cycle features per SNN cluster. [checkpoint:GEX]",
    command = {
      plot <- metadata_w_cell_types_tibble.GEX |>
        dplyr::select(PCA_harmony_SNN_cluster_named, dplyr::any_of(aggregation_continuous_features_vec.GEX)) |>
        tidyr::pivot_longer(cols = -PCA_harmony_SNN_cluster_named, names_to = "feature", values_to = "value") |>
        dplyr::filter(
          value >= stats::quantile(value, probs = 0.02, na.rm = TRUE),
          value <= stats::quantile(value, probs = 0.98, na.rm = TRUE),
          .by = c(PCA_harmony_SNN_cluster_named, feature),
        ) %>%
        ggplot2::ggplot(ggplot2::aes(x = PCA_harmony_SNN_cluster_named, y = value, fill = PCA_harmony_SNN_cluster_named)) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), axis.title.x = ggplot2::element_blank(), legend.position = "none") +
        ggplot2::geom_violin(scale = "width") +
        ggplot2::facet_wrap(~feature, scales = "free")
      save_plots_structured(plot)
    }
  ),
  tarchetypes::tar_file(
    name = module_scores_dot_plot.GEX,
    description = "Dot plot of marker module scores per GEX cell type. [checkpoint:GEX]",
    command = plot_module_scores_dot_for_metadata(
      metadata_tibble = metadata_w_cell_types_tibble.GEX,
      marker_genes_list = UCell_GEX_marker_genes_list,
      cluster_by = "PCA_harmony_SNN_cluster_cell_type"
    ) |>
      save_plots_structured()
  ),
  tarchetypes::tar_file(
    name = categorical_bars_plots.GEX,
    description = "Bar plots of categorical metadata composition per cell type. [checkpoint:GEX]",
    command = plot_categorical_bars_plot(
      metadata_tibble = metadata_w_cell_types_tibble.GEX,
      metadata_cols = aggregation_GEX_categorical_vars,
      cluster_col = "PCA_harmony_SNN_cluster_cell_type"
    ) |>
      save_plots_structured()
  ),
  tarchetypes::tar_file(
    name = marker_volcano_plots.GEX,
    description = "Facetted volcano plot of BPCells marker genes per SNN cluster. [checkpoint:GEX]",
    command = SNN_cluster_marker_tibbles.GEX |>
      plot_markers_volcano_simple() |>
      save_plots_structured()
  ),
  tarchetypes::tar_file(
    name = cell_type_marker_volcano_plots.GEX,
    description = "Facetted volcano plot of BPCells marker genes per GEX cell type. [checkpoint:GEX]",
    command = cell_type_marker_tibbles.GEX |>
      plot_markers_volcano_simple() |>
      save_plots_structured()
  ),
  targets::tar_target(
    name = GEX_Seurat_object,
    description = "Build a GEX-only Seurat compatibility object for manual review before ATAC peak calling. [checkpoint:GEX]",
    command = build_seurat_signac_convenience_object(
      GEX_counts_matrix = aggregated_counts_BPCells_matrix.GEX,
      metadata_tibble = metadata_w_cell_types_tibble.GEX,
      PCA_results = PCA_BPCells.GEX,
      GEX_harmony_embeddings = harmony_embeddings_matrix.GEX,
      GEX_UMAP_embeddings_tibble = UMAP_embeddings_tibble.GEX,
      GEX_non_harmony_UMAP_embeddings_tibble = UMAP_embeddings_tibble.GEX_non_harmony,
      gene_features_df = gene_features_df,
      GEX_dims = aggregation_GEX_data_PCs,
      data_nNNs = aggregation_data_nNNs,
      graph_threads = 15
    ),
    packages = w_def(c("Seurat", "SeuratObject", "BPCells"))
  ),
)
