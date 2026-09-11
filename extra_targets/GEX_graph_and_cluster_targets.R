rlang::list2(
  targets::tar_target(
    name = cell_retention_tibble.GEX,
    description = "Accumulate per-well retention counts through GEX clustering and doublet filtering. [checkpoint:3_GEX-QC]",
    command = summarize_QC_cell_retention(
      before_metadata = GEX_cellranger_kept_metadata_tibble,
      after_metadata = metadata_w_cell_types_tibble.GEX,
      GEM_well_IDs = aggregation_GEM_well_IDs,
      stage = "checkpoint:3_GEX-QC",
      discarded_barcodes = c(
        list("Small GEX clusters" = setdiff(GEX_cellranger_kept_metadata_tibble$barcode_w_prefix,
          metadata_w_clusters_tibble.GEX$barcode_w_prefix)),
        get_doublet_QC_discarded_barcodes(metadata_w_clusters_tibble.GEX,
          metadata_w_cell_types_tibble.GEX, scDblFinder_results_df.GEX,
          class_col = "scDblFinder.class_GEX",
          remove_called_doublets = aggregation_scDblFinder_GEX_remove_called_doublets)
      ),
      previous_stages = cell_retention_tibble.GEX_input
    )
  ),
  tarchetypes::tar_file(
    name = cell_retention_flow_plot.3_GEX_QC,
    description = "Plot cumulative nuclei retention through GEX QC. [checkpoint:3_GEX-QC]",
    command = save_QC_cell_retention_plot(cell_retention_tibble.GEX)
  ),
  targets::tar_target(
    name = PCA_clusters.GEX,
    description = "Leiden clusters from a BPCells nearest-neighbour graph on Harmony-corrected GEX PCA embeddings [part_of_graph:GEX] [part_of_graph:seurat_export]",
    command = cluster_embedding_matrix_BPCells(
      embedding_matrix = harmony_embeddings_matrix.GEX,
      dims = aggregation_GEX_data_PCs,
      k = aggregation_data_nNNs,
      resolution = aggregation_GEX_cluster_res,
      dim_prefix = "PCA_",
      threads = 6,
      min_barcodes = aggregation_cluster_min_barcodes
    ),
    resources = get_tar_resources(cores_req = 6, RAM_GB_req = 32)
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
    name = cluster_UCell_controls.GEX,
    description = "Freeze abundance/detection-matched random signatures from a GEM-well-stratified GEX reference for all modalities",
    command = prepare_cluster_UCell_controls(
      counts_matrix = aggregated_counts_BPCells_matrix.GEX,
      metadata_tibble = metadata_w_clusters_tibble.GEX,
      markers = UCell_GEX_marker_genes_list
    ),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = cluster_UCell_evidence.GEX,
    description = "Cache GEX adjusted marker scores and diagnostic perturbations independently of the assignment threshold",
    command = prepare_cluster_UCell_evidence(
      counts_matrix = aggregated_counts_BPCells_matrix.GEX,
      metadata_tibble = metadata_w_clusters_tibble.GEX,
      control = cluster_UCell_controls.GEX,
      cluster_column = "PCA_harmony_SNN_cluster",
      include_cell_scores = TRUE,
      workers = 2
    ),
    resources = get_tar_resources(cores_req = 2, RAM_GB_req = 32)
  ),
  targets::tar_target(
    name = cluster_UCell_annotation.3_GEX_QC,
    description = "GEX assignments from minimum adjusted-score advantage; stability and GEM-well agreement are diagnostic only. [checkpoint:3_GEX-QC]",
    command = evaluate_cluster_UCell_evidence(cluster_UCell_evidence.GEX,
      min_advantage = aggregation_cluster_annotation_min_advantage),
    resources = get_tar_resources(RAM_GB_req = 8)
  ),
  tarchetypes::tar_file(
    name = cluster_UCell_advantage_plots.3_GEX_QC,
    description = "Faceted GEX adjusted-score bars with matched background and the minimum-advantage competition cutoff. [checkpoint:3_GEX-QC]",
    command = plot_cluster_UCell_advantages(cluster_UCell_annotation.3_GEX_QC) |>
      save_plots_structured(width = 22, height = 15)
  ),
  targets::tar_target(
    name = metadata_w_cell_types_unfiltered_tibble.GEX,
    description = "Attach supported GEX cluster labels or explicit abstentions and per-cell UCell scores before doublet removal",
    command = add_cluster_UCell_annotations(
      metadata_tibble = metadata_w_clusters_tibble.GEX,
      annotation = cluster_UCell_annotation.3_GEX_QC,
      cluster_column = "PCA_harmony_SNN_cluster"
    ),
    resources = get_tar_resources(cores_req = 6, RAM_GB_req = 8, RAM_GB_per_extra_core = 3)
  ),
  tarchetypes::tar_file(
    name = markers_by_cell_type_dot_plot.3_GEX_QC,
    description = "Pre-doublet-filtering marker expression per assigned GEX cell type. [checkpoint:3_GEX-QC]",
    command = {
      plot <- plot_marker_expression_dot_BPCells(
        feature_matrix = aggregated_counts_BPCells_matrix.GEX,
        metadata_tibble = metadata_w_cell_types_unfiltered_tibble.GEX,
        marker_genes_list = UCell_GEX_marker_genes_list,
        group_col = "PCA_harmony_SNN_cluster_cell_type", group_label = "GEX cell type"
      )
      save_plots_structured(plot, width = max(16, 4 + 0.35 * nlevels(plot$data$marker_feature)),
        height = max(9, 4 + 0.25 * nlevels(plot$data$group)))
    },
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  tarchetypes::tar_file(
    name = markers_by_cluster_dot_plot.3_GEX_QC,
    description = "Pre-doublet-filtering marker expression per named GEX cluster. [checkpoint:3_GEX-QC]",
    command = {
      plot <- plot_marker_expression_dot_BPCells(
        feature_matrix = aggregated_counts_BPCells_matrix.GEX,
        metadata_tibble = metadata_w_cell_types_unfiltered_tibble.GEX,
        marker_genes_list = UCell_GEX_marker_genes_list,
        group_col = "PCA_harmony_SNN_cluster_named",
        cell_type_col = "PCA_harmony_SNN_cluster_cell_type", group_label = "GEX cluster"
      )
      save_plots_structured(plot,
        width = max(16, 4 + 0.35 * nlevels(plot$data$marker_feature)),
        height = max(9, 4 + 0.25 * nlevels(plot$data$group)))
    },
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  tarchetypes::tar_file(
    name = module_scores_by_cell_type_dot_plot.3_GEX_QC,
    description = "Pre-doublet-filtering adjusted UCell evidence averaged per GEX cell type, with distance-ordered marker sets. [checkpoint:3_GEX-QC]",
    command = {
      plot <- plot_UCell_annotation_dot(cluster_UCell_annotation.3_GEX_QC,
        metadata_w_cell_types_unfiltered_tibble.GEX, group_by = "cell_type")
      save_plots_structured(plot, width = max(16, 4 + 0.35 * nlevels(plot$data$module)),
        height = max(9, 4 + 0.25 * nlevels(plot$data$cluster)))
    }
  ),
  tarchetypes::tar_file(
    name = module_scores_by_cluster_dot_plot.3_GEX_QC,
    description = "Cached pre-doublet-filtering cluster-adjusted UCell evidence with distance-ordered marker sets and matching row order. [checkpoint:3_GEX-QC]",
    command = {
      plot <- plot_UCell_annotation_dot(cluster_UCell_annotation.3_GEX_QC,
        metadata_w_cell_types_unfiltered_tibble.GEX, group_by = "cluster")
      save_plots_structured(plot,
        width = max(16, 4 + 0.35 * nlevels(plot$data$module)),
        height = max(9, 4 + 0.25 * nlevels(plot$data$cluster)))
    }
  ),
  tarchetypes::tar_file(
    name = marker_set_UCell_summary.3_GEX_QC,
    description = "Per-marker-set matched-background support, winning clusters and closest competition before GEX doublet filtering. [checkpoint:3_GEX-QC]",
    command = {
      output_dir <- get_structured_output_path(list_output = TRUE)
      fs::dir_create(output_dir)
      summary <- summarize_UCell_marker_sets(cluster_UCell_annotation.3_GEX_QC)
      readr::write_tsv(summary, file.path(output_dir, "marker_sets.tsv"))
      writeLines(c(
        "One row per marker set, using GEX clusters before doublet filtering. Each cluster has equal weight.",
        "Adjusted score = mean UCell minus the label-specific matched-control 95th percentile; original negative scores are retained.",
        "clusters_above_min_advantage counts sufficient background evidence only, not separation from competing labels.",
        "best_advantage = maximum across clusters of this set's adjusted score minus max(0, strongest other set's adjusted score).",
        "closest_competitor is the strongest other set in best_advantage_cluster, not the most correlated set across clusters.",
        "clusters_leading uses the deterministic nominated candidate; exact ties never count as assigned.",
        "No assignments can reflect weak evidence, competition or an absent population; this table does not validate biological identity."
      ), file.path(output_dir, "method.txt"))
      output_dir
    }
  ),
  tarchetypes::tar_file(
    name = cluster_UCell_diagnostics.3_GEX_QC,
    description = "Export GEX cluster decisions, per-label evidence, GEM-well agreement and control matching for review. [checkpoint:3_GEX-QC]",
    command = save_cluster_UCell_diagnostics(cluster_UCell_annotation.3_GEX_QC, cluster_UCell_controls.GEX)
  ),
  targets::tar_target(
    name = scDblFinder_GEM_well_tibble.GEX,
    description = "Prepare per GEM well GEX barcodes and cluster labels for scDblFinder",
    command = prepare_scDblFinder_GEM_well_tibble(
      metadata_tibble = metadata_w_cell_types_unfiltered_tibble.GEX,
      cluster_collapse_list = aggregation_scDblFinder_GEX_cell_type_collapse_list,
      cluster_col = "PCA_harmony_SNN_cluster_scDblFinder_group"
    ),
    iteration = "vector"
  ),
  targets::tar_target(
    name = scDblFinder_results_by_GEM_well_tibble.GEX,
    description = "Run GEX scDblFinder independently for each 10x Genomics GEM well from BPCells count slices",
    command = run_scDblFinder_BPCells_GEM_well(
      feature_matrix = aggregated_counts_BPCells_matrix.GEX,
      scDblFinder_GEM_well_tibble = scDblFinder_GEM_well_tibble.GEX,
      output_suffix = "GEX",
      dbr.sd = 1.0
    ),
    pattern = map(scDblFinder_GEM_well_tibble.GEX),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = scDblFinder_results_df.GEX,
    description = "Combine per GEM well GEX scDblFinder classifications [part_of_graph:GEX] [part_of_graph:seurat_export]",
    command = scDblFinder_results_by_GEM_well_tibble.GEX |>
      dplyr::bind_rows() |>
      dplyr::select(-dplyr::any_of("GEM_well_ID")) |>
      tibble::column_to_rownames("barcode_w_prefix"),
    resources = get_tar_resources(RAM_GB_req = 8)
  ),
  tarchetypes::tar_file(
    name = scDblFinder_score_violins_plot.3_GEX_QC,
    description = "Violin plots of scDblFinder doublet scores per cluster and cell type. [checkpoint:3_GEX-QC]",
    command = {
      scDblFinder_metadata <- metadata_w_cell_types_unfiltered_tibble.GEX |>
        dplyr::left_join(scDblFinder_results_df.GEX |> tibble::rownames_to_column("barcode_w_prefix"), by = "barcode_w_prefix")
      plot <- scDblFinder_metadata |>
        dplyr::select(PCA_harmony_SNN_cluster, PCA_harmony_SNN_cluster_cell_type, scDblFinder.score_GEX, scDblFinder.class_GEX, GEM_well_ID) |>
        tidyr::pivot_longer(cols = dplyr::all_of(c("PCA_harmony_SNN_cluster", "PCA_harmony_SNN_cluster_cell_type")), names_to = "cluster_type", values_to = "cluster_id") %>%
        ggplot2::ggplot(ggplot2::aes(x = cluster_id, y = scDblFinder.score_GEX)) +
        ggplot2::geom_violin(scale = "width") +
        ggplot2::geom_jitter(
          data = \(d) dplyr::filter(d, scDblFinder.class_GEX == "doublet"),
          ggplot2::aes(color = GEM_well_ID),
          size = 0.3,
          alpha = 0.5,
          width = 0.2
        ) +
        ggplot2::facet_wrap(~cluster_type, scales = "free_x") +
        ggplot2::labs(subtitle = "Points are cells classified as doublet by scDblFinder, colored by GEM well.") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
        ggplot2::theme(legend.position = "none")

      save_plots_structured(plot)
    },
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = metadata_w_cell_types_tibble.GEX,
    description = "Annotate GEX metadata with scDblFinder QC and apply configured cell- and cluster-level doublet filters [part_of_graph:ATAC] [part_of_graph:GEX] [part_of_graph:seurat_export] [checkpoint:3_GEX-QC]",
    command = filter_metadata_by_scDblFinder(
      metadata_tibble = metadata_w_cell_types_unfiltered_tibble.GEX,
      scDblFinder_results_df = scDblFinder_results_df.GEX,
      class_col = "scDblFinder.class_GEX",
      cluster_col = "PCA_harmony_SNN_cluster",
      remove_called_doublets = aggregation_scDblFinder_GEX_remove_called_doublets,
      max_doublet_fraction_per_cluster = aggregation_scDblFinder_GEX_max_doublet_fraction_per_cluster
    ),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = metadata_w_cell_types_analysis_tibble.GEX,
    description = "Join configured analysis variables onto processed GEX metadata",
    command = prepare_GEX_metadata_tibble(
      metadata_tibble = metadata_w_cell_types_tibble.GEX,
      barcode_vec = metadata_w_cell_types_tibble.GEX$barcode_w_prefix,
      donor_id_metadata_tibble = donor_id_analysis_metadata_tibble,
      GEM_well_metadata_tibble = GEM_well_analysis_metadata_tibble
    )
  ),
  targets::tar_target(
    name = metadata_w_cell_types_annotation_tibble.GEX,
    description = "Join complete donor and GEM well annotations onto processed GEX metadata",
    command = prepare_GEX_metadata_tibble(
      metadata_tibble = metadata_w_cell_types_tibble.GEX,
      barcode_vec = metadata_w_cell_types_tibble.GEX$barcode_w_prefix,
      donor_id_metadata_tibble = donor_id_metadata_tibble,
      GEM_well_metadata_tibble = GEM_well_annotation_metadata_tibble
    )
  ),
  targets::tar_target(
    name = cluster_marker_groups.GEX,
    description = "Prepare accepted cells and ordered cluster IDs for each multicluster GEX cell type",
    command = make_cluster_marker_groups(metadata_w_cell_types_tibble.GEX),
    iteration = "vector"
  ),
  targets::tar_target(
    name = cluster_marker_tibbles.GEX,
    description = "Test clusters against the remaining clusters within one cell type, with BH correction across genes per contrast",
    command = get_cluster_markers_from_matrix(
      feature_matrix = aggregated_counts_BPCells_matrix.GEX,
      group_record = cluster_marker_groups.GEX
    ),
    pattern = map(cluster_marker_groups.GEX),
    iteration = "list",
    resources = get_tar_resources(RAM_GB_req = 8)
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
  targets::tar_target(
    name = categorical_UMAP_var.GEX,
    description = "Categorical GEX UMAP variables to plot one at a time",
    command = aggregation_GEX_categorical_vars,
    iteration = "vector"
  ),
  targets::tar_target(
    name = harmony_continuous_UMAP_spec.GEX,
    description = "Continuous Harmony GEX UMAP variables to plot one at a time",
    command = {
      feature_vars <- intersect(
        unique(c(GEX_marker_genes_vec, interesting_genes, top_variable_genes.GEX)),
        rownames(aggregated_counts_BPCells_matrix.GEX)
      )
      tibble::tibble(
        variable = c(aggregation_non_peak_based_continuous_vars, feature_vars),
        value_source = c(
          rep("metadata", length(aggregation_non_peak_based_continuous_vars)),
          rep("feature", length(feature_vars))
        )
      )
    },
    iteration = "vector"
  ),
  targets::tar_target(
    name = non_harmony_continuous_UMAP_spec.GEX,
    description = "Continuous non-Harmony GEX UMAP variables to plot one at a time",
    command = {
      feature_vars <- intersect(
        unique(c(GEX_marker_genes_vec, interesting_genes)),
        rownames(aggregated_counts_BPCells_matrix.GEX)
      )
      tibble::tibble(
        variable = c(aggregation_non_peak_based_continuous_vars, feature_vars),
        value_source = c(
          rep("metadata", length(aggregation_non_peak_based_continuous_vars)),
          rep("feature", length(feature_vars))
        )
      )
    },
    iteration = "vector"
  ),
  tarchetypes::tar_file(
    name = cross.UMAPs.3_GEX_QC,
    description = "Compute GEX UMAPs across a sweep of PC counts and neighbour counts. [checkpoint:3_GEX-QC]",
    command = {
      sweep_umap <- run_UMAP_from_embedding_matrix(
        embedding_matrix = harmony_embeddings_matrix.GEX,
        dims = seq_len(UMAP_n_dims_seq.GEX),
        n_neighbors = UMAP_neighbors_seq,
        min_dist = aggregation_UMAP_min_dist,
        dim_prefix = "PCA_",
        col_prefix = "GEX_UMAP",
        threads = 6
      ) |>
        tibble::as_tibble(rownames = "barcode_w_prefix")

      metadata_w_cell_types_tibble.GEX |>
        dplyr::select(-dplyr::any_of(c("GEX_UMAP_1", "GEX_UMAP_2"))) |>
        dplyr::left_join(sweep_umap, by = "barcode_w_prefix") |>
        plot_UMAP_from_metadata(variable = "PCA_harmony_SNN_cluster_cell_type", umap_cols = c("GEX_UMAP_1", "GEX_UMAP_2"))
    } |>
      save_plots_structured(dyn_suffix_in_subdir = TRUE, override_suffix = paste0(UMAP_n_dims_seq.GEX, "_", UMAP_neighbors_seq)),
    pattern = cross(UMAP_n_dims_seq.GEX, UMAP_neighbors_seq),
    resources = get_tar_resources(cores_req = 6, RAM_GB_req = 32)
  ),
  tarchetypes::tar_file(
    name = harmony.categorical.UMAPs.3_GEX_QC,
    description = "UMAPs colored by categorical metadata variables. [checkpoint:3_GEX-QC]",
    command = metadata_w_cell_types_analysis_tibble.GEX |>
      plot_UMAP_from_metadata(variable = categorical_UMAP_var.GEX, umap_cols = c("GEX_UMAP_1", "GEX_UMAP_2")) |>
      save_plots_structured(
        dyn_suffix_in_subdir = TRUE,
        override_suffix = stringr::str_replace_all(categorical_UMAP_var.GEX, "[/\\\\]", "_")
      ),
    pattern = map(categorical_UMAP_var.GEX),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  tarchetypes::tar_file(
    name = harmony.continuous.UMAPs.3_GEX_QC,
    description = "UMAPs colored by continuous QC and gene expression features. [checkpoint:3_GEX-QC]",
    command = plot_UMAP_from_metadata(
      metadata_tibble = metadata_w_cell_types_analysis_tibble.GEX,
      variable = harmony_continuous_UMAP_spec.GEX$variable,
      value_source = harmony_continuous_UMAP_spec.GEX$value_source,
      feature_matrix = aggregated_counts_BPCells_matrix.GEX,
      umap_cols = c("GEX_UMAP_1", "GEX_UMAP_2")
    ) |>
      save_plots_structured(
        dyn_suffix_in_subdir = TRUE,
        override_suffix = stringr::str_replace_all(harmony_continuous_UMAP_spec.GEX$variable, "[/\\\\]", "_")
      ),
    pattern = map(harmony_continuous_UMAP_spec.GEX),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  tarchetypes::tar_file(
    name = non_harmony.categorical.UMAPs.3_GEX_QC,
    description = "UMAPs colored by categorical metadata variables on the uncorrected PCA embedding. [checkpoint:3_GEX-QC]",
    command = metadata_w_cell_types_analysis_tibble.GEX |>
      plot_UMAP_from_metadata(variable = categorical_UMAP_var.GEX, umap_cols = c("GEX_non_harmony_UMAP_1", "GEX_non_harmony_UMAP_2")) |>
      save_plots_structured(
        dyn_suffix_in_subdir = TRUE,
        override_suffix = stringr::str_replace_all(categorical_UMAP_var.GEX, "[/\\\\]", "_")
      ),
    pattern = map(categorical_UMAP_var.GEX),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  tarchetypes::tar_file(
    name = non_harmony.continuous.UMAPs.3_GEX_QC,
    description = "UMAPs colored by continuous features on the uncorrected PCA embedding. [checkpoint:3_GEX-QC]",
    command = plot_UMAP_from_metadata(
      metadata_tibble = metadata_w_cell_types_analysis_tibble.GEX,
      variable = non_harmony_continuous_UMAP_spec.GEX$variable,
      value_source = non_harmony_continuous_UMAP_spec.GEX$value_source,
      feature_matrix = aggregated_counts_BPCells_matrix.GEX,
      umap_cols = c("GEX_non_harmony_UMAP_1", "GEX_non_harmony_UMAP_2")
    ) |>
      save_plots_structured(
        dyn_suffix_in_subdir = TRUE,
        override_suffix = stringr::str_replace_all(non_harmony_continuous_UMAP_spec.GEX$variable, "[/\\\\]", "_")
      ),
    pattern = map(non_harmony_continuous_UMAP_spec.GEX),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  tarchetypes::tar_file(
    name = continuous_by_cell_type_violin_plot.3_GEX_QC,
    description = "Violin plots of continuous QC and cell-cycle features per cell type. [checkpoint:3_GEX-QC]",
    command = plot_QC_metric_violins(
      metadata_tibble = metadata_w_cell_types_analysis_tibble.GEX,
      QC_metric_manifest_tibble = QC_metric_manifest_tibble,
      checkpoints = c("1_pre-aggregation-QC", "2_GEX-PCA-QC", "3_GEX-QC"),
      group_col = "PCA_harmony_SNN_cluster_cell_type"
    ) |>
      save_plots_structured(width = max(10, 4 + 0.35 * dplyr::n_distinct(metadata_w_cell_types_analysis_tibble.GEX$PCA_harmony_SNN_cluster_cell_type)))
  ),
  tarchetypes::tar_file(
    name = continuous_by_cluster_violin_plot.3_GEX_QC,
    description = "Violin plots of continuous QC and cell-cycle features per SNN cluster. [checkpoint:3_GEX-QC]",
    command = plot_QC_metric_violins(
      metadata_tibble = metadata_w_cell_types_analysis_tibble.GEX,
      QC_metric_manifest_tibble = QC_metric_manifest_tibble,
      checkpoints = c("1_pre-aggregation-QC", "2_GEX-PCA-QC", "3_GEX-QC"),
      group_col = "PCA_harmony_SNN_cluster_named"
    ) |>
      save_plots_structured(width = max(10, 4 + 0.35 * dplyr::n_distinct(metadata_w_cell_types_analysis_tibble.GEX$PCA_harmony_SNN_cluster_named)))
  ),
  tarchetypes::tar_file(
    name = categorical_bars_plots.3_GEX_QC,
    description = "Bar plots of categorical metadata composition per cell type. [checkpoint:3_GEX-QC]",
    command = plot_categorical_bars_plot(
      metadata_tibble = metadata_w_cell_types_analysis_tibble.GEX,
      metadata_cols = aggregation_GEX_categorical_vars,
      cluster_col = "PCA_harmony_SNN_cluster_cell_type"
    ) |>
      save_plots_structured()
  ),
  tarchetypes::tar_file(
    name = cluster_marker_volcano_plots.3_GEX_QC,
    description = "One cluster-marker volcano file per multicluster cell type; two clusters share one directed comparison. [checkpoint:3_GEX-QC]",
    command = plot_cluster_marker_volcano(cluster_marker_tibbles.GEX) |>
      save_plots_structured(
        dyn_suffix_in_subdir = TRUE,
        override_suffix = cluster_marker_tibbles.GEX$cell_type
      ),
    pattern = map(cluster_marker_tibbles.GEX)
  ),
  tarchetypes::tar_file(
    name = cell_type_marker_volcano_plots.3_GEX_QC,
    description = "Facetted volcano plot of BPCells marker genes per GEX cell type. [checkpoint:3_GEX-QC]",
    command = cell_type_marker_tibbles.GEX |>
      plot_markers_volcano_simple() |>
      save_plots_structured()
  ),
  targets::tar_target(
    name = GEX_Seurat_object.3_GEX_QC,
    description = "Build a GEX-only Seurat compatibility object for manual review before ATAC peak calling. [checkpoint:3_GEX-QC]",
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
