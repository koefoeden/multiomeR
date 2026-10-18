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
    description = "Leiden clusters from the SNN graph of the Harmony-corrected GEX PCA nearest neighbours [part_of_graph:GEX] [part_of_graph:seurat_export] [resource_observation: peak_RAM_GiB=7.3; n_nuclei=713k; aggregation=unpublished; method=slurm_monitor; date=2026-10-10]",
    command = cluster_neighbors(
      knn = knn.GEX,
      k = aggregation_data_nNNs,
      resolution = aggregation_GEX_cluster_res,
      min_barcodes = aggregation_cluster_min_barcodes
    ),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = metadata_w_clusters_tibble.GEX,
    description = "BPCells-native GEX metadata with PCA clusters and UMAP coordinates",
    command = metadata_tibble.GEX |>
      dplyr::filter(.data$barcode_w_prefix %in% names(PCA_clusters.GEX)) |>
      dplyr::left_join(UMAP_embeddings_tibble.GEX, by = "barcode_w_prefix") |>
      dplyr::left_join(UMAP_embeddings_tibble.GEX_non_harmony, by = "barcode_w_prefix") |>
      dplyr::mutate(PCA_harmony_SNN_cluster = unname(PCA_clusters.GEX[.data$barcode_w_prefix])),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = cluster_UCell_controls.GEX,
    description = "Freeze abundance/detection-matched random signatures from a GEM-well-stratified GEX reference for GEX and ATAC labels",
    command = prepare_cluster_UCell_controls(
      counts_matrix = aggregated_counts_BPCells_matrix.GEX,
      metadata_tibble = metadata_w_clusters_tibble.GEX,
      markers = UCell_GEX_marker_genes_list
    ),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = cluster_UCell_evidence.GEX,
    description = "Cache GEX adjusted marker scores independently of the assignment threshold",
    command = prepare_cluster_UCell_evidence(
      counts_matrix = aggregated_counts_BPCells_matrix.GEX,
      metadata_tibble = metadata_w_clusters_tibble.GEX,
      control = cluster_UCell_controls.GEX,
      cluster_column = "PCA_harmony_SNN_cluster",
      include_cell_scores = TRUE,
      workers = 2,
      native_source_file = UCell_native_source_file
    ),
    resources = get_tar_resources(cores_req = 2, RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = cluster_UCell_annotation.3_GEX_QC,
    description = "GEX assignments from minimum adjusted-score advantage. [checkpoint:3_GEX-QC]",
    command = evaluate_cluster_UCell_evidence(cluster_UCell_evidence.GEX,
      min_advantage = aggregation_GEX_cluster_annotation_min_advantage),
    resources = get_tar_resources(RAM_GB_req = 8)
  ),
  tarchetypes::tar_file(
    name = cluster_UCell_advantage_plots.3_GEX_QC,
    description = "Faceted GEX adjusted-score bars with matched background and the minimum-advantage competition cutoff. [checkpoint:3_GEX-QC]",
    command = plot_cluster_UCell_advantages(cluster_UCell_annotation.3_GEX_QC) |>
      save_plots_structured(width = 22, height = 15, dpi = 150)
  ),
  targets::tar_target(
    name = metadata_w_cell_types_unfiltered_tibble.GEX,
    description = "Attach supported GEX cluster labels or explicit abstentions and per-cell UCell scores before doublet removal [resource_observation: peak_RAM_GiB=1.8; n_nuclei=213k; aggregation=mixed_human_31x; method=slurm_monitor; date=2026-10-10]",
    command = add_cluster_UCell_annotations(
      metadata_tibble = metadata_w_clusters_tibble.GEX,
      annotation = cluster_UCell_annotation.3_GEX_QC,
      cluster_column = "PCA_harmony_SNN_cluster"
    ),
    resources = get_tar_resources(cores_req = 6, RAM_GB_req = 16)
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
      save_plots_structured(
        add_cluster_doublet_bars(plot, metadata_w_cell_types_unfiltered_tibble.GEX,
          scDblFinder_results_df.GEX, aggregation_scDblFinder_GEX_max_doublet_fraction_per_cluster),
        width = max(20, 8 + 0.35 * nlevels(plot$data$marker_feature)),
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
      cluster_col = "PCA_harmony_SNN_cluster_cell_type"
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
    command = make_cluster_marker_groups(metadata_w_cell_types_tibble.GEX)
  ),
  targets::tar_target(
    name = cluster_marker_tibbles.GEX,
    description = "Test clusters against the remaining clusters within each multicluster cell type, with BH correction across genes per contrast",
    # One target rather than branches: an aggregation can have no multicluster cell type.
    command = purrr::map(seq_len(nrow(cluster_marker_groups.GEX)), \(i) get_cluster_markers_from_matrix(
      feature_matrix = aggregated_counts_BPCells_matrix.GEX,
      group_record = cluster_marker_groups.GEX[i, ]
    )),
    resources = get_tar_resources(RAM_GB_req = 8)
  ),
  targets::tar_target(
    name = cell_type_marker_tibbles.GEX,
    description = "Compute BPCells marker genes for each GEX cell type against all other cell types [resource_observation: peak_RAM_GiB=5.4; n_nuclei=213k; aggregation=mixed_human_31x; method=slurm_monitor; date=2026-10-09] [resource_observation: peak_RAM_GiB=14.2; n_nuclei=713k; aggregation=unpublished; method=slurm_monitor; date=2026-10-10]",
    command = get_BPCells_markers_from_matrix(
      feature_matrix = aggregated_counts_BPCells_matrix.GEX,
      metadata_tibble = metadata_w_cell_types_tibble.GEX,
      group_col = "PCA_harmony_SNN_cluster_cell_type"
    ),
    resources = get_tar_resources(RAM_GB_req = 32)
  ),
  targets::tar_target(
    name = UMAP_n_dims_seq.GEX,
    description = "Generate a sequence of counts of leading UMAP PCs for cross-parameter UMAP exploration",
    command = round(seq(5, length(aggregation_GEX_data_PCs), length.out = 3))
  ),
  targets::tar_target(
    name = UMAP_neighbors_seq,
    description = "Generate a sequence of neighbor counts for cross-parameter UMAP exploration",
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
    name = cross_UMAP_embeddings_tibbles.GEX,
    description = "Compute GEX UMAPs across a sweep of PC counts and neighbor counts, with one neighbour search per PC count",
    command = {
      embedding_matrix <- select_embedding_dimensions(harmony_embeddings_matrix.GEX,
        dims = aggregation_GEX_data_PCs[seq_len(UMAP_n_dims_seq.GEX)], dim_prefix = "PCA_")
      # The full PC count is the shared search's own; UMAP_neighbors_seq never exceeds its neighbour count.
      knn <- if (UMAP_n_dims_seq.GEX == length(aggregation_GEX_data_PCs)) knn.GEX else
        find_hnsw_neighbors(embedding_matrix, k = max(UMAP_neighbors_seq), threads = 6)
      purrr::map(purrr::set_names(UMAP_neighbors_seq), \(n_neighbors) {
        run_UMAP_from_neighbors(knn, n_neighbors = n_neighbors, min_dist = aggregation_UMAP_min_dist,
          col_prefix = "GEX_UMAP", init = "spca", embedding_matrix = embedding_matrix, threads = 6) |>
          tibble::as_tibble(rownames = "barcode_w_prefix")
      })
    },
    pattern = map(UMAP_n_dims_seq.GEX),
    resources = get_tar_resources(cores_req = 6, RAM_GB_req = 16)
  ),
  tarchetypes::tar_file(
    name = cross.UMAPs.3_GEX_QC,
    description = "Plot GEX cell types on UMAPs across a sweep of PC counts and neighbor counts. [checkpoint:3_GEX-QC]",
    command = {
      metadata <- metadata_w_cell_types_tibble.GEX |>
        dplyr::select(-dplyr::any_of(c("GEX_UMAP_1", "GEX_UMAP_2")))
      cross_UMAP_embeddings_tibbles.GEX |>
        purrr::imap(\(sweep_umap, n_neighbors) {
          dplyr::left_join(metadata, sweep_umap, by = "barcode_w_prefix") |>
            plot_UMAP_from_metadata(variable = "PCA_harmony_SNN_cluster_cell_type", umap_cols = c("GEX_UMAP_1", "GEX_UMAP_2")) +
            ggplot2::labs(subtitle = sprintf("First %s UMAP PCs; neighbours: %s; min_dist: %s. Compare label stability across settings; island spacing is not biological distance.",
              UMAP_n_dims_seq.GEX, n_neighbors, aggregation_UMAP_min_dist))
        }) |>
        purrr::set_names(\(n_neighbors) paste0(UMAP_n_dims_seq.GEX, "_", n_neighbors)) |>
        save_plots_structured(dyn_suffix_in_subdir = TRUE, override_suffix = UMAP_n_dims_seq.GEX)
    },
    pattern = map(cross_UMAP_embeddings_tibbles.GEX, UMAP_n_dims_seq.GEX),
    resources = get_tar_resources(RAM_GB_req = 16)
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
  targets::tar_target(
    name = continuous_UMAP_inputs.GEX,
    description = "Coordinates, metadata variables and gene counts plotted on continuous Harmony and non-Harmony GEX UMAPs, with genes read in one pass",
    command = prepare_continuous_UMAP_inputs(metadata_w_cell_types_analysis_tibble.GEX,
      # The non-Harmony plots show a subset of these variables.
      harmony_continuous_UMAP_spec.GEX,
      aggregated_counts_BPCells_matrix.GEX,
      umap_cols = c("GEX_UMAP_1", "GEX_UMAP_2", "GEX_non_harmony_UMAP_1", "GEX_non_harmony_UMAP_2")),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  tarchetypes::tar_file(
    name = harmony.continuous.UMAPs.3_GEX_QC,
    description = "UMAPs colored by continuous QC and gene expression features. [checkpoint:3_GEX-QC]",
    command = plot_UMAP_from_metadata(
      metadata_tibble = continuous_UMAP_inputs.GEX$metadata_tibble,
      variable = harmony_continuous_UMAP_spec.GEX$variable,
      value_source = harmony_continuous_UMAP_spec.GEX$value_source,
      feature_matrix = continuous_UMAP_inputs.GEX$feature_matrix,
      umap_cols = c("GEX_UMAP_1", "GEX_UMAP_2")
    ) |>
      save_plots_structured(
        dyn_suffix_in_subdir = TRUE,
        override_suffix = stringr::str_replace_all(harmony_continuous_UMAP_spec.GEX$variable, "[/\\\\]", "_")
      ),
    pattern = map(harmony_continuous_UMAP_spec.GEX)
  ),
  targets::tar_target(
    name = metadata_w_cell_types_unfiltered_analysis_tibble.GEX,
    description = "Join doublet scores and configured analysis variables onto GEX metadata before doublet filtering",
    command = prepare_GEX_metadata_tibble(
      metadata_tibble = metadata_w_cell_types_unfiltered_tibble.GEX |>
        dplyr::left_join(
          scDblFinder_results_df.GEX |> tibble::rownames_to_column("barcode_w_prefix"),
          by = "barcode_w_prefix"),
      barcode_vec = metadata_w_cell_types_unfiltered_tibble.GEX$barcode_w_prefix,
      donor_id_metadata_tibble = donor_id_analysis_metadata_tibble,
      GEM_well_metadata_tibble = GEM_well_analysis_metadata_tibble
    )
  ),
  tarchetypes::tar_file(
    name = continuous_by_cell_type_violin_plot.3_GEX_QC,
    description = "Violin plots of continuous QC and cell-cycle features per cell type before doublet filtering. [checkpoint:3_GEX-QC]",
    command = plot_QC_metric_violins(
      metadata_tibble = metadata_w_cell_types_unfiltered_analysis_tibble.GEX,
      QC_metric_manifest_tibble = QC_metric_manifest_tibble,
      checkpoints = c("1_pre-aggregation-QC", "2_GEX-PCA-QC", "3_GEX-QC"),
      group_col = "PCA_harmony_SNN_cluster_cell_type"
    ) |>
      save_plots_structured(width = max(10, 4 + 0.35 * dplyr::n_distinct(metadata_w_cell_types_unfiltered_analysis_tibble.GEX$PCA_harmony_SNN_cluster_cell_type)))
  ),
  tarchetypes::tar_file(
    name = continuous_by_cluster_violin_plot.3_GEX_QC,
    description = "Violin plots of continuous QC and cell-cycle features per SNN cluster before doublet filtering. [checkpoint:3_GEX-QC]",
    command = plot_QC_metric_violins(
      metadata_tibble = metadata_w_cell_types_unfiltered_analysis_tibble.GEX,
      QC_metric_manifest_tibble = QC_metric_manifest_tibble,
      checkpoints = c("1_pre-aggregation-QC", "2_GEX-PCA-QC", "3_GEX-QC"),
      group_col = "PCA_harmony_SNN_cluster_named",
      group_order = get_marker_cell_type_order(
        metadata_w_cell_types_unfiltered_tibble.GEX$PCA_harmony_SNN_cluster_named,
        names(UCell_GEX_marker_genes_list),
        get_marker_group_cell_types(metadata_w_cell_types_unfiltered_tibble.GEX,
          "PCA_harmony_SNN_cluster_named", "PCA_harmony_SNN_cluster_cell_type"))
    ) |>
      save_plots_structured(width = max(10, 4 + 0.35 * dplyr::n_distinct(metadata_w_cell_types_unfiltered_analysis_tibble.GEX$PCA_harmony_SNN_cluster_named)))
  ),
  tarchetypes::tar_file(
    name = categorical_by_cell_type_bars_plots.3_GEX_QC,
    description = "Bar plots of categorical metadata composition per cell type before doublet filtering. [checkpoint:3_GEX-QC]",
    command = plot_categorical_bars_plot(
      metadata_tibble = metadata_w_cell_types_unfiltered_analysis_tibble.GEX,
      metadata_cols = aggregation_GEX_categorical_vars,
      cluster_col = "PCA_harmony_SNN_cluster_cell_type"
    ) |>
      save_plots_structured()
  ),
  tarchetypes::tar_file(
    name = categorical_by_cluster_bars_plots.3_GEX_QC,
    description = "Bar plots of categorical metadata composition per SNN cluster before doublet filtering. [checkpoint:3_GEX-QC]",
    command = plot_categorical_bars_plot(
      metadata_tibble = metadata_w_cell_types_unfiltered_analysis_tibble.GEX,
      metadata_cols = aggregation_GEX_categorical_vars,
      cluster_col = "PCA_harmony_SNN_cluster_named",
      group_order = get_marker_cell_type_order(
        metadata_w_cell_types_unfiltered_tibble.GEX$PCA_harmony_SNN_cluster_named,
        names(UCell_GEX_marker_genes_list),
        get_marker_group_cell_types(metadata_w_cell_types_unfiltered_tibble.GEX,
          "PCA_harmony_SNN_cluster_named", "PCA_harmony_SNN_cluster_cell_type"))
    ) |>
      save_plots_structured(height = max(9, 4 + 0.25 * dplyr::n_distinct(metadata_w_cell_types_unfiltered_analysis_tibble.GEX$PCA_harmony_SNN_cluster_named)))
  ),
  tarchetypes::tar_file(
    name = cluster_composition_table.3_GEX_QC,
    description = "Tabulate each SNN cluster's top category and categories above 5% per configured categorical variable before doublet filtering. [checkpoint:3_GEX-QC]",
    command = save_cluster_composition_table(metadata_w_cell_types_unfiltered_analysis_tibble.GEX,
      "PCA_harmony_SNN_cluster_named", aggregation_proj_spec_categorical_vars)
  ),
  tarchetypes::tar_file(
    name = cluster_marker_volcano_plots.3_GEX_QC,
    description = "One cluster-marker volcano file per multicluster cell type; two clusters share one directed comparison. [checkpoint:3_GEX-QC]",
    command = purrr::map(cluster_marker_tibbles.GEX, plot_cluster_marker_volcano) |>
      purrr::set_names(cluster_marker_groups.GEX$cell_type) |>
      save_plots_structured()
  ),
  tarchetypes::tar_file(
    name = cell_type_marker_volcano_plots.3_GEX_QC,
    description = "Facetted volcano plot of BPCells marker genes per GEX cell type. [checkpoint:3_GEX-QC]",
    command = cell_type_marker_tibbles.GEX |>
      plot_markers_volcano_simple() |>
      (\(plot) plot + ggplot2::labs(title = "GEX markers by assigned cell type",
        caption = stringr::str_wrap(paste(plot$labels$caption,
          "Wilcoxon comparisons against pooled remaining cells; BH adjustment spans all returned gene-by-group tests. Cell-level tests do not establish donor-level replication."), width = 110)))() |>
      save_plots_structured()
  ),
  targets::tar_target(
    name = GEX_Seurat_object.3_GEX_QC,
    description = "Build a GEX-only Seurat compatibility object for manual review before ATAC peak calling. [checkpoint:3_GEX-QC] [resource_observation: peak_RAM_GiB=6.5; n_nuclei=213k; aggregation=mixed_human_31x; method=slurm_monitor; date=2026-10-10] [resource_observation: peak_RAM_GiB=9.3; n_nuclei=713k; aggregation=unpublished; method=slurm_monitor; date=2026-10-10] [resource_observation: peak_RAM_GiB=45.8; n_nuclei=376k; aggregation=comparison_31x; backend=Seurat_SCT; method=slurm_monitor; date=2026-10-10]",
    command = build_seurat_signac_convenience_object(
      GEX_counts_matrix = aggregated_counts_BPCells_matrix.GEX,
      metadata_tibble = metadata_w_cell_types_tibble.GEX,
      PCA_results = PCA_BPCells.GEX,
      GEX_normalization = GEX_normalization.GEX,
      GEX_harmony_embeddings = harmony_embeddings_matrix.GEX,
      GEX_UMAP_embeddings_tibble = UMAP_embeddings_tibble.GEX,
      GEX_non_harmony_UMAP_embeddings_tibble = UMAP_embeddings_tibble.GEX_non_harmony,
      gene_features_df = gene_features_df,
      GEX_dims = aggregation_GEX_data_PCs,
      data_nNNs = aggregation_data_nNNs,
      graph_threads = 15
    ),
    resources = get_tar_resources(RAM_GB_req = if ("Seurat_SCT" %in% aggregation_tibble$aggregation_GEX_PCA_backend) 60 else 16),
    packages = w_def(c("Seurat", "SeuratObject", "BPCells"))
  ),
)
