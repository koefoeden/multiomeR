rlang::list2(
  targets::tar_target(
    name = cell_retention_tibble.WNN,
    description = "Accumulate per-well retention counts through final WNN filtering, including empty wells. [checkpoint:8_multimodal-QC]",
    command = summarize_QC_cell_retention(
      before_metadata = metadata_w_cell_types_tibble.ATAC,
      after_metadata = metadata_w_cell_types_tibble.WNN,
      GEM_well_IDs = aggregation_GEM_well_IDs,
      stage = "checkpoint:8_multimodal-QC",
      discarded_barcodes = list("Small WNN clusters" = setdiff(
        clusters_tibble_raw.WNN$barcode_w_prefix, clusters_tibble.WNN$barcode_w_prefix)),
      previous_stages = cell_retention_tibble.ATAC
    )
  ),
  tarchetypes::tar_file(
    name = cell_retention_flow_plot.8_multimodal_QC,
    description = "Plot cumulative nuclei retention through final WNN QC. [checkpoint:8_multimodal-QC]",
    command = save_QC_cell_retention_plot(cell_retention_tibble.WNN)
  ),
  tarchetypes::tar_file(
    name = nuclei_per_donor_id_bars.8_multimodal_QC,
    description = "Plot the nuclei per donor that pass all QC into the final WNN object. [checkpoint:8_multimodal-QC]",
    command = metadata_w_cell_types_tibble.WNN |>
      dplyr::count(donor_id, name = "n_nuclei") |>
      plot_nuclei_per_donor_id(
        title = "Nuclei per donor after QC",
        caption = "Nuclei retained in the final WNN object after all per-well and aggregation-level QC."
      ) |>
      save_plots_structured(width = 10,
        height = max(6, 3 + 0.22 * dplyr::n_distinct(metadata_w_cell_types_tibble.WNN$donor_id)))
  ),
  WNN_processing_targets = rlang::list2(
    targets::tar_target(
      name = embedding_matrices.WNN,
      description = "Align GEX and ATAC Harmony embeddings for native WNN construction [part_of_graph:WNN] [part_of_graph:seurat_export]",
      command = get_WNN_embedding_matrices(
        GEX_embedding_matrix = harmony_embeddings_matrix.GEX,
        ATAC_embedding_matrix = harmony_embeddings_matrix.ATAC,
        barcode_vec = metadata_w_cell_types_tibble.ATAC$barcode_w_prefix,
        GEX_dims = aggregation_GEX_data_PCs,
        ATAC_dims = aggregation_ATAC_data_PCs
      ),
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    targets::tar_target(
      name = WNN_results_raw,
      description = "Compute native weighted nearest-neighbour graph from aligned GEX and ATAC embeddings before small-cluster filtering, keeping twice the neighbour count so that filtering needs no new search",
      command = weighted_nearest_neighbors_BPCells(
        embeddings_list = embedding_matrices.WNN,
        k = aggregation_data_nNNs,
        candidate_k = 200,
        threads = 6,
        extra_neighbors = aggregation_data_nNNs,
        native_source_file = WNN_native_source_file
      ),
      resources = get_tar_resources(cores_req = 6, RAM_GB_req = 16)
    ),
    targets::tar_target(
      name = clusters_tibble_raw.WNN,
      description = "Cluster cells from the native WNN SNN graph before small-cluster filtering [resource_observation: peak_RAM_GiB=10.8; n_nuclei=213k; aggregation=mixed_human_31x; method=slurm_monitor; date=2026-10-10] [resource_observation: peak_RAM_GiB=5.8; n_nuclei=713k; aggregation=unpublished; method=slurm_monitor; date=2026-10-10]",
      command = cluster_WNN_graph(
        WNN_results = restrict_WNN_results(WNN_results_raw, k = aggregation_data_nNNs),
        resolution = aggregation_WNN_cluster_res,
        min_barcodes = NULL
      ),
      resources = get_tar_resources(RAM_GB_req = 32)
    ),
    targets::tar_target(
      name = clusters_tibble.WNN,
      description = "Filter raw WNN clusters by the minimum barcode threshold [part_of_graph:WNN] [part_of_graph:seurat_export]",
      command = {
        clusters <- clusters_tibble_raw.WNN$WNN_harmony_SNN_cluster
        names(clusters) <- clusters_tibble_raw.WNN$barcode_w_prefix
        clusters <- filter_clusters_by_min_barcodes(
          clusters,
          min_barcodes = aggregation_cluster_min_barcodes
        )

        tibble::tibble(
          barcode_w_prefix = names(clusters),
          WNN_harmony_SNN_cluster = clusters
        )
      }
    ),
    targets::tar_target(
      name = WNN_results,
      description = "Native weighted nearest-neighbour graph restricted to retained WNN cells: each keeps its first neighbours among them, and cells left with too few are reselected among them [part_of_graph:WNN] [part_of_graph:seurat_export]",
      command = restrict_WNN_results(
        WNN_results_raw,
        barcodes = clusters_tibble.WNN$barcode_w_prefix,
        k = aggregation_data_nNNs,
        embeddings_list = embedding_matrices.WNN
      ),
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    targets::tar_target(
      name = UMAP_embeddings_tibble.WNN,
      description = "Compute WNN UMAP coordinates from the retained native weighted neighbour index [part_of_graph:WNN] [part_of_graph:seurat_export]",
      command = run_WNN_UMAP(
        WNN_results = WNN_results,
        n_neighbors = aggregation_UMAP_nNNs,
        min_dist = aggregation_UMAP_min_dist,
        threads = 6
      ),
      resources = get_tar_resources(cores_req = 6, RAM_GB_req = 16)
    ),
    targets::tar_target(
      name = metadata_w_clusters_tibble.WNN,
      description = "Join native WNN weights and clusters onto metadata",
      command = build_WNN_metadata_tibble(
        metadata_tibble = metadata_w_cell_types_tibble.ATAC,
        WNN_results = WNN_results,
        WNN_clusters_tibble = clusters_tibble.WNN
      ),
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    targets::tar_target(
      name = cluster_UCell_controls.WNN,
      description = "Freeze matched random signatures for the WNN marker sets from the same GEX reference",
      command = prepare_cluster_UCell_controls(
        counts_matrix = aggregated_counts_BPCells_matrix.GEX,
        metadata_tibble = metadata_w_clusters_tibble.GEX,
        markers = UCell_WNN_marker_genes_list
      ),
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    targets::tar_target(
      name = cluster_UCell_evidence.WNN,
      description = "Cache adjusted GEX marker scores for WNN clusters independently of stringency",
      command = prepare_cluster_UCell_evidence(
        counts_matrix = aggregated_counts_BPCells_matrix.GEX,
        metadata_tibble = metadata_w_clusters_tibble.WNN,
        control = cluster_UCell_controls.WNN,
        cluster_column = "WNN_harmony_SNN_cluster",
        workers = 2,
        native_source_file = UCell_native_source_file
      ),
      resources = get_tar_resources(cores_req = 2, RAM_GB_req = 16)
    ),
    targets::tar_target(
      name = cluster_UCell_annotation.8_multimodal_QC,
      description = "WNN assignments from minimum adjusted GEX score advantage. [checkpoint:8_multimodal-QC]",
      command = evaluate_cluster_UCell_evidence(cluster_UCell_evidence.WNN,
        min_advantage = aggregation_WNN_cluster_annotation_min_advantage),
      resources = get_tar_resources(RAM_GB_req = 8)
    ),
    tarchetypes::tar_file(
      name = cluster_UCell_advantage_plots.8_multimodal_QC,
      description = "Faceted WNN adjusted GEX scores with matched background and the competition cutoff. [checkpoint:8_multimodal-QC]",
      command = plot_cluster_UCell_advantages(cluster_UCell_annotation.8_multimodal_QC) |>
        save_plots_structured(width = 22, height = 15, dpi = 150)
    ),
    targets::tar_target(
      name = metadata_w_cell_types_tibble.WNN,
      description = "Attach supported WNN cluster labels or explicit abstentions from GEX UCell evidence [part_of_graph:WNN] [part_of_graph:seurat_export] [checkpoint:8_multimodal-QC]",
      command = add_cluster_UCell_annotations(
        metadata_tibble = metadata_w_clusters_tibble.WNN,
        annotation = cluster_UCell_annotation.8_multimodal_QC,
        cluster_column = "WNN_harmony_SNN_cluster"
      ),
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    tarchetypes::tar_file(
      name = cluster_UCell_diagnostics.8_multimodal_QC,
      description = "Export WNN cluster decisions, GEX marker evidence, GEM-well agreement and control matching. [checkpoint:8_multimodal-QC]",
      command = save_cluster_UCell_diagnostics(cluster_UCell_annotation.8_multimodal_QC, cluster_UCell_controls.WNN)
    ),
    targets::tar_target(
      name = metadata_w_cell_types_analysis_tibble.WNN,
      description = "Join configured analysis variables and the WNN, GEX and ATAC UMAP coordinates onto processed WNN metadata for plots",
      command = prepare_GEX_metadata_tibble(
        metadata_tibble = metadata_w_cell_types_tibble.WNN,
        barcode_vec = metadata_w_cell_types_tibble.WNN$barcode_w_prefix,
        donor_id_metadata_tibble = donor_id_analysis_metadata_tibble,
        GEM_well_metadata_tibble = GEM_well_analysis_metadata_tibble
      ) |>
        dplyr::left_join(UMAP_embeddings_tibble.WNN, by = "barcode_w_prefix") |>
        dplyr::left_join(UMAP_embeddings_tibble.GEX, by = "barcode_w_prefix") |>
        dplyr::left_join(UMAP_embeddings_tibble.ATAC, by = "barcode_w_prefix"),
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    targets::tar_target(
      name = metadata_w_cell_types_annotation_tibble.WNN,
      description = "Join complete donor and GEM well annotations onto processed WNN metadata",
      command = prepare_GEX_metadata_tibble(
        metadata_tibble = metadata_w_cell_types_tibble.WNN,
        barcode_vec = metadata_w_cell_types_tibble.WNN$barcode_w_prefix,
        donor_id_metadata_tibble = donor_id_metadata_tibble,
        GEM_well_metadata_tibble = GEM_well_annotation_metadata_tibble
      ),
      resources = get_tar_resources(RAM_GB_req = 16)
    )
  ),
  WNN_plot_targets = rlang::list2(
    targets::tar_target(
      name = categorical_UMAP_var.WNN,
      description = "Categorical WNN UMAP variables to plot one at a time",
      command = aggregation_WNN_categorical_vars,
      iteration = "vector"
    ),
    targets::tar_target(
      name = continuous_UMAP_spec.WNN,
      description = "Continuous WNN UMAP variables to plot one at a time",
      command = {
        feature_vars <- intersect(GEX_marker_genes_vec, rownames(aggregated_counts_BPCells_matrix.GEX))
        tibble::tibble(
          variable = c(aggregation_w_WNN_continuous_vars, feature_vars),
          value_source = c(
            rep("metadata", length(aggregation_w_WNN_continuous_vars)),
            rep("feature", length(feature_vars))
          )
        )
      },
      iteration = "vector"
    ),
    targets::tar_target(
      name = continuous_UMAP_inputs.WNN,
      description = "Coordinates, metadata variables and gene counts plotted on continuous WNN UMAPs, with genes read in one pass",
      command = prepare_continuous_UMAP_inputs(metadata_w_cell_types_analysis_tibble.WNN, continuous_UMAP_spec.WNN,
        aggregated_counts_BPCells_matrix.GEX, umap_cols = c("WNN_UMAP_1", "WNN_UMAP_2")),
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    tarchetypes::tar_file(
      name = categorical.UMAPs.8_multimodal_QC,
      description = "UMAPs colored by categorical metadata variables on the WNN embedding. [checkpoint:8_multimodal-QC]",
      command = metadata_w_cell_types_analysis_tibble.WNN |>
        plot_UMAP_from_metadata(
          variable = categorical_UMAP_var.WNN,
          umap_cols = c("WNN_UMAP_1", "WNN_UMAP_2")
        ) |>
        save_plots_structured(
          dyn_suffix_in_subdir = TRUE,
          override_suffix = stringr::str_replace_all(categorical_UMAP_var.WNN, "[/\\\\]", "_")
        ),
      pattern = map(categorical_UMAP_var.WNN),
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    tarchetypes::tar_file(
      name = continuous.UMAPs.8_multimodal_QC,
      description = "UMAPs colored by continuous QC and gene expression features on the WNN embedding. [checkpoint:8_multimodal-QC]",
      command = plot_UMAP_from_metadata(
        metadata_tibble = continuous_UMAP_inputs.WNN$metadata_tibble,
        variable = continuous_UMAP_spec.WNN$variable,
        value_source = continuous_UMAP_spec.WNN$value_source,
        feature_matrix = continuous_UMAP_inputs.WNN$feature_matrix,
        umap_cols = c("WNN_UMAP_1", "WNN_UMAP_2")
      ) |>
        save_plots_structured(
          dyn_suffix_in_subdir = TRUE,
          override_suffix = stringr::str_replace_all(continuous_UMAP_spec.WNN$variable, "[/\\\\]", "_")
        ),
      pattern = map(continuous_UMAP_spec.WNN)
    ),
    tarchetypes::tar_file(
      name = continuous_by_cell_type_violin_plot.8_multimodal_QC,
      description = "Violin plots of continuous QC and cell-cycle features per cell type among all WNN nuclei. [checkpoint:8_multimodal-QC]",
      command = plot_QC_metric_violins(
        metadata_tibble = metadata_w_cell_types_analysis_tibble.WNN,
        QC_metric_manifest_tibble = QC_metric_manifest_tibble,
        checkpoints = c("1_pre-aggregation-QC", "2_GEX-PCA-QC", "3_GEX-QC", "4_peak-QC", "5_pre-LSI-QC", "6_ATAC-LSI-QC", "7_ATAC-QC", "8_multimodal-QC"),
        group_col = "WNN_harmony_SNN_cluster_cell_type"
      ) |>
        save_plots_structured(width = max(10, 4 + 0.35 * dplyr::n_distinct(metadata_w_cell_types_analysis_tibble.WNN$WNN_harmony_SNN_cluster_cell_type)))
    ),
    tarchetypes::tar_file(
      name = continuous_by_cluster_violin_plot.8_multimodal_QC,
      description = "Violin plots of continuous QC and cell-cycle features per SNN cluster among all WNN nuclei. [checkpoint:8_multimodal-QC]",
      command = plot_QC_metric_violins(
        metadata_tibble = metadata_w_cell_types_analysis_tibble.WNN,
        QC_metric_manifest_tibble = QC_metric_manifest_tibble,
        checkpoints = c("1_pre-aggregation-QC", "2_GEX-PCA-QC", "3_GEX-QC", "4_peak-QC", "5_pre-LSI-QC", "6_ATAC-LSI-QC", "7_ATAC-QC", "8_multimodal-QC"),
        group_col = "WNN_harmony_SNN_cluster_named",
        group_order = get_marker_cell_type_order(
          metadata_w_cell_types_tibble.WNN$WNN_harmony_SNN_cluster_named,
          names(UCell_WNN_marker_genes_list),
          get_marker_group_cell_types(metadata_w_cell_types_tibble.WNN,
            "WNN_harmony_SNN_cluster_named", "WNN_harmony_SNN_cluster_cell_type"))
      ) |>
        save_plots_structured(width = max(10, 4 + 0.35 * dplyr::n_distinct(metadata_w_cell_types_analysis_tibble.WNN$WNN_harmony_SNN_cluster_named)))
    ),
    tarchetypes::tar_file(
      name = categorical_by_cell_type_bars_plots.8_multimodal_QC,
      description = "Bar plots of categorical metadata composition per cell type among all WNN nuclei. [checkpoint:8_multimodal-QC]",
      command = plot_categorical_bars_plot(
        metadata_tibble = metadata_w_cell_types_analysis_tibble.WNN,
        metadata_cols = aggregation_WNN_categorical_vars,
        cluster_col = "WNN_harmony_SNN_cluster_cell_type"
      ) |>
        save_plots_structured()
    ),
    tarchetypes::tar_file(
      name = categorical_by_cluster_bars_plots.8_multimodal_QC,
      description = "Bar plots of categorical metadata composition per SNN cluster among all WNN nuclei. [checkpoint:8_multimodal-QC]",
      command = plot_categorical_bars_plot(
        metadata_tibble = metadata_w_cell_types_analysis_tibble.WNN,
        metadata_cols = aggregation_WNN_categorical_vars,
        cluster_col = "WNN_harmony_SNN_cluster_named",
        group_order = get_marker_cell_type_order(
          metadata_w_cell_types_tibble.WNN$WNN_harmony_SNN_cluster_named,
          names(UCell_WNN_marker_genes_list),
          get_marker_group_cell_types(metadata_w_cell_types_tibble.WNN,
            "WNN_harmony_SNN_cluster_named", "WNN_harmony_SNN_cluster_cell_type"))
      ) |>
        save_plots_structured(height = max(9, 4 + 0.25 * dplyr::n_distinct(metadata_w_cell_types_analysis_tibble.WNN$WNN_harmony_SNN_cluster_named)))
    ),
    tarchetypes::tar_file(
      name = cluster_composition_table.8_multimodal_QC,
      description = "Tabulate each SNN cluster's top category and categories above 5% per configured categorical variable among all WNN nuclei. [checkpoint:8_multimodal-QC]",
      command = save_cluster_composition_table(metadata_w_cell_types_analysis_tibble.WNN,
        "WNN_harmony_SNN_cluster_named", aggregation_proj_spec_categorical_vars)
    ),
    tarchetypes::tar_file(
      name = markers_by_cluster_dot_plot.8_multimodal_QC,
      description = "GEX marker expression per named WNN cluster, with nuclei counts and doublet evidence for retained WNN nuclei. [checkpoint:8_multimodal-QC]",
      command = {
        plot <- plot_marker_expression_dot_BPCells(
          feature_matrix = aggregated_counts_BPCells_matrix.GEX,
          metadata_tibble = metadata_w_cell_types_tibble.WNN,
          marker_genes_list = UCell_WNN_marker_genes_list,
          group_col = "WNN_harmony_SNN_cluster_named",
          cell_type_col = "WNN_harmony_SNN_cluster_cell_type", group_label = "WNN cluster"
        )
        save_plots_structured(
          add_cluster_doublet_bars(plot, metadata_w_cell_types_tibble.WNN,
            scDblFinder_results_df.GEX, max_doublet_fraction = NULL,
            group_col = "WNN_harmony_SNN_cluster_named"),
          width = max(20, 8 + 0.35 * nlevels(plot$data$marker_feature)),
          height = max(9, 4 + 0.25 * nlevels(plot$data$group)))
      },
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    tarchetypes::tar_file(
      name = module_scores_by_cluster_dot_plot.8_multimodal_QC,
      description = "Cached adjusted GEX UCell evidence per WNN cluster, with distance-ordered marker sets. [checkpoint:8_multimodal-QC]",
      command = {
        plot <- plot_UCell_annotation_dot(cluster_UCell_annotation.8_multimodal_QC,
          metadata_w_cell_types_tibble.WNN, group_by = "cluster",
          cluster_column = "WNN_harmony_SNN_cluster", group_label = "WNN cluster",
          cell_scores = prepare_cluster_UCell_evidence(aggregated_counts_BPCells_matrix.GEX,
            metadata_w_cell_types_tibble.WNN, cluster_UCell_controls.WNN, "WNN_harmony_SNN_cluster",
            include_cell_scores = TRUE, native_source_file = UCell_native_source_file)$cell_scores)
        save_plots_structured(plot,
          width = max(16, 4 + 0.35 * nlevels(plot$data$module)),
          height = max(9, 4 + 0.25 * nlevels(plot$data$cluster)))
      },
      resources = get_tar_resources(cores_req = 2, RAM_GB_req = 16)
    ),
    tarchetypes::tar_file(
      name = confusion_matrices_plots.8_multimodal_QC,
      description = "Paired RNA versus WNN and ATAC versus WNN confusion matrices for SNN clusters and cell types. [checkpoint:8_multimodal-QC]",
      command = {
        pairwise_comparison_tibble <- tibble::tribble(
          ~source_label, ~target_label,
          "RNA",         "WNN",
          "ATAC",        "WNN"
        )

        plots <- purrr::pmap(
          pairwise_comparison_tibble,
          \(source_label, target_label) plot_modality_confusion_matrices(
            metadata_tibble = metadata_w_cell_types_tibble.WNN,
            source_label = source_label,
            target_label = target_label
          )
        )
        names(plots) <- paste(pairwise_comparison_tibble$source_label, pairwise_comparison_tibble$target_label, sep = "_vs_")

        save_plots_structured(plots, width = 22, height = 11)
      },
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    targets::tar_target(
      name = WNN_weight_metadata_summary.WNN,
      description = "Summarize retained-cell ATAC-weight associations, support and distributions, pooled and within named WNN clusters",
      command = get_WNN_weight_metadata_summary(
        metadata = metadata_w_cell_types_analysis_tibble.WNN,
        continuous_vars = aggregation_continuous_vars %||% character(),
        categorical_vars = aggregation_categorical_vars %||% character()
      ),
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    tarchetypes::tar_file(
      name = WNN_weight_metadata_associations_plot.8_multimodal_QC,
      description = "Plot continuous and categorical WNN-weight associations pooled and within named WNN clusters. [checkpoint:8_multimodal-QC]",
      command = plot_WNN_weight_metadata_associations(WNN_weight_metadata_summary.WNN) |>
        save_plots_structured(),
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    tarchetypes::tar_file(
      name = cluster_named_dim_tri_plot.UMAPs.8_multimodal_QC,
      description = "3×3 grid of UMAPs with cell-type-named cluster-level identities. [checkpoint:8_multimodal-QC]",
      command = metadata_w_cell_types_analysis_tibble.WNN |>
        plot_3_by_3_clusters_and_reduction_UMAPs_from_metadata(
          cluster_col_suffix = "named"
        ) |>
        save_plots_structured(width = 18, height = 18),
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    tarchetypes::tar_file(
      name = cluster_cell_type_dim_tri_plot.UMAPs.8_multimodal_QC,
      description = "3×3 grid of UMAPs with cell-type identities. [checkpoint:8_multimodal-QC]",
      command = metadata_w_cell_types_analysis_tibble.WNN |>
        plot_3_by_3_clusters_and_reduction_UMAPs_from_metadata(
          cluster_col_suffix = "cell_type"
        ) |>
        save_plots_structured(width = 18, height = 18),
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    targets::tar_target(
      name = cross_UMAP_embeddings_tibble.WNN,
      description = "Compute WNN UMAPs across a sweep of nNN counts [resource_observation: peak_RAM_GiB=6.6; n_nuclei=213k; aggregation=mixed_human_31x; method=slurm_monitor; date=2026-10-10]",
      command = run_WNN_UMAP(
        WNN_results = WNN_results,
        n_neighbors = UMAP_neighbors_seq,
        min_dist = aggregation_UMAP_min_dist,
        threads = 6
      ),
      pattern = map(UMAP_neighbors_seq),
      resources = get_tar_resources(cores_req = 6, RAM_GB_req = 16)
    ),
    tarchetypes::tar_file(
      name = cross.UMAPs.8_multimodal_QC,
      description = "Plot WNN cell types on UMAPs across a sweep of nNN counts. [checkpoint:8_multimodal-QC]",
      command = {
        metadata_w_cell_types_tibble.WNN |>
          dplyr::left_join(cross_UMAP_embeddings_tibble.WNN, by = "barcode_w_prefix") |>
          plot_UMAP_from_metadata(
            variable = "WNN_harmony_SNN_cluster_cell_type",
            umap_cols = c("WNN_UMAP_1", "WNN_UMAP_2")
          ) + ggplot2::labs(subtitle = sprintf("Neighbours: %s. Compare label stability across the WNN UMAP sweep; shared GEX-derived labels are not independent validation.", UMAP_neighbors_seq))
      } |>
        save_plots_structured(
          dyn_suffix_in_subdir = TRUE,
          override_suffix = paste0(UMAP_neighbors_seq)
        ),
      pattern = map(cross_UMAP_embeddings_tibble.WNN, UMAP_neighbors_seq),
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
  ),
)
