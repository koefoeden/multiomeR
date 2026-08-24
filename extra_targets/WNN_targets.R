rlang::list2(
  WNN_processing_targets = rlang::list2(
    targets::tar_target(
      name = embedding_matrices.WNN,
      description = "Align GEX and ATAC Harmony embeddings for native WNN construction [part_of_graph:WNN] [part_of_graph:seurat_export]",
      command = get_WNN_embedding_matrices(
        GEX_embedding_matrix = harmony_embeddings_matrix.GEX,
        ATAC_embedding_matrix = harmony_embeddings_matrix.ATAC,
        barcode_vec = metadata_w_cell_types_tibble.ATAC$barcode_w_prefix,
        GEX_dims = aggregation_GEX_data_PCs,
        ATAC_dims = aggregation_ATAC_data_PCs,
        GEX_dim_prefix = "PCA_",
        ATAC_dim_prefix = "LSI_"
      ),
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    targets::tar_target(
      name = WNN_results_raw,
      description = "Compute native weighted nearest-neighbour graph from aligned GEX and ATAC embeddings before small-cluster filtering",
      command = weighted_nearest_neighbors_BPCells(
        embeddings_list = embedding_matrices.WNN,
        k = aggregation_data_nNNs,
        candidate_k = 200,
        threads = 6
      ),
      resources = get_tar_resources(cores_req = 6, RAM_GB_req = 60)
    ),
    targets::tar_target(
      name = clusters_tibble_raw.WNN,
      description = "Cluster cells from the native WNN SNN graph before small-cluster filtering",
      command = cluster_WNN_graph(
        WNN_results = WNN_results_raw,
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
      description = "Native weighted nearest-neighbour graph restricted to retained WNN cells [part_of_graph:WNN] [part_of_graph:seurat_export]",
      command = {
        retained_barcodes <- rownames(WNN_results_raw$nn_idx)[
          rownames(WNN_results_raw$nn_idx) %in%
            clusters_tibble.WNN$barcode_w_prefix
        ]

        if (length(retained_barcodes) == nrow(WNN_results_raw$nn_idx)) {
          WNN_results_raw
        } else {
          retained_embeddings <- purrr::map(
            embedding_matrices.WNN,
            \(embedding_matrix) {
              embedding_matrix[retained_barcodes, , drop = FALSE]
            }
          )

          weighted_nearest_neighbors_BPCells(
            embeddings_list = retained_embeddings,
            k = aggregation_data_nNNs,
            candidate_k = 200,
            threads = 6
          )
        }
      },
      resources = get_tar_resources(cores_req = 6, RAM_GB_req = 60)
    ),
    targets::tar_target(
      name = UMAP_embeddings_tibble.WNN,
      description = "Compute WNN UMAP coordinates from the retained native weighted neighbour index [part_of_graph:WNN] [part_of_graph:seurat_export]",
      command = run_WNN_UMAP(
        WNN_results = WNN_results,
        n_neighbors = aggregation_UMAP_nNNs,
        min_dist = aggregation_UMAP_min_dist
      ),
      resources = get_tar_resources(cores_req = 6, RAM_GB_req = 32)
    ),
    targets::tar_target(
      name = metadata_w_clusters_tibble.WNN,
      description = "Join native WNN weights, UMAP coordinates, and clusters onto metadata",
      command = build_WNN_metadata_tibble(
        metadata_tibble = metadata_w_cell_types_tibble.ATAC,
        WNN_results = WNN_results,
        UMAP_embeddings_tibble = UMAP_embeddings_tibble.WNN,
        WNN_clusters_tibble = clusters_tibble.WNN,
        GEX_UMAP_embeddings_tibble = UMAP_embeddings_tibble.GEX
      ),
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    targets::tar_target(
      name = metadata_w_cell_types_tibble.WNN,
      description = "Assign cell-type labels to native WNN clusters by GEX UCell marker scores [part_of_graph:WNN] [part_of_graph:seurat_export]",
      command = add_cell_types_to_metadata_from_module_scores(
        metadata_tibble = metadata_w_clusters_tibble.WNN,
        named_marker_genes_list = UCell_GEX_marker_genes_list,
        allow_multiple_cell_types = aggregation_allow_multiple_cell_types,
        cluster_column = "WNN_harmony_SNN_cluster" 
      ),
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    targets::tar_target(
      name = metadata_w_cell_types_analysis_tibble.WNN,
      description = "Join configured analysis variables onto processed WNN metadata",
      command = prepare_GEX_metadata_tibble(
        metadata_tibble = metadata_w_cell_types_tibble.WNN,
        barcode_vec = metadata_w_cell_types_tibble.WNN$barcode_w_prefix,
        donor_id_metadata_tibble = donor_id_analysis_metadata_tibble,
        GEM_well_metadata_tibble = GEM_well_analysis_metadata_tibble
      ),
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    targets::tar_target(
      name = metadata_w_cell_types_subgroup_tibble.WNN,
      description = "Join only configured subgroup model variables onto processed WNN metadata",
      command = prepare_GEX_metadata_tibble(
        metadata_tibble = metadata_w_cell_types_tibble.WNN,
        barcode_vec = metadata_w_cell_types_tibble.WNN$barcode_w_prefix,
        donor_id_metadata_tibble = donor_id_subgroup_metadata_tibble,
        GEM_well_metadata_tibble = GEM_well_subgroup_metadata_tibble
      ),
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
    tarchetypes::tar_file(
      name = categorical.UMAPs.WNN,
      description = "UMAPs colored by categorical metadata variables on the WNN embedding. [checkpoint:multimodal]",
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
      name = continuous.UMAPs.WNN,
      description = "UMAPs colored by continuous QC and gene expression features on the WNN embedding. [checkpoint:multimodal]",
      command = plot_UMAP_from_metadata(
        metadata_tibble = metadata_w_cell_types_analysis_tibble.WNN,
        variable = continuous_UMAP_spec.WNN$variable,
        value_source = continuous_UMAP_spec.WNN$value_source,
        feature_matrix = aggregated_counts_BPCells_matrix.GEX,
        umap_cols = c("WNN_UMAP_1", "WNN_UMAP_2")
      ) |>
        save_plots_structured(
          dyn_suffix_in_subdir = TRUE,
          override_suffix = stringr::str_replace_all(continuous_UMAP_spec.WNN$variable, "[/\\\\]", "_")
        ),
      pattern = map(continuous_UMAP_spec.WNN),
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    tarchetypes::tar_file(
      name = categorical_bars_plots.WNN,
      description = "Bar plots of categorical metadata composition per WNN cell type. [checkpoint:multimodal]",
      command = plot_categorical_bars_plot(
        metadata_tibble = metadata_w_cell_types_analysis_tibble.WNN,
        metadata_cols = aggregation_WNN_categorical_vars,
        cluster_col = "WNN_harmony_SNN_cluster_cell_type"
      ) |>
        save_plots_structured(),
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    tarchetypes::tar_file(
      name = markers_violin_plot.WNN,
      description = "Violin plots of marker gene expression per WNN cell type. [checkpoint:multimodal]",
      command = plot_WNN_marker_expression_violins(
        metadata_tibble = metadata_w_cell_types_tibble.WNN,
        feature_matrix = aggregated_counts_BPCells_matrix.GEX,
        marker_genes = GEX_marker_genes_vec
      ) |>
        save_plots_structured(),
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    tarchetypes::tar_file(
      name = confusion_matrices_plots.WNN,
      description = "Pairwise row-normalized confusion matrices comparing RNA, ATAC, and WNN cluster assignments. [checkpoint:multimodal]",
      command = {
        pairwise_comparison_tibble <- tibble::tribble(
          ~plot_name    , ~source_label , ~target_label , ~source_cluster_col             , ~target_cluster_col             , ~source_cell_type_col               , ~target_cell_type_col               ,
          "RNA_vs_ATAC" , "RNA"         , "ATAC"        , "PCA_harmony_SNN_cluster_named" , "LSI_harmony_SNN_cluster_named" , "PCA_harmony_SNN_cluster_cell_type" , "LSI_harmony_SNN_cluster_cell_type" ,
          "RNA_vs_WNN"  , "RNA"         , "WNN"         , "PCA_harmony_SNN_cluster_named" , "WNN_harmony_SNN_cluster_named" , "PCA_harmony_SNN_cluster_cell_type" , "WNN_harmony_SNN_cluster_cell_type" ,
          "ATAC_vs_WNN" , "ATAC"        , "WNN"         , "LSI_harmony_SNN_cluster_named" , "WNN_harmony_SNN_cluster_named" , "LSI_harmony_SNN_cluster_cell_type" , "WNN_harmony_SNN_cluster_cell_type"
        )

        plots <- purrr::pmap(
          pairwise_comparison_tibble,
          \(
            plot_name,
            source_label,
            target_label,
            source_cluster_col,
            target_cluster_col,
            source_cell_type_col,
            target_cell_type_col
          ) {
            SNN_cluster_plot <- plot_cluster_confusion_matrix(
              metadata_tibble = metadata_w_cell_types_tibble.WNN,
              source_col = source_cluster_col,
              target_col = target_cluster_col,
              source_label = paste(source_label, "SNN clusters"),
              target_label = paste(target_label, "SNN clusters"),
              title = "SNN clusters"
            )
            cell_type_plot <- plot_cluster_confusion_matrix(
              metadata_tibble = metadata_w_cell_types_tibble.WNN,
              source_col = source_cell_type_col,
              target_col = target_cell_type_col,
              source_label = paste(source_label, "cell type"),
              target_label = paste(target_label, "cell type"),
              title = "Cell type"
            )

            patchwork::wrap_plots(
              SNN_cluster_plot,
              cell_type_plot,
              nrow = 1,
              guides = "collect"
            ) +
              patchwork::plot_annotation(
                title = paste(source_label, "vs", target_label)
              )
          }
        )
        names(plots) <- pairwise_comparison_tibble$plot_name

        save_plots_structured(plots, width = 16, height = 7)
      },
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    tarchetypes::tar_file(
      name = ATAC_vs_RNA_weight_boxplots_plot.WNN,
      description = "Boxplots of ATAC vs RNA modality weights per cell type. [checkpoint:multimodal]",
      command = metadata_w_cell_types_tibble.WNN |>
        plot_ATAC_vs_RNA_weight_boxplots() |>
        save_plots_structured(),
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    tarchetypes::tar_file(
      name = cluster_named_dim_tri_plot.WNN,
      description = "3×3 grid of UMAPs with cell-type-named cluster-level identities. [checkpoint:multimodal]",
      command = metadata_w_cell_types_tibble.WNN |>
        plot_3_by_3_clusters_and_reduction_UMAPs_from_metadata(
          cluster_col_suffix = "named"
        ) |>
        save_plots_structured(),
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    tarchetypes::tar_file(
      name = cluster_cell_type_dim_tri_plot.WNN,
      description = "3×3 grid of UMAPs with cell-type identities. [checkpoint:multimodal]",
      command = metadata_w_cell_types_tibble.WNN |>
        plot_3_by_3_clusters_and_reduction_UMAPs_from_metadata(
          cluster_col_suffix = "cell_type"
        ) |>
        save_plots_structured(),
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    tarchetypes::tar_file(
      name = cross.UMAPs.WNN,
      description = "Compute WNN UMAPs across a sweep of nNN counts. [checkpoint:multimodal]",
      command = {
        sweep_umap_tibble <- run_WNN_UMAP(
          WNN_results = WNN_results,
          n_neighbors = UMAP_neighbors_seq,
          min_dist = aggregation_UMAP_min_dist
        )

        metadata_w_cell_types_tibble.WNN |>
          dplyr::select(-dplyr::any_of(c("WNN_UMAP_1", "WNN_UMAP_2"))) |>
          dplyr::left_join(sweep_umap_tibble, by = "barcode_w_prefix") |>
          plot_UMAP_from_metadata(
            variable = "WNN_harmony_SNN_cluster_cell_type",
            umap_cols = c("WNN_UMAP_1", "WNN_UMAP_2")
          )
      } |>
        save_plots_structured(
          dyn_suffix_in_subdir = TRUE,
          override_suffix = paste0(UMAP_neighbors_seq)
        ),
      pattern = map(UMAP_neighbors_seq),
      resources = get_tar_resources(cores_req = 6, RAM_GB_req = 32)
    ),
  ),
)
