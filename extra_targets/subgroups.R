rlang::list2(
  targets::tar_target(
    name = subgroups_to_process_vec,
    description = "Identify cell-type clusters with enough nuclei to process as subgroups [part_of_graph:full_subgroups]",
    # TODO: Generalize this to work with any cluster column - however, seems to be slightly problematic with target-errors regarding unknown object.
    # Probably some tricky interaction with NSE.
    command = {
      filtered_nuclei_tibble <- metadata_w_cell_types_tibble.WNN %>%
        dplyr::group_by(.data[["PCA_harmony_SNN_cluster_cell_type"]]) %>%
        # group_by({{ aggregation_subgroups_col }}) %>%
        dplyr::summarise(n_nuclei = dplyr::n()) %>%
        dplyr::filter(n_nuclei >= aggregation_subgroups_min_nuclei_filter)

      assert_with_info(
        nrow(filtered_nuclei_tibble) > 0,
        glue_info = "No clusters with at least {aggregation_subgroups_min_nuclei_filter} nuclei found. Consider increasing the aggregation_subgroups_min_nuclei_filter parameter."
      )
      filtered_nuclei_tibble %>%
        dplyr::pull(dplyr::all_of("PCA_harmony_SNN_cluster_cell_type")) %>%
        as.character()
    }
  ),
  targets::tar_target(
    name = metadata_tibble.subgroups,
    description = "Subset WNN metadata to cells matching each subgroup cluster value [part_of_graph:full_subgroups]",
    command = metadata_w_cell_types_tibble.WNN |>
      filter_metadata_tibble_by_col_match(
        column_name = aggregation_subgroups_col,
        column_values_pattern = subgroups_to_process_vec
      ),
    pattern = map(subgroups_to_process_vec),
    resources = get_tar_resources(RAM_GB_req = 8)
  ),
  targets::tar_target(
    name = PCA_BPCells.GEX.subgroups,
    description = "Run Seurat SCTransform on BPCells-backed subgroup GEX counts and compute PCA [part_of_graph:full_subgroups]",
    command = run_GEX_PCA_BPCells(
      GEX_counts_matrix = aggregated_counts_BPCells_matrix.GEX,
      metadata_tibble = metadata_tibble.subgroups,
      organism_chr = organism_chr,
      GEX_PCA_backend = aggregation_GEX_PCA_backend,
      SCT_regress_vars = aggregation_subgroups_SCT_regress_vars,
      n_components = utils::tail(aggregation_GEX_data_PCs, 1),
      threads = 6
    ),
    pattern = map(metadata_tibble.subgroups),
    resources = get_tar_resources(cores_req = 6, RAM_GB_req = 32)
  ),
  targets::tar_target(
    name = harmony_embeddings_matrix.GEX.subgroups,
    description = "Harmony-corrected SCTransform subgroup GEX PCA embeddings [part_of_graph:full_subgroups]",
    command = run_harmony_on_embedding_matrix(
      embedding_matrix = PCA_BPCells.GEX.subgroups$cell_embeddings,
      metadata_tibble = metadata_tibble.subgroups,
      harmony_correction_metadata_col_names = c(aggregation_harmony_correction_metadata_col_names, aggregation_subgroups_extra_harmony_covars),
      dims = aggregation_GEX_data_PCs,
      cores = 6,
      dim_prefix = "PCA_"
    ),
    pattern = map(PCA_BPCells.GEX.subgroups, metadata_tibble.subgroups),
    resources = get_tar_resources(cores_req = 6, RAM_GB_req = 32)
  ),
  targets::tar_target(
    name = LSI_BPCells.ATAC.subgroups,
    description = "Run BPCells-native TF-IDF and SVD on subgroup ATAC peak counts [part_of_graph:full_subgroups]",
    command = {
      subgroup_barcodes <- intersect(metadata_tibble.subgroups$barcode_w_prefix, colnames(peak_QC_filtered_BPCells_matrix.ATAC))
      run_ATAC_LSI_BPCells(
        ATAC_peak_BPCells_matrix = peak_QC_filtered_BPCells_matrix.ATAC[, subgroup_barcodes, drop = FALSE],
        n_components = utils::tail(aggregation_ATAC_data_PCs, 1),
        threads = 6
      )
    },
    pattern = map(metadata_tibble.subgroups),
    resources = get_tar_resources(cores_req = 6, RAM_GB_req = 32)
  ),
  targets::tar_target(
    name = harmony_embeddings_matrix.ATAC.subgroups,
    description = "Harmony-corrected BPCells-native subgroup ATAC LSI embeddings [part_of_graph:full_subgroups]",
    command = run_harmony_on_embedding_matrix(
      embedding_matrix = LSI_BPCells.ATAC.subgroups$cell_embeddings,
      metadata_tibble = metadata_tibble.subgroups,
      harmony_correction_metadata_col_names = aggregation_harmony_correction_metadata_col_names,
      dims = aggregation_ATAC_data_PCs,
      cores = 6,
      dim_prefix = "LSI_"
    ),
    pattern = map(LSI_BPCells.ATAC.subgroups, metadata_tibble.subgroups),
    resources = get_tar_resources(cores_req = 6, RAM_GB_req = 32)
  ),
  targets::tar_target(
    name = UMAP_embeddings_tibble.GEX.subgroups,
    description = "Compute subgroup GEX UMAP coordinates from recomputed GEX PCA [part_of_graph:full_subgroups]",
    command = run_UMAP_from_embedding_matrix(
      embedding_matrix = harmony_embeddings_matrix.GEX.subgroups,
      dims = aggregation_UMAP_GEX_PCs,
      n_neighbors = aggregation_UMAP_nNNs,
      min_dist = aggregation_UMAP_min_dist,
      dim_prefix = "PCA_",
      col_prefix = "GEX_sub_UMAP"
    ) |>
      tibble::as_tibble(rownames = "barcode_w_prefix"),
    pattern = map(harmony_embeddings_matrix.GEX.subgroups),
    resources = get_tar_resources(cores_req = 6, RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = PCA_clusters.GEX.subgroups,
    description = "Cluster subgroup cells from recomputed GEX PCA [part_of_graph:full_subgroups]",
    command = cluster_embedding_matrix_BPCells(
      embedding_matrix = harmony_embeddings_matrix.GEX.subgroups,
      dims = aggregation_GEX_data_PCs,
      k = aggregation_data_nNNs,
      resolution = aggregation_subgroups_res,
      dim_prefix = "PCA_",
      threads = 6,
      min_barcodes = aggregation_cluster_min_barcodes
    ),
    pattern = map(harmony_embeddings_matrix.GEX.subgroups),
    resources = get_tar_resources(cores_req = 6, RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = UMAP_embeddings_tibble.ATAC.subgroups,
    description = "Compute subgroup ATAC UMAP coordinates from recomputed ATAC LSI [part_of_graph:full_subgroups]",
    command = run_UMAP_from_embedding_matrix(
      embedding_matrix = harmony_embeddings_matrix.ATAC.subgroups,
      dims = aggregation_UMAP_ATAC_PCs,
      n_neighbors = aggregation_UMAP_nNNs,
      min_dist = aggregation_UMAP_min_dist,
      dim_prefix = "LSI_",
      col_prefix = "ATAC_sub_UMAP"
    ) |>
      tibble::as_tibble(rownames = "barcode_w_prefix"),
    pattern = map(harmony_embeddings_matrix.ATAC.subgroups),
    resources = get_tar_resources(cores_req = 6, RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = LSI_clusters.ATAC.subgroups,
    description = "Cluster subgroup cells from recomputed ATAC LSI [part_of_graph:full_subgroups]",
    command = cluster_embedding_matrix_BPCells(
      embedding_matrix = harmony_embeddings_matrix.ATAC.subgroups,
      dims = aggregation_ATAC_data_PCs,
      k = aggregation_data_nNNs,
      resolution = aggregation_subgroups_res,
      dim_prefix = "LSI_",
      threads = 6,
      min_barcodes = aggregation_cluster_min_barcodes
    ),
    pattern = map(harmony_embeddings_matrix.ATAC.subgroups),
    resources = get_tar_resources(cores_req = 6, RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = embedding_matrices.subgroups,
    description = "Align recomputed subgroup GEX PCA and ATAC LSI embeddings for native WNN construction [part_of_graph:full_subgroups]",
    command = get_WNN_embedding_matrices(
      GEX_embedding_matrix = harmony_embeddings_matrix.GEX.subgroups,
      ATAC_embedding_matrix = harmony_embeddings_matrix.ATAC.subgroups,
      barcode_vec = metadata_tibble.subgroups$barcode_w_prefix,
      GEX_dims = aggregation_GEX_data_PCs,
      ATAC_dims = aggregation_ATAC_data_PCs,
      GEX_dim_prefix = "PCA_",
      ATAC_dim_prefix = "LSI_"
    ),
    pattern = map(harmony_embeddings_matrix.GEX.subgroups, harmony_embeddings_matrix.ATAC.subgroups, metadata_tibble.subgroups),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = WNN_results.subgroups,
    description = "Compute subgroup-native weighted nearest-neighbour graph from aligned GEX and ATAC embeddings [part_of_graph:full_subgroups]",
    command = weighted_nearest_neighbors_BPCells(
      embeddings_list = embedding_matrices.subgroups,
      k = aggregation_data_nNNs,
      candidate_k = 200,
      threads = 6
    ),
    pattern = map(embedding_matrices.subgroups),
    resources = get_tar_resources(cores_req = 6, RAM_GB_req = 32)
  ),
  targets::tar_target(
    name = UMAP_embeddings_tibble.WNN.subgroups,
    description = "Compute subgroup WNN UMAP coordinates from the native weighted neighbour index [part_of_graph:full_subgroups]",
    command = run_WNN_UMAP(
      WNN_results = WNN_results.subgroups,
      n_neighbors = aggregation_UMAP_nNNs,
      min_dist = aggregation_UMAP_min_dist,
      col_prefix = "WNN_sub_UMAP"
    ),
    pattern = map(WNN_results.subgroups),
    resources = get_tar_resources(cores_req = 6, RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = clusters_tibble.WNN.subgroups,
    description = "Cluster subgroup cells from the native WNN SNN graph [part_of_graph:full_subgroups]",
    command = cluster_WNN_graph(
      WNN_results = WNN_results.subgroups,
      resolution = aggregation_subgroups_res / 3,
      cluster_col = "WNN_harmony_SNN_cluster_sub",
      min_barcodes = aggregation_cluster_min_barcodes
    ),
    pattern = map(WNN_results.subgroups),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = metadata_w_clusters_tibble.subgroups,
    description = "Join subgroup-native GEX, ATAC, and WNN embeddings, clusters, and modality weights onto metadata [part_of_graph:full_subgroups]",
    command = build_WNN_metadata_tibble(
      metadata_tibble = metadata_tibble.subgroups,
      WNN_results = WNN_results.subgroups,
      UMAP_embeddings_tibble = UMAP_embeddings_tibble.WNN.subgroups,
      WNN_clusters_tibble = clusters_tibble.WNN.subgroups,
      GEX_UMAP_embeddings_tibble = UMAP_embeddings_tibble.GEX.subgroups
    ) |>
      dplyr::filter(
        .data$barcode_w_prefix %in% names(PCA_clusters.GEX.subgroups),
        .data$barcode_w_prefix %in% names(LSI_clusters.ATAC.subgroups)
      ) |>
      dplyr::left_join(UMAP_embeddings_tibble.ATAC.subgroups, by = "barcode_w_prefix") |>
      dplyr::mutate(
        PCA_harmony_SNN_cluster_sub = PCA_clusters.GEX.subgroups[.data$barcode_w_prefix],
        LSI_harmony_SNN_cluster_sub = LSI_clusters.ATAC.subgroups[.data$barcode_w_prefix]
      ),
    pattern = map(
      metadata_tibble.subgroups,
      WNN_results.subgroups,
      UMAP_embeddings_tibble.WNN.subgroups,
      clusters_tibble.WNN.subgroups,
      UMAP_embeddings_tibble.GEX.subgroups,
      UMAP_embeddings_tibble.ATAC.subgroups,
      PCA_clusters.GEX.subgroups,
      LSI_clusters.ATAC.subgroups
    ),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = metadata_w_clusters_tibble_filtered.subgroups,
    description = "Score module agreement and optionally remove contaminated subgroup clusters [part_of_graph:full_subgroups]",
    command = metadata_w_clusters_tibble.subgroups |>
      add_module_score_agreement_to_metadata(
        marker_genes_list = UCell_GEX_marker_genes_list,
        parent_cluster_col = "PCA_harmony_SNN_cluster_cell_type"
      ) |>
      remove_contaminated_subclusters_from_metadata(
        sub_cluster_col = "PCA_harmony_SNN_cluster_sub",
        parent_cluster_col = "PCA_harmony_SNN_cluster_cell_type",
        prop_threshold = aggregation_prop_threshold_contaminated_cluster_removal
      ),
    pattern = map(metadata_w_clusters_tibble.subgroups),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = marker_tibbles.subgroups,
    description = "Find GEX marker genes for each subgroup GEX sub-cluster using BPCells [part_of_graph:full_subgroups]",
    command = get_BPCells_markers_from_matrix(
      feature_matrix = aggregated_counts_BPCells_matrix.GEX,
      metadata_tibble = metadata_w_clusters_tibble_filtered.subgroups,
      group_col = "PCA_harmony_SNN_cluster_sub"
    ),
    pattern = map(metadata_w_clusters_tibble_filtered.subgroups),
    resources = get_tar_resources(RAM_GB_req = 32)
  ),
  tarchetypes::tar_map(
    values = tibble::tribble(
      ~UMAP_type , ~UMAP_cols                              , ~UMAP_categorical_vars        ,
      "GEX"      , c("GEX_sub_UMAP_1", "GEX_sub_UMAP_2")   , "PCA_harmony_SNN_cluster_sub" ,
      "ATAC"     , c("ATAC_sub_UMAP_1", "ATAC_sub_UMAP_2") , "LSI_harmony_SNN_cluster_sub" ,
      "WNN"      , c("WNN_sub_UMAP_1", "WNN_sub_UMAP_2")   , "WNN_harmony_SNN_cluster_sub"
    ),
    names = UMAP_type,
    delimiter = ".",
    tarchetypes::tar_file(
      name = categorical.UMAPs.subgroups,
      description = "Plot categorical metadata UMAPs for each subgroup and save to file",
      command = {
        categorical_vars <- unique(c(
          UMAP_categorical_vars,
          switch(
            UMAP_type,
            GEX = aggregation_GEX_categorical_vars,
            ATAC = aggregation_ATAC_categorical_vars,
            WNN = aggregation_WNN_categorical_vars
          )
        ))

        plots <- categorical_vars |>
          purrr::set_names() |>
          purrr::map(\(variable) {
            plot_UMAP_from_metadata(
              metadata_tibble = metadata_w_clusters_tibble_filtered.subgroups,
              variable = variable,
              umap_cols = UMAP_cols
            )
          }) |>
          purrr::discard(\(plot) inherits(plot, "empty_plot_list"))

        plots |>
          save_plots_structured(dyn_suffix_in_subdir = TRUE, override_suffix = subgroups_to_process_vec)
      },
      pattern = map(metadata_w_clusters_tibble_filtered.subgroups, subgroups_to_process_vec),
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    tarchetypes::tar_file(
      name = continuous.UMAPs.subgroups,
      description = "Plot continuous feature UMAPs for each subgroup and save to file",
      command = {
        metadata_plots <- aggregation_w_WNN_continuous_vars |>
          purrr::set_names() |>
          purrr::map(\(variable) {
            plot_UMAP_from_metadata(
              metadata_tibble = metadata_w_clusters_tibble_filtered.subgroups,
              variable = variable,
              umap_cols = UMAP_cols
            )
          })
        feature_plots <- intersect(GEX_marker_genes_vec, rownames(aggregated_counts_BPCells_matrix.GEX)) |>
          purrr::set_names() |>
          purrr::map(\(variable) {
            plot_UMAP_from_metadata(
              metadata_tibble = metadata_w_clusters_tibble_filtered.subgroups,
              variable = variable,
              value_source = "feature",
              feature_matrix = aggregated_counts_BPCells_matrix.GEX,
              umap_cols = UMAP_cols
            )
          })

        c(metadata_plots, feature_plots) |>
          purrr::discard(\(plot) inherits(plot, "empty_plot_list"))
      } %>%
        save_plots_structured(
          dyn_suffix_in_subdir = TRUE,
          override_suffix = subgroups_to_process_vec
        ),
      pattern = map(metadata_w_clusters_tibble_filtered.subgroups, subgroups_to_process_vec),
      resources = get_tar_resources(RAM_GB_req = 16)
    ),
    skip_w_dummy_file_if(
      nrow(marker_tibbles.subgroups) == 0,
      tarchetypes::tar_file(
        name = marker_volcano_plots.subgroups,
        description = "Plot volcano plots of subgroup sub-cluster marker genes and save to file",
        command = marker_tibbles.subgroups |>
          plot_markers_volcano_simple() %>%
          save_plots_structured(dyn_suffix_in_subdir = TRUE, override_suffix = subgroups_to_process_vec),
        pattern = map(marker_tibbles.subgroups, subgroups_to_process_vec)
      )
    )
  )
)
