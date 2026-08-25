rlang::list2(
  targets::tar_target(
    name = graph_matrix,
    description = "Build the SNN graph matrix used for SCAVENGE TRS propagation [part_of_graph:genetic_enrichment_single_nucleus]",
    command = {
      if (identical(map_SCAVENGE_graph_input_type, "WNN")) {
        get_SNN_matrix_from_WNN_results(map_SCAVENGE_graph_input)
      } else {
        get_SNN_matrix_from_embedding_matrix(
          embedding_matrix = map_SCAVENGE_graph_input,
          dims = map_SCAVENGE_embedding_dims,
          k = aggregation_data_nNNs,
          dim_prefix = map_SCAVENGE_dim_prefix,
          threads = 6
        )
      }
    },
    resources = get_tar_resources(cores_req = 6, RAM_GB_req = 60)
  ),
  targets::tar_target(
    name = SCAVENGE_result_records,
    description = "Compute cell-level SCAVENGE TRS and cluster-level permutation significance for one GWAS [part_of_graph:genetic_enrichment_single_nucleus]",
    command = get_SCAVENGE_result_from_chromVAR_z_score_record(
      chromVAR_z_score_record = chromVAR_z_score_records.single_nucleus,
      NN_graph = graph_matrix,
      metadata_tibble = metadata_w_cell_types_tibble.WNN,
      graph_name = map_SCAVENGE_graph_name,
      cores = 6,
      permutation_times = genetic_enrichment_SCAVENGE_permutation_times,
      restart_prob = genetic_enrichment_SCAVENGE_restart_prob,
      seed_percent = genetic_enrichment_SCAVENGE_seed_percent,
      native_source_file = SCAVENGE_native_source_file
    ),
    pattern = map(chromVAR_z_score_records.single_nucleus),
    resources = get_tar_resources(RAM_GB_req = 60, cores_req = 6)
  ),
  targets::tar_target(
    name = TRS_tibbles,
    description = "Extract one cell-level SCAVENGE TRS tibble",
    command = purrr::pluck(SCAVENGE_result_records, "TRS_tibble"),
    pattern = map(SCAVENGE_result_records)
  ),
  targets::tar_target(
    name = TRS_summary_tibbles,
    description = "Extract one cluster-level SCAVENGE summary with empirical permutation P-values",
    command = purrr::pluck(SCAVENGE_result_records, "TRS_summary_tibble"),
    pattern = map(SCAVENGE_result_records)
  ),
  targets::tar_target(
    name = TRS_tibble,
    description = "Combine the tibble of SCAVENGE TRS values across GWAS within a specific graph",
    command = bind_rows(TRS_tibbles)
  ),
  tarchetypes::tar_file(
    name = TRS_UMAPs,
    description = "Plot one GWAS SCAVENGE TRS branch on the graph UMAP and save to file. [checkpoint:genetic_enrichment]",
    command = {
      plot <- get_SCAVENGE_TRS_UMAP_plots(
        TRS_tibble = TRS_tibbles,
        metadata_tibble = metadata_w_cell_types_tibble.WNN,
        umap_cols = map_SCAVENGE_umap_cols,
        label_col = paste0(map_SCAVENGE_graph_name, "_cluster_named")
      )
      save_plots_structured(
        plot,
        override_suffix = if (nrow(TRS_tibbles) == 0) NULL else unique(TRS_tibbles$GWAS_ID)[[1]],
        dyn_suffix_in_subdir = TRUE
      )
    },
    pattern = map(TRS_tibbles),
    resources = get_tar_resources(RAM_GB_req = 32)
  )
)
