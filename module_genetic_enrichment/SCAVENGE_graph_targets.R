rlang::list2(
  targets::tar_target(
    name = graph_matrix.SCAVENGE,
    description = "Build the SNN graph matrix used for SCAVENGE TRS propagation",
    command = {
      if (identical(map_SCAVENGE_graph_input_type, "WNN")) {
        get_SNN_matrix_from_WNN_results(map_SCAVENGE_graph_input)
      } else {
        get_SNN_matrix_from_embedding_matrix(
          embedding_matrix = map_SCAVENGE_graph_input,
          dims = map_SCAVENGE_embedding_dims,
          k = aggregation_data_nNNs,
          dim_prefix = map_SCAVENGE_dim_prefix,
          threads = 15
        )
      }
    },
    resources = get_tar_resources(cores_req = 15, RAM_GB_req = 60)
  ),
  targets::tar_target(
    name = TRS_tibbles.SCAVENGE,
    description = "Propagate GWAS_chromVAR Z-scores through the NN graph to compute SCAVENGE TRS for each GWAS",
    command = get_SCAVENGE_TRS_from_ZScore_record(
      ZScore_record = GWAS_chromVAR_ZScore_records.single_cell,
      NN_graph = graph_matrix.SCAVENGE,
      cores = 15,
      permutation_times = genetic_enrichment_SCAVENGE_permutation_times,
      restart_prob = genetic_enrichment_SCAVENGE_restart_prob,
      seed_percent = genetic_enrichment_SCAVENGE_seed_percent
    ),
    pattern = map(GWAS_chromVAR_ZScore_records.single_cell),
    resources = get_tar_resources(RAM_GB_req = 60, cores_req = 15)
  ),
  targets::tar_target(
    name = TRS_summary_tibbles.SCAVENGE,
    description = "Summarize one GWAS SCAVENGE TRS branch across named and cell-type graph clusters",
    command = summarize_SCAVENGE_TRS_by_groups(
      TRS_tibble = TRS_tibbles.SCAVENGE,
      metadata_tibble = metadata_w_cell_types_tibble.WNN,
      graph_name = map_SCAVENGE_graph_name
    ),
    pattern = map(TRS_tibbles.SCAVENGE),
    resources = get_tar_resources(RAM_GB_req = 32)
  ),
  targets::tar_target(
    name = TRS_tibble.SCAVENGE,
    description = "Combine the tibble of SCAVENGE TRS values across GWAS within a specific graph",
    command = bind_rows(TRS_tibbles.SCAVENGE)
  ),
  tarchetypes::tar_file(
    name = SCAVENGE_UMAPs.SCAVENGE,
    description = "Plot one GWAS SCAVENGE TRS branch on the graph UMAP and save to file. [checkpoint:genetic_enrichment]",
    command = {
      plot <- get_SCAVENGE_TRS_UMAP_plots(
        TRS_tibble = TRS_tibbles.SCAVENGE,
        metadata_tibble = metadata_w_cell_types_tibble.WNN,
        umap_cols = map_SCAVENGE_umap_cols
      )
      save_plots_structured(
        plot,
        override_suffix = if (nrow(TRS_tibbles.SCAVENGE) == 0) NULL else unique(TRS_tibbles.SCAVENGE$GWAS_ID)[[1]]
      )
    },
    pattern = map(TRS_tibbles.SCAVENGE),
    resources = get_tar_resources(RAM_GB_req = 32)
  )
)
