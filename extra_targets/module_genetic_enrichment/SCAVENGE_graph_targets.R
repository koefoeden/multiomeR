rlang::list2(
  targets::tar_target(
    name = graph_matrix,
    description = "Build the SNN graph matrix used for SCAVENGE TRS propagation [part_of_graph:genetic_enrichment_single_nucleus] [resource_observation: peak_RAM_GiB=3.8; n_nuclei=213k; aggregation=mixed_human_31x; method=slurm_monitor; date=2026-10-10] [resource_observation: peak_RAM_GiB=24.0; n_nuclei=713k; aggregation=unpublished; method=slurm_monitor; date=2026-10-10] [resource_observation: scaling=about 4.0 GiB more per 100k nuclei from 213k to 713k; date=2026-10-10]",
    command = get_SNN_matrix_from_WNN_results(WNN_results),
    resources = get_tar_resources(cores_req = 6, RAM_GB_req = 32)
  ),
  targets::tar_target(
    name = TRS_tibbles,
    description = "Compute cell-level SCAVENGE TRS for one GWAS [part_of_graph:genetic_enrichment_single_nucleus] [resource_observation: peak_RAM_GiB=4.1; n_nuclei=213k; aggregation=mixed_human_31x; method=slurm_monitor; date=2026-10-10]",
    command = get_SCAVENGE_TRS_tibble(
      chromVAR_z_score_record = chromVAR_z_score_records.single_nucleus,
      NN_graph = graph_matrix,
      restart_prob = genetic_enrichment_SCAVENGE_restart_prob,
      seed_percent = genetic_enrichment_SCAVENGE_seed_percent
    ),
    pattern = map(chromVAR_z_score_records.single_nucleus),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = TRS_tibble,
    description = "Combine the tibble of SCAVENGE TRS values across GWAS within a specific graph",
    command = bind_rows(TRS_tibbles)
  ),
  tarchetypes::tar_file(
    name = TRS_UMAPs,
    description = "Plot one GWAS SCAVENGE TRS branch on the graph UMAP and save to file. [checkpoint:genetic_enrichment] [resource_observation: peak_RAM_GiB=5.9; n_nuclei=213k; aggregation=mixed_human_31x; method=slurm_monitor; date=2026-10-10]",
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
    resources = get_tar_resources(RAM_GB_req = 16)
  )
)
