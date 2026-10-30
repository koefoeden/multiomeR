rlang::list2(
  targets::tar_target(
    name = graph_matrix.SCAVENGE.specificity,
    description = "Build the SNN graph matrix used for SCAVENGE TRS propagation [part_of_graph:genetic_enrichment_single_nucleus] [resource_observation: peak_RAM_GiB=3.8; n_nuclei=213k; aggregation=mixed_human_31x; method=slurm_monitor; date=2026-10-10] [resource_observation: peak_RAM_GiB=24.0; n_nuclei=713k; aggregation=unpublished; method=slurm_monitor; date=2026-10-10] [resource_observation: scaling=about 4.0 GiB more per 100k nuclei from 213k to 713k; date=2026-10-10]",
    command = get_SNN_matrix_from_WNN_results(WNN_results),
    resources = get_tar_resources(cores_req = 6, RAM_GB_req = 32)
  ),
  targets::tar_target(
    name = cell_metadata_tibble.SCAVENGE.specificity,
    description = "Keep each nucleus's graph embedding and cluster labels, so per-GWAS SCAVENGE targets load a small table",
    command = metadata_w_cell_types_tibble.WNN |>
      dplyr::select("barcode_w_prefix", dplyr::all_of(paste0("WNN_harmony_SNN", c("_cluster_named", "_cluster_cell_type")))) |>
      dplyr::left_join(dplyr::select(UMAP_embeddings_tibble.WNN, "barcode_w_prefix", dplyr::all_of(unlist(c("WNN_UMAP_1", "WNN_UMAP_2")))),
        by = "barcode_w_prefix")
  ),
  targets::tar_target(
    name = TRS_tibbles.SCAVENGE.specificity,
    description = "Compute cell-level SCAVENGE TRS for one GWAS [part_of_graph:genetic_enrichment_single_nucleus] [resource_observation: peak_RAM_GiB=4.1; n_nuclei=213k; aggregation=mixed_human_31x; method=slurm_monitor; date=2026-10-10]",
    command = get_SCAVENGE_TRS_tibble(
      chromVAR_z_score_record = z_score_records.chromVAR.specificity,
      NN_graph = graph_matrix.SCAVENGE.specificity,
      restart_prob = genetic_enrichment_SCAVENGE_restart_prob,
      seed_percent = genetic_enrichment_SCAVENGE_seed_percent
    ),
    pattern = map(z_score_records.chromVAR.specificity),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = TRS_tibble.SCAVENGE.specificity,
    description = "Combine cell-level SCAVENGE TRS across GWAS",
    command = bind_rows(TRS_tibbles.SCAVENGE.specificity)
  ),
  tarchetypes::tar_file(
    name = TRS_UMAPs.SCAVENGE.specificity,
    description = "Plot one GWAS SCAVENGE TRS branch on the graph UMAP and save to file. [checkpoint:genetic_enrichment] [resource_observation: peak_RAM_GiB=5.9; n_nuclei=213k; aggregation=mixed_human_31x; method=slurm_monitor; date=2026-10-10]",
    command = {
      plot <- get_SCAVENGE_TRS_UMAP_plots(
        TRS_tibble = TRS_tibbles.SCAVENGE.specificity,
        metadata_tibble = cell_metadata_tibble.SCAVENGE.specificity,
        umap_cols = c("WNN_UMAP_1", "WNN_UMAP_2"),
        label_col = paste0("WNN_harmony_SNN", "_cluster_named")
      )
      save_plots_structured(
        plot,
        override_suffix = if (nrow(TRS_tibbles.SCAVENGE.specificity) == 0) NULL else unique(TRS_tibbles.SCAVENGE.specificity$GWAS_ID)[[1]],
        dyn_suffix_in_subdir = TRUE
      )
    },
    pattern = map(TRS_tibbles.SCAVENGE.specificity),
    resources = get_tar_resources(RAM_GB_req = 16)
  )
)
