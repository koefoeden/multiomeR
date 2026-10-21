rlang::list2(
  targets::tar_target(
    name = TRS_summary_tibble,
    description = "Summarize SCAVENGE TRS per GWAS across named and cell-type graph clusters [part_of_graph:genetic_enrichment_single_nucleus] [resource_observation: peak_RAM_GiB=7.7; n_nuclei=713k; aggregation=unpublished; method=slurm_monitor; date=2026-10-02]",
    command = summarize_SCAVENGE_TRS_by_groups(
      TRS_tibble = TRS_tibble,
      metadata_tibble = SCAVENGE_cell_metadata_tibble,
      graph_name = map_SCAVENGE_graph_name
    ),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = TRS_heatmap_data,
    description = "Prepare SCAVENGE TRS heatmap data from per-GWAS group summaries",
    command = add_GWAS_heatmap_categories(TRS_summary_tibble, GWAS_tibble = GWAS_inputs_tibble),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  tarchetypes::tar_file(
    name = median_by_cell_types_scaled.TRS_heatmap,
    description = "Save scaled median SCAVENGE TRS by WNN cell types. [checkpoint:genetic_enrichment] [resource_observation: peak_RAM_GiB=4.9; n_nuclei=713k; aggregation=unpublished; method=slurm_monitor; date=2026-10-10]",
    command = plot_WNN_TRS_heatmap(
      TRS_heatmap_data,
      GWAS_metadata_tracks_plot = GWAS_metadata_tracks_plot,
      compartments_patterns = genetic_enrichment_compartment_patterns,
      grouping = "cell_types", scaled = TRUE
    ) |>
      save_GWAS_heatmap(TRS_heatmap_data$GWAS_ID),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  tarchetypes::tar_file(
    name = median_by_cell_types_unscaled.TRS_heatmap,
    description = "Save unscaled median SCAVENGE TRS by WNN cell types. [checkpoint:genetic_enrichment]",
    command = plot_WNN_TRS_heatmap(
      TRS_heatmap_data,
      GWAS_metadata_tracks_plot = GWAS_metadata_tracks_plot,
      compartments_patterns = genetic_enrichment_compartment_patterns,
      grouping = "cell_types", scaled = FALSE
    ) |>
      save_GWAS_heatmap(TRS_heatmap_data$GWAS_ID),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  tarchetypes::tar_file(
    name = median_by_clusters_scaled.TRS_heatmap,
    description = "Save scaled median SCAVENGE TRS by WNN clusters. [checkpoint:genetic_enrichment]",
    command = plot_WNN_TRS_heatmap(
      TRS_heatmap_data,
      GWAS_metadata_tracks_plot = GWAS_metadata_tracks_plot,
      compartments_patterns = genetic_enrichment_compartment_patterns,
      grouping = "clusters", scaled = TRUE
    ) |>
      save_GWAS_heatmap(TRS_heatmap_data$GWAS_ID),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  tarchetypes::tar_file(
    name = median_by_clusters_unscaled.TRS_heatmap,
    description = "Save unscaled median SCAVENGE TRS by WNN clusters. [checkpoint:genetic_enrichment]",
    command = plot_WNN_TRS_heatmap(
      TRS_heatmap_data,
      GWAS_metadata_tracks_plot = GWAS_metadata_tracks_plot,
      compartments_patterns = genetic_enrichment_compartment_patterns,
      grouping = "clusters", scaled = FALSE
    ) |>
      save_GWAS_heatmap(TRS_heatmap_data$GWAS_ID),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  tarchetypes::tar_file(
    name = TRS_cluster_summary_plot,
    description = "Plot per-group SCAVENGE TRS medians, interquartile ranges, and full ranges. [checkpoint:genetic_enrichment] [resource_observation: peak_RAM_GiB=5.8; n_nuclei=213k; aggregation=mixed_human_31x; method=slurm_monitor; date=2026-10-10]",
    command = {
      plot <- plot_SCAVENGE_summary_score_intervals(TRS_summary_tibble)
      save_plots_structured(plot)
    },
    resources = get_tar_resources(RAM_GB_req = 16)
  )
)
