rlang::list2(
  targets::tar_target(
    name = TRS_summary_tibble,
    description = "Summarize SCAVENGE TRS per GWAS across named and cell-type graph clusters [part_of_graph:genetic_enrichment_single_nucleus]",
    command = summarize_SCAVENGE_TRS_by_groups(
      TRS_tibble = TRS_tibble,
      metadata_tibble = metadata_w_cell_types_tibble.WNN,
      graph_name = map_SCAVENGE_graph_name
    ),
    resources = get_tar_resources(RAM_GB_req = 32)
  ),
  targets::tar_target(
    name = TRS_heatmap_data,
    description = "Prepare SCAVENGE TRS heatmap data from per-GWAS group summaries",
    command = TRS_summary_tibble |>
      add_GWAS_heatmap_categories(GWAS_tibble = GWAS_inputs_tibble),
    resources = get_tar_resources(RAM_GB_req = 32)
  ),
  tarchetypes::tar_file(
    name = median_by_cell_types_scaled.TRS_heatmap,
    description = "Save scaled median SCAVENGE TRS by WNN cell types. [checkpoint:genetic_enrichment]",
    command = plot_WNN_TRS_heatmap(
      TRS_heatmap_data,
      GWAS_metadata_tracks_plot = GWAS_metadata_tracks_plot,
      compartments_patterns = genetic_enrichment_compartment_patterns,
      grouping = "cell_types", scaled = TRUE
    ) |>
      save_plots_structured(filetype = "png", width = 17,
        height = max(5.5, 0.3 * dplyr::n_distinct(TRS_heatmap_data$GWAS_ID) + 3.5)),
    resources = get_tar_resources(RAM_GB_req = 32)
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
      save_plots_structured(filetype = "png", width = 17,
        height = max(5.5, 0.3 * dplyr::n_distinct(TRS_heatmap_data$GWAS_ID) + 3.5)),
    resources = get_tar_resources(RAM_GB_req = 32)
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
      save_plots_structured(filetype = "png", width = 17,
        height = max(5.5, 0.3 * dplyr::n_distinct(TRS_heatmap_data$GWAS_ID) + 3.5)),
    resources = get_tar_resources(RAM_GB_req = 32)
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
      save_plots_structured(filetype = "png", width = 17,
        height = max(5.5, 0.3 * dplyr::n_distinct(TRS_heatmap_data$GWAS_ID) + 3.5)),
    resources = get_tar_resources(RAM_GB_req = 32)
  ),
  tarchetypes::tar_file(
    name = TRS_cluster_summary_plot,
    description = "Plot per-group SCAVENGE TRS medians, interquartile ranges, and full ranges. [checkpoint:genetic_enrichment]",
    command = {
      plot <- plot_SCAVENGE_summary_score_intervals(TRS_summary_tibble)
      save_plots_structured(plot)
    },
    resources = get_tar_resources(RAM_GB_req = 32)
  )
)
