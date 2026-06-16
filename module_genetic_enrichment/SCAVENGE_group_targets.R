rlang::list2(
  targets::tar_target(
    name = TRS_summary_tibble.SCAVENGE,
    description = "Combine per-GWAS SCAVENGE TRS group summaries across named and cell-type graph clusters",
    command = dplyr::bind_rows(TRS_summary_tibbles.SCAVENGE),
    resources = get_tar_resources(RAM_GB_req = 32)
  ),
  targets::tar_target(
    name = TRS_dotplot_data.SCAVENGE,
    description = "Prepare SCAVENGE TRS dotplot data from per-GWAS group summaries",
    command = TRS_summary_tibble.SCAVENGE |>
      get_SCAVENGE_dotplot_data(GWAS_tibble = GWAS_inputs_tibble),
    resources = get_tar_resources(RAM_GB_req = 32)
  ),
  targets::tar_target(
    name = TRS_dotplot_scaled_data.SCAVENGE,
    description = "Rescale SCAVENGE TRS dotplot data per GWAS and grouping column to the 0-1 range",
    command = scale_GWAS_dotplot_scores(TRS_dotplot_data.SCAVENGE, group_cols = "grouping_col"),
    resources = get_tar_resources(RAM_GB_req = 32)
  ),
  tarchetypes::tar_file(
    name = TRS_dotplot.SCAVENGE,
    description = "Save SCAVENGE TRS dotplots with Open Targets GWAS metadata annotations. [checkpoint:genetic_enrichment]",
    command = {
      plots <- plot_grouped_GWAS_by_cluster_dotplots(
        TRS_dotplot_scaled_data.SCAVENGE,
        split_col = "grouping_col",
        name_suffix = "scaled",
        GWAS_metadata_tracks_plot = GWAS_metadata_tracks_plot,
        compartments_patterns = genetic_enrichment_compartment_patterns,
        scaled = TRUE
      )
      save_plots_structured(
        plots,
        filetype = "png",
        width = 20,
        height = max(7, 0.24 * dplyr::n_distinct(TRS_dotplot_scaled_data.SCAVENGE$GWAS_ID) + 3)
      )
    },
    resources = get_tar_resources(RAM_GB_req = 32)
  ),
  tarchetypes::tar_file(
    name = sig_prop_bars.SCAVENGE,
    description = "Plot proportion of cells with significant SCAVENGE TRS per group and GWAS and save to file. [checkpoint:genetic_enrichment]",
    command = {
      plot <- plot_SCAVENGE_summary_sig_proportion(TRS_summary_tibble.SCAVENGE)
      save_plots_structured(plot)
    },
    resources = get_tar_resources(RAM_GB_req = 32)
  ),
  tarchetypes::tar_file(
    name = TRS_by_cluster_boxplot.SCAVENGE,
    description = "Plot per-group SCAVENGE TRS summary intervals and save to file. [checkpoint:genetic_enrichment]",
    command = {
      plot <- plot_SCAVENGE_summary_score_intervals(TRS_summary_tibble.SCAVENGE)
      save_plots_structured(plot)
    },
    resources = get_tar_resources(RAM_GB_req = 32)
  )
)
