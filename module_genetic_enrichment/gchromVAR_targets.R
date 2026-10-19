rlang::list2(
  targets::tar_target(
    name = posterior_probability_weighting_function,
    description = "Load or define the posteriorProbability weighting function for GWAS_chromVAR",
    command = if (
      is.null(genetic_enrichment_posterior_probability_weighting_function_name) | isFALSE(genetic_enrichment_posterior_probability_weighting_function_name)
    ) {
      NULL
    } else {
      eval(rlang::parse_expr(genetic_enrichment_posterior_probability_weighting_function_name))
    }
  ),
  targets::tar_target(
    name = genetic_enrichment_peak_ranges,
    description = "Extract the ordered ATAC peak ranges used for GWAS_chromVAR weighting",
    command = SummarizedExperiment::rowRanges(chromVAR_obj.ATAC)
  ),
  targets::tar_target(
    name = GWAS_chromVAR_input_records,
    description = "Convert one configured GWAS input to one record for dynamic GWAS_chromVAR branching",
    command = get_GWAS_chromVAR_input_record(
      GWAS_input_tibble = GWAS_analysis_inputs_tibble,
      open_targets_credible_set_dataset_path = open_targets_credible_set_dataset_path
    ),
    pattern = map(GWAS_analysis_inputs_tibble),
    iteration = "list",
    resources = get_tar_resources(RAM_GB_req = 40)
  ),
  targets::tar_target(
    name = GWAS_chromVAR_peak_weight_records,
    description = "Build capped peak posterior-probability weights for one GWAS",
    command = get_GWAS_chromVAR_peak_weight_record(
      GWAS_input_record = GWAS_chromVAR_input_records,
      peak_ranges = genetic_enrichment_peak_ranges,
      posterior_probability_cutoff = genetic_enrichment_posterior_probability_cutoff,
      posterior_probability_weighting_function = posterior_probability_weighting_function
    ),
    pattern = map(GWAS_chromVAR_input_records),
    resources = get_tar_resources(RAM_GB_req = 40)
  ),
  tarchetypes::tar_file(
    name = GWAS_chromVAR_peak_weights_barplot,
    description = "Plot summary bar plot of capped summed SNP PPs within peaks and save to file. [checkpoint:genetic_enrichment]",
    command = {
      plot <- plot_GWAS_chromVAR_peak_weights_summary(GWAS_chromVAR_peak_weight_records)
      save_plots_structured(plot)
    }
  ),
  targets::tar_target(
    name = GWAS_chromVAR_ZScore_chunk_records,
    description = "Compute single-cell GWAS_chromVAR Z-scores for one GWAS and one reusable ATAC chunk",
    command = get_GWAS_chromVAR_ZScore_chunk_record(
      peak_weight_record = GWAS_chromVAR_peak_weight_records,
      RSE_ATAC = chromVAR_obj.ATAC,
      chunk_context_record = chromVAR_chunk_context_records.ATAC
    ),
    pattern = cross(GWAS_chromVAR_peak_weight_records, chromVAR_chunk_context_records.ATAC),
    iteration = "vector",
    resources = get_tar_resources(RAM_GB_req = 60)
  ),
  targets::tar_target(
    name = GWAS_chromVAR_ZScore_chunk_records_by_GWAS,
    description = "Group chunk-level GWAS_chromVAR Z-score records by GWAS for branch-stable recombination",
    command = GWAS_chromVAR_ZScore_chunk_records |>
      dplyr::group_by(GWAS_ID) |>
      targets::tar_group(),
    iteration = "group",
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = GWAS_chromVAR_ZScore_records.single_cell,
    description = "Combine chunk-level single-cell GWAS_chromVAR Z-scores for one GWAS",
    command = combine_GWAS_chromVAR_ZScore_chunk_records(GWAS_chromVAR_ZScore_chunk_records_by_GWAS),
    pattern = map(GWAS_chromVAR_ZScore_chunk_records_by_GWAS)
  ),
  targets::tar_target(
    name = GWAS_chromVAR_ZScores_tibble.single_cell,
    description = "Combine per-GWAS single-cell GWAS_chromVAR Z-score branches into one long tibble",
    command = combine_GWAS_chromVAR_ZScore_records_tibble(GWAS_chromVAR_ZScore_records.single_cell)
  ),
  targets::tar_target(
    name = GWAS_chromVAR_summarized_ZScores_tibble,
    description = "Aggregate single-cell GWAS_chromVAR Z-scores to cluster-donor pseudobulk summaries",
    command = summarize_GWAS_chromVAR_ZScores(
      ZScore_tibble = GWAS_chromVAR_ZScores_tibble.single_cell,
      metadata_tibble = metadata_w_cell_types_tibble.WNN
    )
  ),
  targets::tar_target(
    name = GWAS_chromVAR_summarized_Zscores_dotplot_data,
    description = "Prepare dotplot data from pseudobulk GWAS_chromVAR Z-scores grouped by GWAS and cluster",
    command = GWAS_chromVAR_summarized_ZScores_tibble |>
      get_GWAS_dotplot_data(GWAS_tibble = GWAS_inputs_tibble)
  ),
  targets::tar_target(
    name = GWAS_chromVAR_summarized_Zscores_dotplot_scaled_data,
    description = "Rescale pseudobulk GWAS_chromVAR dotplot data per GWAS to the 0-1 range",
    command = scale_GWAS_dotplot_scores(GWAS_chromVAR_summarized_Zscores_dotplot_data)
  ),
  tarchetypes::tar_file(
    name = GWAS_chromVAR_summarized_Zscores_dotplot,
    description = "Save pseudobulk GWAS_chromVAR Z-score dotplot with Open Targets GWAS metadata annotations. [checkpoint:genetic_enrichment]",
    command = {
      plot <- plot_GWAS_by_cluster_dotplot(
        GWAS_chromVAR_summarized_Zscores_dotplot_data,
        GWAS_metadata_tracks_plot = GWAS_metadata_tracks_plot,
        compartments_patterns = genetic_enrichment_compartment_patterns
      )
      save_plots_structured(
        plot,
        filetype = "png",
        width = 20,
        height = max(7, 0.24 * dplyr::n_distinct(GWAS_chromVAR_summarized_Zscores_dotplot_data$GWAS_ID) + 3)
      )
    }
  ),
  tarchetypes::tar_file(
    name = GWAS_chromVAR_summarized_Zscores_dotplot_scaled,
    description = "Save rescaled pseudobulk GWAS_chromVAR Z-score dotplot with Open Targets GWAS metadata annotations. [checkpoint:genetic_enrichment]",
    command = {
      plot <- plot_GWAS_by_cluster_dotplot(
        GWAS_chromVAR_summarized_Zscores_dotplot_scaled_data,
        GWAS_metadata_tracks_plot = GWAS_metadata_tracks_plot,
        compartments_patterns = genetic_enrichment_compartment_patterns,
        scaled = TRUE
      )
      save_plots_structured(
        plot,
        filetype = "png",
        width = 20,
        height = max(7, 0.24 * dplyr::n_distinct(GWAS_chromVAR_summarized_Zscores_dotplot_scaled_data$GWAS_ID) + 3)
      )
    }
  ),
  tarchetypes::tar_file(
    name = GWAS_chromVAR_summarized_Zscores_boxplot,
    description = "Plot pseudobulk GWAS_chromVAR Z-scores per cluster as boxplots and save to file. [checkpoint:genetic_enrichment]",
    command = {
      plot <- plot_GWAS_by_group(GWAS_chromVAR_summarized_ZScores_tibble)
      save_plots_structured(plot)
    }
  )
)
