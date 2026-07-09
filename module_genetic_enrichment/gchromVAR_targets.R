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
    name = GWAS_input_records,
    description = "Convert one configured GWAS input to one record for dynamic trait-level chromVAR branching",
    command = get_GWAS_chromVAR_input_record(
      GWAS_input_tibble = GWAS_analysis_inputs_tibble,
      open_targets_credible_set_dataset_path = open_targets_credible_set_dataset_path
    ),
    pattern = map(GWAS_analysis_inputs_tibble),
    iteration = "list",
    resources = get_tar_resources(RAM_GB_req = 40)
  ),
  targets::tar_target(
    name = GWAS_peak_weight_records,
    description = "Build capped peak posterior-probability weights for one GWAS [part_of_graph:genetic_enrichment_single_nucleus]",
    command = get_GWAS_chromVAR_peak_weight_record(
      GWAS_input_record = GWAS_input_records,
      peak_ranges = genetic_enrichment_peak_ranges,
      posterior_probability_cutoff = genetic_enrichment_posterior_probability_cutoff,
      posterior_probability_weighting_function = posterior_probability_weighting_function
    ),
    pattern = map(GWAS_input_records),
    resources = get_tar_resources(RAM_GB_req = 40)
  ),
  tarchetypes::tar_file(
    name = GWAS_peak_weights_barplot,
    description = "Plot summary bar plot of capped summed SNP PPs within peaks and save to file. [checkpoint:genetic_enrichment]",
    command = {
      plot <- plot_GWAS_chromVAR_peak_weights_summary(GWAS_peak_weight_records)
      save_plots_structured(plot)
    }
  ),
  targets::tar_target(
    name = chromVAR_z_score_chunk_records.single_nucleus,
    description = "Compute single-nucleus chromVAR z-scores for one GWAS and one reusable ATAC chunk",
    command = get_GWAS_chromVAR_z_score_chunk_record(
      peak_weight_record = GWAS_peak_weight_records,
      RSE_ATAC = chromVAR_obj.ATAC,
      chunk_context_record = chromVAR_chunk_context_records.ATAC
    ),
    pattern = cross(GWAS_peak_weight_records, chromVAR_chunk_context_records.ATAC),
    iteration = "vector",
    resources = get_tar_resources(RAM_GB_req = 60)
  ),
  targets::tar_target(
    name = chromVAR_z_score_chunk_records_by_trait.single_nucleus,
    description = "Group chunk-level single-nucleus chromVAR z-score records by GWAS for branch-stable recombination",
    command = chromVAR_z_score_chunk_records.single_nucleus |>
      dplyr::group_by(GWAS_ID) |>
      targets::tar_group(),
    iteration = "group",
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = chromVAR_z_score_records.single_nucleus,
    description = "Combine chunk-level single-nucleus chromVAR z-scores for SCAVENGE [part_of_graph:genetic_enrichment_single_nucleus]",
    command = combine_GWAS_chromVAR_z_score_chunk_records(chromVAR_z_score_chunk_records_by_trait.single_nucleus),
    pattern = map(chromVAR_z_score_chunk_records_by_trait.single_nucleus)
  )
)
