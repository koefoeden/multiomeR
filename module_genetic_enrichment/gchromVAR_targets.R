rlang::list2(
  targets::tar_target(
    name = genetic_enrichment_peak_ranges,
    description = "Extract the ordered ATAC peak ranges used for GWAS_chromVAR weighting",
    command = SummarizedExperiment::rowRanges(chromVAR_obj.ATAC)
  ),
  targets::tar_target(
    name = GWAS_peak_variant_weight_records,
    description = "Allocate capped peak posterior-probability weights to the overlapping credible-set variants of one GWAS",
    command = get_GWAS_chromVAR_peak_variant_weight_tibble(
      GWAS_ID = GWAS_input_records$GWAS_ID,
      variant_GRanges = S4Vectors::subset(
        GWAS_input_records$credible_set_GRanges,
        posteriorProbability > genetic_enrichment_posterior_probability_cutoff
      ),
      peak_ranges = genetic_enrichment_peak_ranges
    ),
    pattern = map(GWAS_input_records),
    resources = get_tar_resources(RAM_GB_req = 40)
  ),
  targets::tar_target(
    name = GWAS_peak_weight_records,
    description = "Build capped peak posterior-probability weights for one GWAS [part_of_graph:genetic_enrichment_single_nucleus]",
    command = get_GWAS_chromVAR_peak_weight_record(
      peak_variant_weight_tibble = GWAS_peak_variant_weight_records,
      peak_ranges = genetic_enrichment_peak_ranges
    ),
    pattern = map(GWAS_peak_variant_weight_records),
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
