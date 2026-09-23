rlang::list2(
  targets::tar_target(
    name = GWAS_absolute_effect_weighting_status_tibble,
    description = "Infer variant- or locus-level PIP x absolute-beta support for every GWAS",
    command = get_GWAS_absolute_effect_weighting_status_tibble(
      GWAS_input_records = GWAS_input_records,
      posterior_probability_cutoff = genetic_enrichment_posterior_probability_cutoff
    )
  ),
  targets::tar_target(
    name = GWAS_absolute_effect_peak_variant_weight_records,
    description = "Allocate PIP-mass-calibrated and peak-capped PIP x absolute-beta weights to the variants of eligible GWAS inputs",
    command = get_GWAS_absolute_effect_peak_variant_weight_records(
      GWAS_input_records = GWAS_input_records,
      weighting_status_tibble = GWAS_absolute_effect_weighting_status_tibble,
      peak_ranges = genetic_enrichment_peak_ranges,
      posterior_probability_cutoff = genetic_enrichment_posterior_probability_cutoff
    ),
    resources = get_tar_resources(RAM_GB_req = 40)
  ),
  targets::tar_target(
    name = GWAS_absolute_effect_peak_weight_records,
    description = "Build PIP-mass-calibrated and peak-capped PIP x absolute-beta weights for eligible GWAS inputs [part_of_graph:genetic_enrichment_cell_type_absolute_effect]",
    command = purrr::map(
      GWAS_absolute_effect_peak_variant_weight_records,
      get_GWAS_chromVAR_peak_weight_record,
      peak_ranges = genetic_enrichment_peak_ranges
    ),
    resources = get_tar_resources(RAM_GB_req = 40)
  ),
  targets::tar_target(
    name = GWAS_absolute_effect_inputs_tibble,
    description = "Subset GWAS metadata to automatically eligible absolute-effect-weighted inputs",
    command = GWAS_inputs_tibble |>
      dplyr::inner_join(
        GWAS_absolute_effect_weighting_status_tibble |>
          dplyr::filter(eligible) |>
          dplyr::select(GWAS_ID, effect_weighting_route),
        by = "GWAS_ID",
        relationship = "one-to-one"
      ) |>
      dplyr::arrange(config_order)
  ),
  targets::tar_target(
    name = GWAS_absolute_effect_peak_weight_matrix,
    description = "Combine eligible PIP x absolute-beta peak weights into one annotation matrix [part_of_graph:genetic_enrichment_cell_type_absolute_effect]",
    command = get_GWAS_chromVAR_peak_weight_matrix(
      peak_weight_records = GWAS_absolute_effect_peak_weight_records
    ),
    resources = get_tar_resources(RAM_GB_req = 40)
  ),
  targets::tar_target(
    name = GWAS_absolute_effect_metadata_tracks_plot,
    description = "Build metadata tracks for eligible absolute-effect-weighted GWAS inputs",
    command = GWAS_absolute_effect_inputs_tibble |>
      dplyr::arrange(Category, GWAS_ID) |>
      dplyr::mutate(GWAS_ID = factor(GWAS_ID, levels = rev(GWAS_ID))) |>
      plot_GWAS_metadata_tracks()
  ),
  targets::tar_target(
    name = chromVAR_absolute_effect_deviation_tibble.cell_type_pseudobulk,
    description = "Summarize absolute-effect peak contributions into the cell-type heatmap table [part_of_graph:genetic_enrichment_cell_type_absolute_effect]",
    command = summarize_GWAS_chromVAR_peak_contributions(
      peak_contribution_tibble = chromVAR_absolute_effect_peak_contribution_tibble.cell_type_pseudobulk,
      cell_type_support_tibble = cell_type_pseudobulk_support_tibble.ATAC
    ) |>
      dplyr::left_join(
        dplyr::select(GWAS_absolute_effect_inputs_tibble, GWAS_ID, effect_weighting_route),
        by = "GWAS_ID", relationship = "many-to-one") |>
      dplyr::relocate(effect_weighting_route, .after = Category),
    resources = get_tar_resources(RAM_GB_req = 8)
  ),
  tarchetypes::tar_file(
    name = beta_weighted_unscaled.chromVAR_deviation_heatmaps.cell_type_pseudobulk,
    description = "Save absolute-effect-weighted GWAS deviations without within-GWAS scaling. [checkpoint:genetic_enrichment]",
    command = plot_GWAS_chromVAR_deviation_heatmap(
      chromVAR_absolute_effect_deviation_tibble.cell_type_pseudobulk,
      GWAS_metadata_tracks_plot = GWAS_absolute_effect_metadata_tracks_plot,
      compartments_patterns = genetic_enrichment_compartment_patterns,
      standardize = FALSE,
      beta_weighted = TRUE
    ) |>
      save_GWAS_heatmap(chromVAR_absolute_effect_deviation_tibble.cell_type_pseudobulk$GWAS_ID)
  ),
  tarchetypes::tar_file(
    name = beta_weighted_scaled.chromVAR_deviation_heatmaps.cell_type_pseudobulk,
    description = "Save absolute-effect-weighted GWAS deviations standardized within each GWAS. [checkpoint:genetic_enrichment]",
    command = plot_GWAS_chromVAR_deviation_heatmap(
      chromVAR_absolute_effect_deviation_tibble.cell_type_pseudobulk,
      GWAS_metadata_tracks_plot = GWAS_absolute_effect_metadata_tracks_plot,
      compartments_patterns = genetic_enrichment_compartment_patterns,
      standardize = TRUE,
      beta_weighted = TRUE
    ) |>
      save_GWAS_heatmap(chromVAR_absolute_effect_deviation_tibble.cell_type_pseudobulk$GWAS_ID)
  ),
  targets::tar_target(
    name = GWAS_absolute_effect_peak_variant_weight_tibble,
    description = "Combine the absolute-effect peak-to-variant weight allocations of eligible GWAS inputs",
    command = dplyr::bind_rows(GWAS_absolute_effect_peak_variant_weight_records)
  ),
  targets::tar_target(
    name = chromVAR_absolute_effect_peak_contribution_tibble.cell_type_pseudobulk,
    description = "Decompose the existing absolute-effect annotation matrix into exact peak contributions",
    command = if (nrow(GWAS_absolute_effect_inputs_tibble) == 0L) {
      chromVAR_peak_contribution_tibble.cell_type_pseudobulk[0, ]
    } else get_GWAS_chromVAR_peak_contribution_tibble(
      chromVAR_background_record = chromVAR_background_record.cell_type_pseudobulk,
      pseudobulk_ATAC_data_matrix = cell_type_pseudobulk_counts_matrix.ATAC,
      chromVAR_obj = chromVAR_obj.ATAC,
      annotation_matrix = GWAS_absolute_effect_peak_weight_matrix,
      GWAS_inputs_tibble = GWAS_absolute_effect_inputs_tibble),
    resources = get_tar_resources(RAM_GB_req = 60)
  ),
  targets::tar_target(
    name = chromVAR_absolute_effect_variant_contribution_tibble.cell_type_pseudobulk,
    description = "Allocate exact absolute-effect peak contributions to credible-set variants",
    command = if (nrow(chromVAR_absolute_effect_peak_contribution_tibble.cell_type_pseudobulk) == 0L) {
      chromVAR_variant_contribution_tibble.cell_type_pseudobulk[0, ]
    } else get_GWAS_chromVAR_variant_contribution_tibble(
      chromVAR_absolute_effect_peak_contribution_tibble.cell_type_pseudobulk,
      GWAS_absolute_effect_peak_variant_weight_tibble),
    resources = get_tar_resources(RAM_GB_req = 40)
  ),
  targets::tar_target(
    name = chromVAR_absolute_effect_locus_contribution_tibble.cell_type_pseudobulk,
    description = "Sum absolute-effect variant contributions into loci with L2G predictions",
    command = if (nrow(chromVAR_absolute_effect_variant_contribution_tibble.cell_type_pseudobulk) == 0L) {
      chromVAR_locus_contribution_tibble.cell_type_pseudobulk[0, ]
    } else get_GWAS_chromVAR_locus_contribution_tibble(
      chromVAR_absolute_effect_variant_contribution_tibble.cell_type_pseudobulk,
      GWAS_locus_to_gene_tibble),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  tarchetypes::tar_file(
    name = chromVAR_absolute_effect_locus_contribution_per_GWAS_faceted_bars_plots.cell_type_pseudobulk,
    description = "Save one faceted locus contribution figure per eligible GWAS using the absolute-effect heatmap weights. [checkpoint:genetic_enrichment]",
    command = chromVAR_absolute_effect_locus_contribution_tibble.cell_type_pseudobulk |>
      plot_GWAS_absolute_effect_locus_bars() |>
      save_plots_structured(width = 18, height = 2.5 + 2.3 * ceiling(dplyr::n_distinct(
        chromVAR_absolute_effect_locus_contribution_tibble.cell_type_pseudobulk$cluster) / 3)),
    resources = get_tar_resources(RAM_GB_req = 8)
  )

)
