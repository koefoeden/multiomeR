rlang::list2(
  targets::tar_target(
    name = GWAS_absolute_effect_weighting_status_tibble,
    description = "Infer variant- or locus-level PIP x absolute-beta support for every raw-PIP GWAS",
    command = get_GWAS_absolute_effect_weighting_status_tibble(
      GWAS_input_records = GWAS_input_records,
      posterior_probability_cutoff = genetic_enrichment_posterior_probability_cutoff
    )
  ),
  targets::tar_target(
    name = GWAS_absolute_effect_peak_weight_records,
    description = "Build PIP-mass-calibrated and peak-capped PIP x absolute-beta weights for eligible GWAS inputs [part_of_graph:genetic_enrichment_cell_type_absolute_effect]",
    command = get_GWAS_absolute_effect_peak_weight_records(
      GWAS_input_records = GWAS_input_records,
      weighting_status_tibble = GWAS_absolute_effect_weighting_status_tibble,
      peak_ranges = genetic_enrichment_peak_ranges,
      posterior_probability_cutoff = genetic_enrichment_posterior_probability_cutoff
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
      dplyr::mutate(variant_weighting_mode = effect_weighting_route) |>
      dplyr::arrange(config_order)
  ),
  targets::tar_target(
    name = GWAS_absolute_effect_peak_weight_matrix,
    description = "Combine eligible PIP x absolute-beta peak weights into one annotation matrix [part_of_graph:genetic_enrichment_cell_type_absolute_effect]",
    command = get_GWAS_chromVAR_peak_weight_matrix(
      peak_weight_records = GWAS_absolute_effect_peak_weight_records,
      RSE_ATAC = chromVAR_obj.ATAC
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
    description = "Compute the automatically eligible PIP x absolute-beta cell-type chromVAR heatmap table [part_of_graph:genetic_enrichment_cell_type_absolute_effect]",
    command = get_GWAS_absolute_effect_chromVAR_deviation_tibble(
      psbulk_ATAC_data_matrix = cell_type_pseudobulk_counts_matrix.ATAC,
      chromVAR_obj = chromVAR_obj.ATAC,
      annotation_matrix = GWAS_absolute_effect_peak_weight_matrix,
      background_record = chromVAR_background_record.cell_type_pseudobulk,
      GWAS_inputs_tibble = GWAS_absolute_effect_inputs_tibble,
      cell_type_support_tibble = cell_type_pseudobulk_support_tibble.ATAC
    ),
    resources = get_tar_resources(RAM_GB_req = 60)
  ),
  tarchetypes::tar_file(
    name = chromVAR_absolute_effect_deviation_heatmap.cell_type_pseudobulk,
    description = "Save the automatically eligible PIP x absolute-beta cell-type chromVAR relative-deviation heatmap. [checkpoint:genetic_enrichment]",
    command = {
      plot <- plot_GWAS_by_cluster_heatmap(
        chromVAR_absolute_effect_deviation_tibble.cell_type_pseudobulk,
        GWAS_metadata_tracks_plot = GWAS_absolute_effect_metadata_tracks_plot,
        compartments_patterns = genetic_enrichment_compartment_patterns,
        fill_col = "relative_deviation",
        fill_label = "PIP x |beta| relative deviation",
        support_label_col = "support_label"
      )
      if (!inherits(plot, "empty_plot_list")) plot <- plot + patchwork::plot_annotation(
        title = "Effect-magnitude-weighted GWAS accessibility by cell type",
        subtitle = stringr::str_wrap("Compare cell types within a trait; absolute effect weighting removes the direction of genetic effects.", width = 100),
        caption = stringr::str_wrap("Eligible GWAS variants are weighted by PIP x absolute beta. Colour is deviation standardized across cell types within each trait. Stars denote raw chromVAR z >= 2 (*) or z >= 3 (**), not adjusted-p significance. Pooled cell-type counts provide a descriptive comparison, not donor-level inference.", width = 110))
      save_plots_structured(
        plot,
        filetype = "png",
        width = 17,
        height = max(
          5.5,
          0.3 * dplyr::n_distinct(
            chromVAR_absolute_effect_deviation_tibble.cell_type_pseudobulk$GWAS_ID
          ) + 3.5
        )
      )
    }
  ),
  targets::tar_target(
    name = GWAS_absolute_effect_peak_variant_weight_tibble,
    description = "Allocate the absolute-effect heatmap's calibrated, capped peak weights to variants",
    command = get_GWAS_absolute_effect_peak_variant_weights(
      GWAS_input_records, GWAS_absolute_effect_weighting_status_tibble,
      genetic_enrichment_peak_ranges, genetic_enrichment_posterior_probability_cutoff),
    resources = get_tar_resources(RAM_GB_req = 40)
  ),
  targets::tar_target(
    name = chromVAR_absolute_effect_peak_contribution_tibble.cell_type_pseudobulk,
    description = "Decompose the existing absolute-effect annotation matrix into exact peak contributions",
    command = if (nrow(GWAS_absolute_effect_inputs_tibble) == 0L) {
      chromVAR_peak_contribution_tibble.cell_type_pseudobulk[0, ]
    } else get_GWAS_chromVAR_peak_contribution_tibble(
      chromVAR_background_record = chromVAR_background_record.cell_type_pseudobulk,
      psbulk_ATAC_data_matrix = cell_type_pseudobulk_counts_matrix.ATAC,
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
