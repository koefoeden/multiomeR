rlang::list2(
  targets::tar_target(
    name = peak_variant_weight_records.trait_level.attribution,
    description = "Allocate capped trait-level GWAS peak weights to overlapping credible-set variants",
    command = get_GWAS_chromVAR_peak_variant_weight_tibble(
      GWAS_input_record = input_records.trait_level,
      peak_ranges = genetic_enrichment_peak_ranges,
      posterior_probability_cutoff = genetic_enrichment_posterior_probability_cutoff,
      posterior_probability_weighting_function = posterior_probability_weighting_function
    ),
    pattern = map(input_records.trait_level),
    iteration = "list",
    resources = get_tar_resources(RAM_GB_req = 40)
  ),
  targets::tar_target(
    name = peak_variant_weight_tibble.trait_level.attribution,
    description = "Combine peak-to-credible-variant weight allocations across all enabled GWAS",
    command = dplyr::bind_rows(peak_variant_weight_records.trait_level.attribution)
  ),
  targets::tar_target(
    name = peak_attribution_tibble.trait_level.cell_type,
    description = "Decompose each cell-type GWAS chromVAR heatmap value into exact peak contributions [part_of_graph:genetic_enrichment_cell_type_attribution]",
    command = get_GWAS_chromVAR_peak_attribution_tibble(
      chromVAR_background_record = chromVAR_background_record.trait_level.cell_type_pseudobulk,
      psbulk_ATAC_data_matrix = cell_type_pseudobulk_counts_matrix.ATAC,
      chromVAR_obj = chromVAR_obj.ATAC,
      annotation_matrix = peak_weight_matrix.trait_level.pseudobulk,
      GWAS_inputs_tibble = GWAS_inputs_tibble
    ),
    resources = get_tar_resources(RAM_GB_req = 60)
  ),
  targets::tar_target(
    name = chromVAR_deviation_tibble.trait_level.pseudobulk,
    description = "Sum peak contributions into the cell-type GWAS chromVAR heatmap table with z-score support labels [part_of_graph:genetic_enrichment_cell_type_attribution]",
    command = summarize_GWAS_chromVAR_peak_attribution(
      peak_attribution_tibble = peak_attribution_tibble.trait_level.cell_type,
      cell_type_support_tibble = cell_type_pseudobulk_support_tibble.ATAC
    )
  ),
  targets::tar_target(
    name = variant_attribution_tibble.trait_level.cell_type,
    description = "Allocate exact peak-level GWAS chromVAR contributions to credible-set variants [part_of_graph:genetic_enrichment_cell_type_attribution]",
    command = get_GWAS_chromVAR_variant_attribution_tibble(
      peak_attribution_tibble = peak_attribution_tibble.trait_level.cell_type,
      peak_variant_weight_tibble = peak_variant_weight_tibble.trait_level.attribution
    ),
    resources = get_tar_resources(RAM_GB_req = 40)
  ),
  targets::tar_target(
    name = locus_attribution_tibble.trait_level.cell_type,
    description = "Sum exact credible-set variant contributions into cell-type-by-GWAS locus attribution [part_of_graph:genetic_enrichment_cell_type_attribution]",
    command = get_GWAS_chromVAR_locus_attribution_tibble(
      variant_attribution_tibble.trait_level.cell_type
    )
  ),
  targets::tar_target(
    name = attribution_QC_tibble.trait_level.cell_type,
    description = "Reconcile peak, variant, and locus attribution sums with the original chromVAR heatmap values",
    command = get_GWAS_chromVAR_attribution_QC_tibble(
      peak_attribution_tibble = peak_attribution_tibble.trait_level.cell_type,
      variant_attribution_tibble = variant_attribution_tibble.trait_level.cell_type,
      locus_attribution_tibble = locus_attribution_tibble.trait_level.cell_type
    )
  ),
  tarchetypes::tar_file(
    name = locus_attribution_heatmaps.trait_level.cell_type,
    description = "Save one cell-type-by-locus contribution heatmap per enabled GWAS. [checkpoint:genetic_enrichment]",
    command = plot_GWAS_locus_attribution_heatmaps(
      locus_attribution_tibble = locus_attribution_tibble.trait_level.cell_type,
      chromVAR_deviation_tibble = chromVAR_deviation_tibble.trait_level.pseudobulk
    ) |>
      save_plots_structured(width = 14, height = 8)
  ),
  tarchetypes::tar_file(
    name = locus_attribution_waterfalls.trait_level.cell_type,
    description = "Save locus-attribution waterfalls for every enabled GWAS and cell type. [checkpoint:genetic_enrichment]",
    command = locus_attribution_tibble.trait_level.cell_type |>
      plot_GWAS_locus_attribution_waterfalls() |>
      save_plots_structured(width = 12, height = 7)
  ),
  tarchetypes::tar_file(
    name = variant_attribution_detail_plots.trait_level.cell_type,
    description = "Save variant contribution, ATAC coverage, and consensus-peak tracks for top attributed loci. [checkpoint:genetic_enrichment]",
    command = plot_GWAS_variant_attribution_details(
      variant_attribution_tibble = variant_attribution_tibble.trait_level.cell_type,
      locus_attribution_tibble = locus_attribution_tibble.trait_level.cell_type,
      consensus_peak_GRanges = consensus_peak_GRanges.ATAC,
      fragments = combined_BPCells_fragment_obj.ATAC,
      metadata_tibble = metadata_w_cell_types_tibble.WNN
    ) |>
      save_plots_structured(width = 13, height = 11),
    resources = get_tar_resources(RAM_GB_req = 40)
  )
)
