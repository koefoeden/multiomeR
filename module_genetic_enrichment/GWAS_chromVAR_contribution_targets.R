rlang::list2(
  targets::tar_target(
    name = GWAS_peak_variant_weight_records,
    description = "Allocate capped trait-level GWAS peak weights to overlapping credible-set variants",
    command = get_GWAS_chromVAR_peak_variant_weight_tibble(
      GWAS_input_record = GWAS_input_records,
      peak_ranges = genetic_enrichment_peak_ranges,
      posterior_probability_cutoff = genetic_enrichment_posterior_probability_cutoff,
      posterior_probability_weighting_function = posterior_probability_weighting_function
    ),
    pattern = map(GWAS_input_records),
    iteration = "list",
    resources = get_tar_resources(RAM_GB_req = 40)
  ),
  targets::tar_target(
    name = GWAS_peak_variant_weight_tibble,
    description = "Combine peak-to-credible-variant weight allocations across all enabled GWAS",
    command = dplyr::bind_rows(GWAS_peak_variant_weight_records)
  ),
  targets::tar_target(
    name = chromVAR_peak_contribution_tibble.cell_type_pseudobulk,
    description = "Decompose each cell-type GWAS chromVAR heatmap value into exact peak contributions [part_of_graph:genetic_enrichment_cell_type_contributions]",
    command = get_GWAS_chromVAR_peak_contribution_tibble(
      chromVAR_background_record = chromVAR_background_record.cell_type_pseudobulk,
      psbulk_ATAC_data_matrix = cell_type_pseudobulk_counts_matrix.ATAC,
      chromVAR_obj = chromVAR_obj.ATAC,
      annotation_matrix = GWAS_peak_weight_matrix,
      GWAS_inputs_tibble = GWAS_inputs_tibble
    ),
    resources = get_tar_resources(RAM_GB_req = 60)
  ),
  targets::tar_target(
    name = chromVAR_deviation_tibble.cell_type_pseudobulk,
    description = "Sum peak contributions into the cell-type GWAS chromVAR heatmap table with z-score support labels [part_of_graph:genetic_enrichment_cell_type_contributions]",
    command = summarize_GWAS_chromVAR_peak_contributions(
      peak_contribution_tibble = chromVAR_peak_contribution_tibble.cell_type_pseudobulk,
      cell_type_support_tibble = cell_type_pseudobulk_support_tibble.ATAC
    )
  ),
  targets::tar_target(
    name = chromVAR_variant_contribution_tibble.cell_type_pseudobulk,
    description = "Allocate exact peak-level GWAS chromVAR contributions to credible-set variants [part_of_graph:genetic_enrichment_cell_type_contributions]",
    command = get_GWAS_chromVAR_variant_contribution_tibble(
      peak_contribution_tibble = chromVAR_peak_contribution_tibble.cell_type_pseudobulk,
      peak_variant_weight_tibble = GWAS_peak_variant_weight_tibble
    ),
    resources = get_tar_resources(RAM_GB_req = 40)
  ),
  targets::tar_target(
    name = chromVAR_locus_contribution_tibble.cell_type_pseudobulk,
    description = "Sum exact credible-set variant contributions into cell-type-by-GWAS locus contribution [part_of_graph:genetic_enrichment_cell_type_contributions]",
    command = get_GWAS_chromVAR_locus_contribution_tibble(
      chromVAR_variant_contribution_tibble.cell_type_pseudobulk
    )
  ),
  targets::tar_target(
    name = chromVAR_contribution_reconciliation_tibble.cell_type_pseudobulk,
    description = "Reconcile peak, variant, and locus contribution sums with the original chromVAR heatmap values",
    command = get_GWAS_chromVAR_contribution_reconciliation_tibble(
      peak_contribution_tibble = chromVAR_peak_contribution_tibble.cell_type_pseudobulk,
      variant_contribution_tibble = chromVAR_variant_contribution_tibble.cell_type_pseudobulk,
      locus_contribution_tibble = chromVAR_locus_contribution_tibble.cell_type_pseudobulk
    )
  ),
  tarchetypes::tar_file(
    name = chromVAR_locus_contribution_heatmaps.cell_type_pseudobulk,
    description = "Save one cell-type-by-locus contribution heatmap per enabled GWAS. [checkpoint:genetic_enrichment]",
    command = plot_GWAS_locus_contribution_heatmaps(
      locus_contribution_tibble = chromVAR_locus_contribution_tibble.cell_type_pseudobulk,
      chromVAR_deviation_tibble = chromVAR_deviation_tibble.cell_type_pseudobulk
    ) |>
      save_plots_structured(width = 14, height = 8)
  ),
  tarchetypes::tar_file(
    name = chromVAR_locus_contribution_waterfalls.cell_type_pseudobulk,
    description = "Save locus-contribution waterfalls for every enabled GWAS and cell type. [checkpoint:genetic_enrichment]",
    command = chromVAR_locus_contribution_tibble.cell_type_pseudobulk |>
      plot_GWAS_locus_contribution_waterfalls() |>
      save_plots_structured(width = 12, height = 7)
  ),
  tarchetypes::tar_file(
    name = chromVAR_variant_contribution_detail_plots.cell_type_pseudobulk,
    description = "Save variant contribution, ATAC coverage, and consensus-peak tracks for top-contributing loci. [checkpoint:genetic_enrichment]",
    command = plot_GWAS_variant_contribution_details(
      variant_contribution_tibble = chromVAR_variant_contribution_tibble.cell_type_pseudobulk,
      locus_contribution_tibble = chromVAR_locus_contribution_tibble.cell_type_pseudobulk,
      consensus_peak_GRanges = consensus_peak_GRanges.ATAC,
      fragments = combined_BPCells_fragment_obj.ATAC,
      metadata_tibble = metadata_w_cell_types_tibble.WNN
    ) |>
      save_plots_structured(width = 13, height = 11),
    resources = get_tar_resources(RAM_GB_req = 40)
  )
)
