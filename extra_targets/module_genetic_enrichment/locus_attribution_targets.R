# Locus attribution, an experimental route for aggregations that set genetic_enrichment_locus_attribution. The module's
# targets.R sources it after building its mapping values.
locus_attribution_tibble <- genetic_enrichment_tibble |>
  dplyr::filter(purrr::map_lgl(genetic_enrichment_locus_attribution, isTRUE)) |>
  add_aggregation_target_syms(paste0(c(
    "GWAS_peak_variant_weight_records", "GWAS_input_records", "GWAS_inputs_tibble", "GWAS_locus_to_gene_tibble",
    "GWAS_metadata_tracks_plot", "cell_type_pseudobulk_counts_matrix.ATAC", "cell_type_pseudobulk_support_tibble.ATAC",
    "chromVAR_background_record.cell_type_pseudobulk"
  ), ".genetic_enrichment")) |>
  dplyr::rename_with(\(name) sub("[.]genetic_enrichment$", "", name))

if (nrow(locus_attribution_tibble) == 0L) list() else tarchetypes::tar_map(
  values = locus_attribution_tibble,
  names = genetic_enrichment_target_suffix,
  descriptions = NULL,
  delimiter = ".",
  targets::tar_target(
    name = locus_attribution_records,
    description = "Split each credible set's PIP in peaks between cell types by their accessibility, and fit the GWAS's cell-type weights over all its loci",
    command = get_locus_attribution_record(
      peak_variant_weight_tibble = GWAS_peak_variant_weight_records,
      GWAS_input_record = GWAS_input_records,
      pseudobulk_ATAC_data_matrix = cell_type_pseudobulk_counts_matrix.ATAC,
      chromVAR_background_record = chromVAR_background_record.cell_type_pseudobulk,
      support_tibble = cell_type_pseudobulk_support_tibble.ATAC
    ),
    pattern = map(GWAS_peak_variant_weight_records, GWAS_input_records),
    iteration = "list",
    resources = get_tar_resources(RAM_GB_req = 4)
  ),
  targets::tar_target(
    name = locus_attribution_trait_tibble,
    description = "Combine each GWAS's cell-type weights and chance shares",
    command = dplyr::bind_rows(purrr::map(locus_attribution_records, "trait_tibble"))
  ),
  targets::tar_target(
    name = locus_attribution_share_tibble,
    description = "Combine each locus's split of its PIP in peaks between cell types",
    command = dplyr::bind_rows(purrr::map(locus_attribution_records, "share_tibble"))
  ),
  targets::tar_target(
    name = locus_attribution_locus_tibble,
    description = "Summarise each locus's association, likely genes, top cell types and top compartment with their folds",
    command = summarise_locus_attribution(
      locus_tibble = dplyr::bind_rows(purrr::map(locus_attribution_records, "locus_tibble")),
      share_tibble = locus_attribution_share_tibble,
      trait_tibble = locus_attribution_trait_tibble,
      locus_to_gene_tibble = GWAS_locus_to_gene_tibble,
      GWAS_inputs_tibble = GWAS_inputs_tibble,
      compartments_patterns = genetic_enrichment_compartment_patterns
    )
  ),
  tarchetypes::tar_file(
    name = trait_heatmap.locus_attribution,
    description = "Save the GWAS-by-cell-type heatmap of each GWAS's cell-type weights over all its loci against chance. [checkpoint:genetic_enrichment]",
    command = plot_locus_attribution_trait_heatmap(
      locus_attribution_trait_tibble,
      GWAS_metadata_tracks_plot = GWAS_metadata_tracks_plot,
      support_tibble = cell_type_pseudobulk_support_tibble.ATAC,
      compartments_patterns = genetic_enrichment_compartment_patterns
    ) |>
      save_GWAS_heatmap(locus_attribution_trait_tibble$GWAS_ID)
  ),
  tarchetypes::tar_file(
    name = locus_heatmaps.locus_attribution,
    description = "Save one heatmap per GWAS of the cell types whose accessibility its loci with the most PIP in peaks favour. [checkpoint:genetic_enrichment]",
    command = plot_locus_attribution_heatmaps(
      locus_attribution_locus_tibble,
      share_tibble = locus_attribution_share_tibble,
      trait_tibble = locus_attribution_trait_tibble,
      compartments_patterns = genetic_enrichment_compartment_patterns
    ) |>
      save_plots_structured(width = 5 + 0.17 * dplyr::n_distinct(locus_attribution_share_tibble$cluster), height = 11)
  ),
  tarchetypes::tar_file(
    name = locus_table.locus_attribution,
    description = "Export each locus's likely genes, strongest association, PIP in peaks, top cell types and top compartment. [checkpoint:genetic_enrichment]",
    command = save_locus_attribution_table(locus_attribution_locus_tibble)
  ),
  tarchetypes::tar_file(
    name = compartment_table.locus_attribution,
    description = "Export each locus's share and fold over an even share per cell-type compartment, beside the compartment's weight in the trait. [checkpoint:genetic_enrichment]",
    command = summarise_locus_attribution_compartments(
      locus_attribution_locus_tibble,
      share_tibble = locus_attribution_share_tibble,
      trait_tibble = locus_attribution_trait_tibble,
      compartments_patterns = genetic_enrichment_compartment_patterns
    ) |>
      save_locus_attribution_table()
  )
)
