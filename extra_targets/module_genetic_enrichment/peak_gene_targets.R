# Peak-gene strategies against GWAS signal, for aggregations with both the
# genetic_enrichment and peak_gene_correlation modules. Sourced after both
# modules, whose mapping values it joins: the activity of activity by contact
# comes from the peak-gene module's ATAC pseudobulks, which its GEM-well
# restriction may narrow.
SLDSC_peak_gene_tibble <- genetic_enrichment_tibble |>
  dplyr::select("aggregation", "SLDSC_sumstats_index") |>
  dplyr::mutate(genetic_validation_target_suffix = paste("genetic_validation.peak_gene_correlation", aggregation, sep = ".")) |>
  dplyr::inner_join(
    dplyr::select(peak_gene_correlation_restriction$values, "aggregation", "pseudobulk_counts_matrix.ATAC"),
    by = "aggregation"
  ) |>
  add_aggregation_target_syms(c(
    paste0(c(
      "GWAS_input_records", "GWAS_inputs_tibble", "cell_type_pseudobulk_counts_matrix.ATAC",
      "genetic_enrichment_peak_ranges", "tibble.SLDSC.specificity"
    ), ".genetic_enrichment"),
    paste0(c(
      "peak_gene_ABC_scores_tibble.WNN", "peak_gene_ABC_cell_groups_tibble.WNN",
      "peak_gene_correlation_hierarchical_results_tibble.WNN"
    ), ".peak_gene_correlation")
  )) |>
  dplyr::rename_with(\(name) sub("[.](genetic_enrichment|peak_gene_correlation)$", "", name))

# Each ranking orders the peaks of its territories; its annotations keep the top
# shares of each territory's common SNPs, and the territory ranking all of them.
SLDSC_peak_gene_rankings <- tibble::tibble(
  ranking = c("territory", "ABC", "ABC_conditional", "activity", "share", "distance",
    "correlation", "correlation_conditional", "correlation_ABC"),
  shares = c(list(1), rep(list(c(0.1, 0.25, 0.5)), 8L)),
  tiers = c(list(c("candidate", "tested"), c("candidate", "tested"), "candidate"),
    rep(list(c("candidate", "tested")), 3L), rep(list("tested"), 3L))
)

if (nrow(SLDSC_peak_gene_tibble) == 0L) list() else tarchetypes::tar_map(
  values = SLDSC_peak_gene_tibble,
  names = genetic_validation_target_suffix,
  descriptions = NULL,
  delimiter = ".",
  source("extra_targets/module_genetic_enrichment/SLDSC_peak_gene_targets.R")$value
)
