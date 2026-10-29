rlang::list2(
  targets::tar_target(
    name = score_tibble,
    description = "Score each cell type's peaks by activity by contact, peak-gene correlation and both, with their activity, accessibility share, distance and covariate stratum, among activity-by-contact candidate pairs and correlation-tested pairs [part_of_graph:genetic_enrichment_SLDSC] [resource_observation: peak_RAM_GiB=6.8; walltime_s=65; n_rows=7.1M; n_cell_groups=13; aggregation=unpublished; method=VmHWM; date=2026-10-23]",
    command = get_SLDSC_peak_gene_score_tibble(
      ABC_scores_tibble = peak_gene_ABC_scores_tibble.WNN,
      cell_groups_tibble = peak_gene_ABC_cell_groups_tibble.WNN,
      hierarchical_results_tibble = peak_gene_correlation_hierarchical_results_tibble.WNN,
      ATAC_pseudobulk_matrix = pseudobulk_counts_matrix.ATAC,
      cell_type_pseudobulk_counts_matrix = cell_type_pseudobulk_counts_matrix.ATAC,
      peak_ranges = genetic_enrichment_peak_ranges
    ),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  tarchetypes::tar_map(
    values = SLDSC_peak_gene_rankings,
    names = ranking,
    descriptions = NULL,
    delimiter = ".",
    targets::tar_target(
      name = annotation_record,
      description = "Binarise one peak-gene ranking to its top shares of each territory's common SNPs and sum the baseline-LD annotations inside [part_of_graph:genetic_enrichment_SLDSC]",
      command = get_SLDSC_matched_annotation_record(
        get_SLDSC_peak_gene_score_matrix(score_tibble, ranking, tiers, genetic_enrichment_peak_ranges),
        genetic_enrichment_peak_ranges, SLDSC_reference_dir, shares, within_territory = TRUE
      ),
      resources = get_tar_resources(RAM_GB_req = 60)
    ),
    targets::tar_target(
      name = LD_score_record,
      description = "Compute LD scores at the S-LDSC regression SNPs of one peak-gene ranking's annotations [part_of_graph:genetic_enrichment_SLDSC]",
      command = get_SLDSC_LD_score_record(
        annotation_matrix = annotation_record$annotation_matrix,
        annotation_ranges = annotation_record$ranges,
        reference_dir = SLDSC_reference_dir,
        regression_SNP_tibble = SLDSC_regression_SNP_tibble,
        threads = 6L
      ),
      resources = get_tar_resources(cores_req = 6, RAM_GB_req = 32)
    ),
    targets::tar_target(
      name = enrichment_records,
      description = "Estimate one GWAS's S-LDSC heritability and fine-mapped PIP enrichment in one peak-gene ranking's annotations",
      command = get_SLDSC_matched_enrichment_tibble(
        GWAS_input_record = GWAS_input_records,
        sumstats_record = SLDSC_sumstats_records,
        baseline_record = SLDSC_baseline_record,
        annotation_record = annotation_record,
        LD_score_record = LD_score_record
      ) |>
        dplyr::mutate(ranking = ranking, .before = 1L),
      pattern = map(GWAS_input_records, slice(SLDSC_sumstats_records, index = SLDSC_sumstats_index)),
      iteration = "vector",
      resources = get_tar_resources(RAM_GB_req = 32)
    )
  ),
  targets::tar_target(
    name = enrichment_tibble,
    description = "Combine the enrichment of all peak-gene rankings with GWAS categories and whether the single-modality S-LDSC models enrich each GWAS-cell-type pair at q < 0.05 [part_of_graph:genetic_enrichment_SLDSC]",
    command = dplyr::bind_rows(!!!rlang::syms(paste0("enrichment_records.", SLDSC_peak_gene_rankings$ranking))) |>
      dplyr::select(-h2_share_delete_values, -PIP_bootstrap_shares) |>
      tidyr::separate_wider_delim("cluster", "/", names = c("tier", "cluster")) |>
      dplyr::left_join(dplyr::select(GWAS_inputs_tibble, GWAS_ID, Category), by = "GWAS_ID", relationship = "many-to-one") |>
      dplyr::left_join(get_SLDSC_enriched_pair_tibble(SLDSC_cell_type_tibble.cell_type_pseudobulk), by = c("GWAS_ID", "cluster"), relationship = "many-to-one") |>
      dplyr::relocate(tier, ranking, GWAS_ID, Category, cluster, share, enriched_in_models)
  ),
  targets::tar_target(
    name = enrichment_summary_tibble,
    description = "Contrast each peak-gene ranking's enrichment with its territory and baselines in enriched and other GWAS-cell-type pairs, with GWAS-bootstrap intervals [part_of_graph:genetic_enrichment_SLDSC]",
    command = summarise_SLDSC_peak_gene_enrichment(enrichment_tibble),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  tarchetypes::tar_file(
    name = uplift_curves_plot,
    description = "Save the GWAS signal of each peak-gene ranking's top peaks relative to its territory against annotation size. [checkpoint:genetic_enrichment]",
    command = plot_SLDSC_peak_gene_uplift_curves(enrichment_summary_tibble) |>
      save_plots_structured(width = 16, height = 8)
  ),
  tarchetypes::tar_file(
    name = contrast_plot,
    description = "Save the prespecified contrasts of peak-gene strategies with their territories and baselines. [checkpoint:genetic_enrichment]",
    command = plot_SLDSC_peak_gene_contrasts(enrichment_summary_tibble) |>
      save_plots_structured(width = 13, height = 10)
  ),
  tarchetypes::tar_file(
    name = cell_type_plot,
    description = "Save each cell type's GWAS signal in the top peaks of each peak-gene ranking relative to its territory. [checkpoint:genetic_enrichment]",
    command = plot_SLDSC_peak_gene_cell_types(enrichment_summary_tibble) |>
      save_plots_structured(width = 16, height = 4 + 0.35 * dplyr::n_distinct(enrichment_tibble$cluster))
  ),
  tarchetypes::tar_file(
    name = design_plot,
    description = "Save the territory sizes and covariate balance of the peak-gene strategy comparison. [checkpoint:genetic_enrichment]",
    command = plot_SLDSC_peak_gene_design(score_tibble, enrichment_tibble, SLDSC_peak_gene_rankings) |>
      save_plots_structured(width = 14, height = 8)
  )
)
