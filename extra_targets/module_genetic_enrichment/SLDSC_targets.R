rlang::list2(
  targets::tar_target(
    name = LD_score_record.SLDSC.specificity,
    description = "Compute LD scores at the S-LDSC regression SNPs of an all-peaks annotation and, per WNN cell type, each peak's share of depth-normalized pseudobulk accessibility [part_of_graph:genetic_enrichment_SLDSC] [resource_observation: peak_RAM_GiB=5.0; n_nuclei=213k; aggregation=mixed_human_31x; method=slurm_monitor; date=2026-10-12] [resource_observation: peak_RAM_GiB=6.7; n_nuclei=740k; aggregation=unpublished; method=slurm_monitor; date=2026-10-16]",
    command = get_SLDSC_LD_score_record(
      annotation_matrix = get_SLDSC_peak_annotation_matrix(cell_type_pseudobulk_counts_matrix.ATAC),
      annotation_ranges = genetic_enrichment_peak_ranges,
      reference_dir = reference_dir.SLDSC,
      regression_SNP_tibble = regression_SNP_tibble.SLDSC,
      threads = 4L
    ),
    resources = get_tar_resources(cores_req = 4, RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = gene_specificity_record.SLDSC,
    description = "Score the windows (gene body +/- 100 kb) of expressed autosomal protein-coding genes by each gene's share of depth-normalized pseudobulk expression in each WNN cell type [part_of_graph:genetic_enrichment_SLDSC]",
    command = get_SLDSC_gene_specificity_record(pseudobulk_counts_matrix.GEX, marker_validated_Ensembl_annotations_GRanges_list,
      colnames(cell_type_pseudobulk_counts_matrix.ATAC))
  ),
  targets::tar_target(
    name = gene_annotation_record.SLDSC.specificity,
    description = "Mark the windows of all scored genes and, per WNN cell type, of its most specifically expressed tenth of them [part_of_graph:genetic_enrichment_SLDSC]",
    command = get_SLDSC_gene_annotation_record(gene_specificity_record.SLDSC)
  ),
  targets::tar_target(
    name = gene_LD_score_record.SLDSC.specificity,
    description = "Compute LD scores at the S-LDSC regression SNPs of the all-genes and cell-type gene-window annotations [part_of_graph:genetic_enrichment_SLDSC] [resource_observation: peak_RAM_GiB=11.9; n_nuclei=740k; aggregation=unpublished; method=slurm_monitor; date=2026-10-16]",
    command = get_SLDSC_LD_score_record(
      annotation_matrix = gene_annotation_record.SLDSC.specificity$annotation_matrix,
      annotation_ranges = gene_annotation_record.SLDSC.specificity$ranges,
      reference_dir = reference_dir.SLDSC,
      regression_SNP_tibble = regression_SNP_tibble.SLDSC,
      threads = 6L
    ),
    resources = get_tar_resources(cores_req = 6, RAM_GB_req = 32)
  ),
  targets::tar_target(
    name = records.SLDSC.specificity,
    description = "Fit S-LDSC for each cell type of one GWAS: its peak annotation with all peaks, its gene annotation with all genes, and both jointly, each with the baseline-LD model [resource_observation: peak_RAM_GiB=4.3; n_nuclei=213k; aggregation=mixed_human_31x; method=slurm_monitor; date=2026-10-12] [resource_observation: peak_RAM_GiB=10.5; n_nuclei=740k; aggregation=unpublished; method=slurm_monitor; date=2026-10-16] [resource_observation: peak_RAM_GiB=30 in most branches, one crashed twice at its 48 GB limit; n_nuclei=1.16M; aggregation=sizing_91x; method=slurm_monitor; date=2026-10-20]",
    command = get_SLDSC_cell_type_tibble(
      GWAS_input_record = GWAS_input_records,
      sumstats_record = sumstats_records.SLDSC,
      baseline_record = baseline_record.SLDSC,
      peak_LD_score_record = LD_score_record.SLDSC.specificity,
      gene_LD_score_record = gene_LD_score_record.SLDSC.specificity
    ),
    # Select the shared summary-statistics branches of this aggregation's GWAS.
    pattern = map(GWAS_input_records, slice(sumstats_records.SLDSC, index = SLDSC_sumstats_index)),
    iteration = "vector",
    resources = get_tar_resources(RAM_GB_req = 60)
  ),
  targets::tar_target(
    name = tibble.SLDSC.specificity,
    description = "Combine cell-type S-LDSC results with GWAS categories, BH q-values within each model and annotation, each GWAS's top set of cell types by tau* around its top significant one from paired jackknife differences, robustness to deleting influential jackknife blocks at q < 0.05 and cell-type support [part_of_graph:genetic_enrichment_SLDSC]",
    command = records.SLDSC.specificity |>
      dplyr::left_join(dplyr::select(GWAS_inputs_tibble, GWAS_ID, Category), by = "GWAS_ID", relationship = "many-to-one") |>
      dplyr::mutate(q = stats::p.adjust(p, method = "BH"), .by = c(model, annotation)) |>
      # tau* scales each coefficient by a constant, and so do its delete-one-block estimates.
      dplyr::mutate(
        get_top_set_tibble(tau_star, do.call(rbind, tau_delete_values) * tau_star / tau, q < 0.05, jackknife = TRUE),
        .by = c(GWAS_ID, model, annotation)
      ) |>
      dplyr::select(-tau_delete_values) |>
      dplyr::mutate(dplyr::across(dplyr::starts_with(c("share_in_top_", "z_without_top_")), \(x) dplyr::if_else(q < 0.05, x, NA_real_))) |>
      dplyr::left_join(cell_type_pseudobulk_support_tibble.ATAC, by = "cluster", relationship = "many-to-one") |>
      dplyr::relocate(GWAS_ID, Category, model, annotation, cluster, tau_star, z, p, q, in_top_set, z_below_anchor)
  ),
  # Matched-size annotation strategies: each ranks a cell type's territory by its score.
  tarchetypes::tar_map(
    values = tibble::tibble(
      strategy = c("accessible_peaks", "specific_peaks", "specific_genes"),
      score_matrix = rlang::exprs(
        get_SLDSC_measured_matrix(cell_type_pseudobulk_counts_matrix.ATAC),
        get_SLDSC_peak_annotation_matrix(cell_type_pseudobulk_counts_matrix.ATAC)[, -1L],
        gene_specificity_record.SLDSC$specificity_matrix
      ),
      score_ranges = rlang::exprs(genetic_enrichment_peak_ranges, genetic_enrichment_peak_ranges, gene_specificity_record.SLDSC$ranges)
    ),
    names = strategy,
    descriptions = NULL,
    delimiter = ".",
    targets::tar_target(
      name = annotation_record.enrichment,
      description = "Binarise each WNN cell type's annotation score, and ten random ones, to the top 0.5%, 1%, 2% and 5% of common reference SNPs and sum the baseline-LD annotations inside [part_of_graph:genetic_enrichment_SLDSC] [resource_observation: peak_RAM_GiB=4.4; n_nuclei=740k; aggregation=unpublished; method=slurm_monitor; date=2026-10-16] [resource_observation: peak_RAM_GiB=43.3 before a crash at its 48 GB limit; n_nuclei=1.16M; aggregation=sizing_91x; method=slurm_monitor; date=2026-10-19]",
      command = get_SLDSC_matched_annotation_record(score_matrix, score_ranges, reference_dir.SLDSC),
      resources = get_tar_resources(RAM_GB_req = 60)
    ),
    targets::tar_target(
      name = LD_score_record.enrichment,
      description = "Compute LD scores at the S-LDSC regression SNPs of the matched-size cell-type annotations [part_of_graph:genetic_enrichment_SLDSC] [resource_observation: peak_RAM_GiB=13.5; n_nuclei=740k; aggregation=unpublished; method=slurm_monitor; date=2026-10-16]",
      command = get_SLDSC_LD_score_record(
        annotation_matrix = annotation_record.enrichment$annotation_matrix,
        annotation_ranges = annotation_record.enrichment$ranges,
        reference_dir = reference_dir.SLDSC,
        regression_SNP_tibble = regression_SNP_tibble.SLDSC,
        threads = 6L
      ),
      resources = get_tar_resources(cores_req = 6, RAM_GB_req = 32)
    ),
    targets::tar_target(
      name = records.enrichment,
      description = "Estimate one GWAS's S-LDSC heritability enrichment, fitted with the baseline-LD model, and fine-mapped PIP enrichment in each matched-size cell-type annotation and the mean of the random ones [resource_observation: peak_RAM_GiB=12.4; n_nuclei=740k; aggregation=unpublished; method=slurm_monitor; date=2026-10-16]",
      command = get_SLDSC_matched_enrichment_tibble(
        GWAS_input_record = GWAS_input_records,
        sumstats_record = sumstats_records.SLDSC,
        baseline_record = baseline_record.SLDSC,
        annotation_record = annotation_record.enrichment,
        LD_score_record = LD_score_record.enrichment
      ) |>
        dplyr::mutate(strategy = strategy, .before = 1L),
      pattern = map(GWAS_input_records, slice(sumstats_records.SLDSC, index = SLDSC_sumstats_index)),
      iteration = "vector",
      resources = get_tar_resources(RAM_GB_req = 32)
    )
  ),
  targets::tar_target(
    name = tibble.enrichment,
    description = "Combine matched-size enrichment of all annotation strategies, test S-LDSC enrichment above 1 with BH q-values over GWAS and cell types within each strategy and size, find each GWAS's top set of cell types by enrichment from paired jackknife values, test fine-mapped PIP enrichment likewise from paired locus-bootstrap values, and add GWAS categories and whether the single-modality S-LDSC models enrich each GWAS-cell-type pair at q < 0.05 [part_of_graph:genetic_enrichment_SLDSC]",
    command = dplyr::bind_rows(
      records.enrichment.accessible_peaks,
      records.enrichment.specific_peaks,
      records.enrichment.specific_genes
    ) |>
      test_SLDSC_matched_enrichment() |>
      dplyr::left_join(dplyr::select(GWAS_inputs_tibble, GWAS_ID, Category, study_group), by = "GWAS_ID", relationship = "many-to-one") |>
      dplyr::left_join(get_SLDSC_enriched_pair_tibble(tibble.SLDSC.specificity), by = c("GWAS_ID", "cluster"), relationship = "many-to-one") |>
      dplyr::relocate(strategy, GWAS_ID, Category, cluster, share, enriched_in_models)
  ),
  targets::tar_target(
    name = summary_tibble.enrichment,
    description = "Summarise matched-size enrichment over GWAS per strategy, size, statistic and cell type, for enriched and all GWAS-cell-type pairs, with each GWAS's specificity for its enriched cell types [part_of_graph:genetic_enrichment_SLDSC]",
    command = summarise_SLDSC_matched_enrichment(tibble.enrichment)
  ),
  targets::tar_target(
    name = dropped_tibble.enrichment,
    description = "List the matched-size annotations of each strategy that cell types cannot fill",
    command = dplyr::bind_rows(
      accessible_peaks = annotation_record.enrichment.accessible_peaks$dropped_tibble,
      specific_peaks = annotation_record.enrichment.specific_peaks$dropped_tibble,
      specific_genes = annotation_record.enrichment.specific_genes$dropped_tibble,
      .id = "strategy"
    )
  ),
  tarchetypes::tar_file(
    name = curves_per_category_plots.SLDSC.enrichment,
    description = "Save each cell type's S-LDSC heritability enrichment in its most accessible peaks against annotation size, one file per GWAS category with one panel per GWAS and a shared y scale, with each GWAS's top sets. [checkpoint:genetic_enrichment]",
    command = {
      plots <- plot_SLDSC_enrichment_curves(
        tibble.enrichment,
        genetic_enrichment_compartment_patterns,
        dropped_tibble.enrichment
      )
      # Four panels per row; files of one or two GWAS keep room for the caption.
      n_panels <- purrr::map_int(plots, \(plot) dplyr::n_distinct(plot$data$GWAS_ID))
      save_plots_structured(plots, width = pmax(9, 0.8 + 3.3 * pmin(n_panels, 4)), height = 3.6 + 2.9 * ceiling(n_panels / 4))
    }
  ),
  tarchetypes::tar_file(
    name = strategy_curves_plot.SLDSC.enrichment,
    description = "Save S-LDSC heritability and fine-mapped PIP enrichment of peak and gene annotations against annotation size. [checkpoint:genetic_enrichment]",
    command = plot_SLDSC_matched_enrichment_curves(summary_tibble.enrichment, dropped_tibble.enrichment) |>
      save_plots_structured(width = 13, height = 7)
  ),
  tarchetypes::tar_file(
    name = strategy_size_comparison_plot.SLDSC.enrichment,
    description = "Save each cell type's enrichment in peak and gene annotations of matched size. [checkpoint:genetic_enrichment]",
    command = plot_SLDSC_matched_size_comparison(
      summary_tibble.enrichment,
      tibble.enrichment,
      dropped_tibble.enrichment
    ) |>
      save_plots_structured(width = 15, height = 3 + 0.17 * dplyr::n_distinct(tibble.enrichment$cluster))
  ),
  tarchetypes::tar_file(
    name = agreement_plot.comparison.enrichment,
    description = "Save the agreement of S-LDSC heritability and fine-mapped PIP enrichment in matched-size annotations. [checkpoint:genetic_enrichment]",
    command = plot_SLDSC_matched_enrichment_agreement(tibble.enrichment, dropped_tibble.enrichment) |>
      save_plots_structured(width = 14, height = 6.5)
  ),
  tarchetypes::tar_file(
    name = heatmap.SLDSC.specificity,
    description = "Save the cell-type-by-GWAS heatmap of S-LDSC heritability enrichment in cell-type-specific accessibility (tau*) with each GWAS's top set of cell types. [checkpoint:genetic_enrichment]",
    command = plot_SLDSC_specificity_heatmap(
      dplyr::filter(tibble.SLDSC.specificity, model == "peaks"),
      GWAS_inputs_tibble = GWAS_inputs_tibble,
      compartments_patterns = genetic_enrichment_compartment_patterns
    ) |>
      save_GWAS_heatmap(tibble.SLDSC.specificity$GWAS_ID)
  ),
  tarchetypes::tar_file(
    name = heatmaps.SLDSC.enrichment,
    description = "Save one cell-type-by-GWAS heatmap per size of S-LDSC heritability enrichment in each cell type's most accessible peaks, on the enrichment curves' colour limits, with each size's top sets. [checkpoint:genetic_enrichment]",
    command = plot_SLDSC_enrichment_heatmaps(
      tibble.enrichment,
      dplyr::filter(tibble.SLDSC.specificity, model == "peaks"),
      GWAS_inputs_tibble,
      genetic_enrichment_compartment_patterns,
      dropped_tibble.enrichment
    ) |>
      save_GWAS_heatmap(tibble.SLDSC.specificity$GWAS_ID)
  ),
  tarchetypes::tar_file(
    name = heatmaps.fine_mapping.enrichment,
    description = "Save one cell-type-by-GWAS heatmap per size of fine-mapped PIP enrichment in each cell type's most accessible peaks, with each size's FDR stars and top sets. [checkpoint:genetic_enrichment]",
    command = plot_SLDSC_enrichment_heatmaps(
      tibble.enrichment,
      dplyr::filter(tibble.SLDSC.specificity, model == "peaks"),
      GWAS_inputs_tibble,
      genetic_enrichment_compartment_patterns,
      dropped_tibble.enrichment,
      statistic = "PIP"
    ) |>
      save_GWAS_heatmap(tibble.SLDSC.specificity$GWAS_ID)
  ),
  targets::tar_target(
    name = locus_contribution_tibble.fine_mapping.enrichment,
    description = "Decompose one GWAS's fine-mapped PIP enrichment in each cell type's most accessible peaks into credible-set contributions",
    command = get_SLDSC_PIP_locus_contribution_tibble(
      GWAS_input_record = GWAS_input_records,
      annotation_record = annotation_record.enrichment.accessible_peaks,
      baseline_record = baseline_record.SLDSC
    ),
    pattern = map(GWAS_input_records),
    iteration = "vector",
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  tarchetypes::tar_file(
    name = locus_contribution_heatmaps.fine_mapping.enrichment,
    description = "Save one heatmap per GWAS of the credible sets that contribute most to each cell type's fine-mapped PIP enrichment at each size. [checkpoint:genetic_enrichment]",
    command = plot_SLDSC_PIP_locus_contribution_heatmaps(
      locus_contribution_tibble.fine_mapping.enrichment,
      tibble.enrichment,
      GWAS_locus_to_gene_tibble,
      genetic_enrichment_compartment_patterns
    ) |>
      save_plots_structured(width = 17, height = 12)
  ),
  tarchetypes::tar_file(
    name = z_score_plot.comparison.specificity,
    description = "Save per-GWAS scatter plots comparing S-LDSC and pseudobulk chromVAR cell-type z-scores. [checkpoint:genetic_enrichment]",
    command = plot_SLDSC_chromVAR_comparison(
      dplyr::filter(tibble.SLDSC.specificity, model == "peaks"),
      deviation_tibble.chromVAR.specificity
    ) |>
      save_plots_structured(width = 18, discrete_legend_n_distinct_max = Inf, height = 3 + 3 * ceiling(dplyr::n_distinct(
        tibble.SLDSC.specificity$GWAS_ID[!is.na(tibble.SLDSC.specificity$z)]) / 5))
  ),
  tarchetypes::tar_file(
    name = agreement_plot.comparison.specificity,
    description = "Save an overview of S-LDSC and pseudobulk chromVAR agreement per GWAS and per cell type against GWAS and cell-type predictors. [checkpoint:genetic_enrichment]",
    command = plot_SLDSC_chromVAR_agreement(
      dplyr::filter(tibble.SLDSC.specificity, model == "peaks"),
      deviation_tibble.chromVAR.specificity,
      GWAS_inputs_tibble
    ) |>
      save_plots_structured(width = 18, height = 14)
  ),
  tarchetypes::tar_file(
    name = modality_comparison_plot.SLDSC.specificity,
    description = "Save which enriched GWAS-cell-type pairs cell-type-specific accessibility and expression explain, alone and jointly, in S-LDSC. [checkpoint:genetic_enrichment]",
    command = plot_SLDSC_modality_comparison(tibble.SLDSC.specificity) |>
      save_plots_structured(width = 16, height = 5 + 0.2 * dplyr::n_distinct(tibble.SLDSC.specificity$GWAS_ID))
  ),
  tarchetypes::tar_file(
    name = block_influence_plot.SLDSC.specificity,
    description = "Save how much of each S-LDSC coefficient at q < 0.05 its most influential jackknife blocks carry, and whether it holds without them. [checkpoint:genetic_enrichment]",
    command = plot_SLDSC_block_influence(tibble.SLDSC.specificity) |>
      save_plots_structured(width = 11, height = 5)
  )
)
