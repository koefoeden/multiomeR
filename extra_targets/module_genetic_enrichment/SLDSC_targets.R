rlang::list2(
  # A pseudobulk of very few nuclei, such as a single unassigned nucleus, has no stable profile and,
  # with every cell type weighted equally, would distort the specificity of all the others.
  targets::tar_target(
    name = SLDSC_cell_types.cell_type_pseudobulk,
    description = "Select the WNN cell types with at least 100 nuclei for S-LDSC [part_of_graph:genetic_enrichment_SLDSC]",
    command = cell_type_pseudobulk_support_tibble.ATAC$cluster[cell_type_pseudobulk_support_tibble.ATAC$n_cells >= 100]
  ),
  targets::tar_target(
    name = SLDSC_LD_score_record.cell_type_pseudobulk,
    description = "Compute LD scores at the S-LDSC regression SNPs of an all-peaks annotation and, per WNN cell type, each peak's share of depth-normalized pseudobulk accessibility [part_of_graph:genetic_enrichment_SLDSC] [resource_observation: peak_RAM_GiB=5.0; n_nuclei=213k; aggregation=mixed_human_31x; method=slurm_monitor; date=2026-10-12] [resource_observation: peak_RAM_GiB=6.7; n_nuclei=740k; aggregation=unpublished; method=slurm_monitor; date=2026-10-16]",
    command = get_SLDSC_LD_score_record(
      annotation_matrix = get_SLDSC_peak_annotation_matrix(cell_type_pseudobulk_counts_matrix.ATAC, SLDSC_cell_types.cell_type_pseudobulk),
      annotation_ranges = genetic_enrichment_peak_ranges,
      reference_dir = SLDSC_reference_dir,
      regression_SNP_tibble = SLDSC_regression_SNP_tibble,
      threads = 4L
    ),
    resources = get_tar_resources(cores_req = 4, RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = SLDSC_gene_specificity_record.cell_type_pseudobulk,
    description = "Score the windows (gene body +/- 100 kb) of expressed autosomal protein-coding genes by each gene's share of depth-normalized pseudobulk expression in each WNN cell type [part_of_graph:genetic_enrichment_SLDSC]",
    command = get_SLDSC_gene_specificity_record(pseudobulk_counts_matrix.GEX, marker_validated_Ensembl_annotations_GRanges_list,
      SLDSC_cell_types.cell_type_pseudobulk)
  ),
  targets::tar_target(
    name = SLDSC_gene_annotation_record.cell_type_pseudobulk,
    description = "Mark the windows of all scored genes and, per WNN cell type, of its most specifically expressed tenth of them [part_of_graph:genetic_enrichment_SLDSC]",
    command = get_SLDSC_gene_annotation_record(SLDSC_gene_specificity_record.cell_type_pseudobulk)
  ),
  targets::tar_target(
    name = SLDSC_gene_LD_score_record.cell_type_pseudobulk,
    description = "Compute LD scores at the S-LDSC regression SNPs of the all-genes and cell-type gene-window annotations [part_of_graph:genetic_enrichment_SLDSC] [resource_observation: peak_RAM_GiB=11.9; n_nuclei=740k; aggregation=unpublished; method=slurm_monitor; date=2026-10-16]",
    command = get_SLDSC_LD_score_record(
      annotation_matrix = SLDSC_gene_annotation_record.cell_type_pseudobulk$annotation_matrix,
      annotation_ranges = SLDSC_gene_annotation_record.cell_type_pseudobulk$ranges,
      reference_dir = SLDSC_reference_dir,
      regression_SNP_tibble = SLDSC_regression_SNP_tibble,
      threads = 6L
    ),
    resources = get_tar_resources(cores_req = 6, RAM_GB_req = 32)
  ),
  targets::tar_target(
    name = SLDSC_cell_type_records.cell_type_pseudobulk,
    description = "Fit S-LDSC for each cell type of one GWAS: its peak annotation with all peaks, its gene annotation with all genes, and both jointly, each with the baseline-LD model [resource_observation: peak_RAM_GiB=4.3; n_nuclei=213k; aggregation=mixed_human_31x; method=slurm_monitor; date=2026-10-12] [resource_observation: peak_RAM_GiB=10.5; n_nuclei=740k; aggregation=unpublished; method=slurm_monitor; date=2026-10-16] [resource_observation: peak_RAM_GiB=30 in most branches, one crashed twice at its 48 GB limit; n_nuclei=1.16M; aggregation=sizing_91x; method=slurm_monitor; date=2026-10-20]",
    command = get_SLDSC_cell_type_tibble(
      GWAS_input_record = GWAS_input_records,
      sumstats_record = SLDSC_sumstats_records,
      baseline_record = SLDSC_baseline_record,
      peak_LD_score_record = SLDSC_LD_score_record.cell_type_pseudobulk,
      gene_LD_score_record = SLDSC_gene_LD_score_record.cell_type_pseudobulk
    ),
    # Select the shared summary-statistics branches of this aggregation's GWAS.
    pattern = map(GWAS_input_records, slice(SLDSC_sumstats_records, index = SLDSC_sumstats_index)),
    iteration = "vector",
    resources = get_tar_resources(RAM_GB_req = 60)
  ),
  targets::tar_target(
    name = SLDSC_cell_type_tibble.cell_type_pseudobulk,
    description = "Combine cell-type S-LDSC results with GWAS categories, BH q-values within each model and annotation, each GWAS's top set of cell types by tau* around its top significant one from paired jackknife differences, robustness to deleting influential jackknife blocks at q < 0.05 and cell-type support [part_of_graph:genetic_enrichment_SLDSC]",
    command = SLDSC_cell_type_records.cell_type_pseudobulk |>
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
      strategy = c("peaks", "genes"),
      score_matrix = rlang::exprs(
        get_SLDSC_peak_annotation_matrix(cell_type_pseudobulk_counts_matrix.ATAC, SLDSC_cell_types.cell_type_pseudobulk)[, -1L],
        SLDSC_gene_specificity_record.cell_type_pseudobulk$specificity_matrix
      ),
      score_ranges = rlang::exprs(genetic_enrichment_peak_ranges, SLDSC_gene_specificity_record.cell_type_pseudobulk$ranges)
    ),
    names = strategy,
    descriptions = NULL,
    delimiter = ".",
    targets::tar_target(
      name = SLDSC_matched_annotation_record.cell_type_pseudobulk,
      description = "Binarise each WNN cell type's annotation score, and a random one, to the top 0.5%, 1%, 2% and 5% of common reference SNPs and sum the baseline-LD annotations inside [part_of_graph:genetic_enrichment_SLDSC] [resource_observation: peak_RAM_GiB=4.4; n_nuclei=740k; aggregation=unpublished; method=slurm_monitor; date=2026-10-16] [resource_observation: peak_RAM_GiB=43.3 before a crash at its 48 GB limit; n_nuclei=1.16M; aggregation=sizing_91x; method=slurm_monitor; date=2026-10-19]",
      command = get_SLDSC_matched_annotation_record(score_matrix, score_ranges, SLDSC_reference_dir),
      resources = get_tar_resources(RAM_GB_req = 60)
    ),
    targets::tar_target(
      name = SLDSC_matched_LD_score_record.cell_type_pseudobulk,
      description = "Compute LD scores at the S-LDSC regression SNPs of the matched-size cell-type annotations [part_of_graph:genetic_enrichment_SLDSC] [resource_observation: peak_RAM_GiB=13.5; n_nuclei=740k; aggregation=unpublished; method=slurm_monitor; date=2026-10-16]",
      command = get_SLDSC_LD_score_record(
        annotation_matrix = SLDSC_matched_annotation_record.cell_type_pseudobulk$annotation_matrix,
        annotation_ranges = SLDSC_matched_annotation_record.cell_type_pseudobulk$ranges,
        reference_dir = SLDSC_reference_dir,
        regression_SNP_tibble = SLDSC_regression_SNP_tibble,
        threads = 6L
      ),
      resources = get_tar_resources(cores_req = 6, RAM_GB_req = 32)
    ),
    targets::tar_target(
      name = SLDSC_matched_enrichment_records.cell_type_pseudobulk,
      description = "Estimate one GWAS's S-LDSC heritability enrichment, fitted with the baseline-LD model, and fine-mapped PIP enrichment in each matched-size cell-type annotation [resource_observation: peak_RAM_GiB=12.4; n_nuclei=740k; aggregation=unpublished; method=slurm_monitor; date=2026-10-16]",
      command = get_SLDSC_matched_enrichment_tibble(
        GWAS_input_record = GWAS_input_records,
        sumstats_record = SLDSC_sumstats_records,
        baseline_record = SLDSC_baseline_record,
        annotation_record = SLDSC_matched_annotation_record.cell_type_pseudobulk,
        LD_score_record = SLDSC_matched_LD_score_record.cell_type_pseudobulk
      ) |>
        dplyr::mutate(strategy = strategy, .before = 1L),
      pattern = map(GWAS_input_records, slice(SLDSC_sumstats_records, index = SLDSC_sumstats_index)),
      iteration = "vector",
      resources = get_tar_resources(RAM_GB_req = 32)
    )
  ),
  targets::tar_target(
    name = SLDSC_matched_enrichment_tibble.cell_type_pseudobulk,
    description = "Combine matched-size enrichment of all annotation strategies with GWAS categories and whether the single-modality S-LDSC models enrich each GWAS-cell-type pair at q < 0.05 [part_of_graph:genetic_enrichment_SLDSC]",
    command = dplyr::bind_rows(
      SLDSC_matched_enrichment_records.cell_type_pseudobulk.peaks,
      SLDSC_matched_enrichment_records.cell_type_pseudobulk.genes
    ) |>
      dplyr::left_join(dplyr::select(GWAS_inputs_tibble, GWAS_ID, Category), by = "GWAS_ID", relationship = "many-to-one") |>
      dplyr::left_join(
        SLDSC_cell_type_tibble.cell_type_pseudobulk |>
          dplyr::filter(model != "joint") |>
          dplyr::summarise(enriched_in_models = any(q < 0.05), .by = c(GWAS_ID, cluster)),
        by = c("GWAS_ID", "cluster"),
        relationship = "many-to-one"
      ) |>
      dplyr::relocate(strategy, GWAS_ID, Category, cluster, share, enriched_in_models)
  ),
  targets::tar_target(
    name = SLDSC_matched_enrichment_summary_tibble.cell_type_pseudobulk,
    description = "Summarise matched-size enrichment over GWAS per strategy, size, statistic and cell type, for enriched and all GWAS-cell-type pairs, with each GWAS's specificity for its enriched cell types [part_of_graph:genetic_enrichment_SLDSC]",
    command = summarise_SLDSC_matched_enrichment(SLDSC_matched_enrichment_tibble.cell_type_pseudobulk)
  ),
  tarchetypes::tar_file(
    name = SLDSC_matched_enrichment_curves_plot.cell_type_pseudobulk,
    description = "Save S-LDSC heritability and fine-mapped PIP enrichment of peak and gene annotations against annotation size. [checkpoint:genetic_enrichment]",
    command = plot_SLDSC_matched_enrichment_curves(SLDSC_matched_enrichment_summary_tibble.cell_type_pseudobulk) |>
      save_plots_structured(width = 13, height = 7)
  ),
  tarchetypes::tar_file(
    name = SLDSC_matched_size_comparison_plot.cell_type_pseudobulk,
    description = "Save each cell type's enrichment in peak and gene annotations of matched size. [checkpoint:genetic_enrichment]",
    command = plot_SLDSC_matched_size_comparison(
      SLDSC_matched_enrichment_summary_tibble.cell_type_pseudobulk,
      SLDSC_matched_enrichment_tibble.cell_type_pseudobulk
    ) |>
      save_plots_structured(width = 15, height = 3 + 0.17 * dplyr::n_distinct(SLDSC_matched_enrichment_tibble.cell_type_pseudobulk$cluster))
  ),
  tarchetypes::tar_file(
    name = SLDSC_matched_enrichment_agreement_plot.cell_type_pseudobulk,
    description = "Save the agreement of S-LDSC heritability and fine-mapped PIP enrichment in matched-size annotations. [checkpoint:genetic_enrichment]",
    command = plot_SLDSC_matched_enrichment_agreement(SLDSC_matched_enrichment_tibble.cell_type_pseudobulk) |>
      save_plots_structured(width = 14, height = 6.5)
  ),
  tarchetypes::tar_file(
    name = SLDSC_heatmap.cell_type_pseudobulk,
    description = "Save the cell-type-by-GWAS S-LDSC heritability enrichment heatmap with each GWAS's top set of cell types. [checkpoint:genetic_enrichment]",
    command = plot_SLDSC_cell_type_heatmap(
      dplyr::filter(SLDSC_cell_type_tibble.cell_type_pseudobulk, model == "peaks"),
      GWAS_inputs_tibble = GWAS_inputs_tibble,
      compartments_patterns = genetic_enrichment_compartment_patterns
    ) |>
      save_GWAS_heatmap(SLDSC_cell_type_tibble.cell_type_pseudobulk$GWAS_ID)
  ),
  tarchetypes::tar_file(
    name = SLDSC_chromVAR_comparison_plot.cell_type_pseudobulk,
    description = "Save per-GWAS scatter plots comparing S-LDSC and pseudobulk chromVAR cell-type z-scores. [checkpoint:genetic_enrichment]",
    command = plot_SLDSC_chromVAR_comparison(
      dplyr::filter(SLDSC_cell_type_tibble.cell_type_pseudobulk, model == "peaks"),
      chromVAR_deviation_tibble.cell_type_pseudobulk
    ) |>
      save_plots_structured(width = 18, discrete_legend_n_distinct_max = Inf, height = 3 + 3 * ceiling(dplyr::n_distinct(
        SLDSC_cell_type_tibble.cell_type_pseudobulk$GWAS_ID[!is.na(SLDSC_cell_type_tibble.cell_type_pseudobulk$z)]) / 5))
  ),
  tarchetypes::tar_file(
    name = SLDSC_chromVAR_agreement_plot.cell_type_pseudobulk,
    description = "Save an overview of S-LDSC and pseudobulk chromVAR agreement per GWAS and per cell type against GWAS and cell-type predictors. [checkpoint:genetic_enrichment]",
    command = plot_SLDSC_chromVAR_agreement(
      dplyr::filter(SLDSC_cell_type_tibble.cell_type_pseudobulk, model == "peaks"),
      chromVAR_deviation_tibble.cell_type_pseudobulk,
      GWAS_inputs_tibble
    ) |>
      save_plots_structured(width = 18, height = 14)
  ),
  tarchetypes::tar_file(
    name = SLDSC_modality_comparison_plot.cell_type_pseudobulk,
    description = "Save which enriched GWAS-cell-type pairs cell-type-specific accessibility and expression explain, alone and jointly, in S-LDSC. [checkpoint:genetic_enrichment]",
    command = plot_SLDSC_modality_comparison(SLDSC_cell_type_tibble.cell_type_pseudobulk) |>
      save_plots_structured(width = 16, height = 5 + 0.2 * dplyr::n_distinct(SLDSC_cell_type_tibble.cell_type_pseudobulk$GWAS_ID))
  ),
  tarchetypes::tar_file(
    name = SLDSC_block_influence_plot.cell_type_pseudobulk,
    description = "Save how much of each S-LDSC coefficient at q < 0.05 its most influential jackknife blocks carry, and whether it holds without them. [checkpoint:genetic_enrichment]",
    command = plot_SLDSC_block_influence(SLDSC_cell_type_tibble.cell_type_pseudobulk) |>
      save_plots_structured(width = 11, height = 5)
  )
)
