rlang::list2(
  targets::tar_target(
    name = inputs_tibble.gene_level.pseudobulk,
    description = "Restrict GWAS inputs to the configured gene-level pseudobulk GWAS chromVAR traits [part_of_graph:genetic_enrichment_pseudobulk]",
    command = get_GWAS_gene_chromVAR_inputs_tibble(
      GWAS_analysis_inputs_tibble = GWAS_analysis_inputs_tibble,
      selected_GWAS_IDs = genetic_enrichment_psbulk_GWAS_gene_chromVAR_GWAS_IDs
    ),
    iteration = "vector"
  ),
  targets::tar_target(
    name = credible_set_variants_tibble.gene_level.pseudobulk,
    description = "Expose Open Targets 95 percent credible-set variant rows for configured gene-level GWAS chromVAR traits [part_of_graph:genetic_enrichment_pseudobulk]",
    command = get_open_targets_credible_set_variants_tibble(
      GWAS_inputs_tibble = inputs_tibble.gene_level.pseudobulk,
      open_targets_credible_set_dataset_path = open_targets_credible_set_dataset_path
    ) |>
      dplyr::select(-is95CredibleSet, -is99CredibleSet),
    resources = get_tar_resources(RAM_GB_req = 40)
  ),
  targets::tar_target(
    name = L2G_tibble.gene_level.pseudobulk,
    description = "Read Open Targets L2G evidence for the gene-level GWAS chromVAR credible-set loci [part_of_graph:genetic_enrichment_pseudobulk]",
    command = get_GWAS_gene_chromVAR_L2G_tibble(
      GWAS_gene_chromVAR_credible_set_variants_tibble = credible_set_variants_tibble.gene_level.pseudobulk,
      open_targets_gwas_credible_sets_evidence_dataset_path = open_targets_gwas_credible_sets_evidence_dataset_path
    ),
    resources = get_tar_resources(RAM_GB_req = 40)
  ),
  targets::tar_target(
    name = peak_weight_records.gene_level.pseudobulk,
    description = "Build sparse peak-by-target L2G-weighted GWAS chromVAR annotations for one configured GWAS",
    command = get_GWAS_gene_chromVAR_peak_weight_record(
      GWAS_input_tibble = inputs_tibble.gene_level.pseudobulk,
      GWAS_gene_chromVAR_credible_set_variants_tibble = credible_set_variants_tibble.gene_level.pseudobulk,
      GWAS_gene_chromVAR_L2G_tibble = L2G_tibble.gene_level.pseudobulk,
      peak_ranges = genetic_enrichment_peak_ranges,
      posterior_probability_cutoff = genetic_enrichment_posterior_probability_cutoff,
      posterior_probability_weighting_function = posterior_probability_weighting_function
    ),
    pattern = map(inputs_tibble.gene_level.pseudobulk),
    iteration = "list",
    resources = get_tar_resources(RAM_GB_req = 60)
  ),
  targets::tar_target(
    name = feature_metadata_tibble.gene_level.pseudobulk,
    description = "Summarize support metrics for each GWAS-gene chromVAR annotation feature",
    command = get_GWAS_gene_chromVAR_feature_metadata_tibble(
      GWAS_gene_chromVAR_peak_weight_records = peak_weight_records.gene_level.pseudobulk,
      min_n_peaks = 3
    ) |>
      add_open_targets_target_metadata(
        open_targets_target_dataset_path = open_targets_target_dataset_path
      )
  ),
  targets::tar_target(
    name = peak_weight_matrix.gene_level.pseudobulk,
    description = "Combine supported peak-by-GWAS-gene annotations into one sparse peak-by-feature matrix [part_of_graph:genetic_enrichment_pseudobulk]",
    command = get_GWAS_gene_chromVAR_peak_weight_matrix(
      GWAS_gene_chromVAR_peak_weight_records = peak_weight_records.gene_level.pseudobulk,
      GWAS_gene_chromVAR_feature_metadata_tibble = feature_metadata_tibble.gene_level.pseudobulk,
      RSE_ATAC = chromVAR_obj.ATAC
    ),
    resources = get_tar_resources(RAM_GB_req = 60)
  ),
  targets::tar_target(
    name = activity_matrix.gene_level.pseudobulk,
    description = "Compute pseudobulk GWAS-gene chromVAR activity scores from BPCells-backed pseudobulk ATAC counts [part_of_graph:genetic_enrichment_pseudobulk]",
    command = if (ncol(peak_weight_matrix.gene_level.pseudobulk) == 0) {
      matrix(
        nrow = 0,
        ncol = ncol(pseudobulk_counts_matrix.ATAC),
        dimnames = list(character(), colnames(pseudobulk_counts_matrix.ATAC))
      )
    } else {
      get_pseudobulk_chromVAR_activity_matrix(
        psbulk_ATAC_data_matrix = pseudobulk_counts_matrix.ATAC,
        chromVAR_obj = chromVAR_obj.ATAC,
        annotation_matrix = peak_weight_matrix.gene_level.pseudobulk,
        normalize = FALSE
      )
    },
    resources = get_tar_resources(RAM_GB_req = 80)
  ),
  targets::tar_target(
    name = model_tibble.gene_level.pseudobulk,
    description = "Build a dynamic tibble with one row per pseudobulk GWAS-gene chromVAR model",
    command = if (nrow(activity_matrix.gene_level.pseudobulk) == 0) {
      model_tibble.trait_level.pseudobulk[0, ]
    } else {
      model_tibble.trait_level.pseudobulk
    },
    iteration = "vector"
  ),
  targets::tar_target(
    name = filtered_matrix.gene_level.pseudobulk,
    description = "Filter the pseudobulk GWAS-gene chromVAR activity matrix for each model",
    command = filter_psbulk_data_matrix(
      psbulk_data_matrix = activity_matrix.gene_level.pseudobulk,
      psbulk_feature_dynamic_tibble = model_tibble.gene_level.pseudobulk,
      extended_donor_id_metadata_tibble = genetic_enrichment_donor_id_metadata_tibble.extended,
      sample_depth_tibble = pseudobulk_depth_tibble.ATAC,
      min_sample_counts = genetic_enrichment_psbulk_GWAS_chromVAR_min_ATAC_counts
    ),
    pattern = map(model_tibble.gene_level.pseudobulk),
    resources = get_tar_resources(RAM_GB_req = 60)
  ),
  targets::tar_target(
    name = fit.gene_level.pseudobulk,
    description = "Fit limma models to pseudobulk GWAS-gene chromVAR activity matrices for each model",
    command = fit_psbulk_feature_matrix_model(
      psbulk_feature_matrix = filtered_matrix.gene_level.pseudobulk,
      extended_donor_id_metadata_tibble = genetic_enrichment_donor_id_metadata_tibble.extended,
      psbulk_feature_dynamic_tibble = model_tibble.gene_level.pseudobulk
    ),
    pattern = map(filtered_matrix.gene_level.pseudobulk, model_tibble.gene_level.pseudobulk),
    resources = get_tar_resources(RAM_GB_req = 80)
  ),
  targets::tar_target(
    name = results_per_model_tibble.gene_level.pseudobulk,
    description = "Extract GWAS-gene chromVAR differential activity results for each model",
    command = get_psbulk_feature_model_results(
      psbulk_feature_matrix_fit = fit.gene_level.pseudobulk,
      psbulk_feature_dynamic_tibble = model_tibble.gene_level.pseudobulk
    ),
    pattern = map(fit.gene_level.pseudobulk, model_tibble.gene_level.pseudobulk)
  ),
  targets::tar_target(
    name = results_tibble.gene_level.pseudobulk,
    description = "Combine GWAS-gene chromVAR differential activity results and add within-GWAS FDR [part_of_graph:genetic_enrichment_pseudobulk]",
    command = results_per_model_tibble.gene_level.pseudobulk |>
      format_psbulk_GWAS_gene_chromVAR_results(
        GWAS_gene_chromVAR_feature_metadata_tibble = feature_metadata_tibble.gene_level.pseudobulk
      )
  ),
  targets::tar_target(
    name = volcano_plots.gene_level.pseudobulk,
    description = "Plot volcano plots of pseudobulk GWAS-gene chromVAR results per contrast",
    command = plot_psbulk_GWAS_gene_chromVAR_volcanoes(
      psbulk_GWAS_gene_chromVAR_results_tibble = results_tibble.gene_level.pseudobulk,
      model_name = model_tibble.gene_level.pseudobulk$model_name
    ),
    pattern = map(model_tibble.gene_level.pseudobulk)
  ),
  tarchetypes::tar_file(
    name = volcano_plot_files.gene_level.pseudobulk,
    description = "Save GWAS-gene chromVAR volcano plots per model and contrast to file. [checkpoint:genetic_enrichment]",
    command = save_plots_structured(
      plots = volcano_plots.gene_level.pseudobulk,
      override_suffix = model_tibble.gene_level.pseudobulk$model_name,
      dyn_suffix_in_subdir = TRUE
    ),
    pattern = map(volcano_plots.gene_level.pseudobulk, model_tibble.gene_level.pseudobulk)
  ),
  tarchetypes::tar_file(
    name = QC_plots.gene_level.pseudobulk,
    description = "Save QC plots for GWAS-gene chromVAR support and model results. [checkpoint:genetic_enrichment]",
    command = plot_psbulk_GWAS_gene_chromVAR_QC(
      psbulk_GWAS_gene_chromVAR_results_tibble = results_tibble.gene_level.pseudobulk,
      GWAS_gene_chromVAR_feature_metadata_tibble = feature_metadata_tibble.gene_level.pseudobulk
    ) |>
      save_plots_structured()
  ),
  targets::tar_target(
    name = locus_check_tibble.gene_level.pseudobulk,
    description = "Select top significant GWAS-gene chromVAR effects for locus sanity-check tracks",
    command = get_top_GWAS_gene_chromVAR_locus_check_tibble(
      psbulk_GWAS_gene_chromVAR_results_tibble = results_tibble.gene_level.pseudobulk,
      GWAS_gene_chromVAR_credible_set_variants_tibble = credible_set_variants_tibble.gene_level.pseudobulk,
      GWAS_gene_chromVAR_L2G_tibble = L2G_tibble.gene_level.pseudobulk,
      open_targets_target_dataset_path = open_targets_target_dataset_path,
      n_top = genetic_enrichment_psbulk_GWAS_gene_chromVAR_locus_plots_n_top,
      FDR_threshold = genetic_enrichment_psbulk_GWAS_gene_chromVAR_locus_plots_FDR_threshold,
      flank = genetic_enrichment_psbulk_GWAS_gene_chromVAR_locus_plots_flank
    )
  ),
  targets::tar_target(
    name = locus_tracks_plots.gene_level.pseudobulk,
    description = "Render locus sanity-check tracks for top significant GWAS-gene chromVAR effects",
    command = plot_GWAS_gene_chromVAR_locus_tracks(
      locus_check_tibble = locus_check_tibble.gene_level.pseudobulk,
      GWAS_gene_chromVAR_credible_set_variants_tibble = credible_set_variants_tibble.gene_level.pseudobulk,
      GWAS_gene_chromVAR_L2G_tibble = L2G_tibble.gene_level.pseudobulk,
      consensus_peak_GRanges = consensus_peak_GRanges.ATAC,
      fragments = combined_BPCells_fragment_obj.ATAC,
      metadata_tibble = metadata_w_cell_types_tibble.WNN
    ),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  tarchetypes::tar_file(
    name = locus_tracks_plot_files.gene_level.pseudobulk,
    description = "Save locus sanity-check tracks for top significant GWAS-gene chromVAR effects. [checkpoint:genetic_enrichment]",
    command = locus_tracks_plots.gene_level.pseudobulk |>
      save_plots_structured(
        override_suffix = "top_locus_checks",
        width = 12,
        height = 14
      )
  )
)
