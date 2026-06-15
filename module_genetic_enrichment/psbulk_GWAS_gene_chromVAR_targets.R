rlang::list2(
  targets::tar_target(
    name = GWAS_gene_chromVAR_inputs_tibble,
    description = "Restrict GWAS inputs to the configured gene-level pseudobulk GWAS chromVAR traits",
    command = get_GWAS_gene_chromVAR_inputs_tibble(
      GWAS_analysis_inputs_tibble = GWAS_analysis_inputs_tibble,
      selected_GWAS_IDs = genetic_enrichment_psbulk_GWAS_gene_chromVAR_GWAS_IDs
    ),
    iteration = "vector"
  ),
  targets::tar_target(
    name = GWAS_gene_chromVAR_credible_set_variants_tibble,
    description = "Expose Open Targets 95 percent credible-set variant rows for configured gene-level GWAS chromVAR traits",
    command = get_open_targets_credible_set_variants_tibble(
      GWAS_inputs_tibble = GWAS_gene_chromVAR_inputs_tibble,
      open_targets_credible_set_dataset_path = open_targets_credible_set_dataset_path
    ) |>
      dplyr::select(-is95CredibleSet, -is99CredibleSet),
    resources = get_tar_resources(RAM_GB_req = 40)
  ),
  targets::tar_target(
    name = GWAS_gene_chromVAR_L2G_tibble,
    description = "Read Open Targets L2G evidence for the gene-level GWAS chromVAR credible-set loci",
    command = get_GWAS_gene_chromVAR_L2G_tibble(
      GWAS_gene_chromVAR_credible_set_variants_tibble = GWAS_gene_chromVAR_credible_set_variants_tibble,
      open_targets_gwas_credible_sets_evidence_dataset_path = open_targets_gwas_credible_sets_evidence_dataset_path
    ),
    resources = get_tar_resources(RAM_GB_req = 40)
  ),
  targets::tar_target(
    name = GWAS_gene_chromVAR_peak_weight_records,
    description = "Build sparse peak-by-target L2G-weighted GWAS chromVAR annotations for one configured GWAS",
    command = get_GWAS_gene_chromVAR_peak_weight_record(
      GWAS_input_tibble = GWAS_gene_chromVAR_inputs_tibble,
      GWAS_gene_chromVAR_credible_set_variants_tibble = GWAS_gene_chromVAR_credible_set_variants_tibble,
      GWAS_gene_chromVAR_L2G_tibble = GWAS_gene_chromVAR_L2G_tibble,
      peak_ranges = genetic_enrichment_peak_ranges,
      posterior_probability_cutoff = genetic_enrichment_posterior_probability_cutoff,
      posterior_probability_weighting_function = posterior_probability_weighting_function
    ),
    pattern = map(GWAS_gene_chromVAR_inputs_tibble),
    iteration = "list",
    resources = get_tar_resources(RAM_GB_req = 60)
  ),
  targets::tar_target(
    name = GWAS_gene_chromVAR_feature_metadata_tibble,
    description = "Summarize support metrics for each GWAS-gene chromVAR annotation feature",
    command = get_GWAS_gene_chromVAR_feature_metadata_tibble(
      GWAS_gene_chromVAR_peak_weight_records = GWAS_gene_chromVAR_peak_weight_records,
      min_n_peaks = 3
    ) |>
      add_open_targets_target_metadata(
        open_targets_target_dataset_path = open_targets_target_dataset_path
      )
  ),
  targets::tar_target(
    name = GWAS_gene_chromVAR_peak_weight_matrix,
    description = "Combine supported peak-by-GWAS-gene annotations into one sparse peak-by-feature matrix",
    command = get_GWAS_gene_chromVAR_peak_weight_matrix(
      GWAS_gene_chromVAR_peak_weight_records = GWAS_gene_chromVAR_peak_weight_records,
      GWAS_gene_chromVAR_feature_metadata_tibble = GWAS_gene_chromVAR_feature_metadata_tibble,
      RSE_ATAC = chromVAR_obj.ATAC
    ),
    resources = get_tar_resources(RAM_GB_req = 60)
  ),
  targets::tar_target(
    name = psbulk_GWAS_gene_chromVAR_activity_matrix,
    description = "Compute pseudobulk GWAS-gene chromVAR activity scores from BPCells-backed pseudobulk ATAC counts",
    command = if (ncol(GWAS_gene_chromVAR_peak_weight_matrix) == 0) {
      matrix(
        nrow = 0,
        ncol = ncol(pseudobulk_counts_matrix.ATAC),
        dimnames = list(character(), colnames(pseudobulk_counts_matrix.ATAC))
      )
    } else {
      get_pseudobulk_chromVAR_activity_matrix(
        psbulk_ATAC_data_matrix = pseudobulk_counts_matrix.ATAC,
        chromVAR_obj = chromVAR_obj.ATAC,
        annotation_matrix = GWAS_gene_chromVAR_peak_weight_matrix,
        normalize = FALSE
      )
    },
    resources = get_tar_resources(RAM_GB_req = 80)
  ),
  targets::tar_target(
    name = psbulk_GWAS_gene_chromVAR_dynamic_tibble,
    description = "Build a dynamic tibble with one row per pseudobulk GWAS-gene chromVAR model",
    command = if (nrow(psbulk_GWAS_gene_chromVAR_activity_matrix) == 0) {
      psbulk_GWAS_chromVAR_dynamic_tibble[0, ]
    } else {
      psbulk_GWAS_chromVAR_dynamic_tibble
    },
    iteration = "vector"
  ),
  targets::tar_target(
    name = psbulk_GWAS_gene_chromVAR_filtered_mat_per_model,
    description = "Filter the pseudobulk GWAS-gene chromVAR activity matrix for each model",
    command = filter_psbulk_data_matrix(
      psbulk_data_matrix = psbulk_GWAS_gene_chromVAR_activity_matrix,
      psbulk_feature_dynamic_tibble = psbulk_GWAS_gene_chromVAR_dynamic_tibble,
      extended_donor_id_metadata_tibble = genetic_enrichment_donor_id_metadata_tibble.extended,
      sample_depth_tibble = pseudobulk_depth_tibble.ATAC,
      min_sample_counts = genetic_enrichment_psbulk_GWAS_chromVAR_min_ATAC_counts
    ),
    pattern = map(psbulk_GWAS_gene_chromVAR_dynamic_tibble),
    resources = get_tar_resources(cores_req = 1, RAM_GB_req = 60)
  ),
  targets::tar_target(
    name = psbulk_GWAS_gene_chromVAR_fit,
    description = "Fit limma models to pseudobulk GWAS-gene chromVAR activity matrices for each model",
    command = fit_psbulk_feature_matrix_model(
      psbulk_feature_matrix = psbulk_GWAS_gene_chromVAR_filtered_mat_per_model,
      extended_donor_id_metadata_tibble = genetic_enrichment_donor_id_metadata_tibble.extended,
      psbulk_feature_dynamic_tibble = psbulk_GWAS_gene_chromVAR_dynamic_tibble
    ),
    pattern = map(psbulk_GWAS_gene_chromVAR_filtered_mat_per_model, psbulk_GWAS_gene_chromVAR_dynamic_tibble),
    resources = get_tar_resources(RAM_GB_req = 80)
  ),
  targets::tar_target(
    name = psbulk_GWAS_gene_chromVAR_results_per_model_tibble,
    description = "Extract GWAS-gene chromVAR differential activity results for each model",
    command = get_psbulk_feature_model_results(
      psbulk_feature_matrix_fit = psbulk_GWAS_gene_chromVAR_fit,
      psbulk_feature_dynamic_tibble = psbulk_GWAS_gene_chromVAR_dynamic_tibble
    ),
    pattern = map(psbulk_GWAS_gene_chromVAR_fit, psbulk_GWAS_gene_chromVAR_dynamic_tibble)
  ),
  targets::tar_target(
    name = psbulk_GWAS_gene_chromVAR_results_tibble,
    description = "Combine GWAS-gene chromVAR differential activity results and add within-GWAS FDR",
    command = psbulk_GWAS_gene_chromVAR_results_per_model_tibble |>
      format_psbulk_GWAS_gene_chromVAR_results(
        GWAS_gene_chromVAR_feature_metadata_tibble = GWAS_gene_chromVAR_feature_metadata_tibble
      )
  ),
  targets::tar_target(
    name = psbulk_GWAS_gene_chromVAR_volcano_plots,
    description = "Plot volcano plots of pseudobulk GWAS-gene chromVAR results per contrast",
    command = plot_psbulk_GWAS_gene_chromVAR_volcanoes(
      psbulk_GWAS_gene_chromVAR_results_tibble = psbulk_GWAS_gene_chromVAR_results_tibble,
      model_name = psbulk_GWAS_gene_chromVAR_dynamic_tibble$model_name
    ),
    pattern = map(psbulk_GWAS_gene_chromVAR_dynamic_tibble)
  ),
  tarchetypes::tar_file(
    name = psbulk_GWAS_gene_chromVAR_volcano_plots_files,
    description = "Save GWAS-gene chromVAR volcano plots per model and contrast to file. [checkpoint:genetic_enrichment]",
    command = save_plots_structured(
      plots = psbulk_GWAS_gene_chromVAR_volcano_plots,
      override_suffix = psbulk_GWAS_gene_chromVAR_dynamic_tibble$model_name,
      dyn_suffix_in_subdir = TRUE
    ),
    pattern = map(psbulk_GWAS_gene_chromVAR_volcano_plots, psbulk_GWAS_gene_chromVAR_dynamic_tibble)
  ),
  tarchetypes::tar_file(
    name = psbulk_GWAS_gene_chromVAR_QC_plots,
    description = "Save QC plots for GWAS-gene chromVAR support and model results. [checkpoint:genetic_enrichment]",
    command = plot_psbulk_GWAS_gene_chromVAR_QC(
      psbulk_GWAS_gene_chromVAR_results_tibble = psbulk_GWAS_gene_chromVAR_results_tibble,
      GWAS_gene_chromVAR_feature_metadata_tibble = GWAS_gene_chromVAR_feature_metadata_tibble
    ) |>
      save_plots_structured()
  ),
  targets::tar_target(
    name = psbulk_GWAS_gene_chromVAR_KCNJ12_ATAC_accessibility_track,
    description = "Build the hardcoded KCNJ12 BPCells ATAC coverage track for locus-effect sanity-check plotting",
    command = {
      hardcoded_GWAS_ID <- "T2D_Suzuki2024_locus_effect"
      hardcoded_target_id <- "ENSG00000184185"
      hardcoded_contrast <- "fibers_vs_non_myogenic"

      locus_check_tibble <- get_GWAS_gene_chromVAR_locus_check_tibble(
        psbulk_GWAS_gene_chromVAR_results_tibble = psbulk_GWAS_gene_chromVAR_results_tibble,
        GWAS_gene_chromVAR_credible_set_variants_tibble = GWAS_gene_chromVAR_credible_set_variants_tibble,
        GWAS_gene_chromVAR_L2G_tibble = GWAS_gene_chromVAR_L2G_tibble,
        open_targets_target_dataset_path = open_targets_target_dataset_path,
        GWAS_ID = hardcoded_GWAS_ID,
        target_id = hardcoded_target_id,
        contrast = hardcoded_contrast
      )
      region <- GenomicRanges::GRanges(
        seqnames = locus_check_tibble$chr[[1]],
        ranges = IRanges::IRanges(
          start = locus_check_tibble$locus_start[[1]],
          end = locus_check_tibble$locus_end[[1]]
        )
      )

      metadata_tibble <- metadata_w_cell_types_tibble.WNN
      fragments <- BPCells::select_cells(
        combined_BPCells_fragment_obj.ATAC,
        metadata_tibble$barcode_w_prefix
      )
      fragment_cell_names <- BPCells::cellNames(fragments)
      metadata <- metadata_tibble |>
        dplyr::distinct(.data$barcode_w_prefix, .keep_all = TRUE) |>
        dplyr::filter(
          .data$barcode_w_prefix %in% fragment_cell_names,
          !is.na(.data[["PCA_harmony_SNN_cluster_cell_type"]])
        ) |>
        dplyr::arrange(match(.data$barcode_w_prefix, fragment_cell_names))
      fragments <- BPCells::select_cells(fragments, metadata$barcode_w_prefix)

      groups <- metadata[["PCA_harmony_SNN_cluster_cell_type"]]
      cell_read_counts <- if ("atac_fragments" %in% colnames(metadata)) {
        metadata$atac_fragments
      } else {
        metadata$nCount_ATAC
      }
      coverage_tibble <- BPCells::trackplot_coverage(
        fragments = fragments,
        region = region,
        groups = groups,
        cell_read_counts = cell_read_counts,
        group_order = gtools::mixedsort(unique(as.character(groups))),
        bins = 500,
        return_data = TRUE
      )

      make_BPCells_ATAC_coverage_track_from_tibble(coverage_tibble, region)
    },
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  tarchetypes::tar_file(
    name = psbulk_GWAS_gene_chromVAR_KCNJ12_locus_tracks_plot,
    description = "Save a hardcoded KCNJ12 locus-effect sanity-check track plot for the T2D Suzuki 2024 fiber contrast. [checkpoint:genetic_enrichment]",
    command = {
      # First-pass sanity check plot. If this remains useful, generalize this
      # block to selected volcano hits instead of keeping one hardcoded gene.
      hardcoded_GWAS_ID <- "T2D_Suzuki2024_locus_effect"
      hardcoded_target_id <- "ENSG00000184185"
      hardcoded_contrast <- "fibers_vs_non_myogenic"

      locus_check_tibble <- get_GWAS_gene_chromVAR_locus_check_tibble(
        psbulk_GWAS_gene_chromVAR_results_tibble = psbulk_GWAS_gene_chromVAR_results_tibble,
        GWAS_gene_chromVAR_credible_set_variants_tibble = GWAS_gene_chromVAR_credible_set_variants_tibble,
        GWAS_gene_chromVAR_L2G_tibble = GWAS_gene_chromVAR_L2G_tibble,
        open_targets_target_dataset_path = open_targets_target_dataset_path,
        GWAS_ID = hardcoded_GWAS_ID,
        target_id = hardcoded_target_id,
        contrast = hardcoded_contrast
      )
      region <- GenomicRanges::GRanges(
        seqnames = locus_check_tibble$chr[[1]],
        ranges = IRanges::IRanges(
          start = locus_check_tibble$locus_start[[1]],
          end = locus_check_tibble$locus_end[[1]]
        )
      )

      L2G_tibble <- GWAS_gene_chromVAR_L2G_tibble |>
        dplyr::filter(
          .data$GWAS_ID == hardcoded_GWAS_ID,
          .data$targetId == hardcoded_target_id
        )
      variant_tibble <- GWAS_gene_chromVAR_credible_set_variants_tibble |>
        dplyr::filter(.data$GWAS_ID == hardcoded_GWAS_ID) |>
        dplyr::semi_join(
          L2G_tibble |> dplyr::select(studyLocusId),
          by = "studyLocusId"
        )

      plot <- BPCells::trackplot_combine(
        tracks = c(
          list(
            make_open_targets_target_locus_track(locus_check_tibble, region),
            psbulk_GWAS_gene_chromVAR_KCNJ12_ATAC_accessibility_track,
            make_consensus_peak_locus_track(consensus_peak_GRanges.ATAC, region)
          ),
          make_open_targets_variant_PIP_tracks(variant_tibble, region)
        ),
        title = stringr::str_glue(
          "{hardcoded_GWAS_ID} {locus_check_tibble$gene_symbol}: {hardcoded_contrast}; ",
          "logFC={round(locus_check_tibble$logFC, 3)}, ",
          "FDR={format(locus_check_tibble$FDR, scientific = TRUE, digits = 3)}, ",
          "L2G={round(locus_check_tibble$L2G_score, 3)}"
        )
      )

      plot |>
        save_plots_structured(
          override_suffix = "KCNJ12_T2D_Suzuki2024_locus_effect_fibers_vs_non_myogenic",
          width = 12,
          height = 14
        )
    }
  ),
)
