rlang::list2(
  targets::tar_target(
    name = dynamic_tibble,
    description = "Build a dynamic tibble with one row per pseudobulk differential model",
    command = tibble::enframe(models, name = "model_name", value = "model") |>
      dplyr::mutate(cell_type_subset = purrr::map_chr(model, ~ get_pseudobulk_cell_type_subset_label(.x$cell_type_subset))) |>
      dplyr::arrange(model_name),
    iteration = "vector"
  ),
  targets::tar_target(
    name = filtered_mat_per_model,
    description = "Filter the pseudobulk data matrix for each model",
    command = filter_pseudobulk_data_matrix(
      pseudobulk_data_matrix = map_pseudobulk_data_matrix,
      pseudobulk_feature_dynamic_tibble = dynamic_tibble,
      extended_donor_id_metadata_tibble = donor_id_metadata_tibble.analysis,
      sample_depth_tibble = if (identical(map_analysis_suffix, "motif_family_accessibility")) pseudobulk_depth_tibble.ATAC else NULL,
      min_sample_counts = if (identical(map_analysis_suffix, "motif_family_accessibility")) differential_analyses_motif_family_accessibility_min_ATAC_counts else NULL
    ),
    pattern = map(dynamic_tibble),
    resources = get_tar_resources(RAM_GB_req = 60)
  ),
  tarchetypes::tar_file(
    name = cohort_tsv,
    description = "Export model-specific donor eligibility, sample counts, contributing wells and exclusions. [checkpoint:differential_analyses]",
    command = get_feature_model_cohort(map_pseudobulk_data_matrix, filtered_mat_per_model, feature_matrix_fit,
      donor_id_metadata_tibble.analysis, dynamic_tibble$model[[1]], metadata_w_cell_types_tibble.WNN) |>
      save_differential_model_table(dynamic_tibble$model_name),
    pattern = map(filtered_mat_per_model, feature_matrix_fit, dynamic_tibble)
  ),
  targets::tar_target(
    name = feature_matrix_fit,
    description = "Fit edgeR/limma models to pseudobulk feature matrices for each model [part_of_graph:differential_analyses]",
    command = fit_pseudobulk_feature_matrix_model(
      pseudobulk_feature_matrix = filtered_mat_per_model,
      extended_donor_id_metadata_tibble = donor_id_metadata_tibble.analysis,
      pseudobulk_feature_dynamic_tibble = dynamic_tibble
    ),
    pattern = map(filtered_mat_per_model, dynamic_tibble),
    resources = get_tar_resources(
      cores_req = 6,
      RAM_GB_req = 60
    )
  ),
  targets::tar_target(
    name = results_tibble,
    description = "Extract differential testing results tibble from each fitted model [part_of_graph:differential_analyses]",
    command = get_pseudobulk_feature_model_results(
      pseudobulk_feature_matrix_fit = feature_matrix_fit,
      pseudobulk_feature_dynamic_tibble = dynamic_tibble
    ),
    pattern = map(feature_matrix_fit, dynamic_tibble)
  ),
  targets::tar_target(
    name = top_features_tibble,
    description = "Collect the top 40 most significant features per model and contrast",
    command = results_tibble |>
      tibble::as_tibble() |>
      dplyr::arrange(PValue) |>
      dplyr::group_by(model, contrast) |>
      dplyr::slice_head(n = 40) |>
      dplyr::ungroup(),
    pattern = map(results_tibble)
  ),
  targets::tar_target(
    name = significant_elements_tibble,
    description = "Count FDR-significant pseudobulk differential elements per model, contrast, and direction",
    command = get_pseudobulk_differential_significant_elements_tibble(results_tibble)
  ),
  tarchetypes::tar_file(
    name = significant_elements_plot,
    description = "Plot signed counts of FDR-significant pseudobulk differential elements per model and contrast. [checkpoint:differential_analyses]",
    command = significant_elements_tibble |>
      plot_pseudobulk_differential_significant_elements(modality = switch(map_analysis_suffix,
        gene_expression = "gene expression", chromatin_accessibility = "chromatin accessibility", motif_family_accessibility = "motif-family accessibility",
        transcription_factor_activity = "TF activity")) |>
      save_plots_structured()
  ),
  targets::tar_target(
    name = top_feature_open_targets_evidence_tibble,
    description = "Optionally annotate top pseudobulk differential features with Open Targets GWAS evidence",
    command = {
      top_features_vec <- top_features_tibble |>
        dplyr::pull(feature_id) |>
        unique()

      if (!identical(map_analysis_suffix, "gene_expression")) {
        return(tibble::tibble(feature_id = top_features_vec, ensembl_id = NA_character_, OT_GWAS_evidence = NA))
      }

      get_OT_GWAS_gene_evidence_tibble(
        gene_symbols = top_features_vec,
        gene_features_df = gene_features_df,
        efo_id = differential_analyses_pseudobulk_OT_GWAS_efo_id
      ) |>
        dplyr::rename(feature_id = gene)
    },
    pattern = map(top_features_tibble)
  ),
  tarchetypes::tar_file(
    name = volcano_plots,
    description = "Plot volcano plots of pseudobulk differential results per contrast. [checkpoint:differential_analyses]",
    command = {
      if (stringr::str_detect(map_analysis_suffix, "chromatin_accessibility")) {
        plots_list <- plot_pseudobulk_chromatin_accessibility_volcanoes(
          pseudobulk_differential_results_tibble = results_tibble,
          ATAC_consensus_peak_GRanges = consensus_peak_annotated_GRanges.ATAC
        )
      } else {
        plots_list <- plot_pseudobulk_feature_volcanoes(
          pseudobulk_differential_results_tibble = results_tibble,
          pseudobulk_differential_top_features_tibble = top_features_tibble,
          pseudobulk_differential_top_feature_open_targets_evidence_tibble = top_feature_open_targets_evidence_tibble,
          feature_labels = if (identical(map_analysis_suffix, "motif_family_accessibility")) JASPAR_motif_family_labels else character()
        )
      }
      save_plots_structured(
        plots = plots_list,
        override_suffix = dynamic_tibble$model_name,
        dyn_suffix_in_subdir = TRUE
      )
    },
    pattern = map(results_tibble, top_features_tibble, dynamic_tibble, top_feature_open_targets_evidence_tibble)
  ),
  tarchetypes::tar_file(
    name = p_value_distribution_plot,
    description = "Save P-value density plots across all pseudobulk differential models. [checkpoint:differential_analyses]",
    command = results_tibble |>
      dplyr::bind_rows() |>
      plot_pseudobulk_p_value_distribution() |>
      save_plots_structured()
  )
)
