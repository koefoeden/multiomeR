rlang::list2(
  targets::tar_target(
    name = GWAS_chromVAR_peak_weight_matrix,
    description = "Combine per-GWAS peak posterior-probability weights into one peak-by-GWAS annotation matrix",
    command = get_GWAS_chromVAR_peak_weight_matrix(
      peak_weight_records = GWAS_chromVAR_peak_weight_records,
      RSE_ATAC = chromVAR_obj.ATAC
    ),
    resources = get_tar_resources(RAM_GB_req = 40)
  ),
  targets::tar_target(
    name = psbulk_GWAS_chromVAR_activity_matrix,
    description = "Compute pseudobulk GWAS chromVAR activity scores from BPCells-backed pseudobulk ATAC counts",
    command = get_pseudobulk_chromVAR_activity_matrix(
      psbulk_ATAC_data_matrix = pseudobulk_counts_matrix.ATAC,
      chromVAR_obj = chromVAR_obj.ATAC,
      annotation_matrix = GWAS_chromVAR_peak_weight_matrix,
      normalize = TRUE
    ),
    resources = get_tar_resources(RAM_GB_req = 60)
  ),
  targets::tar_target(
    name = psbulk_GWAS_chromVAR_models,
    description = "Build named list of configured pseudobulk GWAS chromVAR model specifications",
    command = normalize_psbulk_feature_models(genetic_enrichment_psbulk_GWAS_chromVAR_models)
  ),
  targets::tar_target(
    name = psbulk_GWAS_chromVAR_dynamic_tibble,
    description = "Build a dynamic tibble with one row per pseudobulk GWAS chromVAR model",
    command = tibble::enframe(psbulk_GWAS_chromVAR_models, name = "model_name", value = "model") |>
      dplyr::mutate(cell_type_subset = purrr::map_chr(model, ~ get_psbulk_cell_type_subset_label(.x$cell_type_subset))) |>
      dplyr::arrange(model_name),
    iteration = "vector"
  ),
  targets::tar_target(
    name = psbulk_GWAS_chromVAR_filtered_mat_per_model,
    description = "Filter the pseudobulk GWAS chromVAR activity matrix for each model",
    command = filter_psbulk_data_matrix(
      psbulk_data_matrix = psbulk_GWAS_chromVAR_activity_matrix,
      psbulk_feature_dynamic_tibble = psbulk_GWAS_chromVAR_dynamic_tibble,
      extended_donor_id_metadata_tibble = genetic_enrichment_donor_id_metadata_tibble.extended,
      sample_depth_tibble = pseudobulk_depth_tibble.ATAC,
      min_sample_counts = genetic_enrichment_psbulk_GWAS_chromVAR_min_ATAC_counts
    ),
    pattern = map(psbulk_GWAS_chromVAR_dynamic_tibble),
    resources = get_tar_resources(cores_req = 1, RAM_GB_req = 60)
  ),
  targets::tar_target(
    name = psbulk_GWAS_chromVAR_fit,
    description = "Fit limma models to pseudobulk GWAS chromVAR activity matrices for each model",
    command = fit_psbulk_feature_matrix_model(
      psbulk_feature_matrix = psbulk_GWAS_chromVAR_filtered_mat_per_model,
      extended_donor_id_metadata_tibble = genetic_enrichment_donor_id_metadata_tibble.extended,
      psbulk_feature_dynamic_tibble = psbulk_GWAS_chromVAR_dynamic_tibble
    ),
    pattern = map(psbulk_GWAS_chromVAR_filtered_mat_per_model, psbulk_GWAS_chromVAR_dynamic_tibble)
  ),
  targets::tar_target(
    name = psbulk_GWAS_chromVAR_results_per_model_tibble,
    description = "Extract GWAS chromVAR differential activity results for each model",
    command = get_psbulk_feature_model_results(
      psbulk_feature_matrix_fit = psbulk_GWAS_chromVAR_fit,
      psbulk_feature_dynamic_tibble = psbulk_GWAS_chromVAR_dynamic_tibble
    ),
    pattern = map(psbulk_GWAS_chromVAR_fit, psbulk_GWAS_chromVAR_dynamic_tibble)
  ),
  targets::tar_target(
    name = psbulk_GWAS_chromVAR_results_tibble,
    description = "Combine GWAS chromVAR differential activity results and add GWAS-category FDR",
    command = psbulk_GWAS_chromVAR_results_per_model_tibble |>
      format_psbulk_GWAS_chromVAR_results(GWAS_inputs_tibble = GWAS_inputs_tibble)
  ),
  targets::tar_target(
    name = psbulk_GWAS_chromVAR_significant_elements_tibble,
    description = "Count FDR-significant pseudobulk GWAS chromVAR elements per model, contrast, and direction",
    command = get_psbulk_DX_significant_elements_tibble(
      psbulk_GWAS_chromVAR_results_tibble,
      FDR_col = "FDR_within_GWAS_category"
    )
  ),
  tarchetypes::tar_file(
    name = psbulk_GWAS_chromVAR_dotplot,
    description = "Save GWAS-by-contrast dot plot of pseudobulk GWAS chromVAR differential activity. [checkpoint:genetic_enrichment]",
    command = {
      score_plots <- psbulk_GWAS_chromVAR_results_tibble |>
        plot_psbulk_GWAS_chromVAR_dotplot(
          contrast_compartment_patterns = genetic_enrichment_compartment_patterns
        )

      if (!inherits(score_plots, "empty_plot_list")) {
        score_plots <- score_plots |>
          purrr::imap(\(score_plot, model_name) {
            patchwork::wrap_plots(
              GWAS_metadata_tracks_plot,
              score_plot,
              nrow = 1,
              widths = c(5.8, 9),
              guides = "collect"
            )
          })
      }

      score_plots |>
      save_plots_structured(
        width = 20,
        height = max(7, 0.24 * dplyr::n_distinct(psbulk_GWAS_chromVAR_results_tibble$GWAS_ID) + 3)
      )
    }
  ),
  tarchetypes::tar_file(
    name = psbulk_GWAS_chromVAR_QC_plots,
    description = "Save QC plots for pseudobulk GWAS chromVAR model support and P-value distributions. [checkpoint:genetic_enrichment]",
    command = psbulk_GWAS_chromVAR_results_tibble |>
      plot_psbulk_GWAS_chromVAR_QC() |>
      save_plots_structured()
  )
)
