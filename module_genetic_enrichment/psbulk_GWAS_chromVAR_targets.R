rlang::list2(
  targets::tar_target(
    name = peak_weight_matrix.trait_level.pseudobulk,
    description = "Combine per-GWAS peak posterior-probability weights into one peak-by-GWAS annotation matrix",
    command = get_GWAS_chromVAR_peak_weight_matrix(
      peak_weight_records = peak_weight_records.trait_level,
      RSE_ATAC = chromVAR_obj.ATAC
    ),
    resources = get_tar_resources(RAM_GB_req = 40)
  ),
  targets::tar_target(
    name = activity_matrix.trait_level.pseudobulk,
    description = "Compute trait-level pseudobulk chromVAR activity scores from BPCells-backed pseudobulk ATAC counts",
    command = get_pseudobulk_chromVAR_activity_matrix(
      psbulk_ATAC_data_matrix = pseudobulk_counts_matrix.ATAC,
      chromVAR_obj = chromVAR_obj.ATAC,
      annotation_matrix = peak_weight_matrix.trait_level.pseudobulk,
      normalize = FALSE
    ),
    resources = get_tar_resources(RAM_GB_req = 60)
  ),
  targets::tar_target(
    name = models.trait_level.pseudobulk,
    description = "Build named list of configured trait-level pseudobulk chromVAR model specifications",
    command = normalize_psbulk_feature_models(genetic_enrichment_psbulk_GWAS_chromVAR_models)
  ),
  targets::tar_target(
    name = model_tibble.trait_level.pseudobulk,
    description = "Build a dynamic tibble with one row per trait-level pseudobulk chromVAR model",
    command = tibble::enframe(models.trait_level.pseudobulk, name = "model_name", value = "model") |>
      dplyr::mutate(cell_type_subset = purrr::map_chr(model, ~ get_psbulk_cell_type_subset_label(.x$cell_type_subset))) |>
      dplyr::arrange(model_name),
    iteration = "vector"
  ),
  targets::tar_target(
    name = filtered_matrix.trait_level.pseudobulk,
    description = "Filter the trait-level pseudobulk chromVAR activity matrix for each model",
    command = filter_psbulk_data_matrix(
      psbulk_data_matrix = activity_matrix.trait_level.pseudobulk,
      psbulk_feature_dynamic_tibble = model_tibble.trait_level.pseudobulk,
      extended_donor_id_metadata_tibble = genetic_enrichment_donor_id_metadata_tibble.extended,
      sample_depth_tibble = pseudobulk_depth_tibble.ATAC,
      min_sample_counts = genetic_enrichment_psbulk_GWAS_chromVAR_min_ATAC_counts
    ),
    pattern = map(model_tibble.trait_level.pseudobulk),
    resources = get_tar_resources(cores_req = 1, RAM_GB_req = 60)
  ),
  targets::tar_target(
    name = fit.trait_level.pseudobulk,
    description = "Fit limma models to trait-level pseudobulk chromVAR activity matrices for each model",
    command = fit_psbulk_feature_matrix_model(
      psbulk_feature_matrix = filtered_matrix.trait_level.pseudobulk,
      extended_donor_id_metadata_tibble = genetic_enrichment_donor_id_metadata_tibble.extended,
      psbulk_feature_dynamic_tibble = model_tibble.trait_level.pseudobulk
    ),
    pattern = map(filtered_matrix.trait_level.pseudobulk, model_tibble.trait_level.pseudobulk)
  ),
  targets::tar_target(
    name = results_per_model_tibble.trait_level.pseudobulk,
    description = "Extract trait-level pseudobulk chromVAR differential activity results for each model",
    command = get_psbulk_feature_model_results(
      psbulk_feature_matrix_fit = fit.trait_level.pseudobulk,
      psbulk_feature_dynamic_tibble = model_tibble.trait_level.pseudobulk
    ),
    pattern = map(fit.trait_level.pseudobulk, model_tibble.trait_level.pseudobulk)
  ),
  targets::tar_target(
    name = results_tibble.trait_level.pseudobulk,
    description = "Combine trait-level pseudobulk chromVAR differential activity results and add GWAS-category FDR",
    command = results_per_model_tibble.trait_level.pseudobulk |>
      format_psbulk_GWAS_chromVAR_results(GWAS_inputs_tibble = GWAS_inputs_tibble)
  ),
  targets::tar_target(
    name = significant_elements_tibble.trait_level.pseudobulk,
    description = "Count FDR-significant trait-level pseudobulk chromVAR elements per model, contrast, and direction",
    command = get_psbulk_DX_significant_elements_tibble(
      results_tibble.trait_level.pseudobulk,
      FDR_col = "FDR_within_GWAS_category"
    )
  ),
  tarchetypes::tar_file(
    name = dotplot.trait_level.pseudobulk,
    description = "Save GWAS-by-contrast dot plot of trait-level pseudobulk chromVAR differential activity. [checkpoint:genetic_enrichment]",
    command = {
      score_plots <- results_tibble.trait_level.pseudobulk |>
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
        height = max(7, 0.24 * dplyr::n_distinct(results_tibble.trait_level.pseudobulk$GWAS_ID) + 3)
      )
    }
  ),
  tarchetypes::tar_file(
    name = QC_plots.trait_level.pseudobulk,
    description = "Save QC plots for trait-level pseudobulk chromVAR model support and P-value distributions. [checkpoint:genetic_enrichment]",
    command = results_tibble.trait_level.pseudobulk |>
      plot_psbulk_GWAS_chromVAR_QC() |>
      save_plots_structured()
  )
)
