rlang::list2(
  targets::tar_target(
    name = dynamic_tibble,
    description = "Build a dynamic tibble with one row per pseudobulk DX model",
    command = tibble::enframe(models, name = "model_name", value = "model") |>
      dplyr::mutate(cell_type_subset = purrr::map_chr(model, ~ get_psbulk_cell_type_subset_label(.x$cell_type_subset))) |>
      dplyr::arrange(model_name),
    iteration = "vector"
  ),
  targets::tar_target(
    name = filtered_mat_per_model,
    description = "Filter the pseudobulk data matrix for each model",
    command = filter_psbulk_data_matrix(
      psbulk_data_matrix = map_psbulk_data_matrix,
      psbulk_feature_dynamic_tibble = dynamic_tibble,
      extended_donor_id_metadata_tibble = donor_id_metadata_tibble.extended,
      sample_depth_tibble = if (identical(map_psbulk_DX_tar_suffix, "DTFA")) pseudobulk_depth_tibble.ATAC else NULL,
      min_sample_counts = if (identical(map_psbulk_DX_tar_suffix, "DTFA")) differential_analyses_DTFA_min_ATAC_counts else NULL
    ),
    pattern = map(dynamic_tibble),
    resources = get_tar_resources(RAM_GB_req = 60)
  ),
  targets::tar_target(
    name = feature_matrix_fit,
    description = "Fit edgeR/limma models to pseudobulk feature matrices for each model [part_of_graph:differential_analyses]",
    command = fit_psbulk_feature_matrix_model(
      psbulk_feature_matrix = filtered_mat_per_model,
      extended_donor_id_metadata_tibble = donor_id_metadata_tibble.extended,
      psbulk_feature_dynamic_tibble = dynamic_tibble
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
    command = get_psbulk_feature_model_results(
      psbulk_feature_matrix_fit = feature_matrix_fit,
      psbulk_feature_dynamic_tibble = dynamic_tibble
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
    description = "Count FDR-significant pseudobulk DX elements per model, contrast, and direction",
    command = get_psbulk_DX_significant_elements_tibble(results_tibble)
  ),
  tarchetypes::tar_file(
    name = significant_elements_plot,
    description = "Plot signed counts of FDR-significant pseudobulk DX elements per model and contrast. [checkpoint:differential_analyses]",
    command = significant_elements_tibble |>
      plot_psbulk_DX_significant_elements() |>
      save_plots_structured()
  ),
  targets::tar_target(
    name = top_feature_OT_evidence_tibble,
    description = "Optionally annotate top pseudobulk DX features with Open Targets GWAS evidence",
    command = {
      top_features_vec <- top_features_tibble |>
        dplyr::pull(feature_id) |>
        unique()

      if (stringr::str_detect(map_psbulk_DX_tar_suffix, "DCA")) {
        return(tibble::tibble(feature_id = top_features_vec, ensembl_id = NA_character_, OT_GWAS_evidence = NA))
      }

      get_OT_GWAS_gene_evidence_tibble(
        gene_symbols = top_features_vec,
        gene_features_df = gene_features_df,
        efo_id = differential_analyses_psbulk_DX_OT_GWAS_efo_id
      ) |>
        dplyr::rename(feature_id = gene)
    },
    pattern = map(top_features_tibble)
  ),
  tarchetypes::tar_file(
    name = volcano_plots_files,
    description = "Plot volcano plots of pseudobulk DX results per contrast. [checkpoint:differential_analyses]",
    command = {
      if (stringr::str_detect(map_psbulk_DX_tar_suffix, "DCA")) {
        plots_list <- plot_psbulk_DCA_volcanoes(
          psbulk_DX_results_tibble = results_tibble,
          ATAC_consensus_peak_GRanges = consensus_peak_annotated_GRanges.ATAC
        )
      } else {
        plots_list <- plot_psbulk_DGE_volcanoes(
          psbulk_DX_results_tibble = results_tibble,
          psbulk_DX_top_features_tibble = top_features_tibble,
          psbulk_DX_top_feature_OT_evidence_tibble = top_feature_OT_evidence_tibble
        )
      }
      save_plots_structured(
        plots = plots_list,
        override_suffix = dynamic_tibble$model_name,
        dyn_suffix_in_subdir = TRUE
      )
    },
    pattern = map(results_tibble, top_features_tibble, dynamic_tibble, top_feature_OT_evidence_tibble)
  ),
  # tarchetypes::tar_file(
  #   name = volcano_plots_files,
  #   description = "Save volcano plots of pseudobulk DX results per model and contrast to file",
  #   command = save_plots_structured(
  #     plots = volcano_plots,
  #     override_suffix = dynamic_tibble$model_name,
  #     dyn_suffix_in_subdir = TRUE
  #   ),
  #   pattern = map(volcano_plots, dynamic_tibble)
  # ),
  targets::tar_target(
    name = PValue_density_plot,
    description = "Plot P-value density histograms across all pseudobulk DX models",
    command = results_tibble |>
      dplyr::bind_rows() |>
      plot_psbulk_DX_PValue_density()
  ),
  tarchetypes::tar_file(
    name = PValue_density_plot_file,
    description = "Save P-value density histograms across all pseudobulk DX models to file. [checkpoint:differential_analyses]",
    command = save_plots_structured(PValue_density_plot)
  )
)
