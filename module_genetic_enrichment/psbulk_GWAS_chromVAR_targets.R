rlang::list2(
  targets::tar_target(
    name = peak_weight_matrix.trait_level.pseudobulk,
    description = "Combine per-GWAS peak posterior-probability weights into one peak-by-GWAS annotation matrix [part_of_graph:genetic_enrichment_pseudobulk]",
    command = get_GWAS_chromVAR_peak_weight_matrix(
      peak_weight_records = peak_weight_records.trait_level,
      RSE_ATAC = chromVAR_obj.ATAC
    ),
    resources = get_tar_resources(RAM_GB_req = 40)
  ),
  tarchetypes::tar_file(
    name = cell_type_pseudobulk_counts_BPCells_matrix_dir.ATAC,
    description = "Write ATAC counts summed per WNN cell type to BPCells for descriptive GWAS chromVAR plots [part_of_graph:genetic_enrichment_pseudobulk]",
    command = {
      out_dir <- get_structured_file_path()
      get_BPCells_group_pseudobulk_matrix(
        feature_matrix = peak_QC_filtered_BPCells_matrix.ATAC,
        metadata_tibble = metadata_w_cell_types_tibble.WNN,
        group_col = "PCA_harmony_SNN_cluster_cell_type",
        threads = 6
      ) |>
        methods::as("dgCMatrix") |>
        BPCells::write_matrix_dir(out_dir, overwrite = TRUE)
      out_dir
    },
    resources = get_tar_resources(cores_req = 6, RAM_GB_req = 60)
  ),
  targets::tar_target(
    name = cell_type_pseudobulk_counts_matrix.ATAC,
    description = "Open BPCells-backed ATAC counts summed per WNN cell type",
    command = BPCells::open_matrix_dir(cell_type_pseudobulk_counts_BPCells_matrix_dir.ATAC)
  ),
  targets::tar_target(
    name = cell_type_pseudobulk_support_tibble.ATAC,
    description = "Compute cell counts and ATAC depth support per WNN cell type",
    command = get_group_pseudobulk_depth_tibble(cell_type_pseudobulk_counts_matrix.ATAC) |>
      dplyr::left_join(
        get_group_cell_count_tibble(
          metadata_tibble = metadata_w_cell_types_tibble.WNN,
          group_col = "PCA_harmony_SNN_cluster_cell_type"
        ),
        by = "cluster"
      ) |>
      dplyr::mutate(modality = "ATAC")
  ),
  targets::tar_target(
    name = activity_matrix.trait_level.pseudobulk,
    description = "Compute trait-level pseudobulk chromVAR activity scores from BPCells-backed pseudobulk ATAC counts [part_of_graph:genetic_enrichment_pseudobulk]",
    command = get_pseudobulk_chromVAR_activity_matrix(
      psbulk_ATAC_data_matrix = pseudobulk_counts_matrix.ATAC,
      chromVAR_obj = chromVAR_obj.ATAC,
      annotation_matrix = peak_weight_matrix.trait_level.pseudobulk,
      normalize = FALSE
    ),
    resources = get_tar_resources(RAM_GB_req = 60)
  ),
  targets::tar_target(
    name = chromVAR_deviation_record.trait_level.cell_type_pseudobulk,
    description = "Compute trait-level GWAS chromVAR deviations and retain their exact cell-type background model [part_of_graph:genetic_enrichment_pseudobulk]",
    command = get_pseudobulk_chromVAR_deviation_record(
      psbulk_ATAC_data_matrix = cell_type_pseudobulk_counts_matrix.ATAC,
      chromVAR_obj = chromVAR_obj.ATAC,
      annotation_matrix = peak_weight_matrix.trait_level.pseudobulk
    ),
    resources = get_tar_resources(RAM_GB_req = 60)
  ),
  targets::tar_target(
    name = chromVAR_deviation_tibble.trait_level.pseudobulk,
    description = "Format cell-type pseudobulk GWAS chromVAR deviations and z-score support labels for descriptive plots",
    command = chromVAR_deviation_record.trait_level.cell_type_pseudobulk$deviation_SE |>
      format_cell_type_GWAS_chromVAR_deviations(
        GWAS_inputs_tibble = GWAS_inputs_tibble,
        cell_type_support_tibble = cell_type_pseudobulk_support_tibble.ATAC
      )
  ),
  targets::tar_target(
    name = chromVAR_deviation_heatmap_data.trait_level.pseudobulk,
    description = "Prepare heatmap data from cell-type pseudobulk GWAS chromVAR relative deviations",
    command = chromVAR_deviation_tibble.trait_level.pseudobulk
  ),
  targets::tar_target(
    name = models.trait_level.pseudobulk,
    description = "Build named list of configured trait-level pseudobulk chromVAR model specifications [part_of_graph:differential_analyses] [part_of_graph:genetic_enrichment_pseudobulk]",
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
    resources = get_tar_resources(RAM_GB_req = 60)
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
    description = "Combine trait-level pseudobulk chromVAR differential activity results and add GWAS-category FDR [part_of_graph:genetic_enrichment_pseudobulk]",
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
    name = coefficient_ranges.trait_level.pseudobulk,
    description = "Save GWAS-by-contrast coefficient range plots of trait-level pseudobulk chromVAR differential activity. [checkpoint:genetic_enrichment]",
    command = {
      coefficient_plots <- results_tibble.trait_level.pseudobulk |>
        plot_psbulk_GWAS_chromVAR_coefficient_ranges(
          contrast_compartment_patterns = genetic_enrichment_compartment_patterns
        )

      if (!inherits(coefficient_plots, "empty_plot_list")) {
        coefficient_plots <- coefficient_plots |>
          purrr::imap(\(coefficient_plot, model_name) {
            patchwork::wrap_plots(
              GWAS_metadata_tracks_plot,
              coefficient_plot,
              nrow = 1,
              widths = c(5.8, 9),
              guides = "collect"
            )
          })
      }

      coefficient_plots |>
        save_plots_structured(
          width = max(20, 0.6 * dplyr::n_distinct(results_tibble.trait_level.pseudobulk$contrast) + 13),
          height = max(9, 0.55 * dplyr::n_distinct(results_tibble.trait_level.pseudobulk$GWAS_ID) + 3)
        )
    }
  ),
  tarchetypes::tar_file(
    name = chromVAR_deviation_heatmap.trait_level.pseudobulk,
    description = "Save cell-type pseudobulk GWAS chromVAR relative-deviation heatmap with z-score support labels and nuclei counts. [checkpoint:genetic_enrichment]",
    command = {
      plot <- plot_GWAS_by_cluster_heatmap(
        chromVAR_deviation_heatmap_data.trait_level.pseudobulk,
        GWAS_metadata_tracks_plot = GWAS_metadata_tracks_plot,
        compartments_patterns = genetic_enrichment_compartment_patterns,
        fill_col = "median_score",
        fill_label = "Relative deviation",
        support_label_col = "support_label"
      )
      save_plots_structured(
        plot,
        filetype = "png",
        width = 20,
        height = max(7, 0.24 * dplyr::n_distinct(chromVAR_deviation_heatmap_data.trait_level.pseudobulk$GWAS_ID) + 3)
      )
    }
  ),
  tarchetypes::tar_file(
    name = chromVAR_deviation_heatmap_unscaled.trait_level.pseudobulk,
    description = "Save cell-type pseudobulk GWAS chromVAR raw-deviation heatmap with z-score support labels and nuclei counts. [checkpoint:genetic_enrichment]",
    command = {
      plot <- plot_GWAS_by_cluster_heatmap(
        chromVAR_deviation_heatmap_data.trait_level.pseudobulk,
        GWAS_metadata_tracks_plot = GWAS_metadata_tracks_plot,
        compartments_patterns = genetic_enrichment_compartment_patterns,
        fill_col = "deviation",
        fill_label = "Deviation",
        support_label_col = "support_label"
      )
      save_plots_structured(
        plot,
        filetype = "png",
        width = 20,
        height = max(7, 0.24 * dplyr::n_distinct(chromVAR_deviation_heatmap_data.trait_level.pseudobulk$GWAS_ID) + 3)
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
