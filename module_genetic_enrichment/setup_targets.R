rlang::list2(
  targets::tar_target(
    name = GWAS_analysis_inputs_config_tibble,
    description = "Drop plot-only GWAS fields and stabilize row order before source-specific input branching",
    command = GWAS_inputs_config_tibble |>
      dplyr::select(-Category, -config_order) |>
      dplyr::arrange(GWAS_ID),
    iteration = "vector"
  ),
  targets::tar_target(
    name = GWAS_source_files,
    description = "Track one local fine-mapped GWAS Parquet file per configured local source",
    command = resolve_GWAS_source_file(
      sourceId = GWAS_analysis_inputs_config_tibble$sourceId,
      sourceType = GWAS_analysis_inputs_config_tibble$sourceType
    ),
    pattern = map(GWAS_analysis_inputs_config_tibble),
    format = "file"
  ),
  targets::tar_target(
    name = GWAS_input_records,
    description = "Normalize one Open Targets study or local fine-mapped GWAS file for trait-level chromVAR branching",
    command = get_GWAS_chromVAR_input_record(
      GWAS_config_tibble = GWAS_analysis_inputs_config_tibble,
      local_source_file = GWAS_source_files,
      open_targets_study_dataset_path = open_targets_study_dataset_path,
      open_targets_credible_set_dataset_path = open_targets_credible_set_dataset_path
    ),
    pattern = map(GWAS_analysis_inputs_config_tibble, GWAS_source_files),
    iteration = "list",
    resources = get_tar_resources(RAM_GB_req = 40)
  ),
  targets::tar_target(
    name = GWAS_inputs_tibble,
    description = "Combine source-neutral GWAS metadata for this aggregation [part_of_graph:genetic_enrichment_single_nucleus] [part_of_graph:genetic_enrichment_cell_type_contributions]",
    command = purrr::map_dfr(GWAS_input_records, "metadata_tibble") |>
      dplyr::left_join(
        GWAS_inputs_config_tibble |>
          dplyr::select(GWAS_ID, Category, config_order),
        by = "GWAS_ID",
        relationship = "one-to-one"
      ) |>
      dplyr::arrange(config_order) |>
      dplyr::relocate(Category, .after = GWAS_ID),
    iteration = "vector"
  ),
  targets::tar_target(
    name = GWAS_metadata_tracks_plot,
    description = "Build the shared source-neutral GWAS metadata track plot for genetic enrichment score plots",
    command = GWAS_inputs_tibble |>
      dplyr::arrange(Category, GWAS_ID) |>
      dplyr::mutate(GWAS_ID = factor(GWAS_ID, levels = rev(GWAS_ID))) |>
      plot_GWAS_metadata_tracks()
  )
)
