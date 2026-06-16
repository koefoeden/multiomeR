rlang::list2(
  tarchetypes::tar_file(
    name = genetic_enrichment_donor_id_metadata_tsv.extended,
    description = "Locate the extended donor ID metadata TSV for genetic enrichment pseudobulk analyses",
    command = genetic_enrichment_extended_donor_id_metadata_tsv %||%
      aggregation_donor_id_metadata_tsv,
    deployment = "main"
  ),
  targets::tar_target(
    name = genetic_enrichment_donor_id_metadata_tibble.extended,
    description = "Read the genetic enrichment extended donor ID metadata TSV into a tibble",
    command = read_keyed_metadata_tibble(genetic_enrichment_donor_id_metadata_tsv.extended, "donor_id")
  ),
  targets::tar_target(
    name = GWAS_inputs_tibble,
    description = "Resolve configured Open Targets GWAS metadata and fine-mapping methods for this aggregation [part_of_graph:genetic_enrichment_single_nucleus] [part_of_graph:genetic_enrichment_pseudobulk]",
    command = {
      study_ids <- unique(GWAS_inputs_config_tibble$studyId)
      available_methods_tibble <- arrow::open_dataset(open_targets_credible_set_dataset_path) |>
        dplyr::filter(studyType == "gwas", studyId %in% study_ids) |>
        dplyr::select(studyId, finemappingMethod, confidence) |>
        dplyr::distinct() |>
        dplyr::collect()

      missing_studies <- setdiff(study_ids, available_methods_tibble$studyId)
      if (length(missing_studies) > 0) {
        stop("Configured Open Targets studyId(s) missing from credible_set: ", paste(missing_studies, collapse = ", "))
      }

      open_targets_release <- basename(dirname(open_targets_credible_set_dataset_path))
      method_priority <- c("SuSie", "SuSiE-inf", "PICS")
      GWAS_inputs_config_tibble |>
        dplyr::mutate(
          finemappingMethod = purrr::pmap_chr(
            dplyr::pick(GWAS_ID, studyId, requested_finemappingMethod),
            \(GWAS_ID, studyId, requested_finemappingMethod) {
              study_id <- studyId
              available <- available_methods_tibble |>
                dplyr::filter(.data$studyId == .env$study_id) |>
                dplyr::pull(finemappingMethod) |>
                unique()

              if (requested_finemappingMethod != "auto") {
                if (!requested_finemappingMethod %in% available) {
                  stop(
                    "Requested Open Targets finemappingMethod is unavailable for GWAS_ID=",
                    GWAS_ID,
                    ", studyId=",
                    studyId,
                    ". Requested: ",
                    requested_finemappingMethod,
                    ". Available: ",
                    paste(available, collapse = ", ")
                  )
                }
                return(requested_finemappingMethod)
              }

              selected <- method_priority[method_priority %in% available]
              if (length(selected) == 0) {
                stop(
                  "No supported Open Targets finemappingMethod found for GWAS_ID=",
                  GWAS_ID,
                  ", studyId=",
                  studyId,
                  ". Available: ",
                  paste(available, collapse = ", ")
                )
              }
              selected[[1]]
            }
          ),
          open_targets_release = open_targets_release
        ) |>
        dplyr::select(Category, GWAS_ID, studyId, finemappingMethod, variant_weighting_mode, open_targets_release)
    },
    iteration = "vector"
  ),
  targets::tar_target(
    name = GWAS_analysis_inputs_tibble,
    description = "Drop plot-only GWAS metadata before expensive GWAS chromVAR branching",
    command = GWAS_inputs_tibble |>
      dplyr::select(GWAS_ID, studyId, finemappingMethod, variant_weighting_mode, open_targets_release),
    iteration = "vector"
  ),
  targets::tar_target(
    name = GWAS_open_targets_metadata_tibble,
    description = "Summarize Open Targets study, ancestry, sample-size, and credible-set metadata for configured GWAS heatmaps",
    command = get_open_targets_GWAS_metadata_tibble(
      GWAS_inputs_tibble = GWAS_inputs_tibble,
      open_targets_study_dataset_path = open_targets_study_dataset_path,
      open_targets_credible_set_dataset_path = open_targets_credible_set_dataset_path
    ),
    iteration = "vector"
  ),
  targets::tar_target(
    name = GWAS_metadata_tracks_plot,
    description = "Build the shared Open Targets GWAS metadata track plot for genetic enrichment score plots",
    command = GWAS_open_targets_metadata_tibble |>
      dplyr::arrange(Category, GWAS_ID) |>
      dplyr::mutate(GWAS_ID = factor(GWAS_ID, levels = rev(GWAS_ID))) |>
      plot_GWAS_metadata_tracks()
  )
)
