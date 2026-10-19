genetic_enrichment_tibble <- build_module_tibble("genetic_enrichment", aggregation_tibble, aggregation_tibble_all_from_yaml) |>
  dplyr::mutate(genetic_enrichment_target_suffix = stringr::str_c("genetic_enrichment", aggregation, sep = "."))

genetic_enrichment_GWAS_config_tibble <- if (nrow(genetic_enrichment_tibble) == 0) {
  tibble::tibble(aggregation = character(), GWAS_inputs_config_tibble = list(), SLDSC_sumstats_index = list())
} else {
  GWAS_config_tibble <- purrr::map_dfr(seq_len(nrow(genetic_enrichment_tibble)), \(i) {
    row <- genetic_enrichment_tibble[i, , drop = FALSE]
    tibble::enframe(row$genetic_enrichment_GWAS_studies[[1]], name = "GWAS_ID", value = "GWAS_config") |>
      dplyr::mutate(
        aggregation = row$aggregation,
        config_order = dplyr::row_number(),
        Category = purrr::map2_chr(GWAS_config, GWAS_ID, \(config, GWAS_ID) {
          if (is.null(config$Category)) {
            stop("Missing Category for GWAS_ID: ", GWAS_ID)
          }
          as.character(config$Category)
        }),
        sourceId = purrr::map2_chr(GWAS_config, GWAS_ID, \(config, GWAS_ID) {
          if (is.null(config$sourceId)) {
            stop("Missing sourceId for GWAS_ID: ", GWAS_ID)
          }
          as.character(config$sourceId)
        }),
        sourceType = purrr::map_chr(sourceId, classify_GWAS_source),
        summaryStatistics = purrr::map_chr(GWAS_config, \(config) as.character(config$summaryStatistics %||% NA)),
        requested_finemappingMethod = purrr::pmap_chr(
          dplyr::pick(GWAS_config, sourceType, GWAS_ID),
          \(GWAS_config, sourceType, GWAS_ID) {
            finemappingMethod <- GWAS_config$finemappingMethod
            if (identical(sourceType, "local_file")) {
              if (!is.null(finemappingMethod)) {
                stop("finemappingMethod must come from the local file for GWAS_ID: ", GWAS_ID)
              }
              return(NA_character_)
            }
            if (is.null(finemappingMethod) || is.na(finemappingMethod) || identical(finemappingMethod, "")) {
              return("auto")
            }
            method_lower <- stringr::str_to_lower(finemappingMethod)
            dplyr::case_when(
              method_lower == "auto" ~ "auto",
              method_lower == "susie" ~ "SuSie",
              method_lower == "susie-inf" ~ "SuSiE-inf",
              method_lower == "pics" ~ "PICS",
              TRUE ~ as.character(finemappingMethod)
            )
          }
        )
      ) |>
      dplyr::select(
        aggregation,
        config_order,
        Category,
        GWAS_ID,
        sourceId,
        sourceType,
        requested_finemappingMethod,
        summaryStatistics
      )
  })

  # S-LDSC fetches each distinct source's summary statistics once for all
  # aggregations. Each aggregation slices these shared branches in the
  # GWAS_ID order of its GWAS_input_records.
  SLDSC_sumstats_source_tibble <- GWAS_config_tibble |>
    dplyr::distinct(sourceId, summaryStatistics) |>
    dplyr::arrange(sourceId, summaryStatistics)
  SLDSC_sumstats_index_tibble <- GWAS_config_tibble |>
    dplyr::left_join(
      dplyr::mutate(SLDSC_sumstats_source_tibble, SLDSC_sumstats_index = dplyr::row_number()),
      by = c("sourceId", "summaryStatistics"),
      relationship = "many-to-one"
    ) |>
    dplyr::arrange(GWAS_ID) |>
    dplyr::summarise(SLDSC_sumstats_index = list(SLDSC_sumstats_index), .by = aggregation)

  GWAS_config_tibble |>
    dplyr::select(-summaryStatistics) |>
    dplyr::group_by(aggregation) |>
    tidyr::nest(
      GWAS_inputs_config_tibble = c(
        config_order,
        Category,
        GWAS_ID,
        sourceId,
        sourceType,
        requested_finemappingMethod
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::left_join(SLDSC_sumstats_index_tibble, by = "aggregation", relationship = "one-to-one")
}

genetic_enrichment_tibble <- genetic_enrichment_tibble |>
  dplyr::left_join(genetic_enrichment_GWAS_config_tibble, by = "aggregation") |>
  add_aggregation_target_syms(c(
    "marker_validated_Ensembl_annotations_GRanges_list",
    "pseudobulk_counts_matrix.GEX",
    "metadata_w_cell_types_tibble.WNN",
    "combined_BPCells_fragment_obj.ATAC",
    "consensus_peak_GRanges.ATAC",
    "peak_QC_filtered_BPCells_matrix.ATAC",
    "chromVAR_obj.ATAC",
    "chromVAR_peak_expectation.ATAC",
    "chromVAR_background_bins.ATAC",
    "chromVAR_chunk_context_records.ATAC",
    "harmony_embeddings_matrix.GEX",
    "harmony_embeddings_matrix.ATAC",
    "WNN_results",
    "UMAP_embeddings_tibble.WNN"
  ))

rlang::list2(
  if (nrow(genetic_enrichment_tibble) > 0L) source("extra_targets/module_genetic_enrichment/shared_targets.R")$value,
  tarchetypes::tar_map(
    values = genetic_enrichment_tibble,
    names = genetic_enrichment_target_suffix,
    descriptions = NULL,
    delimiter = ".",
    source("extra_targets/module_genetic_enrichment/setup_targets.R")$value,
    source("extra_targets/module_genetic_enrichment/gchromVAR_targets.R")$value,
    source("extra_targets/module_genetic_enrichment/GWAS_chromVAR_cell_type_targets.R")$value,
    source("extra_targets/module_genetic_enrichment/GWAS_chromVAR_absolute_effect_targets.R")$value,
    source("extra_targets/module_genetic_enrichment/GWAS_chromVAR_contribution_targets.R")$value,
    source("extra_targets/module_genetic_enrichment/SLDSC_targets.R")$value,
    tarchetypes::tar_map(
      values = tibble::tibble(
        map_SCAVENGE_tar_suffix = "WNN_harmony_SNN.SCAVENGE.single_nucleus",
        map_SCAVENGE_graph_name = "WNN_harmony_SNN",
        map_SCAVENGE_umap_cols = list(c("WNN_UMAP_1", "WNN_UMAP_2"))
      ),
      names = map_SCAVENGE_tar_suffix,
      descriptions = NULL,
      delimiter = ".",
      source("extra_targets/module_genetic_enrichment/SCAVENGE_graph_targets.R")$value,
      source("extra_targets/module_genetic_enrichment/SCAVENGE_group_targets.R")$value
    )
  )
)
