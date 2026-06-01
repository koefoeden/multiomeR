add_genetic_enrichment_module_target_syms <- function(module_tibble, target_names) {
  target_sym_cols <- target_names |>
    purrr::set_names() |>
    purrr::map(\(target_name) {
      purrr::map(
        module_tibble$genetic_enrichment_target_suffix,
        \(suffix) rlang::sym(stringr::str_c(target_name, suffix, sep = "."))
      )
    })

  module_tibble |>
    dplyr::mutate(!!!target_sym_cols)
}

genetic_enrichment_aggregation_tibble <- aggregation_tibble |>
  dplyr::filter(aggregation_has_module(modules, "genetic_enrichment"))

genetic_enrichment_config_tibble <- read_module_config_tibble(
  config_file = "module_genetic_enrichment/cfg.yaml"
)

validate_module_config(
  module_name = "genetic_enrichment",
  module_config_tibble = genetic_enrichment_config_tibble,
  module_aggregation_tibble = genetic_enrichment_aggregation_tibble,
  aggregation_tibble = aggregation_tibble_all_from_yaml
)

genetic_enrichment_tibble <- genetic_enrichment_aggregation_tibble |>
  dplyr::left_join(genetic_enrichment_config_tibble, by = "aggregation") |>
  dplyr::mutate(genetic_enrichment_target_suffix = stringr::str_c("genetic_enrichment", aggregation, sep = "."))

genetic_enrichment_GWAS_config_tibble <- if (nrow(genetic_enrichment_tibble) == 0) {
  tibble::tibble(aggregation = character(), GWAS_inputs_config_tibble = list())
} else {
  GWAS_config_tibble <- purrr::map_dfr(seq_len(nrow(genetic_enrichment_tibble)), \(i) {
    row <- genetic_enrichment_tibble[i, , drop = FALSE]
    open_targets_studies <- row$genetic_enrichment_open_targets_studies[[1]]
    if (is.null(open_targets_studies)) {
      stop("genetic_enrichment_open_targets_studies must be configured for aggregation: ", row$aggregation)
    }

    tibble::enframe(open_targets_studies, name = "GWAS_ID", value = "open_targets_config") |>
      dplyr::mutate(
        aggregation = row$aggregation,
        Category = purrr::map2_chr(open_targets_config, GWAS_ID, \(config, GWAS_ID) {
          if (is.null(config$Category)) {
            stop("Missing Open Targets Category for GWAS_ID: ", GWAS_ID)
          }
          as.character(config$Category)
        }),
        studyId = purrr::map2_chr(open_targets_config, GWAS_ID, \(config, GWAS_ID) {
          if (is.null(config$studyId)) {
            stop("Missing Open Targets studyId for GWAS_ID: ", GWAS_ID)
          }
          as.character(config$studyId)
        }),
        requested_finemappingMethod = purrr::map_chr(open_targets_config, \(config) {
          finemappingMethod <- config$finemappingMethod
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
        }),
        variant_weighting_mode = purrr::map2_chr(open_targets_config, GWAS_ID, \(config, GWAS_ID) {
          normalize_open_targets_variant_weighting_mode(config$variant_weighting_mode, GWAS_ID = GWAS_ID)
        })
      ) |>
      dplyr::select(aggregation, Category, GWAS_ID, studyId, requested_finemappingMethod, variant_weighting_mode)
  })

  duplicated_GWAS <- GWAS_config_tibble |>
    dplyr::count(aggregation, GWAS_ID) |>
    dplyr::filter(n > 1)
  if (nrow(duplicated_GWAS) > 0) {
    stop("Each configured trait must map to one GWAS per aggregation. Duplicated GWAS_ID(s): ", paste(duplicated_GWAS$GWAS_ID, collapse = ", "))
  }

  GWAS_config_tibble |>
    dplyr::group_by(aggregation) |>
    tidyr::nest(GWAS_inputs_config_tibble = c(Category, GWAS_ID, studyId, requested_finemappingMethod, variant_weighting_mode)) |>
    dplyr::ungroup()
}

genetic_enrichment_tibble <- genetic_enrichment_tibble |>
  dplyr::left_join(genetic_enrichment_GWAS_config_tibble, by = "aggregation") |>
  add_aggregation_target_syms(c(
    "metadata_w_cell_types_tibble.WNN",
    "combined_BPCells_fragment_obj.ATAC",
    "consensus_peak_GRanges.ATAC",
    "pseudobulk_counts_matrix.ATAC",
    "pseudobulk_depth_tibble.ATAC",
    "chromVAR_obj.ATAC",
    "chromVAR_peak_expectation.ATAC",
    "chromVAR_background_bins.ATAC",
    "chromVAR_chunk_context_records.ATAC",
    "harmony_embeddings_matrix.GEX",
    "harmony_embeddings_matrix.ATAC",
    "WNN_results"
  ))

genetic_enrichment_gene_chromVAR_tibble <- genetic_enrichment_tibble |>
  dplyr::filter(purrr::map_lgl(
    .data$genetic_enrichment_psbulk_GWAS_gene_chromVAR_GWAS_IDs,
    \(x) !is.null(x) && length(x) > 0
  )) |>
  add_genetic_enrichment_module_target_syms(c(
    "genetic_enrichment_donor_id_metadata_tibble.extended",
    "GWAS_analysis_inputs_tibble",
    "genetic_enrichment_peak_ranges",
    "posterior_probability_weighting_function",
    "psbulk_GWAS_chromVAR_dynamic_tibble"
  ))

rlang::list2(
  tarchetypes::tar_map(
    values = genetic_enrichment_tibble,
    names = genetic_enrichment_target_suffix,
    descriptions = NULL,
    delimiter = ".",
    source("module_genetic_enrichment/setup_targets.R")$value,
    source("module_genetic_enrichment/gchromVAR_targets.R")$value,
    source("module_genetic_enrichment/psbulk_GWAS_chromVAR_targets.R")$value,
    tarchetypes::tar_map(
      values = tibble::tibble(
        map_SCAVENGE_tar_suffix = c("SCAVENGE.PCA_harmony_SNN", "SCAVENGE.LSI_harmony_SNN", "SCAVENGE.WNN_harmony_SNN"),
        map_SCAVENGE_graph_name = c("PCA_harmony_SNN", "LSI_harmony_SNN", "WNN_harmony_SNN"),
        map_SCAVENGE_graph_input_type = c("embedding", "embedding", "WNN"),
        map_SCAVENGE_graph_input = rlang::syms(c(
          "harmony_embeddings_matrix.GEX",
          "harmony_embeddings_matrix.ATAC",
          "WNN_results"
        )),
        map_SCAVENGE_embedding_dims = rlang::syms(c(
          "aggregation_GEX_data_PCs",
          "aggregation_ATAC_data_PCs",
          "aggregation_GEX_data_PCs"
        )),
        map_SCAVENGE_dim_prefix = c("PCA_", "LSI_", NA_character_),
        map_SCAVENGE_umap_cols = list(
          c("GEX_UMAP_1", "GEX_UMAP_2"),
          c("LSI_UMAP_1", "LSI_UMAP_2"),
          c("WNN_UMAP_1", "WNN_UMAP_2")
        )
      ),
      names = map_SCAVENGE_tar_suffix,
      descriptions = NULL,
      delimiter = ".",
      source("module_genetic_enrichment/SCAVENGE_graph_targets.R")$value,
      source("module_genetic_enrichment/SCAVENGE_group_targets.R")$value
    )
  ),
  tarchetypes::tar_map(
    values = genetic_enrichment_gene_chromVAR_tibble,
    names = genetic_enrichment_target_suffix,
    descriptions = NULL,
    delimiter = ".",
    source("module_genetic_enrichment/psbulk_GWAS_gene_chromVAR_targets.R")$value
  )
)
