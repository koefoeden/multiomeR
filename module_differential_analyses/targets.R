differential_analyses_tibble <- build_module_tibble("differential_analyses", aggregation_tibble, aggregation_tibble_all_from_yaml) |>
  dplyr::mutate(differential_analyses_target_suffix = stringr::str_c("differential_analyses", aggregation, sep = ".")) |>
  add_aggregation_target_syms(c(
    "metadata_w_cell_types_tibble.WNN",
    "pseudobulk_counts_matrix.GEX",
    "pseudobulk_counts_matrix.ATAC",
    "pseudobulk_depth_tibble.GEX",
    "pseudobulk_depth_tibble.ATAC",
    "pseudobulk_motif_family_accessibility_matrix.ATAC",
    "consensus_peak_annotated_GRanges.ATAC",
    "gene_features_df",
    "organism_chr"
  ))

rlang::list2(
  if (nrow(differential_analyses_tibble) > 0L) source("module_differential_analyses/shared_targets.R")$value,
  tarchetypes::tar_map(
    values = differential_analyses_tibble,
    names = differential_analyses_target_suffix,
    descriptions = NULL,
    delimiter = ".",
    source("module_differential_analyses/setup_and_cell_type_composition_targets.R")$value,
    targets::tar_target(
      name = pseudobulk_CollecTRI_TF_activity_matrix.GEX,
      description = "Infer signed CollecTRI ULM TF activities from GEX pseudobulks [part_of_graph:differential_analyses]",
      command = get_pseudobulk_CollecTRI_TF_activity_matrix(
        pseudobulk_GEX_counts_matrix = pseudobulk_counts_matrix.GEX,
        CollecTRI_network_tibble = CollecTRI_human_network_tibble
      ),
      packages = w_def("decoupleR"),
      resources = get_tar_resources(RAM_GB_req = 32)
    ),
    tarchetypes::tar_map(
      values = tibble::tribble(
        ~map_analysis_suffix, ~map_pseudobulk_data_matrix,
        "gene_expression", rlang::sym("pseudobulk_counts_matrix.GEX"),
        "chromatin_accessibility", rlang::sym("pseudobulk_counts_matrix.ATAC"),
        "motif_family_accessibility", rlang::sym("pseudobulk_motif_family_accessibility_matrix.ATAC"),
        "transcription_factor_activity", rlang::sym("pseudobulk_CollecTRI_TF_activity_matrix.GEX")
      ),
      names = map_analysis_suffix,
      descriptions = NULL,
      delimiter = ".",
      source("module_differential_analyses/pseudobulk_differential_targets.R")$value
    ),
    source("module_differential_analyses/cross_modality_targets.R")$value,
    tarchetypes::tar_map(
      values = tibble::tribble(
        ~map_gene_set_enrichment_suffix, ~map_MSigDB_collection, ~map_MSigDB_subcollection,
        "Hallmark.gene_set_enrichment.gene_expression", "H", NA_character_,
        "Reactome.gene_set_enrichment.gene_expression", "C2", "CP:REACTOME"
      ),
      names = map_gene_set_enrichment_suffix,
      descriptions = NULL,
      delimiter = ".",
      source("module_differential_analyses/gene_set_enrichment_targets.R")$value
    )
  )
)
