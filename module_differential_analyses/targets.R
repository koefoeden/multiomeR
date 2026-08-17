differential_analyses_aggregation_tibble <- aggregation_tibble |>
  dplyr::filter(aggregation_has_module(modules, "differential_analyses"))

differential_analyses_config_tibble <- read_module_config_tibble(
  config_file = "module_differential_analyses/cfg.yaml",
  module_name = "differential_analyses",
  module_aggregation_tibble = differential_analyses_aggregation_tibble,
  aggregation_tibble = aggregation_tibble_all_from_yaml
)

differential_analyses_tibble <- differential_analyses_aggregation_tibble |>
  dplyr::left_join(differential_analyses_config_tibble, by = "aggregation") |>
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
  tarchetypes::tar_map(
    values = differential_analyses_tibble,
    names = differential_analyses_target_suffix,
    descriptions = NULL,
    delimiter = ".",
    source("module_differential_analyses/setup_and_DCTC_targets.R")$value,
    source("module_differential_analyses/cross_modality_targets.R")$value,
    tarchetypes::tar_map(
      values = tibble::tribble(
        ~map_psbulk_DX_tar_suffix , ~map_psbulk_data_matrix                      ,
        "DGE"                     , rlang::sym("pseudobulk_counts_matrix.GEX")   ,
        "DCA"                     , rlang::sym("pseudobulk_counts_matrix.ATAC")  ,
        "DTFA"                    , rlang::sym("pseudobulk_motif_family_accessibility_matrix.ATAC")
      ),
      names = map_psbulk_DX_tar_suffix,
      descriptions = NULL,
      delimiter = ".",
      source("module_differential_analyses/psbulk_DX_targets.R")$value
    ),
    tarchetypes::tar_map(
      values = tibble::tribble(
        ~map_psbulk_DX_GSEA_tar_suffix, ~map_MSigDB_collection, ~map_MSigDB_subcollection,
        "H.GSEA.DGE",                   "H",                   NA_character_,
        "CP_REACTOME.GSEA.DGE",         "C2",                  "CP:REACTOME"
      ),
      names = map_psbulk_DX_GSEA_tar_suffix,
      descriptions = NULL,
      delimiter = ".",
      source("module_differential_analyses/GSEA_targets.R")$value
    )
  )
)
