differential_analyses_restriction <- build_module_tibble("differential_analyses", aggregation_tibble, aggregation_tibble_all_from_yaml) |>
  dplyr::mutate(differential_analyses_target_suffix = stringr::str_c("differential_analyses", aggregation, sep = ".")) |>
  add_aggregation_target_syms(c(
    "metadata_w_cell_types_tibble.WNN",
    "pseudobulk_counts_matrix.GEX",
    "pseudobulk_counts_matrix.ATAC",
    "pseudobulk_depth_tibble.GEX",
    "pseudobulk_depth_tibble.ATAC",
    "pseudobulk_motif_family_accessibility_matrix.ATAC",
    "chromVAR_obj.ATAC",
    "chromVAR_TF_motif_family_matrix.ATAC",
    "aggregated_cellranger_ref_list",
    "consensus_peak_annotated_GRanges.ATAC",
    "gene_features_df",
    "organism_chr"
  )) |>
  restrict_module_to_GEM_well_datasets("differential_analyses")

rlang::list2(
  if (nrow(differential_analyses_restriction$values) > 0L) source("extra_targets/module_differential_analyses/shared_targets.R")$value,
  differential_analyses_restriction$targets,
  tarchetypes::tar_map(
    values = differential_analyses_restriction$values,
    names = differential_analyses_target_suffix,
    descriptions = NULL,
    delimiter = ".",
    source("extra_targets/module_differential_analyses/setup_and_cell_type_composition_targets.R")$value,
    tarchetypes::tar_map(
      values = tibble::tribble(
        ~map_analysis_suffix, ~map_pseudobulk_data_matrix,
        "gene_expression", rlang::sym("pseudobulk_counts_matrix.GEX"),
        "chromatin_accessibility", rlang::sym("pseudobulk_counts_matrix.ATAC"),
        "motif_family_accessibility", rlang::sym("pseudobulk_motif_family_accessibility_matrix.ATAC"),
        "TF_motif_accessibility", rlang::sym("pseudobulk_TF_motif_accessibility_matrix")
      ),
      names = map_analysis_suffix,
      descriptions = NULL,
      delimiter = ".",
      source("extra_targets/module_differential_analyses/pseudobulk_differential_targets.R")$value
    ),
    source("extra_targets/module_differential_analyses/regulon_target_enrichment_targets.R")$value,
    source("extra_targets/module_differential_analyses/cross_modality_targets.R")$value,
    tarchetypes::tar_map(
      values = tibble::tribble(
        ~map_gene_set_enrichment_suffix, ~map_MSigDB_collection, ~map_MSigDB_subcollection,
        "Hallmark.gene_set_enrichment.gene_expression", "H", NA_character_,
        "Reactome.gene_set_enrichment.gene_expression", "C2", "CP:REACTOME"
      ),
      names = map_gene_set_enrichment_suffix,
      descriptions = NULL,
      delimiter = ".",
      source("extra_targets/module_differential_analyses/gene_set_enrichment_targets.R")$value
    )
  )
)
