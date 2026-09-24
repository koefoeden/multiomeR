assert_genome_packages_installed()

# Mapping tibbles ---------------------------------------------------------------
GEM_well_config_file <- configuration_path("cfg_GEM_wells.tsv")
GEM_well_tibble_all <- build_GEM_well_tibble(GEM_well_config_file)

aggregation_tibble_all_from_yaml <- read_aggregation_config_tibble(config_file = configuration_path("cfg_aggregations.yaml"))
if (length(validation_aggregations())) {
  message("Validation aggregations: ", paste(validation_aggregations(), collapse = ", "))
}
aggregation_tibble <- build_aggregation_tibble(
  aggregation_tibble_all_from_yaml = aggregation_tibble_all_from_yaml,
  GEM_well_tibble = GEM_well_tibble_all
)
peak_calling_aggregation_tibble <- aggregation_tibble |>
  add_aggregation_target_syms(c(
    "combined_BPCells_fragment_obj.ATAC",
    "peak_calling_cluster_discovery_tibble.ATAC",
    "aggregated_cellranger_ref_list"
  ))
GEM_well_tibble <- build_active_GEM_well_tibble(GEM_well_tibble_all)
roadmap_EDACC_names <- get_roadmap_EDACC_names(aggregation_tibble = aggregation_tibble)

pipeline <- rlang::list2(
  source("extra_targets/setup_targets.R")$value,
  tarchetypes::tar_map(
    values = GEM_well_tibble,
    names = GEM_well_ID,
    descriptions = NULL,
    delimiter = ".",
    source("extra_targets/per_GEM_well_targets.R")$value
  ),
  tarchetypes::tar_map(
    values = add_aggregation_target_syms(aggregation_tibble, "peaks_per_cluster_narrowPeaks.peaks.ATAC"),
    names = aggregation,
    descriptions = NULL,
    delimiter = ".",
    source("extra_targets/general_aggregation_targets.R")$value,
    source("extra_targets/GEX_merge_and_dim_reduc_targets.R")$value,
    source("extra_targets/GEX_graph_and_cluster_targets.R")$value,
    source("extra_targets/ATAC_targets.R")$value,
    source("extra_targets/WNN_targets.R")$value,
    source("extra_targets/Seurat_Signac_export_targets.R")$value
  ),
  tarchetypes::tar_map(
    values = dplyr::filter(peak_calling_aggregation_tibble, aggregation_ATAC_peak_calling_method == "macs3"),
    names = aggregation,
    descriptions = NULL,
    delimiter = ".",
    source("extra_targets/ATAC_MACS3_targets.R")$value
  ),
  tarchetypes::tar_map(
    values = dplyr::filter(peak_calling_aggregation_tibble, aggregation_ATAC_peak_calling_method == "bpcells_tile"),
    names = aggregation,
    descriptions = NULL,
    delimiter = ".",
    source("extra_targets/ATAC_tile_targets.R")$value
  ),
  source("module_differential_analyses/targets.R")$value,
  source("module_genetic_enrichment/targets.R")$value,
  source("module_peak_gene_correlation/targets.R")$value
)

pipeline
