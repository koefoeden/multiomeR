# Mapping tibbles ---------------------------------------------------------------
GEM_well_tibble_all <- build_GEM_well_tibble()

aggregation_tibble_all_from_yaml <- read_aggregation_config_tibble(config_file = "cfg_aggregations.yaml")
aggregation_tibble <- build_aggregation_tibble(
  aggregation_tibble_all_from_yaml = aggregation_tibble_all_from_yaml,
  GEM_well_tibble = GEM_well_tibble_all
) |>
  dplyr::mutate(
    aggregation_ATAC_peak_calling_input_sym = purrr::map2(
      aggregation,
      aggregation_ATAC_peak_calling_method,
      \(aggregation_name, peak_calling_method) {
        input_target <- if (identical(peak_calling_method, "macs3")) {
          "fragments_per_peak_calling_cluster_discovery.fragments.ATAC"
        } else {
          "peak_calling_cluster_discovery_tibble.ATAC"
        }
        rlang::sym(stringr::str_c(input_target, aggregation_name, sep = "."))
      }
    ),
    aggregation_ATAC_peak_calling_pattern = purrr::pmap(
      list(
        aggregation,
        aggregation_ATAC_peak_calling_method,
        aggregation_ATAC_peak_calling_input_sym
      ),
      \(aggregation_name, peak_calling_method, peak_calling_input_sym) {
        discovery_tibble_sym <- rlang::sym(stringr::str_c(
          "peak_calling_cluster_discovery_tibble.ATAC",
          aggregation_name,
          sep = "."
        ))
        if (identical(peak_calling_method, "macs3")) {
          rlang::call2("map", peak_calling_input_sym, discovery_tibble_sym)
        } else {
          rlang::call2("map", discovery_tibble_sym)
        }
      }
    )
  )
aggregation_MACS3_tibble <- aggregation_tibble |>
  dplyr::filter(aggregation_ATAC_peak_calling_method == "macs3") |>
  add_aggregation_target_syms(c(
    "combined_BPCells_fragment_obj.ATAC",
    "peak_calling_cluster_discovery_tibble.ATAC"
  ))
GEM_well_tibble <- build_active_GEM_well_tibble(GEM_well_tibble_all)
dataset_tibble <- build_dataset_tibble(GEM_well_tibble = GEM_well_tibble)
roadmap_EDACC_names <- get_roadmap_EDACC_names(aggregation_tibble = aggregation_tibble)

known_aggregation_modules <- c(
  "differential_analyses",
  "genetic_enrichment",
  "peak_gene_correlation"
)
validate_aggregation_module_names(
  aggregation_tibble = aggregation_tibble,
  known_modules = known_aggregation_modules
)

peak_gene_correlation_aggregation_tibble <- aggregation_tibble |>
  dplyr::filter(aggregation_has_module(modules, "peak_gene_correlation")) |>
  add_aggregation_target_syms(c(
    "marker_validated_Ensembl_annotations_GRanges_list",
    "aggregated_counts_BPCells_matrix.GEX",
    "metadata_w_cell_types_tibble.WNN",
    "consensus_peak_GRanges.ATAC",
    "consensus_peak_BPCells_matrix.ATAC",
    "harmony_embeddings_matrix.ATAC",
    "combined_BPCells_fragment_obj.ATAC"
  ))

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
    values = dataset_tibble,
    names = dataset,
    descriptions = NULL,
    delimiter = ".",
    source("extra_targets/per_dataset_targets.R")$value
  ),
  tarchetypes::tar_map(
    values = aggregation_tibble,
    names = aggregation,
    descriptions = NULL,
    delimiter = ".",
    source("extra_targets/general_aggregation_targets.R")$value,
    source("extra_targets/GEX_merge_and_dim_reduc_targets.R")$value,
    source("extra_targets/GEX_graph_and_cluster_targets.R")$value,
    source("extra_targets/ATAC_targets.R")$value,
    source("extra_targets/WNN_targets.R")$value,
    source("extra_targets/Seurat_Signac_export_targets.R")$value,
    source("extra_targets/subgroups.R")$value
  ),
  tarchetypes::tar_map(
    values = peak_gene_correlation_aggregation_tibble,
    names = aggregation,
    descriptions = NULL,
    delimiter = ".",
    source("extra_targets/peak_gene_correlation_targets.R")$value
  ),
  tarchetypes::tar_map(
    values = aggregation_MACS3_tibble,
    names = aggregation,
    descriptions = NULL,
    delimiter = ".",
    source("extra_targets/ATAC_MACS3_targets.R")$value
  ),
  source("module_differential_analyses/targets.R")$value,
  source("module_genetic_enrichment/targets.R")$value
)

pipeline
