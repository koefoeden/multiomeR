# Mapping tibbles ---------------------------------------------------------------
dataset_tibble_from_yaml <- read_config_tibble(config_file = "cfg_datasets.yaml", key_col = "dataset")
reaction_tibble <- build_reaction_tibble(dataset_tibble_from_yaml = dataset_tibble_from_yaml)

aggregation_tibble_all_from_yaml <- read_aggregation_config_tibble(config_file = "cfg_aggregations.yaml")
aggregation_tibble <- build_aggregation_tibble(
  aggregation_tibble_all_from_yaml = aggregation_tibble_all_from_yaml,
  reaction_tibble = reaction_tibble
)
dataset_tibble <- build_dataset_tibble(
  reaction_tibble = reaction_tibble,
  dataset_tibble_from_yaml = dataset_tibble_from_yaml
)
roadmap_EDACC_names <- get_roadmap_EDACC_names(aggregation_tibble = aggregation_tibble)

known_aggregation_modules <- c("differential_analyses", "genetic_enrichment")
validate_aggregation_module_names(
  aggregation_tibble = aggregation_tibble,
  known_modules = known_aggregation_modules
)

pipeline <- rlang::list2(
  source("extra_targets/setup_targets.R")$value,
  tarchetypes::tar_map(
    values = reaction_tibble,
    names = reaction_ID,
    descriptions = NULL,
    delimiter = ".",
    source("extra_targets/per_reaction_targets.R")$value
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
    source("extra_targets/peak_gene_correlation_targets.R")$value,
    source("extra_targets/subgroups.R")$value
  ),
  source("module_differential_analyses/targets.R")$value,
  source("module_genetic_enrichment/targets.R")$value
)

pipeline
