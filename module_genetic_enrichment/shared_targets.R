rlang::list2(
  tarchetypes::tar_file(
    name = open_targets_credible_set_dataset_path,
    description = "Download the Open Targets 26.03 credible_set Parquet dataset for shared GWAS consumers [part_of_graph:genetic_enrichment_single_nucleus]",
    command = download_open_targets_dataset(open_targets_dataset_urls()[["credible_set"]])
  ),
  tarchetypes::tar_file(
    name = open_targets_study_dataset_path,
    description = "Download the Open Targets 26.03 study Parquet dataset for shared GWAS metadata consumers",
    command = download_open_targets_dataset(open_targets_dataset_urls()[["study"]])
  ),
  tarchetypes::tar_file(
    name = open_targets_gwas_credible_sets_evidence_dataset_path,
    description = "Download the Open Targets 26.03 GWAS credible-set evidence Parquet dataset for L2G annotation",
    command = download_open_targets_dataset(open_targets_dataset_urls()[["evidence_gwas_credible_sets"]])
  ),
  tarchetypes::tar_file(
    name = open_targets_target_dataset_path,
    description = "Download the Open Targets 26.03 target Parquet dataset for L2G gene metadata",
    command = download_open_targets_dataset(open_targets_dataset_urls()[["target"]])
  )
)
