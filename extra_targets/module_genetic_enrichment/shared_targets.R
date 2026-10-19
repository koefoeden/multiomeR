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
  ),
  tarchetypes::tar_file(
    name = SLDSC_reference_dir,
    description = "Download the GRCh38 1000 Genomes EUR LD reference, baseline-LD v2.2 model and HapMap3 regression weights for S-LDSC (Zenodo 10515792) [part_of_graph:genetic_enrichment_SLDSC]",
    command = download_SLDSC_reference(
      url = SLDSC_reference_source()$url,
      md5 = SLDSC_reference_source()$md5
    )
  ),
  targets::tar_target(
    name = SLDSC_baseline_record,
    description = "Read baseline-LD v2.2 LD scores, regression weights and reference alleles at the S-LDSC regression SNPs [part_of_graph:genetic_enrichment_SLDSC]",
    command = read_SLDSC_baseline_record(SLDSC_reference_dir),
    resources = get_tar_resources(RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = SLDSC_regression_SNP_tibble,
    description = "Extract the S-LDSC regression SNPs with their reference alleles",
    command = SLDSC_baseline_record$SNP_tibble
  ),
  targets::tar_target(
    name = SLDSC_sumstats_sources,
    description = "List the distinct GWAS sources and summary-statistics settings of all genetic-enrichment aggregations for S-LDSC",
    command = tibble::tibble(!!!SLDSC_sumstats_source_tibble),
    iteration = "vector",
    deployment = "main"
  ),
  targets::tar_target(
    name = SLDSC_sumstats_files,
    description = "Track one GWAS source's configured local summary-statistics file and its GWAS-SSF metadata",
    command = resolve_SLDSC_sumstats_files(SLDSC_sumstats_sources$summaryStatistics),
    pattern = map(SLDSC_sumstats_sources),
    format = "file",
    deployment = "main"
  ),
  targets::tar_target(
    name = SLDSC_sumstats_rows,
    description = "Fetch one GWAS source's genome-wide summary statistics and keep their S-LDSC columns at the regression SNPs [part_of_graph:genetic_enrichment_SLDSC]",
    command = fetch_SLDSC_sumstats_rows(SLDSC_sumstats_sources, SLDSC_sumstats_files, SLDSC_regression_SNP_tibble$SNP),
    pattern = map(SLDSC_sumstats_sources, SLDSC_sumstats_files),
    iteration = "list",
    # One source at a time keeps EBI's FTP server from refusing connections.
    deployment = "main"
  ),
  targets::tar_target(
    name = SLDSC_sumstats_records,
    description = "Munge one GWAS source's summary statistics for S-LDSC [part_of_graph:genetic_enrichment_SLDSC]",
    command = get_SLDSC_sumstats_record(SLDSC_sumstats_rows, SLDSC_regression_SNP_tibble),
    pattern = map(SLDSC_sumstats_rows),
    iteration = "list",
    deployment = "main"
  )
)
