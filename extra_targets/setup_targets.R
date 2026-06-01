rlang::list2(
  tarchetypes::tar_file(
    name = open_targets_credible_set_dataset_path,
    description = "Download the Open Targets 26.03 credible_set Parquet dataset for shared GWAS consumers",
    command = download_open_targets_dataset("https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.03/output/credible_set/")
  ),
  tarchetypes::tar_file(
    name = open_targets_study_dataset_path,
    description = "Download the Open Targets 26.03 study Parquet dataset for shared GWAS metadata consumers",
    command = download_open_targets_dataset("https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.03/output/study/")
  ),
  tarchetypes::tar_file(
    name = open_targets_gwas_credible_sets_evidence_dataset_path,
    description = "Download the Open Targets 26.03 GWAS credible-set evidence Parquet dataset for L2G consumers",
    command = download_open_targets_dataset("https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.03/output/evidence_gwas_credible_sets/")
  ),
  tarchetypes::tar_file(
    name = open_targets_target_dataset_path,
    description = "Download the Open Targets 26.03 target Parquet dataset for target gene metadata",
    command = download_open_targets_dataset("https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.03/output/target/")
  ),
  targets::tar_target(
    name = Ensembl_gene_annotation_GRanges_list,
    description = "Download all supported Ensembl gene annotations once for downstream reference-specific lookup",
    command = with_annotation_hub_cache_lock({
      AnnotationHub::setAnnotationHubOption("CACHE", file.path(targets::tar_config_get("store"), "files", "AnnotationHub"))
      annot_hub_interface <- AnnotationHub::AnnotationHub(ask = FALSE)

      c("AH113665", "AH113713", "AH75011", "AH75036") |>
        purrr::set_names() |>
        purrr::map(\(annotation_hub_id) {
          get_gene_annotation_GRanges_from_EnsDb(annot_hub_interface[[annotation_hub_id]])
        })
    }),
    resources = get_tar_resources(cores_req = 1, RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = TF_motif_matrix_list,
    description = "Load vertebrate TF motif matrices",
    command = {
      JASPAR2024_db <- RSQLite::dbConnect(
        RSQLite::SQLite(),
        JASPAR2024::db(JASPAR2024::JASPAR2024())
      )
      on.exit(RSQLite::dbDisconnect(JASPAR2024_db))

      TF_motif_matrix_list <- TFBSTools::getMatrixSet(
        x = JASPAR2024_db,
        list(tax_group = "vertebrates", collection = "CORE", all_versions = FALSE)
      )
      motif_names <- vapply(seq_along(TF_motif_matrix_list), \(idx) TFBSTools::name(TF_motif_matrix_list[[idx]]), character(1))
      names(TF_motif_matrix_list) <- paste0(stringr::str_to_upper(motif_names), "__", names(TF_motif_matrix_list))
      TF_motif_matrix_list
    }
  ),
  targets::tar_target(
    name = chromHMMs_list_general,
    description = "Load and liftover Roadmap Epigenomics chromHMM state annotations for all configured EDACC names",
    command = with_annotation_hub_cache_lock({
      AnnotationHub::setAnnotationHubOption("CACHE", file.path(targets::tar_config_get("store"), "files", "AnnotationHub"))
      annot_hub_interface <- AnnotationHub::AnnotationHub(ask = FALSE)
      get_roadmap_chromHMMs_from_annotation_hub(
        annot_hub_interface = annot_hub_interface,
        roadmap_EDACC_names = roadmap_EDACC_names
      )
    }),
    packages = w_def("AnnotationHub"),
    resources = get_tar_resources(cores_req = 1, RAM_GB_req = 16)
  )
)
