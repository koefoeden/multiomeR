rlang::list2(
  tarchetypes::tar_file(
    name = cellranger_reference_json_files,
    description = "Discover and track supplied Cell Ranger reference metadata",
    command = c(!!list.files("reference_metadata", pattern = "^reference[.]json$", recursive = TRUE, full.names = TRUE)),
    deployment = "main"
  ),
  tarchetypes::tar_file(
    name = GEM_well_config_tsv,
    description = "Track the canonical GEM well processing and metadata TSV",
    command = "cfg_GEM_wells.tsv",
    deployment = "main"
  ),
  tarchetypes::tar_file(
    name = QC_metric_manifest_tsv,
    description = "Track the informational QC metric manifest",
    command = "QC_metric_manifest.tsv",
    deployment = "main"
  ),
  targets::tar_target(
    name = QC_metric_manifest_tibble,
    description = "Read and validate the informational QC metric manifest",
    command = {
      manifest <- readr::read_tsv(
        QC_metric_manifest_tsv,
        col_types = readr::cols(
          metric_id = readr::col_character(),
          display_name = readr::col_character(),
          description = readr::col_character(),
          available_from_checkpoint = readr::col_character(),
          plot_min_q = readr::col_double(),
          plot_max_q = readr::col_double(),
          do_plot = readr::col_logical()
        ),
        show_col_types = FALSE
      )
      expected_columns <- c(
        "metric_id",
        "display_name",
        "description",
        "available_from_checkpoint",
        "plot_min_q",
        "plot_max_q",
        "do_plot"
      )
      required_values <- manifest |>
        dplyr::select(
          metric_id,
          display_name,
          description,
          available_from_checkpoint,
          do_plot
        )
      plotting_quantiles <- c(manifest$plot_min_q, manifest$plot_max_q)
      checkpoint_names <- readr::read_tsv(
        QC_checkpoint_manifest_tsv,
        col_types = readr::cols(.default = readr::col_character())
      )$checkpoint_name
      invalid_quantiles <- any(
        !is.na(plotting_quantiles) &
          (!is.finite(plotting_quantiles) |
            plotting_quantiles < 0 |
            plotting_quantiles > 1)
      )
      invalid_quantile_intervals <- any(
        !is.na(manifest$plot_min_q) &
          !is.na(manifest$plot_max_q) &
          manifest$plot_min_q >= manifest$plot_max_q
      )
      if (
        !identical(colnames(manifest), expected_columns) ||
          anyDuplicated(manifest$metric_id) ||
          any(!stats::complete.cases(required_values)) ||
          any(!manifest$available_from_checkpoint %in% checkpoint_names) ||
          invalid_quantiles ||
          invalid_quantile_intervals
      ) {
        stop(
          "QC_metric_manifest.tsv must contain the expected columns, complete ",
          "required values, unique metric IDs, known checkpoint names, and valid plotting quantiles.",
          call. = FALSE
        )
      }
      manifest
    }
  ),
  tarchetypes::tar_file(
    name = QC_checkpoint_manifest_tsv,
    description = "Track the QC checkpoint names and review guidance",
    command = "QC_checkpoint_manifest.tsv",
    deployment = "main"
  ),
  tarchetypes::tar_file(
    name = amulet_BPCells_native_source_file,
    description = "Track the native BPCells fragment iterator used by AMULET",
    command = "src/amulet_bpcells.cpp"
  ),
  tarchetypes::tar_file(
    name = SCAVENGE_native_source_file,
    description = "Track the shared-memory SCAVENGE permutation random-walk kernel",
    command = "src/scavenge_random_walk.cpp"
  ),
  tarchetypes::tar_file(
    name = WNN_native_source_file,
    description = "Track the native WNN small-SNN bandwidth kernel",
    command = "src/wnn_snn_bandwidth.cpp"
  ),
  tarchetypes::tar_file(
    name = JASPAR2026_vertebrate_familial_root_motifs_tf,
    description = "Track the 233 official JASPAR2026 CORE vertebrate familial root motifs [part_of_graph:ATAC] [part_of_graph:seurat_export] [part_of_graph:differential_analyses]",
    command = "resources/JASPAR2026_vertebrate_familial_root_motifs.tf"
  ),
  tarchetypes::tar_file(
    name = JASPAR2026_vertebrate_motif_families_tsv,
    description = "Track the official JASPAR2026 CORE vertebrate familial motif membership map [part_of_graph:ATAC] [part_of_graph:seurat_export] [part_of_graph:differential_analyses]",
    command = "resources/JASPAR2026_vertebrate_motif_families.tsv"
  ),
  tarchetypes::tar_file(
    name = CollecTRI_human_network_csv,
    description = "Download and checksum the published human CollecTRI signed TF-target network [part_of_graph:differential_analyses]",
    command = download_CollecTRI_human_network(
      network_url = "https://rescued.omnipathdb.org/CollecTRI.csv",
      expected_sha256 = "86c90b30f2cc75c189da1f0a8c353d1547287cd656a9fac1c678634285bcb4e0"
    )
  ),
  tarchetypes::tar_file(
    name = open_targets_credible_set_dataset_path,
    description = "Download the Open Targets 26.03 credible_set Parquet dataset for shared GWAS consumers [part_of_graph:genetic_enrichment_single_nucleus]",
    command = download_open_targets_dataset("https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.03/output/credible_set/")
  ),
  tarchetypes::tar_file(
    name = open_targets_study_dataset_path,
    description = "Download the Open Targets 26.03 study Parquet dataset for shared GWAS metadata consumers",
    command = download_open_targets_dataset("https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.03/output/study/")
  ),
  tarchetypes::tar_file(
    name = open_targets_gwas_credible_sets_evidence_dataset_path,
    description = "Download the Open Targets 26.03 GWAS credible-set evidence Parquet dataset for L2G annotation",
    command = download_open_targets_dataset("https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.03/output/evidence_gwas_credible_sets/")
  ),
  tarchetypes::tar_file(
    name = open_targets_target_dataset_path,
    description = "Download the Open Targets 26.03 target Parquet dataset for L2G gene metadata",
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
    name = JASPAR_familial_root_motif_matrix_list,
    description = "Load the 233 JASPAR2026 CORE vertebrate familial root motifs as sequence-scanning inputs [part_of_graph:ATAC] [part_of_graph:seurat_export] [part_of_graph:differential_analyses]",
    command = {
      root_motifs <- read_JASPAR_familial_root_PFMatrixList(
        JASPAR2026_vertebrate_familial_root_motifs_tf
      )
      if (!identical(names(root_motifs), sprintf("cluster_%03d", seq_len(233L)))) {
        stop("The vendored JASPAR2026 root-motif file must contain cluster_001 through cluster_233 in order.")
      }
      root_motifs
    }
  ),
  targets::tar_target(
    name = JASPAR_motif_family_members_tibble,
    description = "Map individual JASPAR2026 CORE vertebrate motifs to 233 sequence-similarity families [part_of_graph:ATAC] [part_of_graph:seurat_export] [part_of_graph:differential_analyses]",
    command = {
      family_members <- readr::read_tsv(
        JASPAR2026_vertebrate_motif_families_tsv,
        show_col_types = FALSE
      ) |>
        dplyr::mutate(motif_feature = paste0(stringr::str_to_upper(TF_name), "__", motif_id))
      root_motif_ids <- vapply(
        seq_along(JASPAR_familial_root_motif_matrix_list),
        function(index) TFBSTools::ID(JASPAR_familial_root_motif_matrix_list[[index]]),
        character(1)
      )

      validation <- c(
        member_rows = nrow(family_members) == 1019L,
        motif_families = dplyr::n_distinct(family_members$motif_family) == 233L,
        unique_motif_features = !anyDuplicated(family_members$motif_feature),
        root_motif_names = setequal(
          family_members$motif_family,
          root_motif_ids
        )
      )
      if (!all(validation)) {
        stop(
          "The vendored JASPAR2026 motif-family map failed: ",
          paste(names(validation)[!validation], collapse = ", "),
          "."
        )
      }

      family_members
    }
  ),
  tarchetypes::tar_file(
    name = JASPAR_motif_family_annotations_tsv,
    description = "Track readable JASPAR2026 motif-family labels and complete membership lookup",
    command = "resources/JASPAR2026_vertebrate_motif_family_annotations.tsv"
  ),
  targets::tar_target(
    name = JASPAR_motif_family_labels,
    description = "Read unique display labels for all JASPAR2026 motif families without changing scoring identifiers",
    command = {
      annotations <- readr::read_tsv(JASPAR_motif_family_annotations_tsv, show_col_types = FALSE)
      stopifnot(
        !anyDuplicated(annotations$motif_family),
        !anyDuplicated(annotations$display_name),
        !anyNA(annotations$display_name),
        all(nzchar(annotations$display_name)),
        setequal(annotations$motif_family, JASPAR_motif_family_members_tibble$motif_family)
      )
      stats::setNames(annotations$display_name, annotations$motif_family)
    }
  ),
  targets::tar_target(
    name = CollecTRI_human_network_tibble,
    description = "Validate and reformat the checksum-pinned human CollecTRI signed TF-target network [part_of_graph:differential_analyses]",
    command = read_CollecTRI_human_network(CollecTRI_human_network_csv)
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
