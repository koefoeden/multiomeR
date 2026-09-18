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
    command = GEM_well_config_file,
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
  targets::tar_target(
    name = Ensembl_gene_annotation_GRanges_list,
    description = "Download all supported Ensembl gene annotations once for downstream reference-specific lookup",
    command = with_annotation_hub_cache_lock({
      AnnotationHub::setAnnotationHubOption("CACHE", file.path(targets::tar_config_get("store"), "files", "AnnotationHub"))
      annot_hub_interface <- AnnotationHub::AnnotationHub(ask = FALSE)

      annotation_hub_ensembl_ids() |>
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
