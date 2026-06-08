rlang::list2(
  targets::tar_target(
    name = aggregated_cellranger_ref_list,
    description = "Use the validated CellRanger reference metadata shared by this aggregation",
    command = aggregation_cellranger_ref_list_syms[[1]]
  ),
  targets::tar_target(
    name = gene_features_df,
    description = "Use the CellRanger gene metadata from this aggregation reference",
    command = aggregation_gene_features_df_syms[[1]]
  ),
  targets::tar_target(
    name = aggregation_unfiltered_cells_n,
    description = "Sum full QC metadata barcode counts across all reactions in this aggregation",
    command = sum(unlist(aggregation_unfiltered_cells_n_vecs_syms))
  ),
  targets::tar_target(
    name = aggregation_excluded_BCs_list,
    description = "Combine per-reaction QC exclusion lists into a single aggregation-level list",
    command = collapse_duplicate_names(purrr::flatten(
      aggregation_excluded_barcodes_by_type_list_syms
    ))
  ),
  tarchetypes::tar_file(
    name = aggregation_excluded_barcodes_by_type_upset,
    description = "Plot an UpSet plot of aggregation-level QC exclusion overlaps across all reactions and save to file. [checkpoint:GEX]",
    command = plot_upset_from_excluded_BCs_list(
      QC_excluded_BCs_list = aggregation_excluded_BCs_list,
      n_total = aggregation_unfiltered_cells_n
    ) |>
      save_plots_structured(width = 20, height = 7)
  ),
  targets::tar_target(
    name = aggregation_excluded_cellranger_only_BCs_list,
    description = "Combine per-reaction CellRanger-only QC exclusion lists into a single aggregation-level list",
    command = collapse_duplicate_names(purrr::flatten(
      aggregation_excluded_cellranger_only_barcodes_by_type_list_syms
    ))
  ),
  tarchetypes::tar_file(
    name = aggregation_excluded_cellranger_only_barcodes_by_type_upset,
    description = "Plot an UpSet plot of aggregation-level CellRanger-only QC exclusion overlaps across all reactions and save to file. [checkpoint:GEX]",
    command = plot_upset_from_excluded_BCs_list(
      QC_excluded_BCs_list = aggregation_excluded_cellranger_only_BCs_list,
      n_total = nrow(dplyr::bind_rows(aggregation_cellranger_kept_metadata_tibble_syms))
    ) |>
      save_plots_structured(width = 20, height = 7)
  ),
  targets::tar_target(
    name = UCell_GEX_marker_genes_list,
    description = "Validate configured UCell GEX marker genes against Cell Ranger feature names",
    command = {
      marker_genes_requested <- aggregation_GEX_marker_genes |>
        unlist(use.names = FALSE) |>
        stringr::str_remove_all("[+-]$") |>
        unique()
      missing_marker_genes <- setdiff(marker_genes_requested, gene_features_df$name)
      if (length(missing_marker_genes) > 0) {
        stop(
          "Configured marker gene(s) are missing from Cell Ranger gene feature names in ",
          targets::tar_name(),
          ": ",
          paste(missing_marker_genes, collapse = ", "),
          ". Update aggregation_GEX_marker_genes in ",
          "cfg_aggregations.yaml",
          " to the symbols used by this Cell Ranger reference.",
          call. = FALSE
        )
      }
      aggregation_GEX_marker_genes
    }
  ),
  targets::tar_target(
    name = GEX_marker_genes_vec,
    description = "Flatten validated GEX marker genes to plain feature names",
    command = {
      UCell_GEX_marker_genes_list |>
        unlist(use.names = FALSE) |>
        stringr::str_remove_all("[+-]$") |>
        unique()
    }
  ),
  targets::tar_target(
    name = reference_Ensembl_annotations_GRanges_list,
    description = "Download Ensembl gene and transcript annotations for this aggregation reference",
    command = with_annotation_hub_cache_lock({
      AnnotationHub::setAnnotationHubOption("CACHE", file.path(targets::tar_config_get("store"), "files", "AnnotationHub"))
      annot_hub_interface <- AnnotationHub::AnnotationHub(ask = FALSE)
      ensembl_db_interface <- annot_hub_interface[[aggregated_cellranger_ref_list$annot_hub_code]] |>
        GenomeInfoDb::`seqlevelsStyle<-`("UCSC")

      tx_grangeslist <- S4Vectors::SimpleList(
        TSS = ensembl_db_interface |> ensembldb::promoters(upstream = 100, downstream = 100) |> GenomicRanges::trim(),
        promoter = ensembl_db_interface |> ensembldb::promoters(upstream = 1500, downstream = 500) |> GenomicRanges::trim(),
        fiveUTR = ensembl_db_interface |> ensembldb::fiveUTRsByTranscript() |> unlist(),
        threeUTR = ensembl_db_interface |> ensembldb::threeUTRsByTranscript() |> unlist(),
        CDS = ensembl_db_interface |> ensembldb::cdsBy() |> unlist(),
        intron = ensembl_db_interface |> ensembldb::intronsByTranscript() |> unlist(),
        proximal = ensembl_db_interface |> ensembldb::promoters(upstream = 1E4L, downstream = 1E4L) |> GenomicRanges::trim(),
        distal = ensembl_db_interface |> ensembldb::promoters(upstream = 1E5L, downstream = 1E5L) |> GenomicRanges::trim()
      ) |>
        S4Vectors::endoapply(GenomicRanges::granges)

      gene_biotypes_to_keep <- c(
        "protein_coding",
        "lncRNA",
        "processed_pseudogene",
        "IG_C_gene",
        "IG_D_gene",
        "IG_J_gene",
        "IG_V_gene",
        "TR_C_gene",
        "TR_D_gene",
        "TR_J_gene",
        "TR_V_gene"
      )

      genes_granges <- ensembl_db_interface |>
        ensembldb::genes() |>
        GenomeInfoDb::keepStandardChromosomes(pruning.mode = "coarse") |>
        subset(gene_biotype %in% gene_biotypes_to_keep)
      genes_granges <- genes_granges[!duplicated(genes_granges$gene_name)]
      names(genes_granges) <- genes_granges$gene_name

      list(genes = genes_granges, transcripts = tx_grangeslist)
    }),
    resources = get_tar_resources(cores_req = 1, RAM_GB_req = 16)
  ),
  targets::tar_target(
    name = marker_validated_Ensembl_annotations_GRanges_list,
    description = "Validate configured marker genes against Ensembl annotation names",
    command = {
      gene_GRanges <- reference_Ensembl_annotations_GRanges_list$genes
      missing_marker_genes <- setdiff(GEX_marker_genes_vec, gene_GRanges$gene_name)
      if (length(missing_marker_genes) > 0) {
        stop(
          "Configured marker gene(s) are missing from Ensembl annotation names in ",
          targets::tar_name(),
          ": ",
          paste(missing_marker_genes, collapse = ", "),
          ". Update aggregation_GEX_marker_genes in ",
          "cfg_aggregations.yaml",
          " to the symbols used by this reference, or update the reference/annotation mapping.",
          call. = FALSE
        )
      }
      reference_Ensembl_annotations_GRanges_list
    }
  ),
  targets::tar_target(
    name = GEX_cellranger_kept_metadata_tibble,
    description = "Combine CellRanger-kept metadata for cells that pass per-reaction QC in this aggregation",
    command = {
      excluded_barcodes_vec <- unique(unlist(
        aggregation_excluded_cellranger_only_barcodes_by_type_list_syms
      ))
      cellranger_kept_metadata_tibble <- dplyr::bind_rows(
        aggregation_cellranger_kept_metadata_tibble_syms
      ) |>
        dplyr::filter(!.data$barcode_w_prefix %in% excluded_barcodes_vec)

      duplicate_barcodes_vec <- cellranger_kept_metadata_tibble |>
        dplyr::count(.data$barcode_w_prefix) |>
        dplyr::filter(.data$n > 1) |>
        dplyr::pull(.data$barcode_w_prefix)

      if (length(duplicate_barcodes_vec) > 0) {
        stop(
          "Duplicate barcode_w_prefix values in aggregation '",
          aggregation,
          "': ",
          paste(utils::head(duplicate_barcodes_vec, 20), collapse = ", ")
        )
      }

      cellranger_kept_metadata_tibble
    }
  ),
  targets::tar_target(
    name = organism_chr,
    description = "Extract organism string (Homo_sapiens or Mus_musculus) from the CellRanger reference",
    command = aggregated_cellranger_ref_list$organism
  ),
  tarchetypes::tar_file(
    name = donor_id_metadata_tsv,
    description = "Track the donor ID metadata TSV file for change detection",
    command = aggregation_donor_id_metadata_tsv,
    deployment = "main"
  ),
  tarchetypes::tar_file(
    name = reaction_ID_metadata_tsv,
    description = "Track the reaction ID metadata TSV file for change detection",
    command = aggregation_reaction_ID_metadata_tsv,
    deployment = "main"
  ),
  targets::tar_target(
    name = donor_id_metadata_tibble,
    description = "Read the donor ID metadata TSV into a tibble with donor_id coerced to character",
    command = dplyr::mutate(
      readr::read_tsv(donor_id_metadata_tsv),
      donor_id = as.character(donor_id)
    )
  ),
  targets::tar_target(
    name = reaction_ID_metadata_tibble,
    description = "Read the reaction ID metadata TSV into a tibble with TENX_reaction_ID coerced to character",
    command = dplyr::mutate(
      readr::read_tsv(reaction_ID_metadata_tsv),
      TENX_reaction_ID = as.character(TENX_reaction_ID)
    )
  ),
  targets::tar_target(
    name = interesting_genes,
    description = "Collect unique gene names from marker annotations and other configured genes",
    command = unique(unname(c(
      unlist(aggregation_GEX_marker_genes),
      aggregation_other_interesting_genes
    )))
  ),
  tarchetypes::tar_file(
    name = nuclei_per_donor_id_bars,
    description = "Plot bar chart of nuclei counts per donor ID and save to file",
    command = metadata_w_cell_types_tibble.WNN |>
      plot_nuclei_per_donor_id() |>
      save_plots_structured(width = 10)
  ),
  tarchetypes::tar_file(
    name = pseudobulk_counts_BPCells_matrix_dir.GEX,
    description = "Write pseudobulk GEX counts per WNN cell-type-donor combination to BPCells",
    command = {
      out_dir <- get_structured_file_path()
      get_BPCells_pseudobulk_matrix(
        feature_matrix = aggregated_counts_BPCells_matrix.GEX,
        metadata_tibble = metadata_w_cell_types_tibble.WNN,
        cluster_col = "PCA_harmony_SNN_cluster_cell_type",
        threads = 15
      ) |>
        methods::as("dgCMatrix") |>
        BPCells::write_matrix_dir(out_dir, overwrite = TRUE)
      out_dir
    },
    resources = get_tar_resources(cores_req = 15, RAM_GB_req = 60)
  ),
  targets::tar_target(
    name = pseudobulk_counts_matrix.GEX,
    description = "Open BPCells-backed pseudobulk GEX counts per WNN cell-type-donor combination",
    command = BPCells::open_matrix_dir(pseudobulk_counts_BPCells_matrix_dir.GEX)
  ),
  targets::tar_target(
    name = pseudobulk_depth_tibble.GEX,
    description = "Compute pseudobulk GEX count depths per WNN cell-type-donor combination",
    command = get_pseudobulk_depth_tibble(pseudobulk_counts_matrix.GEX) |>
      dplyr::mutate(modality = "GEX")
  ),
  tarchetypes::tar_file(
    name = pseudobulk_counts_BPCells_matrix_dir.ATAC,
    description = "Write pseudobulk ATAC counts per WNN cell-type-donor combination to BPCells",
    command = {
      out_dir <- get_structured_file_path()
      get_BPCells_pseudobulk_matrix(
        feature_matrix = peak_QC_filtered_BPCells_matrix.ATAC,
        metadata_tibble = metadata_w_cell_types_tibble.WNN,
        cluster_col = "PCA_harmony_SNN_cluster_cell_type",
        threads = 15
      ) |>
        methods::as("dgCMatrix") |>
        BPCells::write_matrix_dir(out_dir, overwrite = TRUE)
      out_dir
    },
    resources = get_tar_resources(cores_req = 15, RAM_GB_req = 60)
  ),
  targets::tar_target(
    name = pseudobulk_counts_matrix.ATAC,
    description = "Open BPCells-backed pseudobulk ATAC counts per WNN cell-type-donor combination",
    command = BPCells::open_matrix_dir(pseudobulk_counts_BPCells_matrix_dir.ATAC)
  ),
  targets::tar_target(
    name = pseudobulk_depth_tibble.ATAC,
    description = "Compute pseudobulk ATAC count depths per WNN cell-type-donor combination",
    command = get_pseudobulk_depth_tibble(pseudobulk_counts_matrix.ATAC) |>
      dplyr::mutate(modality = "ATAC")
  ),
  targets::tar_target(
    name = pseudobulk_activity_matrix.TFA,
    description = "Compute pseudobulk TF activity scores from BPCells-backed pseudobulk ATAC counts",
    command = get_pseudobulk_activity_matrix.TFA(
      psbulk_ATAC_data_matrix = pseudobulk_counts_matrix.ATAC,
      chromVAR_obj = chromVAR_obj.ATAC,
      chromVAR_motif_matrix = chromVAR_TF_motif_matrix.ATAC
    ),
    resources = get_tar_resources(RAM_GB_req = 60)
  )
)
