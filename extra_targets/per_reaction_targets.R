rlang::list2(
  tarchetypes::tar_file(
    name = cellranger_summary_file,
    description = "Track the CellRanger summary CSV file as a sentinel for the completed reaction output directory",
    command = file.path(reaction_cellranger_count_dir, "outs", "summary.csv")
  ),
  tarchetypes::tar_file(
    name = cellranger_barcodes_tsv,
    description = "Write CellRanger-native barcodes to a TSV file for use by cellsnp-lite and vireo [part_of_graph:parallel] [part_of_graph:seurat_export]",
    command = {
      cellranger_h5_file <- file.path(dirname(cellranger_summary_file), "filtered_feature_bc_matrix.h5")
      cellranger_h5_file_con <- hdf5r::H5File$new(cellranger_h5_file, mode = "r")
      on.exit(cellranger_h5_file_con$close_all())
      out_file <- get_structured_file_path(filetype = "tsv")
      writeLines(cellranger_h5_file_con[["matrix/barcodes"]][], con = out_file)
      out_file
    }
  ),
  targets::tar_target(
    name = gene_features_df,
    description = "Read and format CellRanger gene feature metadata",
    command = {
      cellranger_h5_file <- file.path(dirname(cellranger_summary_file), "filtered_feature_bc_matrix.h5")
      cellranger_h5_file_con <- hdf5r::H5File$new(cellranger_h5_file, mode = "r")
      on.exit(cellranger_h5_file_con$close_all())

      gene_features_df <- tibble::tibble(
        id = cellranger_h5_file_con[["matrix/features/id"]][],
        name = cellranger_h5_file_con[["matrix/features/name"]][],
        feature_type = cellranger_h5_file_con[["matrix/features/feature_type"]][],
        genome = cellranger_h5_file_con[["matrix/features/genome"]][],
        interval = cellranger_h5_file_con[["matrix/features/interval"]][]
      ) |>
        dplyr::filter(.data$feature_type == "Gene Expression") |>
        dplyr::mutate(
          interval = dplyr::case_when(
            stringr::str_starts(.data$name, "(?i)^MT-") & .data$interval == "NA" ~ "chrM:NA-NA",
            .default = .data$interval
          )
        ) |>
        tidyr::separate_wider_delim(cols = "interval", delim = ":", names = c("seqnames", "interval")) |>
        tidyr::separate_wider_delim(cols = "interval", delim = "-", names = c("start", "end")) |>
        dplyr::mutate(
          start = dplyr::na_if(.data$start, "NA") |> as.numeric(),
          end = dplyr::na_if(.data$end, "NA") |> as.numeric()
        ) |>
        dplyr::transmute(
          gene_name_unique = make.unique(.data$name),
          seqnames,
          start,
          end,
          id,
          name,
          genome,
          feature_type
        ) |>
        as.data.frame()

      rownames(gene_features_df) <- gene_features_df$gene_name_unique
      gene_features_df
    }
  ),
  targets::tar_target(
    name = cellranger_ref_list,
    description = "Load CellRanger reference JSON and annotate with Gencode version and AnnotationHub code",
    command = {
      cellranger_ref_list <- jsonlite::read_json(fs::path(dataset_cellranger_arc_refdata_dir, "reference.json"))
      cellranger_ref_list$gencode_version <- stringr::str_remove(
        cellranger_ref_list[["input_gtf_files"]][[1]],
        ".primary_assembly.*"
      )
      cellranger_ref_list$annot_hub_code <- c(
        "gencode.v44" = "AH113665",
        "gencode.vM33" = "AH113713",
        "gencode.v32" = "AH75011",
        "gencode.vM23" = "AH75036"
      )[[cellranger_ref_list$gencode_version]]
      if (is.null(cellranger_ref_list$annot_hub_code)) {
        stop("Unsupported Gencode version: ", cellranger_ref_list$gencode_version)
      }
      cellranger_ref_list
    },
  ),
  tarchetypes::tar_file(
    name = cellbender_h5_file,
    description = "Track the CellBender output h5 file for this reaction, or no files if not configured [part_of_graph:parallel] [part_of_graph:seurat_export]",
    command = if (is.na(reaction_cellbender_h5_file) || reaction_cellbender_h5_file == "") {
      character(0)
    } else {
      reaction_cellbender_h5_file
    }
  ),
  tarchetypes::tar_file(
    name = GEX_counts_BPcells_matrix_dir,
    description = "Convert the CellBender or CellRanger GEX h5 matrix to a BPCells directory with reaction-prefixed barcodes [part_of_graph:parallel] [part_of_graph:seurat_export]",
    command = {
      out_dir <- get_structured_file_path()

      gene_expression_matrix <-
        if (length(cellbender_h5_file) > 0L) {
          read_cellbender_h5_matrix(cellbender_h5_file, feature_type = "Gene Expression")
        } else {
          cellranger_h5_file <- file.path(dirname(cellranger_summary_file), "filtered_feature_bc_matrix.h5")
          BPCells::open_matrix_10x_hdf5(cellranger_h5_file, feature_type = "Gene Expression")
        }

      source_label <- if (length(cellbender_h5_file) > 0L) "CellBender GEX matrix" else "CellRanger GEX matrix"
      if (nrow(gene_expression_matrix) != nrow(gene_features_df)) {
        stop(
          source_label,
          " row count (",
          nrow(gene_expression_matrix),
          ") does not match gene feature row count (",
          nrow(gene_features_df),
          ")."
        )
      }

      existing_feature_names <- rownames(gene_expression_matrix)
      has_feature_names <- !is.null(existing_feature_names) &&
        length(existing_feature_names) == nrow(gene_expression_matrix) &&
        !all(is.na(existing_feature_names) | existing_feature_names == "")

      if (
        has_feature_names &&
          !identical(existing_feature_names, gene_features_df$gene_name_unique) &&
          !identical(existing_feature_names, gene_features_df$name) &&
          !identical(existing_feature_names, gene_features_df$id)
      ) {
        stop(source_label, " feature names do not match CellRanger gene feature order.")
      }

      rownames(gene_expression_matrix) <- gene_features_df$gene_name_unique
      colnames(gene_expression_matrix) <- paste0(reaction_ID, "_", colnames(gene_expression_matrix))

      BPCells::write_matrix_dir(gene_expression_matrix, out_dir, overwrite = TRUE)
      out_dir
    }
  ),
  targets::tar_target(
    name = GEX_counts_BPCells_matrix,
    description = "Open the matrix from the GEX BPCells directory",
    command = BPCells::open_matrix_dir(GEX_counts_BPcells_matrix_dir)
  ),
  tarchetypes::tar_file(
    name = fragments_w_prefix_bpcells_dir,
    description = "Convert the CellRanger fragment file to a BPCells directory with reaction-prefixed barcodes [part_of_graph:parallel] [part_of_graph:seurat_export]",
    command = {
      out_dir <- get_structured_file_path()

      fragment_file <- file.path(dirname(cellranger_summary_file), "atac_fragments.tsv.gz")
      fragments <- BPCells::open_fragments_10x(fragment_file)
      fragments <- BPCells::prefix_cell_names(fragments, paste0(reaction_ID, "_"))
      BPCells::write_fragments_dir(fragments, out_dir, overwrite = TRUE)

      out_dir
    }
  ),
  targets::tar_target(
    name = fragments_w_prefix_bpcells,
    description = "Convert the CellRanger fragment file to a BPCells directory with reaction-prefixed barcodes",
    command = BPCells::open_fragments_dir(fragments_w_prefix_bpcells_dir)
  ),
  targets::tar_target(
    name = ATAC_qc_metrics_tibble,
    description = "Compute TSS enrichment and nucleosome signal from BPCells ATAC fragments [part_of_graph:parallel] [part_of_graph:seurat_export]",
    command = {
      barcode_prefix <- paste0(reaction_ID, "_")

      ATAC_qc_metrics_tibble <-
        BPCells::qc_scATAC(
          fragments = fragments_w_prefix_bpcells,
          genes = Ensembl_gene_annotation_GRanges_list[[cellranger_ref_list$annot_hub_code]],
          blacklist = GenomicRanges::GRanges()
        ) |>
        tibble::as_tibble()

      if (!all(startsWith(ATAC_qc_metrics_tibble$cellName, barcode_prefix))) {
        stop("All BPCells ATAC QC cell names must start with reaction prefix: ", barcode_prefix)
      }

      ATAC_qc_metrics_tibble |>
        dplyr::transmute(
          barcode = substring(cellName, nchar(barcode_prefix) + 1),
          TSS.enrichment = TSSEnrichment,
          nucleosome_signal = monoNucleosomal / subNucleosomal
        )
    },
  ),
  targets::tar_target(
    name = GEX_basic_metadata_tibble,
    description = "Compute basic per-cell RNA QC metadata from the BPCells GEX matrix",
    command = {
      barcode_prefix <- paste0(reaction_ID, "_")

      if (!all(startsWith(colnames(GEX_counts_BPCells_matrix), barcode_prefix))) {
        stop("All GEX BPCells matrix column names must start with reaction prefix: ", barcode_prefix)
      }

      col_stats <- BPCells::matrix_stats(GEX_counts_BPCells_matrix, col_stats = "mean")$col_stats
      nCount_RNA <- col_stats["mean", ] * nrow(GEX_counts_BPCells_matrix)
      nFeature_RNA <- col_stats["nonzero", ]

      mitochondrial_features <- stringr::str_detect(rownames(GEX_counts_BPCells_matrix), "(?i)^MT-")
      mitochondrial_counts <- if (any(mitochondrial_features)) {
        colSums(GEX_counts_BPCells_matrix[mitochondrial_features, , drop = FALSE])
      } else {
        purrr::set_names(rep(0, ncol(GEX_counts_BPCells_matrix)), colnames(GEX_counts_BPCells_matrix))
      }

      tibble::tibble(
        orig.ident = reaction_ID,
        barcode = substring(colnames(GEX_counts_BPCells_matrix), nchar(barcode_prefix) + 1),
        nCount_RNA = unname(nCount_RNA),
        nFeature_RNA = unname(nFeature_RNA),
        log10_nCount_RNA = log10(nCount_RNA),
        RNA_mito_percent = unname(mitochondrial_counts / nCount_RNA * 100),
        novelty = log10(nFeature_RNA) / log10(nCount_RNA)
      )
    }
  ),
  targets::tar_target(
    name = per_barcode_metrics_tibble,
    description = "Read per-barcode CellRanger metrics and filter to called cells [part_of_graph:seurat_export]",
    command = {
      per_barcode_metrics_file <- file.path(dirname(cellranger_summary_file), "per_barcode_metrics.csv")
      per_barcode_metrics_file |>
        readr::read_csv(col_types = paste0(c(rep(x = "c", 3), rep("d", 28)), collapse = "")) |>
        dplyr::filter(is_cell == 1)
    }
  ),
  targets::tar_target(
    name = amulet_metrics_tibble,
    description = "Run AMULET doublet detection on ATAC fragments and return per-barcode metrics [part_of_graph:parallel] [part_of_graph:seurat_export]",
    command = {
      if (isTRUE(dataset_run_amulet)) {
        fragment_file <- file.path(dirname(cellranger_summary_file), "atac_fragments.tsv.gz")
        fragment_file |>
          scDblFinder::amulet(
            barcodes = cellranger_barcodes_tsv,
            minFrags = 1000,
            BPPARAM = BiocParallel::MulticoreParam(workers = 6)
          ) |>
          dplyr::rename_with(~ paste0("amulet_", .x), .cols = 2:dplyr::last_col()) |>
          tibble::rownames_to_column("barcode")
      } else {
        tibble::tibble(
          barcode = readLines(cellranger_barcodes_tsv),
          amulet_nFrags = NA_integer_,
          amulet_uniqFrags = NA_integer_,
          amulet_nAbove2 = NA_integer_,
          amulet_total.nAbove2 = NA_integer_,
          amulet_p.value = NA_real_,
          amulet_q.value = NA_real_
        )
      }
    },
    resources = get_tar_resources(RAM_GB_req = 60, cores_req = 6) # needs quite a lot of memory - consider decreasing workers if this is a problem.
  ),
  tarchetypes::tar_file(
    name = cellsnp_dir,
    description = "Run cellsnp-lite to genotype each barcode at SNP positions, or return no files if no VCF is configured [part_of_graph:parallel] [part_of_graph:seurat_export]",
    command = if (!is.na(reaction_donors_VCF_file)) {
      ATAC_bam_file <- file.path(dirname(cellranger_summary_file), "atac_possorted_bam.bam")
      get_cellsnp_dir(
        bam_file = ATAC_bam_file,
        cellranger_barcodes_tsv = cellranger_barcodes_tsv,
        reaction_donors_VCF_file = reaction_donors_VCF_file,
        cores = 6
      )
    } else {
      character(0)
    },
    resources = get_tar_resources(cores_req = 6)
  ),
  targets::tar_target(
    name = vireo_donor_ids_tibble,
    description = "Run vireo to demultiplex donors from cellsnp-lite genotypes and return per-barcode donor assignments [part_of_graph:parallel] [part_of_graph:seurat_export]",
    command = get_vireo_donor_ids_tibble(
      cellsnp_dir,
      reaction_donors_VCF_file = reaction_donors_VCF_file,
      reaction_n_donors = reaction_n_donors,
      reaction_donor_id = reaction_donor_id,
      cellranger_barcodes_tsv = cellranger_barcodes_tsv,
      cores = 4 # consider decreasing this if increasing memory from 60 > 120 GB doesn't fix the issue.
    ),
    resources = get_tar_resources(cores_req = 6, RAM_GB_req = 60)
  ),
  targets::tar_target(
    name = unfiltered_cells_n_vecs,
    description = "Count barcodes in the full metadata table used by per-reaction QC exclusion filters",
    command = nrow(full_metadata_tibble)
  ),
  targets::tar_target(
    name = full_metadata_tibble,
    description = "Join GEX metadata with ATAC QC, vireo, per-barcode, and AMULET metrics into a comprehensive per-barcode metadata tibble [part_of_graph:parallel] [part_of_graph:seurat_export]",
    command = GEX_basic_metadata_tibble |>
      dplyr::full_join(ATAC_qc_metrics_tibble, by = "barcode") |>
      dplyr::full_join(vireo_donor_ids_tibble, by = "barcode") |>
      dplyr::full_join(per_barcode_metrics_tibble, by = "barcode") |>
      dplyr::full_join(amulet_metrics_tibble, by = "barcode") |>
      dplyr::mutate(
        discarded_by_cellranger = is.na(is_cell),
        discarded_by_amulet = if (isTRUE(dataset_run_amulet)) is.na(amulet_uniqFrags) else FALSE,
        not_found_in_GEX_matrix = is.na(orig.ident),
        TENX_reaction_ID = reaction_ID,
        dataset = dataset,
        barcode_w_prefix = paste0(reaction_ID, "_", barcode),
        vireo_type = tidyr::replace_na(vireo_type, "discarded_by_cellranger"),
        ATAC_TSS_fragments_frac = atac_TSS_fragments / atac_fragments,
        ATAC_peak_region_frac = atac_peak_region_fragments / atac_fragments,
        ATAC_peak_region_cutsites_frac = atac_peak_region_cutsites / atac_fragments,
        GEX_exonic_to_intronic_umis_frac = gex_exonic_umis / gex_intronic_umis
      ),
    resources = get_tar_resources(RAM_GB_req = 60)
  ),
  targets::tar_target(
    name = excluded_barcodes_by_type_list,
    description = "Identify barcodes failing per-reaction QC thresholds and group them by exclusion reason",
    command = get_excluded_BCs(full_metadata_tibble, QC_exclude_vector = dataset_QC_exclude_list_per_reaction)
  ),
  tarchetypes::tar_file(
    name = excluded_barcodes_by_type_upset,
    description = "Plot an UpSet plot of per-reaction QC exclusion overlaps and save to file",
    command = plot_upset_from_excluded_BCs_list(
      QC_excluded_BCs_list = excluded_barcodes_by_type_list,
      n_total = unfiltered_cells_n_vecs
    ) |>
      save_plots_structured(width = 10, height = 10)
  ),
  targets::tar_target(
    name = excluded_cellranger_only_barcodes_by_type_list,
    description = "Identify CellRanger-only barcodes failing per-reaction QC and group by exclusion reason [part_of_graph:parallel] [part_of_graph:seurat_export]",
    command = full_metadata_tibble |>
      dplyr::filter(!discarded_by_cellranger) |>
      get_excluded_BCs(QC_exclude_vector = dataset_QC_exclude_list_per_reaction)
  ),
  tarchetypes::tar_file(
    name = excluded_cellranger_only_barcodes_by_type_upset,
    description = "Plot an UpSet plot of CellRanger-only QC exclusion overlaps and save to file",
    command = plot_upset_from_excluded_BCs_list(
      QC_excluded_BCs_list = excluded_cellranger_only_barcodes_by_type_list,
      n_total = nrow(cellranger_kept_metadata_tibble)
    ) |>
      save_plots_structured(width = 10, height = 10)
  ),
  targets::tar_target(
    name = cellranger_kept_metadata_tibble,
    description = "Filter metadata to CellRanger-called cells",
    command = full_metadata_tibble |>
      dplyr::filter_out(
        discarded_by_cellranger
      )
  )
)
