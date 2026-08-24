cellranger_arc_subset_feature_metadata <- function(matrix_h5) {
  h5_file <- hdf5r::H5File$new(matrix_h5, mode = "r")
  on.exit(h5_file$close_all())
  feature_group <- h5_file[["matrix/features"]]
  required_datasets <- c("id", "name", "feature_type")
  missing_datasets <- setdiff(required_datasets, names(feature_group))
  if (length(missing_datasets)) {
    stop(
      "10x feature metadata are missing dataset(s): ",
      paste(missing_datasets, collapse = ", "),
      call. = FALSE
    )
  }

  optional_datasets <- setdiff(names(feature_group), c("_all_tag_keys", required_datasets))
  list(
    feature_ids = feature_group[["id"]][],
    feature_names = feature_group[["name"]][],
    feature_types = feature_group[["feature_type"]][],
    feature_metadata = stats::setNames(
      lapply(optional_datasets, \(dataset) feature_group[[dataset]][]),
      optional_datasets
    )
  )
}

cellranger_arc_subset_select_barcodes <- function(count_dir, n_cells) {
  matrix_h5 <- file.path(count_dir, "outs", "filtered_feature_bc_matrix.h5")
  metrics_file <- file.path(count_dir, "outs", "per_barcode_metrics.csv")
  matrix_barcodes <- colnames(BPCells::open_matrix_10x_hdf5(matrix_h5))
  metrics <- readr::read_csv(metrics_file, show_col_types = FALSE)
  required_columns <- c("barcode", "is_cell", "atac_fragments", "gex_umis_count")
  missing_columns <- setdiff(required_columns, names(metrics))
  if (length(missing_columns)) {
    stop(
      "per_barcode_metrics.csv is missing column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  candidates <- metrics |>
    dplyr::filter(.data$is_cell == 1, .data$barcode %in% matrix_barcodes) |>
    dplyr::arrange(
      dplyr::desc(.data$atac_fragments),
      dplyr::desc(.data$gex_umis_count),
      .data$barcode
    )
  if (nrow(candidates) < n_cells) {
    stop(
      "Requested ", n_cells, " cells, but only ", nrow(candidates),
      " filtered matrix barcodes have Cell Ranger cell metrics.",
      call. = FALSE
    )
  }
  candidates$barcode[seq_len(n_cells)]
}

cellranger_arc_subset_write_matrix_h5 <- function(input_h5, output_h5, barcodes) {
  matrix <- BPCells::open_matrix_10x_hdf5(input_h5)
  barcode_indices <- match(barcodes, colnames(matrix))
  if (anyNA(barcode_indices)) {
    stop(
      "Matrix is missing selected barcode(s): ",
      paste(barcodes[is.na(barcode_indices)], collapse = ", "),
      call. = FALSE
    )
  }
  features <- cellranger_arc_subset_feature_metadata(input_h5)
  BPCells::write_matrix_10x_hdf5(
    mat = matrix[, barcode_indices, drop = FALSE],
    path = output_h5,
    barcodes = barcodes,
    feature_ids = features$feature_ids,
    feature_names = features$feature_names,
    feature_types = features$feature_types,
    feature_metadata = features$feature_metadata,
    gzip_level = 1L
  )
  invisible(output_h5)
}

cellranger_arc_subset_write_matrix_market <- function(matrix_h5, output_dir) {
  matrix <- BPCells::open_matrix_10x_hdf5(matrix_h5)
  features <- cellranger_arc_subset_feature_metadata(matrix_h5)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  matrix_file <- file.path(output_dir, "matrix.mtx")
  Matrix::writeMM(methods::as(matrix, "dgCMatrix"), matrix_file)
  R.utils::gzip(matrix_file, overwrite = FALSE, remove = TRUE)

  barcodes_connection <- gzfile(file.path(output_dir, "barcodes.tsv.gz"), open = "wt")
  writeLines(colnames(matrix), barcodes_connection)
  close(barcodes_connection)

  feature_tibble <- tibble::tibble(
    id = features$feature_ids,
    name = features$feature_names,
    feature_type = features$feature_types
  )
  features_connection <- gzfile(file.path(output_dir, "features.tsv.gz"), open = "wt")
  utils::write.table(
    feature_tibble,
    file = features_connection,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  close(features_connection)
  invisible(output_dir)
}

cellranger_arc_subset_write_fragments <- function(input_file, output_file, barcode_file) {
  bgzip <- Sys.which("bgzip")
  if (!nzchar(bgzip)) {
    stop("bgzip is required to write the subset fragment file.", call. = FALSE)
  }
  awk_program <- paste0(
    "BEGIN { FS = OFS = \"\\t\" } ",
    "NR == FNR { keep[$1] = 1; next } ($4 in keep)"
  )
  command <- paste(
    "gzip -dc", shQuote(input_file), "|",
    "awk", shQuote(awk_program), shQuote(barcode_file), "-", "|",
    shQuote(bgzip), "-c", ">", shQuote(output_file)
  )
  status <- system(paste("bash -o pipefail -c", shQuote(command)))
  if (!identical(status, 0L) || !file.exists(output_file) || file.size(output_file) == 0) {
    stop("Failed to create the subset fragment file.", call. = FALSE)
  }
  Rsamtools::indexTabix(output_file, format = "bed")
  invisible(output_file)
}

cellranger_arc_subset_write_bam <- function(input_bam, output_bam, barcodes) {
  input_bam_file <- Rsamtools::BamFile(
    input_bam,
    index = paste0(input_bam, ".bai"),
    yieldSize = 1000000L
  )
  parameter <- Rsamtools::ScanBamParam(tagFilter = list(CB = barcodes))
  Rsamtools::filterBam(
    file = input_bam_file,
    destination = output_bam,
    indexDestination = TRUE,
    param = parameter
  )
  invisible(output_bam)
}

cellranger_arc_subset_manifest <- function(input_dir, include_bams) {
  input_files <- fs::dir_ls(input_dir, recurse = TRUE, type = "file", fail = TRUE)
  relative_path <- fs::path_rel(input_files, start = input_dir)
  subset_paths <- c(
    "outs/filtered_feature_bc_matrix.h5",
    "outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
    "outs/filtered_feature_bc_matrix/features.tsv.gz",
    "outs/filtered_feature_bc_matrix/matrix.mtx.gz",
    "outs/raw_feature_bc_matrix.h5",
    "outs/raw_feature_bc_matrix/barcodes.tsv.gz",
    "outs/raw_feature_bc_matrix/features.tsv.gz",
    "outs/raw_feature_bc_matrix/matrix.mtx.gz",
    "outs/atac_fragments.tsv.gz",
    "outs/atac_fragments.tsv.gz.tbi",
    "outs/per_barcode_metrics.csv"
  )
  copied_paths <- c("outs/summary.csv", "outs/atac_peaks.bed")
  bam_paths <- c(
    "outs/atac_possorted_bam.bam",
    "outs/atac_possorted_bam.bam.bai",
    "outs/gex_possorted_bam.bam",
    "outs/gex_possorted_bam.bam.bai"
  )
  action <- dplyr::case_when(
    relative_path %in% subset_paths ~ "subset",
    relative_path %in% copied_paths ~ "copied_not_recomputed",
    relative_path %in% bam_paths & include_bams ~ "subset",
    relative_path %in% bam_paths ~ "omitted_use_include_bams",
    TRUE ~ "omitted_derived_or_runtime_artifact"
  )
  tibble::tibble(
    relative_path = relative_path,
    input_bytes = as.numeric(fs::file_size(input_files)),
    action = action
  ) |>
    dplyr::arrange(.data$relative_path)
}

create_cellranger_arc_count_subset <- function(
  input_dir,
  output_dir = NULL,
  n_cells = 300L,
  include_bams = FALSE
) {
  input_dir <- fs::path_abs(input_dir)
  if (is.null(output_dir)) {
    output_dir <- paste0(input_dir, "_subset")
  }
  output_dir <- fs::path_abs(output_dir)
  if (!dir.exists(input_dir)) {
    stop("Input Cell Ranger ARC count directory does not exist: ", input_dir, call. = FALSE)
  }
  if (dir.exists(output_dir) || file.exists(output_dir)) {
    stop("Output path already exists: ", output_dir, call. = FALSE)
  }
  if (
    !is.numeric(n_cells) || length(n_cells) != 1L || is.na(n_cells) ||
      n_cells < 2 || n_cells != as.integer(n_cells)
  ) {
    stop("n_cells must be one integer of at least 2.", call. = FALSE)
  }
  n_cells <- as.integer(n_cells)

  required_files <- file.path(
    input_dir,
    c(
      "outs/summary.csv",
      "outs/filtered_feature_bc_matrix.h5",
      "outs/raw_feature_bc_matrix.h5",
      "outs/atac_fragments.tsv.gz",
      "outs/atac_fragments.tsv.gz.tbi",
      "outs/atac_peaks.bed",
      "outs/per_barcode_metrics.csv"
    )
  )
  if (isTRUE(include_bams)) {
    required_files <- c(
      required_files,
      file.path(
        input_dir,
        "outs",
        c(
          "atac_possorted_bam.bam",
          "atac_possorted_bam.bam.bai",
          "gex_possorted_bam.bam",
          "gex_possorted_bam.bam.bai"
        )
      )
    )
  }
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files)) {
    stop("Required input file(s) are missing: ", paste(missing_files, collapse = ", "), call. = FALSE)
  }

  staging_dir <- paste0(output_dir, ".partial-", Sys.getpid())
  dir.create(file.path(staging_dir, "outs"), recursive = TRUE)
  keep_staging <- FALSE
  on.exit({
    if (!keep_staging && dir.exists(staging_dir)) {
      unlink(staging_dir, recursive = TRUE)
    }
  }, add = TRUE)

  barcodes <- cellranger_arc_subset_select_barcodes(input_dir, n_cells)
  barcode_file <- file.path(staging_dir, "selected_barcodes.tsv")
  writeLines(barcodes, barcode_file)

  for (matrix_name in c("filtered_feature_bc_matrix", "raw_feature_bc_matrix")) {
    input_h5 <- file.path(input_dir, "outs", paste0(matrix_name, ".h5"))
    output_h5 <- file.path(staging_dir, "outs", paste0(matrix_name, ".h5"))
    cellranger_arc_subset_write_matrix_h5(input_h5, output_h5, barcodes)
    cellranger_arc_subset_write_matrix_market(
      output_h5,
      file.path(staging_dir, "outs", matrix_name)
    )
  }

  metrics <- readr::read_csv(
    file.path(input_dir, "outs", "per_barcode_metrics.csv"),
    show_col_types = FALSE
  ) |>
    dplyr::filter(.data$barcode %in% barcodes) |>
    dplyr::slice(match(barcodes, .data$barcode))
  if (nrow(metrics) != n_cells || anyNA(metrics$barcode)) {
    stop("Subset per-barcode metrics do not cover every selected barcode.", call. = FALSE)
  }
  readr::write_csv(metrics, file.path(staging_dir, "outs", "per_barcode_metrics.csv"))

  cellranger_arc_subset_write_fragments(
    input_file = file.path(input_dir, "outs", "atac_fragments.tsv.gz"),
    output_file = file.path(staging_dir, "outs", "atac_fragments.tsv.gz"),
    barcode_file = barcode_file
  )

  copied_files <- c("summary.csv", "atac_peaks.bed")
  copied <- file.copy(
    from = file.path(input_dir, "outs", copied_files),
    to = file.path(staging_dir, "outs", copied_files),
    overwrite = FALSE
  )
  if (!all(copied)) {
    stop("Failed to copy Cell Ranger summary or peak definition.", call. = FALSE)
  }

  if (isTRUE(include_bams)) {
    for (bam_name in c("atac_possorted_bam.bam", "gex_possorted_bam.bam")) {
      cellranger_arc_subset_write_bam(
        input_bam = file.path(input_dir, "outs", bam_name),
        output_bam = file.path(staging_dir, "outs", bam_name),
        barcodes = barcodes
      )
    }
  }

  manifest <- cellranger_arc_subset_manifest(input_dir, include_bams = include_bams)
  readr::write_tsv(manifest, file.path(staging_dir, "SUBSET_MANIFEST.tsv"))
  metadata <- tibble::tibble(
    key = c(
      "input_dir",
      "output_dir",
      "created_at",
      "n_cells",
      "selection",
      "include_bams",
      "barcode_xxhash64"
    ),
    value = c(
      input_dir,
      output_dir,
      format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      n_cells,
      "highest Cell Ranger ATAC fragment count; ties by GEX UMI count and barcode",
      include_bams,
      digest::digest(barcodes, algo = "xxhash64")
    )
  )
  readr::write_tsv(metadata, file.path(staging_dir, "SUBSET_METADATA.tsv"))

  validation <- Seurat::Read10X_h5(
    file.path(staging_dir, "outs", "filtered_feature_bc_matrix.h5"),
    use.names = TRUE
  )
  if (
    !identical(sort(names(validation)), sort(c("Gene Expression", "Peaks"))) ||
      any(vapply(validation, ncol, integer(1)) != n_cells) ||
      !file.exists(file.path(staging_dir, "outs", "atac_fragments.tsv.gz.tbi"))
  ) {
    stop("Generated subset failed matrix or fragment-index validation.", call. = FALSE)
  }

  if (!file.rename(staging_dir, output_dir)) {
    stop("Failed to atomically move the generated subset into place.", call. = FALSE)
  }
  keep_staging <- TRUE
  message("Created Cell Ranger ARC count subset: ", output_dir)
  message("Selected cells: ", n_cells)
  invisible(output_dir)
}
