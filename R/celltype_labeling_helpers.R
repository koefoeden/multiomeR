get_chromVAR_genome_obj <- function(genome) {
  switch(
    genome,
    "GRCh38" = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
    "mm10" = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10,
    "GRCm39" = BSgenome.Mmusculus.UCSC.mm39::BSgenome.Mmusculus.UCSC.mm39,
    stop(paste("Unsupported genome:", genome))
  )
}

get_peak_names_from_GRanges <- function(peak_ranges) {
  peak_names <- names(peak_ranges)
  if (is.null(peak_names) || anyNA(peak_names) || any(peak_names == "")) {
    peak_names <- stringr::str_c(
      as.character(GenomicRanges::seqnames(peak_ranges)),
      GenomicRanges::start(peak_ranges),
      GenomicRanges::end(peak_ranges),
      sep = "-"
    )
  }
  peak_names
}

row_sums_for_peak_count_matrix <- function(counts_matrix) {
  if (inherits(counts_matrix, "IterableMatrix")) {
    BPCells::rowSums(counts_matrix)
  } else {
    Matrix::rowSums(counts_matrix)
  }
}

get_motif_matrix_from_peak_ranges <- function(peak_ranges, TF_motif_matrix_list, genome_obj) {
  motif_matches <- motifmatchr::matchMotifs(
    pwms = TF_motif_matrix_list,
    subject = peak_ranges,
    genome = genome_obj
  )
  motif_matrix <- SummarizedExperiment::assay(motif_matches, "motifMatches")
  rownames(motif_matrix) <- get_peak_names_from_GRanges(peak_ranges)
  motif_matrix
}

get_peak_ranges_for_peak_names <- function(peak_names, ATAC_peak_GRanges, genome_obj) {
  peak_range_idx <- match(peak_names, get_peak_names_from_GRanges(ATAC_peak_GRanges))
  if (anyNA(peak_range_idx)) {
    stop("Some peak matrix rows were not found in ATAC_peak_GRanges.")
  }

  peak_ranges <- ATAC_peak_GRanges[peak_range_idx]
  GenomeInfoDb::seqinfo(peak_ranges) <- GenomeInfoDb::seqinfo(genome_obj)[GenomeInfoDb::seqlevels(peak_ranges)]
  peak_ranges <- GenomicRanges::trim(peak_ranges)
  names(peak_ranges) <- peak_names
  peak_ranges
}

#' Get motif matrix from ATAC peak names
#'
#' Extract motif annotations for named ATAC peaks after resolving their ranges.
#'
#' @param ATAC_peak_names Character vector of peak names to extract from
#'   `ATAC_peak_GRanges`.
#' @param ATAC_peak_GRanges GRanges for ATAC peaks, aligned by peak name to the corresponding ATAC matrix rows.
#' @param TF_motif_matrix_list Motif annotation resources used by
#'   `get_motif_matrix_from_peak_ranges()`.
#' @param genome Genome build key used to choose chromosome sizes, blacklist resources, and external-tool parameters.
#' @return Peak-by-motif matrix for the requested peak names.
#' @keywords internal

get_motif_matrix_from_ATAC_peak_names <- function(ATAC_peak_names,
                                                  ATAC_peak_GRanges,
                                                  TF_motif_matrix_list,
                                                  genome = "GRCh38") {
  genome_obj <- get_chromVAR_genome_obj(genome)
  peak_ranges <- get_peak_ranges_for_peak_names(
    peak_names = ATAC_peak_names,
    ATAC_peak_GRanges = ATAC_peak_GRanges,
    genome_obj = genome_obj
  )

  motif_matrix <- get_motif_matrix_from_peak_ranges(peak_ranges, TF_motif_matrix_list, genome_obj)
  motif_matrix[ATAC_peak_names, , drop = FALSE]
}

#' Split cols by nonzeros
#'
#' Split matrix columns into chunks below the dgCMatrix nonzero limit.
#'
#' @param col_nonzeros Numeric vector giving the number of nonzero entries in
#'   each column, in matrix column order.
#' @param max_nonzeros Maximum total nonzero entries allowed in one chunk.
#'   Individual columns at or above this value are rejected because they cannot
#'   fit in any chunk.
#' @return A list of integer column-index vectors, each safe to coerce to a
#'   `dgCMatrix` without exceeding `max_nonzeros`.
#' @keywords internal

split_cols_by_nonzeros <- function(col_nonzeros, max_nonzeros) {
  col_nonzeros <- as.numeric(col_nonzeros)
  if (any(col_nonzeros >= max_nonzeros)) {
    stop("At least one cell column exceeds the dgCMatrix nonzero limit.")
  }

  chunks <- vector("list", length(col_nonzeros))
  n_chunks <- 0L
  chunk_start <- 1L
  current_nonzeros <- 0
  for (idx in seq_along(col_nonzeros)) {
    if (idx > chunk_start && current_nonzeros + col_nonzeros[[idx]] >= max_nonzeros) {
      n_chunks <- n_chunks + 1L
      chunks[[n_chunks]] <- seq.int(chunk_start, idx - 1L)
      chunk_start <- idx
      current_nonzeros <- 0
    }
    current_nonzeros <- current_nonzeros + col_nonzeros[[idx]]
  }
  n_chunks <- n_chunks + 1L
  chunks[[n_chunks]] <- seq.int(chunk_start, length(col_nonzeros))
  chunks[seq_len(n_chunks)]
}

#' Get chromVAR obj from peak matrix
#'
#' Build a GC-bias annotated chromVAR SummarizedExperiment from ATAC peaks.
#'
#' @param ATAC_peak_matrix Peak-by-cell ATAC matrix used for chromVAR, LSI, or Seurat/Signac export.
#' @param ATAC_peak_GRanges GRanges for ATAC peaks, aligned by peak name to the corresponding ATAC matrix rows.
#' @param genome Genome build key used to choose chromosome sizes, blacklist resources, and external-tool parameters.
#' @return chromVAR SummarizedExperiment with nonzero peaks, aligned ranges, and
#'   GC/flanking bias metadata.
#' @keywords internal

get_chromVAR_obj_from_peak_matrix <- function(ATAC_peak_matrix,
                                              ATAC_peak_GRanges,
                                              genome = "GRCh38") {
  genome_obj <- get_chromVAR_genome_obj(genome)
  peak_matrix <- ATAC_peak_matrix
  peak_matrix <- peak_matrix[BPCells::rowSums(peak_matrix) > 0, , drop = FALSE]

  peak_ranges <- get_peak_ranges_for_peak_names(
    peak_names = rownames(peak_matrix),
    ATAC_peak_GRanges = ATAC_peak_GRanges,
    genome_obj = genome_obj
  )
  names(peak_ranges) <- rownames(peak_matrix)

  chromVAR_obj <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = peak_matrix),
    rowRanges = peak_ranges
  ) |>
    betterChromVAR::addGCBias(genome = genome_obj)

  row_data <- data.frame(SummarizedExperiment::rowData(chromVAR_obj))
  row_data[is.na(row_data)] <- 0
  SummarizedExperiment::rowData(chromVAR_obj) <- row_data
  chromVAR_obj
}

get_chromVAR_peak_expectation <- function(chromVAR_obj) {
  as.numeric(row_sums_for_peak_count_matrix(SummarizedExperiment::assay(chromVAR_obj, "counts"))) / ncol(chromVAR_obj)
}

get_chromVAR_background_bins <- function(chromVAR_obj, expectation = NULL) {
  if (is.null(expectation)) {
    expectation <- get_chromVAR_peak_expectation(chromVAR_obj)
  }

  betterChromVAR::getBackgroundBins(
    x = expectation,
    bias = SummarizedExperiment::rowData(chromVAR_obj)$bias,
    flbias = SummarizedExperiment::rowData(chromVAR_obj)$flbias,
    verbose = FALSE
  )
}

#' Get chromVAR cell chunk tibble
#'
#' Split chromVAR cells into chunks below the sparse-matrix nonzero limit.
#'
#' @param chromVAR_obj chromVAR SummarizedExperiment containing deviations, annotations, and background metadata.
#' @param chunk_nonzero_limit Maximum total nonzero matrix entries allowed in one chromVAR computation chunk.
#' @return Chunk tibble with `chunk_id`, list-column `cell_idx`, cell counts, and
#'   nonzero counts.
#' @keywords internal

get_chromVAR_cell_chunk_tibble <- function(chromVAR_obj, chunk_nonzero_limit = 2^27) {
  counts_matrix <- SummarizedExperiment::assay(chromVAR_obj, "counts")
  if (!inherits(counts_matrix, "IterableMatrix")) {
    return(tibble::tibble(
      chunk_id = "chunk_001",
      cell_idx = list(seq_len(ncol(counts_matrix))),
      n_cells = ncol(counts_matrix),
      n_nonzero = as.numeric(Matrix::nnzero(counts_matrix))
    ))
  }

  col_nonzeros <- BPCells::matrix_stats(counts_matrix, col_stats = "nonzero")$col_stats["nonzero", ]
  chunks <- split_cols_by_nonzeros(
    col_nonzeros = col_nonzeros,
    max_nonzeros = chunk_nonzero_limit
  )

  tibble::tibble(
    chunk_id = sprintf("chunk_%03d", seq_along(chunks)),
    cell_idx = chunks,
    n_cells = lengths(chunks),
    n_nonzero = purrr::map_dbl(chunks, \(chunk_idx) sum(col_nonzeros[chunk_idx]))
  )
}

#' Get chromVAR chunk context record
#'
#' Materialize one chromVAR cell chunk with its computed background matrix.
#'
#' @param chunk_record One row from the chromVAR chunk tibble, identifying the cells/features for a dynamic branch.
#' @param chromVAR_obj chromVAR SummarizedExperiment containing deviations, annotations, and background metadata.
#' @param background_bins chromVAR background-bin object aligned to the full peak set.
#' @param expectation chromVAR expectation matrix or vector aligned to peaks and cells.
#' @return A single branch record, usually a list or one-row tibble, carrying all inputs needed by a dynamic target branch.
#' @keywords internal

get_chromVAR_chunk_context_record <- function(chunk_record, chromVAR_obj, background_bins, expectation) {
  counts_matrix <- SummarizedExperiment::assay(chromVAR_obj, "counts")
  chunk_idx <- chunk_record$cell_idx[[1]]
  chunk_counts <- methods::as(counts_matrix[, chunk_idx, drop = FALSE], "dgCMatrix")
  chunk_background <- betterChromVAR::computeBackgrounds(
    object = chunk_counts,
    bins = background_bins,
    expectation = expectation,
    verbose = FALSE
  )

  list(
    chunk_id = chunk_record$chunk_id[[1]],
    cell_names = colnames(chunk_counts),
    counts = chunk_counts,
    background = chunk_background
  )
}

align_chromVAR_annotations_to_obj <- function(annotations, chromVAR_obj) {
  peak_names <- get_peak_names_from_GRanges(SummarizedExperiment::rowRanges(chromVAR_obj))
  missing_peaks <- setdiff(peak_names, rownames(annotations))
  if (length(missing_peaks) > 0) {
    stop("Missing annotation rows for ", length(missing_peaks), " chromVAR peak(s).")
  }
  annotations[peak_names, , drop = FALSE]
}

#' Compute chromVAR annotation chunk result
#'
#' Compute betterChromVAR deviations for one cell chunk and one annotation matrix.
#'
#' @param annotations Peak-by-feature annotation matrix. Rows are matched to
#'   `chromVAR_obj` peak names and converted to sparse matrix form.
#' @param chromVAR_obj Template chromVAR SummarizedExperiment providing peak
#'   ranges and row metadata shared across chunks.
#' @param chunk_context_record List returned by the chunk-preparation helper,
#'   containing `chunk_id`, chunk count matrix, and precomputed background bins.
#' @param compute Character vector passed to `betterChromVAR::computeDeviationsAnalytic()`;
#'   commonly `deviations`, `z`, or both.
#' @return A betterChromVAR result list for the chunk with an added `chunk_id`
#'   element for deterministic recombination.
#' @keywords internal

compute_chromVAR_annotation_chunk_result <- function(annotations, chromVAR_obj, chunk_context_record, compute = c("deviations", "z")) {
  annotations <- align_chromVAR_annotations_to_obj(
    annotations = annotations,
    chromVAR_obj = chromVAR_obj
  )
  annotations <- Matrix::Matrix(annotations, sparse = TRUE)

  chunk_obj <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = chunk_context_record$counts),
    rowRanges = SummarizedExperiment::rowRanges(chromVAR_obj)
  )
  SummarizedExperiment::rowData(chunk_obj) <- SummarizedExperiment::rowData(chromVAR_obj)

  out <- betterChromVAR::computeDeviationsAnalytic(
    object = chunk_obj,
    background = chunk_context_record$background,
    annotations = annotations,
    verbose = FALSE,
    retSE = FALSE,
    compute = compute
  )
  out$chunk_id <- chunk_context_record$chunk_id
  out
}

#' Combine chromVAR chunk results
#'
#' Recombine per-chunk betterChromVAR results into full-cell outputs.
#'
#' @param chunk_results List of chunk results from `compute_chromVAR_annotation_chunk_result()`;
#'   each element must contain `chunk_id`, `z`, `deviations`, and `total`.
#' @param chromVAR_obj Template chromVAR SummarizedExperiment whose cell and
#'   metadata order defines the combined output.
#' @param annotations Peak-by-feature annotation matrix used to compute feature
#'   counts in the combined row metadata.
#' @return A list containing the original `chromVAR_obj`, combined chromVAR
#'   SummarizedExperiment, full z-score matrix, and motif/annotation matrix.
#' @keywords internal

combine_chromVAR_chunk_results <- function(chunk_results, chromVAR_obj, annotations) {
  chunk_results <- chunk_results[order(purrr::map_chr(chunk_results, "chunk_id"))]
  chromVAR_z_scores <- do.call(cbind, purrr::map(chunk_results, "z"))
  chromVAR_deviations <- do.call(cbind, purrr::map(chunk_results, "deviations"))
  total <- Reduce(`+`, purrr::map(chunk_results, "total"))

  sd_deviations <- matrixStats::rowSds(chromVAR_z_scores, na.rm = TRUE)
  p_sd <- stats::pchisq(
    (ncol(chromVAR_obj) - 1) * (sd_deviations^2),
    df = ncol(chromVAR_obj) - 1,
    lower.tail = FALSE
  )
  dev_metadata <- data.frame(
    N = Matrix::colSums(annotations),
    total = total,
    variability = sd_deviations,
    var.pval = p_sd,
    var.adjPval = stats::p.adjust(p = p_sd, method = "BH")
  )

  chromVAR_dev <- SummarizedExperiment::SummarizedExperiment(
    assays = list(deviations = chromVAR_deviations, z = chromVAR_z_scores),
    colData = SummarizedExperiment::colData(chromVAR_obj),
    rowData = dev_metadata,
    metadata = S4Vectors::metadata(chromVAR_obj)
  )

  list(
    chromVAR_obj = chromVAR_obj,
    chromVAR_dev = chromVAR_dev,
    chromVAR_z_scores = chromVAR_z_scores,
    motifs = annotations
  )
}

#' Run betterChromVAR from peak inputs
#'
#' Build chromVAR objects and motif deviations directly from peak matrices.
#'
#' @param peak_matrix Peak-by-cell accessibility matrix. Row names must align
#'   with `peak_ranges`; columns are cells.
#' @param peak_ranges GRanges of consensus peaks; names must match peak rows used in peak-weight or accessibility matrices.
#' @param motif_matrix Peak-by-motif annotation matrix with rows in peak order.
#' @param genome_obj Genome object passed to chromVAR/betterChromVAR for sequence-aware background and motif calculations.
#' @param chunk_nonzero_limit Maximum nonzero entries allowed per cell chunk
#'   before coercing to `dgCMatrix`.
#' @return A list containing the source chromVAR object, combined deviations
#'   SummarizedExperiment, z-score matrix, and motif annotation matrix.
#' @keywords internal

run_betterChromVAR_from_peak_inputs <- function(peak_matrix,
                                                peak_ranges,
                                                motif_matrix,
                                                genome_obj,
                                                chunk_nonzero_limit = 2^31 - 2) {
  names(peak_ranges) <- rownames(peak_matrix)

  chromVAR_obj <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = peak_matrix),
    rowRanges = peak_ranges
  ) |>
    betterChromVAR::addGCBias(genome = genome_obj)

  row.data <- data.frame(SummarizedExperiment::rowData(x = chromVAR_obj))
  row.data[is.na(x = row.data)] <- 0
  SummarizedExperiment::rowData(x = chromVAR_obj) <- row.data

  expectation <- as.numeric(BPCells::rowSums(peak_matrix)) / ncol(peak_matrix)
  bg <- betterChromVAR::getBackgroundBins(
    x = expectation,
    bias = SummarizedExperiment::rowData(chromVAR_obj)$bias,
    flbias = SummarizedExperiment::rowData(chromVAR_obj)$flbias,
    verbose = FALSE
  )

  chunks <- split_cols_by_nonzeros(
    col_nonzeros = BPCells::matrix_stats(peak_matrix, col_stats = "nonzero")$col_stats["nonzero", ],
    max_nonzeros = chunk_nonzero_limit
  )

  # If this becomes runtime-limited, promote these column chunks to dynamic targets.
  # Keep the expectation/background model global so chunk outputs stay comparable.
  chunk_results <- purrr::map(chunks, \(chunk_idx) {
    chunk_counts <- methods::as(peak_matrix[, chunk_idx, drop = FALSE], "dgCMatrix")
    chunk_obj <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = chunk_counts),
      rowRanges = peak_ranges
    )
    SummarizedExperiment::rowData(chunk_obj) <- SummarizedExperiment::rowData(chromVAR_obj)

    chunk_bg <- betterChromVAR::computeBackgrounds(
      object = chunk_counts,
      bins = bg,
      expectation = expectation,
      verbose = FALSE
    )
    betterChromVAR::computeDeviationsAnalytic(
      object = chunk_obj,
      background = chunk_bg,
      annotations = motif_matrix,
      verbose = FALSE,
      retSE = FALSE,
      compute = c("deviations", "z")
    )
  })

  chromVAR_z_scores <- do.call(cbind, purrr::map(chunk_results, "z"))
  chromVAR_deviations <- do.call(cbind, purrr::map(chunk_results, "deviations"))
  total <- Reduce(`+`, purrr::map(chunk_results, "total"))

  sd_deviations <- matrixStats::rowSds(chromVAR_z_scores, na.rm = TRUE)
  p_sd <- stats::pchisq(
    (ncol(peak_matrix) - 1) * (sd_deviations^2),
    df = ncol(peak_matrix) - 1,
    lower.tail = FALSE
  )
  dev_metadata <- data.frame(
    N = Matrix::colSums(motif_matrix),
    total = total,
    variability = sd_deviations,
    var.pval = p_sd,
    var.adjPval = stats::p.adjust(p = p_sd, method = "BH")
  )

  chromVAR_dev <- SummarizedExperiment::SummarizedExperiment(
    assays = list(deviations = chromVAR_deviations, z = chromVAR_z_scores),
    colData = SummarizedExperiment::colData(chromVAR_obj),
    rowData = dev_metadata,
    metadata = S4Vectors::metadata(chromVAR_obj)
  )

  list(
    chromVAR_obj = chromVAR_obj,
    chromVAR_dev = chromVAR_dev,
    chromVAR_z_scores = chromVAR_z_scores,
    motifs = motif_matrix
  )
}

get_marker_TF_activities_from_chromVAR_BPCells_z_scores <- function(chromVAR_z_scores_BPCells_matrix, metadata_tibble, group_col) {
  if (!inherits(chromVAR_z_scores_BPCells_matrix, "IterableMatrix")) {
    stop("chromVAR_z_scores_BPCells_matrix must be a BPCells IterableMatrix.")
  }

  keep_barcodes <- intersect(colnames(chromVAR_z_scores_BPCells_matrix), metadata_tibble$barcode_w_prefix)
  marker_mat <- chromVAR_z_scores_BPCells_matrix[, keep_barcodes, drop = FALSE]
  groups <- stats::setNames(metadata_tibble[[group_col]], metadata_tibble$barcode_w_prefix)[keep_barcodes]

  BPCells::marker_features(marker_mat, groups, method = "wilcoxon") |>
    dplyr::mutate(
      avg_diff = foreground_mean - background_mean,
      p_val = p_val_raw,
      p_val_adj = stats::p.adjust(p_val_raw, method = "BH"),
      cluster = foreground,
      gene = feature
    )
}

#' Add cell types to metadata from module scores
#'
#' Assign cluster-level cell-type labels from marker module score columns.
#'
#' @param metadata_tibble Cell metadata tibble containing the clustering column
#'   and one numeric module-score column per marker set name.
#' @param named_marker_genes_list Named list of marker gene vectors. Only names
#'   are used here; they define which score columns to compare.
#' @param allow_multiple_cell_types Logical; when `TRUE`, the top two marker
#'   labels are joined with `&` if their cluster mean scores differ by less than
#'   `std_threshold`.
#' @param cluster_column Metadata column containing the cluster labels to name.
#' @param max_q_cap Upper quantile used to cap scaled module scores before
#'   computing cluster means, reducing the influence of extreme cells.
#' @param std_threshold Minimum standardized-score gap required to call only the
#'   top marker label when multiple labels are allowed.
#' @return `metadata_tibble` with `<cluster_column>_cell_type` and
#'   `<cluster_column>_named` columns joined back by cluster.
#' @keywords internal

add_cell_types_to_metadata_from_module_scores <- function(
  metadata_tibble,
  named_marker_genes_list,
  allow_multiple_cell_types,
  cluster_column = "LSI_harmony_SNN_cluster",
  max_q_cap = 0.95,
  std_threshold = 0.5
) {
  cell_types <- names(named_marker_genes_list)
  named_cluster_column <- paste0(cluster_column, "_named")
  cell_type_column <- paste0(cluster_column, "_cell_type")

  missing_score_cols <- setdiff(cell_types, colnames(metadata_tibble))
  if (length(missing_score_cols) > 0) {
    warning("Dropping marker modules without metadata score columns: ", paste(missing_score_cols, collapse = ", "))
    cell_types <- intersect(cell_types, colnames(metadata_tibble))
  }
  if (length(cell_types) == 0) {
    stop("No configured marker module score columns were found in metadata.")
  }

  mean_score_per_cell_type_and_cluster_capped <- metadata_tibble |>
    dplyr::mutate(dplyr::across(dplyr::any_of(cell_types), ~ scale(.x, center = TRUE, scale = TRUE))) |>
    dplyr::group_by(.data[[cluster_column]]) |>
    dplyr::summarise(dplyr::across(
      dplyr::any_of(cell_types),
      ~ {
        hi <- stats::quantile(.x, probs = max_q_cap, na.rm = TRUE)
        vals <- .x[.x <= hi]
        mean(vals, na.rm = TRUE)
      }
    ))

  cell_type_cluster_mapping_table <- mean_score_per_cell_type_and_cluster_capped |>
    dplyr::rowwise() |>
    dplyr::mutate(
      !!cell_type_column := {
        score_vec <- dplyr::c_across(dplyr::all_of(cell_types))
        names(score_vec) <- cell_types
        sorted_scores_vec <- score_vec |> sort(decreasing = TRUE)
        max_score <- sorted_scores_vec[1]
        second_score <- sorted_scores_vec[2]

        if (!is.na(second_score) && (max_score - second_score) < std_threshold && allow_multiple_cell_types) {
          paste0(names(max_score), "&", names(second_score))
        } else {
          names(max_score)
        }
      }
    ) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      !!cluster_column := as.character(.data[[cluster_column]]) |> get_mixsorted_factor(),
      !!cell_type_column := .data[[cell_type_column]] |> get_mixsorted_factor(),
      !!named_cluster_column := stringr::str_c(.data[[cluster_column]], "-", .data[[cell_type_column]]) |> get_mixsorted_factor()
    )

  metadata_tibble |>
    dplyr::left_join(cell_type_cluster_mapping_table, by = cluster_column)
}
