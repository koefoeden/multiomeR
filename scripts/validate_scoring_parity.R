#!/usr/bin/env Rscript

if (!exists("calculate_BPCells_UCell_scores_from_matrix", mode = "function")) {
  source("R/bootstrap_helpers.R")
  load_project_runtime()
}

fail <- function(...) {
  stop(paste0(...), call. = FALSE)
}

expect_identical_scores <- function(observed, expected, label) {
  observed <- as.matrix(observed)
  expected <- as.matrix(expected)

  if (!identical(dim(observed), dim(expected))) {
    fail(label, ": dimension mismatch")
  }
  if (!identical(dimnames(observed), dimnames(expected))) {
    fail(label, ": dimname mismatch")
  }
  if (!identical(observed, expected)) {
    max_delta <- max(abs(observed - expected))
    fail(label, ": values are not identical; max delta = ", format(max_delta, digits = 16))
  }

  invisible(TRUE)
}

make_counts_matrix <- function() {
  set.seed(219)
  feature_names <- sprintf("gene%03d", seq_len(500))
  cell_names <- sprintf("cell%02d", seq_len(37))
  counts <- matrix(
    stats::rpois(length(feature_names) * length(cell_names), lambda = 2),
    nrow = length(feature_names),
    dimnames = list(feature_names, cell_names)
  )
  counts[sample(length(counts), 1200)] <- 0L
  counts[seq(1, 500, by = 17), ] <- counts[seq(1, 500, by = 17), ] %% 3L
  storage.mode(counts) <- "integer"
  counts
}

make_bpcells_matrix <- function(counts) {
  matrix_dir <- tempfile("bpcells_ucell_parity_")
  sparse_counts <- Matrix::Matrix(counts, sparse = TRUE)
  sparse_counts <- BPCells::convert_matrix_type(sparse_counts, type = "uint32_t")
  invisible(BPCells::write_matrix_dir(sparse_counts, matrix_dir))
  BPCells::open_matrix_dir(matrix_dir)
}

compare_score_signatures_ucell <- function(counts, bpcells_counts, marker_genes, missing_genes) {
  expected <- suppressWarnings(UCell::ScoreSignatures_UCell(
    matrix = Matrix::Matrix(counts, sparse = TRUE),
    features = marker_genes,
    maxRank = 80,
    w_neg = 0.75,
    name = "",
    chunk.size = 7,
    missing_genes = missing_genes,
    BPPARAM = BiocParallel::SerialParam(),
    ncores = 1,
    ties.method = "average"
  ))

  observed <- calculate_BPCells_UCell_scores_from_matrix(
    counts_matrix = bpcells_counts,
    features = marker_genes,
    max_rank = 80,
    chunk_size = 7,
    w_neg = 0.75,
    ties_method = "average",
    missing_genes = missing_genes
  )

  expect_identical_scores(
    observed = observed,
    expected = expected[, colnames(observed), drop = FALSE],
    label = paste0("UCell::ScoreSignatures_UCell missing_genes=", missing_genes)
  )
}

compare_add_module_score_ucell <- function(counts, bpcells_counts, marker_genes) {
  seurat_obj <- SeuratObject::CreateSeuratObject(counts = Matrix::Matrix(counts, sparse = TRUE))
  seurat_obj <- suppressWarnings(UCell::AddModuleScore_UCell(
    obj = seurat_obj,
    features = marker_genes,
    maxRank = 80,
    chunk.size = 7,
    BPPARAM = BiocParallel::SerialParam(),
    ncores = 1,
    storeRanks = FALSE,
    w_neg = 0.75,
    assay = "RNA",
    slot = "counts",
    ties.method = "average",
    missing_genes = "impute",
    force.gc = FALSE,
    name = ""
  ))

  observed <- calculate_BPCells_UCell_scores_from_matrix(
    counts_matrix = bpcells_counts,
    features = marker_genes,
    max_rank = 80,
    chunk_size = 7,
    w_neg = 0.75,
    ties_method = "average",
    missing_genes = "impute"
  )
  expected <- seurat_obj@meta.data[, colnames(observed), drop = FALSE]

  expect_identical_scores(
    observed = observed,
    expected = expected,
    label = "UCell::AddModuleScore_UCell"
  )
}

reference_seurat_module_scores <- function(normalized_matrix, features, nbin, ctrl, seed) {
  assay <- suppressWarnings(SeuratObject::CreateAssay5Object(data = normalized_matrix))
  set.seed(seed)
  scores <- Seurat::AddModuleScore(
    object = assay,
    features = features,
    nbin = nbin,
    ctrl = ctrl,
    name = "Module",
    slot = "data"
  )[]
  colnames(scores) <- names(features)
  scores
}

compare_seurat_module_scores <- function(counts, bpcells_counts) {
  features <- list(
    S.Score = c("gene005", "gene019", "gene033", "gene047"),
    G2M.Score = c("gene061", "gene075", "gene089", "gene103")
  )

  expected <- reference_seurat_module_scores(
    normalized_matrix = counts,
    features = features,
    nbin = 10,
    ctrl = 4,
    seed = 11
  )
  observed <- calculate_BPCells_module_scores_from_matrix(
    normalized_data = bpcells_counts,
    features = features,
    nbin = 10,
    ctrl = 4,
    seed = 11
  )

  expect_identical_scores(
    observed = observed,
    expected = expected,
    label = "Seurat::AddModuleScore"
  )
}

compare_cell_cycle_scores <- function(counts, bpcells_counts) {
  s_features <- c("gene005", "gene019", "gene033", "gene047")
  g2m_features <- c("gene061", "gene075", "gene089", "gene103")
  features <- list(S.Score = s_features, G2M.Score = g2m_features)

  expected_scores <- reference_seurat_module_scores(
    normalized_matrix = counts,
    features = features,
    nbin = 10,
    ctrl = 4,
    seed = 11
  )
  expected_phase <- apply(expected_scores, 1, function(scores) {
    if (all(scores < 0)) {
      return("G1")
    }
    if (sum(scores == max(scores)) > 1) {
      return("Undecided")
    }
    c("S", "G2M")[which(scores == max(scores))]
  })
  expected <- data.frame(
    S.Score = expected_scores[, "S.Score"],
    G2M.Score = expected_scores[, "G2M.Score"],
    Phase = expected_phase,
    row.names = rownames(expected_scores)
  )

  observed <- calculate_BPCells_cell_cycle_scores_from_matrix(
    normalized_data = bpcells_counts,
    s.features = s_features,
    g2m.features = g2m_features,
    nbin = 10,
    seed = 11
  )

  expect_identical_scores(
    observed = observed[, c("S.Score", "G2M.Score"), drop = FALSE],
    expected = expected[, c("S.Score", "G2M.Score"), drop = FALSE],
    label = "cell-cycle Seurat::AddModuleScore scores"
  )
  if (!identical(observed$Phase, expected$Phase)) {
    fail("cell-cycle Seurat::AddModuleScore phase mismatch")
  }
}

compare_metadata_join <- function(bpcells_counts, marker_genes) {
  metadata <- tibble::tibble(
    barcode_w_prefix = c("cell03", "cell01", "missing_cell", "cell37"),
    batch = c("a", "a", "b", "b")
  )

  scored_metadata <- add_GEX_UCell_scores_to_metadata(
    metadata_tibble = metadata,
    named_marker_genes_list = marker_genes,
    GEX_counts_matrix = bpcells_counts,
    max_rank = 80,
    chunk_size = 7,
    w_neg = 0.75,
    missing_genes = "impute"
  )

  if (!identical(scored_metadata$barcode_w_prefix, metadata$barcode_w_prefix)) {
    fail("add_GEX_UCell_scores_to_metadata: metadata row order changed")
  }
  missing_row <- scored_metadata$barcode_w_prefix == "missing_cell"
  if (!all(is.na(scored_metadata[missing_row, names(marker_genes)]))) {
    fail("add_GEX_UCell_scores_to_metadata: unmatched metadata barcode should keep NA scores")
  }
  if (anyNA(scored_metadata[!missing_row, names(marker_genes)])) {
    fail("add_GEX_UCell_scores_to_metadata: matched metadata barcodes received NA scores")
  }

  invisible(TRUE)
}

counts <- make_counts_matrix()
bpcells_counts <- make_bpcells_matrix(counts)
marker_genes <- list(
  alpha = c("gene003+", "gene017", "gene029-", "missing_alpha+"),
  beta = c("gene041", "gene053+", "gene067-", "gene079-", "missing_beta-"),
  gamma = c("gene101+", "gene113", "gene127-")
)

compare_score_signatures_ucell(counts, bpcells_counts, marker_genes, missing_genes = "impute")
compare_score_signatures_ucell(counts, bpcells_counts, marker_genes, missing_genes = "skip")
compare_add_module_score_ucell(counts, bpcells_counts, marker_genes)
compare_seurat_module_scores(counts, bpcells_counts)
compare_cell_cycle_scores(counts, bpcells_counts)
compare_metadata_join(bpcells_counts, marker_genes)

cat("scoring parity ok\n")
