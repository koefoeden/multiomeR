load_scoring_test_runtime <- function() {
  if (!exists("calculate_BPCells_UCell_scores_from_matrix", mode = "function")) {
    load_project_test_runtime()
  }
}

expect_identical_scores <- function(observed, expected, label) {
  observed <- as.matrix(observed)
  expected <- as.matrix(expected)

  testthat::expect_identical(observed, expected, info = label)
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

run_ucell_reference <- function(counts, marker_genes, method, missing_genes = "impute") {
  require_reference_version("UCell", "2.14.0")

  # UCell 2.14.0 can emit non-fatal R stack-imbalance warnings for the
  # missing-gene fixture under R 4.5. Run only the reference implementation in
  # a disposable process so those package-level warnings cannot corrupt this
  # validation session; callr still returns the exact matrix being compared.
  callr::r(
    function(counts, marker_genes, method, missing_genes) {
      counts <- Matrix::Matrix(counts, sparse = TRUE)
      if (identical(method, "score_signatures")) {
        return(UCell::ScoreSignatures_UCell(
          matrix = counts,
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
      }

      seurat_obj <- SeuratObject::CreateSeuratObject(counts = counts)
      seurat_obj <- UCell::AddModuleScore_UCell(
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
      )
      seurat_obj@meta.data
    },
    args = list(
      counts = counts,
      marker_genes = marker_genes,
      method = method,
      missing_genes = missing_genes
    ),
    show = FALSE
  )
}

make_scoring_fixture <- function() {
  counts <- make_counts_matrix()
  list(
    counts = counts,
    bpcells_counts = make_bpcells_matrix(counts),
    marker_genes = list(
      alpha = c("gene003+", "gene017", "gene029-", "missing_alpha+"),
      beta = c("gene041", "gene053+", "gene067-", "gene079-", "missing_beta-"),
      gamma = c("gene101+", "gene113", "gene127-")
    )
  )
}

compare_score_signatures_ucell <- function(counts, bpcells_counts, marker_genes, missing_genes) {
  expected <- run_ucell_reference(
    counts = counts,
    marker_genes = marker_genes,
    method = "score_signatures",
    missing_genes = missing_genes
  )

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
  expected <- run_ucell_reference(
    counts = counts,
    marker_genes = marker_genes,
    method = "add_module_score"
  )

  observed <- calculate_BPCells_UCell_scores_from_matrix(
    counts_matrix = bpcells_counts,
    features = marker_genes,
    max_rank = 80,
    chunk_size = 7,
    w_neg = 0.75,
    ties_method = "average",
    missing_genes = "impute"
  )
  expected <- expected[, colnames(observed), drop = FALSE]

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
  testthat::expect_identical(observed$Phase, expected$Phase)
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

  testthat::expect_identical(
    scored_metadata$barcode_w_prefix,
    metadata$barcode_w_prefix
  )
  missing_row <- scored_metadata$barcode_w_prefix == "missing_cell"
  testthat::expect_true(
    all(is.na(scored_metadata[missing_row, names(marker_genes)]))
  )
  testthat::expect_false(
    anyNA(scored_metadata[!missing_row, names(marker_genes)])
  )
}

testthat::test_that("integration: BPCells UCell scores match imputed reference scores", {
  load_scoring_test_runtime()
  fixture <- make_scoring_fixture()
  compare_score_signatures_ucell(
    fixture$counts,
    fixture$bpcells_counts,
    fixture$marker_genes,
    missing_genes = "impute"
  )
})

testthat::test_that("integration: BPCells UCell scores match skipped reference scores", {
  load_scoring_test_runtime()
  fixture <- make_scoring_fixture()
  compare_score_signatures_ucell(
    fixture$counts,
    fixture$bpcells_counts,
    fixture$marker_genes,
    missing_genes = "skip"
  )
})

testthat::test_that("integration: BPCells UCell scores match AddModuleScore_UCell", {
  load_scoring_test_runtime()
  fixture <- make_scoring_fixture()
  compare_add_module_score_ucell(
    fixture$counts,
    fixture$bpcells_counts,
    fixture$marker_genes
  )
})

testthat::test_that("integration: BPCells module scores match Seurat", {
  load_scoring_test_runtime()
  fixture <- make_scoring_fixture()
  compare_seurat_module_scores(fixture$counts, fixture$bpcells_counts)
})

testthat::test_that("integration: BPCells cell-cycle scores and phases match Seurat", {
  load_scoring_test_runtime()
  fixture <- make_scoring_fixture()
  compare_cell_cycle_scores(fixture$counts, fixture$bpcells_counts)
})

testthat::test_that("integration: UCell scores join metadata without changing row order", {
  load_scoring_test_runtime()
  fixture <- make_scoring_fixture()
  compare_metadata_join(fixture$bpcells_counts, fixture$marker_genes)
})
