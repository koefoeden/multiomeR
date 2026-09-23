load_scoring_test_runtime <- function() {
  if (!exists("summarize_cluster_UCell_counts", mode = "function")) {
    load_project_test_runtime()
  }
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
  matrix_dir <- tempfile("bpcells_scoring_parity_")
  sparse_counts <- Matrix::Matrix(counts, sparse = TRUE)
  sparse_counts <- BPCells::convert_matrix_type(sparse_counts, type = "uint32_t")
  invisible(BPCells::write_matrix_dir(sparse_counts, matrix_dir))
  BPCells::open_matrix_dir(matrix_dir)
}

# UCell 2.14.0 can emit non-fatal R stack-imbalance warnings under R 4.5. Run all
# reference calls in one disposable process so those package-level warnings
# cannot corrupt this session; callr still returns the exact matrices compared.
run_ucell_references <- function(counts, calls) {
  require_reference_version("UCell", "2.14.0")
  callr::r(
    function(counts, calls) {
      counts <- Matrix::Matrix(counts, sparse = TRUE)
      lapply(calls, function(arguments) as.matrix(do.call(UCell::ScoreSignatures_UCell, c(
        list(matrix = counts, name = "", ncores = 1, ties.method = "average",
          BPPARAM = BiocParallel::SerialParam()),
        arguments
      ))))
    },
    args = list(counts = counts, calls = calls),
    show = FALSE
  )
}

testthat::test_that("integration: BPCells UCell scores match ScoreSignatures_UCell", {
  load_scoring_test_runtime()
  counts <- make_counts_matrix()
  bpcells_counts <- make_bpcells_matrix(counts)
  markers <- list(
    alpha = c("gene003+", "gene017", "gene029-", "missing_alpha+"),
    beta = c("gene041", "gene053+", "gene067-", "gene079-", "missing_beta-"),
    gamma = c("gene101+", "gene113", "gene127-")
  )
  modes <- c(impute = "impute", skip = "skip")
  expected <- run_ucell_references(counts, lapply(modes, function(mode) list(
    features = markers, maxRank = 80, w_neg = 0.75, chunk.size = 7, missing_genes = mode
  )))

  for (mode in modes) {
    observed <- calculate_BPCells_UCell_scores_from_matrix(
      counts_matrix = bpcells_counts,
      features = markers,
      max_rank = 80,
      chunk_size = 7,
      w_neg = 0.75,
      ties_method = "average",
      missing_genes = mode
    )
    testthat::expect_identical(
      as.matrix(observed),
      expected[[mode]][, colnames(observed), drop = FALSE],
      info = paste0("missing_genes = ", mode)
    )
  }
})

testthat::test_that("integration: cluster annotation evidence matches UCell averaged within clusters", {
  load_scoring_test_runtime()
  counts <- make_counts_matrix()
  bpcells_counts <- make_bpcells_matrix(counts)
  markers <- list(
    positive = c("gene001", "gene003", "gene007"),
    positive_pair = c("gene010", "gene020"),
    signed = c("gene041+", "gene053", "gene067-"),
    signed_pair = c("gene079", "gene089-"),
    negative_only = "gene101-"
  )
  metadata <- data.frame(
    barcode_w_prefix = colnames(counts),
    cluster = rep(c("1", "2", "3"), length.out = ncol(counts)),
    GEM_well_ID = c("singleton", rep(c("well_A", "well_B"), 18L))
  )
  control <- prepare_cluster_UCell_controls(bpcells_counts, metadata, markers)
  control$n_controls <- 9L
  control$draws <- control$draws[, seq_len(control$n_controls), drop = FALSE]
  control$max_rank <- 80L
  score <- function(chunk_size, workers) {
    summaries <- summarize_cluster_UCell_counts(bpcells_counts, metadata, control, "cluster",
      chunk_size = chunk_size, workers = workers, include_cell_scores = TRUE)
    score_cluster_UCell_summaries(summaries, control)
  }
  scored <- score(chunk_size = 7L, workers = 1L)

  # Score every observed and matched-control signature per cell with UCell.
  signature_set <- function(genes) c(list(genes), lapply(seq_len(control$n_controls), function(index) {
    paste0(control$reference_genes[control$draws[sub("-$", "", genes), index]],
      ifelse(grepl("-$", genes), "-", ""))
  }))
  variants <- lapply(control$markers, list)
  signatures <- unlist(lapply(names(variants), function(label) {
    unlist(lapply(seq_along(variants[[label]]), function(variant) {
      stats::setNames(signature_set(variants[[label]][[variant]]),
        sprintf("%s.%d.%d", label, variant, seq_len(control$n_controls + 1L)))
    }), recursive = FALSE)
  }), recursive = FALSE)
  reference <- run_ucell_references(counts, list(list(features = signatures, maxRank = 80, w_neg = 1)))[[1]]
  cluster_means <- apply(reference, 2L, function(values) tapply(values, metadata$cluster, mean))
  reference_evidence <- function(label, variant) {
    columns <- sprintf("%s.%d.%d", label, variant, seq_len(control$n_controls + 1L))
    observed <- cluster_means[, columns[1L]]
    null <- cluster_means[, columns[-1L], drop = FALSE]
    q95 <- apply(null, 1L, stats::quantile, probs = 0.95, names = FALSE)
    data.frame(cluster = rownames(cluster_means), label, mean_score = unname(observed),
      control_mean = unname(rowMeans(null)), control_q95 = q95, excess = unname(observed) - q95)
  }

  testthat::expect_equal(
    unname(scored$cell_scores[metadata$barcode_w_prefix, names(markers)]),
    unname(reference[metadata$barcode_w_prefix, sprintf("%s.1.1", names(markers))]),
    tolerance = 1e-12
  )
  expected <- do.call(rbind, lapply(names(markers), reference_evidence, variant = 1L))
  observed <- scored$evidence[order(match(scored$evidence$label, names(markers)), scored$evidence$cluster), names(expected)]
  rownames(observed) <- rownames(expected) <- NULL
  testthat::expect_equal(observed, expected, tolerance = 1e-12)
  # Chunking and fork workers must not change any aggregated statistic.
  rescored <- score(chunk_size = 11L, workers = 2L)
  testthat::expect_equal(rescored$evidence, scored$evidence, tolerance = 1e-14)
  testthat::expect_equal(rescored$cell_scores[rownames(scored$cell_scores), ], scored$cell_scores,
    tolerance = 1e-14)
})

testthat::test_that("integration: cell-cycle scores and phases match Seurat::CellCycleScoring", {
  load_scoring_test_runtime()
  require_reference_version("Seurat", "5.5.0")
  cell_cycle_genes <- Seurat::cc.genes.updated.2019
  genes <- c(cell_cycle_genes$s.genes, cell_cycle_genes$g2m.genes, sprintf("gene%04d", seq_len(2900L)))
  cells <- sprintf("cell%02d", seq_len(60L))
  set.seed(223)
  expression <- outer(exp(stats::rnorm(length(genes))), rep(1, length(cells)))
  cycling <- list(cell_cycle_genes$s.genes, cell_cycle_genes$g2m.genes)
  for (phase in 1:2) {
    selected <- (phase - 1L) * 20L + seq_len(20L)
    expression[genes %in% cycling[[phase]], selected] <- 3 * expression[genes %in% cycling[[phase]], selected]
  }
  counts <- matrix(stats::rpois(length(expression), expression), length(genes), dimnames = list(genes, cells))

  observed <- add_cell_cycle_scores_to_cell_attr(
    make_bpcells_matrix(counts), data.frame(row.names = cells)
  )
  object <- Seurat::NormalizeData(
    SeuratObject::CreateSeuratObject(counts = Matrix::Matrix(counts, sparse = TRUE)),
    verbose = FALSE
  )
  expected <- Seurat::CellCycleScoring(
    object,
    s.features = cell_cycle_genes$s.genes,
    g2m.features = cell_cycle_genes$g2m.genes
  )[[]]

  # Production normalizes in BPCells, whose lower-precision arithmetic moves scores by
  # about 1e-8 relative to Seurat's double-precision NormalizeData(); bins and phases are unchanged.
  testthat::expect_equal(observed[cells, c("S.Score", "G2M.Score")], expected[cells, c("S.Score", "G2M.Score")],
    tolerance = 1e-6)
  testthat::expect_identical(observed[cells, "Phase"], expected[cells, "Phase"])
  testthat::expect_setequal(observed$Phase, c("G1", "S", "G2M"))
})
