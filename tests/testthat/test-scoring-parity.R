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

testthat::test_that("cluster UCell summaries preserve per-cell scores and group means", {
  load_scoring_test_runtime()
  counts <- make_counts_matrix()
  markers <- list(A = c("gene001", "gene003", "gene007"), B = c("gene010", "gene020"))
  metadata <- data.frame(barcode_w_prefix = colnames(counts),
    cluster = rep(c("1", "2"), length.out = ncol(counts)),
    GEM_well_ID = c("singleton", rep(c("well_A", "well_B"), 18L)))
  control <- prepare_cluster_UCell_controls(counts, metadata, markers)
  testthat::expect_true(colnames(counts)[1] %in% control$reference_barcodes)
  testthat::expect_false(any(control$reference_genes[control$draws] %in% unlist(markers)))
  testthat::expect_true(all(apply(control$draws, 2L, anyDuplicated) == 0L))
  control$max_rank <- 80L
  summaries <- summarize_cluster_UCell_counts(counts, metadata, control, "cluster",
    chunk_size = 7L, workers = 1L, include_cell_scores = TRUE)
  expected <- calculate_BPCells_UCell_scores_from_matrix(counts, markers, max_rank = 80L)
  testthat::expect_equal(summaries$cell_scores[rownames(expected), ], as.matrix(expected), tolerance = 1e-14)
  for (index in which(summaries$groups$kind == "cluster")) {
    selected <- metadata$cluster == summaries$groups$cluster[index]
    observed <- vapply(markers, function(genes) score_UCell_rank_means(
      summaries$rank_means[, summaries$groups$key[index], drop = FALSE], genes, 80L), numeric(1L))
    testthat::expect_equal(observed, colMeans(expected[selected, ]), tolerance = 1e-14)
  }
  parallel_summary <- summarize_cluster_UCell_counts(counts, metadata, control, "cluster",
    chunk_size = 7L, workers = 2L, include_cell_scores = TRUE)
  testthat::expect_equal(parallel_summary$rank_means, summaries$rank_means, tolerance = 1e-14)
  annotation <- evaluate_cluster_UCell_evidence(score_cluster_UCell_summaries(summaries, control))
  joined <- add_cluster_UCell_annotations(metadata, annotation, "cluster")
  testthat::expect_identical(joined$barcode_w_prefix, metadata$barcode_w_prefix)
  testthat::expect_false(anyNA(joined$cluster_cell_type))
  testthat::expect_equal(unname(as.matrix(joined[, names(markers)])), unname(as.matrix(expected)), tolerance = 1e-14)
})

testthat::test_that("cluster evidence abstains on unsupported and competing signatures", {
  load_scoring_test_runtime()
  genes <- sprintf("gene%03d", seq_len(500L))
  markers <- list(A = genes[1:3], B = genes[4:6])
  reference <- data.frame(gene = genes, abundance = rep(0.002, 500L), detection = rep(0.2, 500L))
  control <- build_UCell_controls(reference, markers)
  control$max_rank <- 500L
  means <- matrix(0.01, 500L, 3L, dimnames = list(genes, c("clear", "tie", "absent")))
  means[markers$A, "clear"] <- 0.9
  means[unlist(markers), "tie"] <- 0.9
  detection <- matrix(0.5, 500L, 3L, dimnames = dimnames(means))
  evidence <- score_UCell_group_evidence(means, control, detection)
  result <- list(decisions = assign_UCell_cluster_evidence(evidence))
  testthat::expect_identical(result$decisions$status, c("Assigned", "Unassigned", "Unassigned"))
  testthat::expect_identical(result$decisions$label, c("A", NA_character_, NA_character_))
  metadata <- data.frame(barcode_w_prefix = c("x", "y", "z"), cluster = c("clear", "tie", "absent"))
  annotated <- add_cluster_UCell_annotations(metadata, result, "cluster")
  testthat::expect_identical(annotated$cluster_scDblFinder_group,
    c("A", "Unassigned_cluster_tie", "Unassigned_cluster_absent"))
  # Unavailable cell-stability diagnostics do not veto a sufficient advantage.
  singleton_means <- cbind(means, means)
  singleton_detection <- cbind(detection, detection)
  colnames(singleton_means) <- colnames(singleton_detection) <- as.character(seq_len(6L))
  singleton_groups <- data.frame(key = as.character(seq_len(6L)),
    kind = rep(c("cluster", "block"), each = 3L),
    cluster = rep(colnames(means), 2L), subgroup = rep(c("", "1"), each = 3L), cells = 1L)
  singleton <- evaluate_cluster_UCell_evidence(score_cluster_UCell_summaries(list(rank_means = singleton_means,
    detection = singleton_detection, groups = singleton_groups, n_blocks = 10L), control))
  testthat::expect_identical(singleton$decisions$status, c("Assigned", "Unassigned", "Unassigned"))
  testthat::expect_true(all(is.na(singleton$decisions$cell_stability)))
  permissive <- assign_UCell_cluster_evidence(evidence, min_advantage = 0)
  strict <- assign_UCell_cluster_evidence(evidence, min_advantage = 0.95)
  testthat::expect_identical(permissive$status, c("Assigned", "Unassigned", "Unassigned"))
  testthat::expect_identical(strict$candidate, permissive$candidate)
  testthat::expect_true(all(strict$status == "Unassigned"))
  boundary <- result$decisions$advantage[1]
  testthat::expect_identical(assign_UCell_cluster_evidence(evidence, boundary)$status[1], "Assigned")
  testthat::expect_identical(assign_UCell_cluster_evidence(evidence, boundary + 1e-8)$status[1], "Unassigned")
  testthat::expect_error(assign_UCell_cluster_evidence(evidence, -0.1))
  testthat::expect_equal(normalize_marker_panel(list(A = "gene001-", B = "gene002+"), genes),
    list(A = "gene001-", B = "gene002"))
})

testthat::test_that("signed cluster scores and controls match per-cell UCell before averaging", {
  load_scoring_test_runtime()
  counts <- make_counts_matrix()
  markers <- list(A = c("gene001+", "gene003", "gene007-"),
    B = c("gene010", "gene020-"), negative_only = "gene030-")
  metadata <- data.frame(barcode_w_prefix = colnames(counts),
    cluster = rep(c("1", "2"), length.out = ncol(counts)), GEM_well_ID = "well")
  control <- prepare_cluster_UCell_controls(counts, metadata, markers)
  control$n_controls <- 9L
  control$draws <- control$draws[, seq_len(control$n_controls), drop = FALSE]
  control$max_rank <- 80L
  summaries <- summarize_cluster_UCell_counts(counts, metadata, control, "cluster",
    chunk_size = 7L, workers = 1L, include_cell_scores = TRUE)
  reference <- function(signatures) callr::r(function(counts, signatures) {
    as.matrix(UCell::ScoreSignatures_UCell(Matrix::Matrix(counts, sparse = TRUE),
      features = signatures, maxRank = 80L, w_neg = 1, name = "", ncores = 1L))
  }, args = list(counts = counts, signatures = signatures))
  expected <- reference(markers)
  testthat::expect_equal(summaries$cell_scores[rownames(expected), ], expected, tolerance = 1e-14)
  for (label in names(control$markers)) {
    genes <- control$markers[[label]]
    variants <- c(list(genes), if (length(genes) > 1L) lapply(genes, function(gene) setdiff(genes, gene)))
    for (variant in seq_along(variants)) {
      signature <- variants[[variant]]
      controls <- lapply(seq_len(control$n_controls), function(index) {
        selected <- control$reference_genes[control$draws[sub("-$", "", signature), index]]
        paste0(selected, ifelse(grepl("-$", signature), "-", ""))
      })
      signatures <- stats::setNames(c(list(signature), controls), paste0("signature", seq_len(control$n_controls + 1L)))
      scores <- reference(signatures)
      for (group in which(summaries$groups$kind == "cluster")) {
        selected <- metadata$cluster == summaries$groups$cluster[group]
        testthat::expect_equal(unname(summaries$signed_means[[label]][[variant]][, group]),
          unname(colMeans(scores[selected, , drop = FALSE])), tolerance = 1e-14)
      }
    }
  }
  other <- summarize_cluster_UCell_counts(counts, metadata, control, "cluster",
    chunk_size = 11L, workers = 2L, include_cell_scores = TRUE)
  testthat::expect_equal(other$signed_means, summaries$signed_means, tolerance = 1e-14)
  scored <- score_cluster_UCell_summaries(summaries, control)
  testthat::expect_true(all(is.finite(scored$evidence$excess)))
  testthat::expect_true(all(is.finite(scored$marker_deletions$excess)))
  # A signed score clipped after averaging is generally a different quantity.
  gap <- matrix(c(0.9, 0.1, 0.1, 0.9), 2L, dimnames = list(c("gene001", "gene007"), NULL))
  exact <- mean(pmax(0, gap[1, ] - gap[2, ]))
  shortcut <- max(0, mean(gap[1, ]) - mean(gap[2, ]))
  testthat::expect_gt(exact, shortcut)
})
