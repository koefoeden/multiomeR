source_project_file("R/peak_gene_correlation_helpers.R")
source_project_file("R/peak_gene_filter_helpers.R")

testthat::test_that("support presets filter hypotheses without requiring co-detection", {
  ids <- paste0("a", seq_len(24))
  depth <- data.frame(aggregate_id = ids, donor_id = rep(letters[1:4], each = 6),
    GEX_depth = rep(100, 24), ATAC_depth = rep(100, 24))
  RNA <- rbind(high = rep(c(10, 10, 10, 0, 0, 0), 4),
    low = rep(5, 24), one_donor = c(rep(10, 6), rep(0, 18)))
  ATAC <- rbind(peak = rep(c(0, 0, 0, 5, 5, 5), 4))
  colnames(RNA) <- colnames(ATAC) <- ids
  raw <- list(cell_group = "test", chr = "chr1", GEX_counts = RNA, ATAC_counts = ATAC,
    aggregate_depth_tibble = depth)
  norm <- normalize_peak_gene_correlation_aggregate_matrices(raw)
  pairs <- data.frame(chr = "chr1", peak = "peak", TargetGeneID = rownames(RNA),
    gene_matrix_feature = rownames(RNA))
  result <- lapply(c("lenient", "moderate", "strict"), function(level)
    filter_peak_gene_candidate_pairs(raw, norm, pairs, level))
  testthat::expect_equal(result[[1]]$candidate_pairs$TargetGeneID, c("high", "low"))
  testthat::expect_equal(result[[2]]$candidate_pairs$TargetGeneID, "high")
  testthat::expect_equal(result[[3]]$candidate_pairs$TargetGeneID, "high")
  testthat::expect_equal(result[[1]]$diagnostics$n_retained, c(2L, 1L, 1L))
  testthat::expect_true(all(result[[1]]$diagnostics$n_low_donor_support >= 1))
  testthat::expect_error(filter_peak_gene_candidate_pairs(raw, norm, pairs, "invalid"),
    "lenient, moderate, or strict")
  empty <- filter_peak_gene_candidate_pairs(raw, norm, pairs[0, ], "lenient")
  testthat::expect_equal(nrow(empty$candidate_pairs), 0L)
  testthat::expect_equal(empty$diagnostics$n_retained, c(0L, 0L, 0L))
})
