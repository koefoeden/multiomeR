source_project_file("R/peak_gene_finemapping_helpers.R")
source_project_file("R/peak_gene_hierarchical_helpers.R")

testthat::test_that("the fine-mapping reference restores the FDR of every pair of a linked gene", {
  set.seed(42)
  results <- tidyr::expand_grid(cell_group = c("A", "B"), chr = c("chr1", "chr2"), gene = 1:5, peak = paste0("p", 1:20)) |>
    dplyr::mutate(
      TargetGeneID = paste0(.data$chr, "_g", .data$gene), link_class = "distal", isSelfPromoter = FALSE,
      hierarchical_status = "estimable", hierarchical_coefficient = stats::rnorm(dplyr::n()),
      hierarchical_gene_z = stats::rnorm(dplyr::n()), hierarchical_FDR = stats::runif(dplyr::n())^3
    )
  linked_genes <- unique(paste(results$cell_group, results$TargetGeneID)[is_peak_gene_link(results)])
  testthat::expect_gt(length(linked_genes), 0L)
  reference <- make_peak_gene_finemapping_reference(results)
  for (branch in dplyr::group_split(results, .data$cell_group, .data$chr)) {
    # Branch results also hold background rows, which the restored result drops.
    input <- dplyr::bind_rows(dplyr::select(branch, -"hierarchical_FDR"), dplyr::mutate(branch[1, ], link_class = "background", hierarchical_FDR = NULL))
    observed <- restore_peak_gene_correlation_FDR(input, reference)
    linked <- paste(branch$cell_group, branch$TargetGeneID) %in% linked_genes
    testthat::expect_identical(observed$hierarchical_FDR, dplyr::if_else(linked, branch$hierarchical_FDR, NA_real_))
  }
})
