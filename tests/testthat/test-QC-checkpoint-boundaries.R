testthat::test_that("numbered checkpoint selections stop before later analysis stages", {
  withr::local_dir(multiomeR_project_root)
  load_project_test_runtime()
  manifest <- targets::tar_manifest(
    fields = c(name, command, description), callr_function = NULL,
    envir = globalenv()
  )
  dependencies <- stats::setNames(lapply(manifest$command, function(command) {
    intersect(all.names(str2lang(command)), manifest$name)
  }), manifest$name)
  ancestors <- function(selected) {
    previous <- character()
    while (!setequal(previous, selected)) {
      previous <- selected
      selected <- union(selected, unlist(dependencies[selected], use.names = FALSE))
    }
    selected
  }
  stages <- c(
    "1_pre-aggregation-QC", "2_GEX-PCA-QC", "3_GEX-QC", "4_peak-QC",
    "5_pre-LSI-QC", "6_ATAC-LSI-QC", "7_ATAC-QC", "8_multimodal-QC"
  )
  forbidden <- c(
    "^(PCA_BPCells|LSI_BPCells|WNN_results)",
    "^(PCA_clusters|LSI_BPCells|WNN_results)",
    "^(consensus_peak_BPCells_matrix|LSI_BPCells|WNN_results)",
    "^(QC_filtered_BCs[.]ATAC|LSI_BPCells|WNN_results)",
    "^(LSI_BPCells|WNN_results)",
    "^(LSI_clusters|scDblFinder_results_df[.]ATAC|motif_family_chromVAR|WNN_results)",
    "^WNN_results", "[.]subgroups[.]"
  )
  aggregation_suffixes <- sub("^aggregated_cellranger_ref_list", "", grep(
    "^aggregated_cellranger_ref_list[.]", manifest$name, value = TRUE
  ))
  testthat::expect_gt(length(aggregation_suffixes), 0L)
  for (suffix in aggregation_suffixes) {
    for (index in seq_along(stages)) {
      selected <- manifest$name[
        endsWith(manifest$name, suffix) &
          grepl(paste0("[checkpoint:", stages[index], "]"), manifest$description, fixed = TRUE)
      ]
      testthat::expect_true(length(selected) > 0L, info = paste(stages[index], suffix))
      leaked <- grep(forbidden[index], ancestors(selected), value = TRUE)
      testthat::expect_true(length(leaked) == 0L, info = paste(stages[index], suffix, paste(leaked, collapse = ", ")))
    }
  }
})
