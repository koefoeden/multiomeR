source_project_file("R/cellranger_reference_helpers.R")

testthat::test_that("aggregation input validation rejects incompatible gene definitions", {
  features <- data.frame(
    id = c("ENSG1", "ENSG2"), name = c("GENE1", "GENE2"),
    seqnames = "chr1", start = c(10, 100), end = c(20, 200)
  )
  resolve <- function(other) resolve_aggregation_gene_features(
    list(features, other), c("reference_well", "other_well"), "mixed"
  )
  testthat::expect_identical(resolve(features), features)
  testthat::expect_error(resolve(features[2:1, ]), "other_well")
  testthat::expect_error(resolve(features[1, ]), "incompatible ordered gene definitions")
  renamed <- features
  renamed$name[1] <- "RENAMED"
  testthat::expect_error(resolve(renamed), "other_well")
  moved <- features
  moved$start[1] <- 11
  testthat::expect_error(resolve(moved), "other_well")
})
