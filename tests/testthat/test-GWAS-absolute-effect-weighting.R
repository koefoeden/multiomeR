source_project_file("R/celltype_labeling_helpers.R")
source_project_file("R/GWAS_chromVAR_helpers.R")
source_project_file("R/GWAS_chromVAR_absolute_effect_helpers.R")

make_test_record <- function(GWAS_ID, beta) {
  variants <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 4),
    ranges = IRanges::IRanges(start = c(100, 200, 300, 400), width = 1)
  )
  S4Vectors::mcols(variants) <- S4Vectors::DataFrame(
    studyLocusId = c("locus_1", "locus_1", "locus_2", "locus_2"),
    posteriorProbability = c(0.6, 0.3, 0.7, 0.2),
    posteriorProbability_raw = c(0.6, 0.3, 0.7, 0.2),
    beta = beta
  )
  list(
    GWAS_ID = GWAS_ID,
    studyId = GWAS_ID,
    finemappingMethod = "test",
    variant_weighting_mode = "raw_PIP",
    credible_set_GRanges = variants
  )
}

make_absolute_effect_fixture <- function() {
  records <- list(
    make_test_record("variant", c(2, 1, 4, 2)),
    make_test_record("locus", c(2, NA, 4, NA)),
    make_test_record("incomplete", c(2, NA, NA, NA))
  )
  list(
    records = records,
    status = get_GWAS_absolute_effect_weighting_status_tibble(
      GWAS_input_records = records,
      posterior_probability_cutoff = 0.001
    )
  )
}

testthat::test_that("absolute-effect weighting routes reflect beta coverage", {
  fixture <- make_absolute_effect_fixture()
  status <- fixture$status

  testthat::expect_identical(status$effect_weighting_route, c(
    "variant_absolute_effect",
    "locus_absolute_effect",
    NA_character_
  ))
  testthat::expect_identical(status$eligible, c(TRUE, TRUE, FALSE))
  testthat::expect_identical(
    status$status[[3]],
    "skipped_incomplete_locus_beta_coverage"
  )
})

testthat::test_that("variant and locus absolute-effect weights are normalized", {
  fixture <- make_absolute_effect_fixture()
  variant_ranges <- get_GWAS_absolute_effect_variant_GRanges(
    GWAS_input_record = fixture$records[[1]],
    effect_weighting_route = "variant_absolute_effect",
    posterior_probability_cutoff = 0.001
  )
  locus_ranges <- get_GWAS_absolute_effect_variant_GRanges(
    GWAS_input_record = fixture$records[[2]],
    effect_weighting_route = "locus_absolute_effect",
    posterior_probability_cutoff = 0.001
  )

  testthat::expect_equal(
    S4Vectors::mcols(variant_ranges)$absolute_effect_weight,
    c(1.2, 0.3, 2.8, 0.4) / (4.7 / 1.8)
  )
  testthat::expect_equal(
    S4Vectors::mcols(locus_ranges)$absolute_effect_weight,
    c(1.2, 0.6, 2.8, 0.8) / (5.4 / 1.8)
  )
})

testthat::test_that("eligible absolute-effect weights map to peaks and remain bounded", {
  fixture <- make_absolute_effect_fixture()
  peak_ranges <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 4),
    ranges = IRanges::IRanges(start = c(90, 190, 290, 390), width = 21)
  )
  peak_records <- get_GWAS_absolute_effect_peak_weight_records(
    GWAS_input_records = fixture$records,
    weighting_status_tibble = fixture$status,
    peak_ranges = peak_ranges,
    posterior_probability_cutoff = 0.001
  )

  testthat::expect_identical(
    purrr::map_chr(peak_records, "GWAS_ID"),
    c("variant", "locus")
  )
  testthat::expect_equal(
    unname(peak_records[[1]]$peak_weights_vec),
    pmin(c(1.2, 0.3, 2.8, 0.4) / (4.7 / 1.8), 1)
  )
  testthat::expect_equal(
    unname(peak_records[[2]]$peak_weights_vec),
    pmin(c(1.2, 0.6, 2.8, 0.8) / (5.4 / 1.8), 1)
  )
})

source_project_file("R/GWAS_chromVAR_contribution_helpers.R")

testthat::test_that("variant allocations exactly reproduce absolute-effect heatmap peak weights", {
  fixture <- make_absolute_effect_fixture()
  # Include a peak shared by variants and loci, exercising the cap and allocation.
  peaks <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(90, 190, 390), c(310, 210, 410)))
  weights <- get_GWAS_absolute_effect_peak_weight_records(fixture$records, fixture$status, peaks, 0.001)
  allocations <- get_GWAS_absolute_effect_peak_variant_weights(fixture$records, fixture$status, peaks, 0.001)
  for (record in weights) {
    observed <- allocations |>
      dplyr::filter(GWAS_ID == record$GWAS_ID) |>
      dplyr::summarise(weight = sum(peak_variant_weight), .by = peak_name)
    testthat::expect_equal(observed$weight[match(get_peak_names_from_GRanges(peaks), observed$peak_name)],
      unname(record$peak_weights_vec))
  }
  testthat::expect_false("incomplete" %in% allocations$GWAS_ID)
  testthat::expect_equal(nrow(get_GWAS_absolute_effect_peak_variant_weights(
    fixture$records, dplyr::mutate(fixture$status, eligible = FALSE), peaks)), 0L)
})

testthat::test_that("detail loci combine both rankings without duplicates or changing focal cells", {
  ordinary <- tibble::tibble(GWAS_ID = "trait", cluster = "focal", studyLocusId = letters[1:4],
    relative_deviation = 2, relative_deviation_contribution = c(4, 3, 2, 1))
  weighted <- dplyr::mutate(ordinary, relative_deviation_contribution = c(1, 4, 2, 3))
  selected <- select_GWAS_detail_loci(ordinary, weighted, n_top_loci = 2L)
  testthat::expect_setequal(selected$studyLocusId, c("a", "b", "d"))
  testthat::expect_equal(selected$selection[selected$studyLocusId == "b"], "Ordinary + Absolute effect")
  testthat::expect_false(anyDuplicated(selected$studyLocusId) > 0L)
  testthat::expect_setequal(select_GWAS_detail_loci(ordinary, weighted[0, ], n_top_loci = 2L)$studyLocusId, c("a", "b"))
  testthat::expect_s3_class(plot_GWAS_absolute_effect_locus_bars(weighted[0, ]), "empty_plot_list")
})
