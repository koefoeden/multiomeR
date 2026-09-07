testthat::test_that("blacklist filtering preserves peaks with no hits and removes every hit", {
  load_project_test_runtime()
  peaks <- tempfile(fileext = ".narrowPeak")
  on.exit(unlink(peaks))
  writeLines(c(
    "chr1\t100\t200\tpeak1\t10\t.\t2\t3\t4\t50",
    "chr1\t500\t600\tpeak2\t20\t.\t2\t3\t4\t50",
    "chr1\t900\t1000\tpeak3\t30\t.\t2\t3\t4\t50"
  ), peaks)
  read_peaks <- function(blacklist) {
    get_peak_GRanges_w_fixed_width(peaks, extend_summits = 50, blacklist_GRanges = blacklist)
  }
  unfiltered <- read_peaks(GenomicRanges::GRanges("chr1:2000-2100"))
  testthat::expect_equal(as.character(unfiltered), c("chr1:101-200", "chr1:501-600", "chr1:901-1000"))
  testthat::expect_equal(S4Vectors::mcols(unfiltered)$name, c("peak1", "peak2", "peak3"))
  testthat::expect_identical(read_peaks(GenomicRanges::GRanges()), unfiltered)
  testthat::expect_identical(
    read_peaks(GenomicRanges::GRanges(c("chr1:520-530", "chr1:540-550"))),
    unfiltered[c(1L, 3L)]
  )
  testthat::expect_identical(read_peaks(GenomicRanges::GRanges("chr1:1-1100")), unfiltered[integer()])
})

testthat::test_that("empty peak files are distinct from missing peak files", {
  load_project_test_runtime()
  peaks <- tempfile(fileext = ".narrowPeak")
  testthat::expect_error(
    get_peak_GRanges_w_fixed_width(peaks, blacklist_GRanges = GenomicRanges::GRanges()),
    "Peak file does not exist"
  )
  file.create(peaks)
  on.exit(unlink(peaks))
  testthat::expect_identical(
    get_peak_GRanges_w_fixed_width(peaks, blacklist_GRanges = GenomicRanges::GRanges()),
    empty_peak_GRanges()
  )
})
