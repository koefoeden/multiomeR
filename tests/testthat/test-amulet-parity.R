load_amulet_test_runtime <- function() {
  if (!exists("calculate_amulet_metrics_BPCells", mode = "function")) {
    load_project_test_runtime()
  }
}

# scDblFinder 1.24.0 emits stack-imbalance warnings under R 4.5, so every reference
# result is computed in one disposable process and compared in this session.
run_reference_amulet <- function(calls) {
  require_reference_version("scDblFinder", "1.24.0")
  callr::r(
    function(calls) lapply(calls, function(call) do.call(
      getExportedValue("scDblFinder", call$fun), c(list(x = call$x, verbose = FALSE), call$arguments)
    )),
    args = list(calls = calls),
    show = FALSE
  )
}

make_multichromosome_fragments <- function() {
  set.seed(20260712)
  n_fragments <- 12000L
  cell_ids <- sprintf("cell%03d", sample.int(40L, n_fragments, replace = TRUE))
  chromosomes <- sample(c("chr1", "chr2", "chrM"), n_fragments, replace = TRUE, prob = c(0.47, 0.47, 0.06))
  starts <- sample.int(200000L, n_fragments, replace = TRUE)
  hotspot_indices <- seq_len(800L)
  cell_ids[hotspot_indices] <- sample(sprintf("cell%03d", seq_len(8L)), length(hotspot_indices), replace = TRUE)
  chromosomes[hotspot_indices] <- "chr1"
  starts[hotspot_indices] <- sample(10000:10100, length(hotspot_indices), replace = TRUE)
  GenomicRanges::GRanges(
    seqnames = chromosomes,
    ranges = IRanges::IRanges(start = starts, width = sample.int(500L, n_fragments, replace = TRUE)),
    name = cell_ids,
    cell_id = cell_ids
  )
}

testthat::test_that("integration: AMULET loci and metrics match scDblFinder", {
  load_amulet_test_runtime()
  fragment_file <- system.file("extdata", "example_fragments.tsv.gz", package = "scDblFinder")
  fragments <- BPCells::write_fragments_dir(
    BPCells::open_fragments_10x(fragment_file), tempfile("amulet_bpcells_reference_"), overwrite = TRUE
  )
  barcodes <- BPCells::cellNames(fragments)
  synthetic_GRanges <- make_multichromosome_fragments()
  synthetic_fragments <- BPCells::write_fragments_dir(
    BPCells::convert_to_fragments(synthetic_GRanges), tempfile("amulet_bpcells_multichromosome_"), overwrite = TRUE
  )
  expected <- run_reference_amulet(list(
    bundled = list(fun = "amulet", x = fragment_file, arguments = list()),
    production = list(fun = "amulet", x = fragment_file, arguments = list(barcodes = barcodes)),
    loci = list(fun = "getFragmentOverlaps", x = synthetic_GRanges, arguments = list(
      regionsToExclude = NULL, minFrags = 0L, removeHighOverlapSites = FALSE, ret = "loci")),
    synthetic = list(fun = "amulet", x = synthetic_GRanges, arguments = list(
      regionsToExclude = NULL, minFrags = 0L, removeHighOverlapSites = TRUE))
  ))

  testthat::expect_identical(
    calculate_amulet_metrics_BPCells(fragments = fragments, verbose = FALSE),
    expected$bundled,
    info = "Cell Ranger fragment-file metrics"
  )
  # Production scores the prefixed Cell Ranger-called barcodes without a fragment minimum.
  prefix <- "GEM_well_1_"
  production <- calculate_amulet_metrics_BPCells(
    fragments = BPCells::prefix_cell_names(fragments, prefix),
    barcodes = paste0(prefix, barcodes),
    verbose = FALSE
  )
  rownames(production) <- base::substring(rownames(production), base::nchar(prefix) + 1L)
  testthat::expect_identical(production, expected$production, info = "prefixed production settings")
  testthat::expect_identical(
    get_amulet_fragment_overlaps_BPCells(
      fragments = synthetic_fragments, regions_to_exclude = NULL, min_fragments = 0L,
      remove_high_overlap_sites = FALSE, cellranger_end_inclusive = FALSE, return_type = "loci", verbose = FALSE
    ),
    expected$loci,
    info = "multi-chromosome loci and ordering"
  )
  testthat::expect_identical(
    calculate_amulet_metrics_BPCells(
      fragments = synthetic_fragments, regions_to_exclude = NULL, min_fragments = 0L,
      remove_high_overlap_sites = TRUE, cellranger_end_inclusive = FALSE, verbose = FALSE
    ),
    expected$synthetic,
    info = "multi-chromosome AMULET metrics"
  )
})

testthat::test_that("integration: PCR-duplicate expansion fails explicitly", {
  load_amulet_test_runtime()
  fragments <- BPCells::open_fragments_10x(
    system.file("extdata", "example_fragments.tsv.gz", package = "scDblFinder")
  )
  testthat::expect_error(
    get_amulet_fragment_overlaps_BPCells(
      fragments = fragments,
      unique_fragments = FALSE,
      verbose = FALSE
    ),
    "PCR-duplicate counts"
  )
})
