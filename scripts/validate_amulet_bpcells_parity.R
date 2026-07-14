#!/usr/bin/env Rscript

if (!exists("calculate_amulet_metrics_BPCells", mode = "function")) {
  source("R/bootstrap_helpers.R")
  load_project_runtime()
}

fail <- function(...) {
  stop(paste0(...), call. = FALSE)
}

expect_identical_result <- function(observed, expected, label) {
  if (!identical(observed, expected)) {
    comparison <- paste(utils::capture.output(all.equal(observed, expected)), collapse = "\n")
    fail(label, " is not identical to scDblFinder 1.24.0:\n", comparison)
  }
  message("PASS: ", label)
  invisible(TRUE)
}

scDblFinder_reference_version <- as.character(utils::packageVersion("scDblFinder"))
if (!identical(scDblFinder_reference_version, "1.24.0")) {
  fail(
    "scDblFinder reference version changed from 1.24.0 to ",
    scDblFinder_reference_version,
    "; review exact AMULET parity before updating the pin"
  )
}
validate_amulet_BPCells_native_abi()

run_reference_amulet <- function(fragment_input, arguments = list()) {
  callr::r(
    function(fragment_input, arguments) {
      do.call(
        scDblFinder::amulet,
        c(list(x = fragment_input, verbose = FALSE), arguments)
      )
    },
    args = list(fragment_input = fragment_input, arguments = arguments),
    show = FALSE
  )
}

run_reference_fragment_overlaps <- function(fragment_GRanges, arguments = list()) {
  callr::r(
    function(fragment_GRanges, arguments) {
      do.call(
        scDblFinder::getFragmentOverlaps,
        c(list(x = fragment_GRanges, verbose = FALSE), arguments)
      )
    },
    args = list(fragment_GRanges = fragment_GRanges, arguments = arguments),
    show = FALSE
  )
}

fragment_file <- system.file(
  "extdata",
  "example_fragments.tsv.gz",
  package = "scDblFinder"
)
fragment_dir <- tempfile("amulet_bpcells_reference_")
fragments <- BPCells::open_fragments_10x(fragment_file)
fragments <- BPCells::write_fragments_dir(
  fragments,
  fragment_dir,
  overwrite = TRUE
)
expected_metrics <- run_reference_amulet(fragment_file)
observed_metrics <- calculate_amulet_metrics_BPCells(
  fragments = fragments,
  verbose = FALSE
)
expect_identical_result(
  observed = observed_metrics,
  expected = expected_metrics,
  label = "Cell Ranger fragment-file metrics"
)

prefixed_fragments <- BPCells::prefix_cell_names(fragments, "GEM_well_1_")
prefixed_barcodes <- paste0("GEM_well_1_", rownames(expected_metrics))
observed_prefixed_metrics <- calculate_amulet_metrics_BPCells(
  fragments = prefixed_fragments,
  barcodes = prefixed_barcodes,
  verbose = FALSE
)
observed_pipeline_tibble <- observed_prefixed_metrics |>
  dplyr::rename_with(
    ~ base::paste0("amulet_", .x),
    .cols = 2:dplyr::last_col()
  ) |>
  tibble::rownames_to_column("barcode") |>
  dplyr::mutate(
    barcode = base::substring(.data$barcode, base::nchar("GEM_well_1_") + 1L)
  )
expected_pipeline_tibble <- expected_metrics |>
  dplyr::rename_with(
    ~ base::paste0("amulet_", .x),
    .cols = 2:dplyr::last_col()
  ) |>
  tibble::rownames_to_column("barcode")
expect_identical_result(
  observed = observed_pipeline_tibble,
  expected = expected_pipeline_tibble,
  label = "amulet_metrics_tibble output contract"
)
rownames(observed_prefixed_metrics) <- base::substring(
  rownames(observed_prefixed_metrics),
  base::nchar("GEM_well_1_") + 1L
)
expect_identical_result(
  observed = observed_prefixed_metrics,
  expected = expected_metrics,
  label = "prefixed pipeline-fragment metrics"
)

set.seed(20260712)
n_fragments <- 12000L
cell_ids <- sprintf("cell%03d", sample.int(40L, n_fragments, replace = TRUE))
chromosomes <- sample(
  c("chr1", "chr2", "chrM"),
  n_fragments,
  replace = TRUE,
  prob = c(0.47, 0.47, 0.06)
)
starts <- sample.int(200000L, n_fragments, replace = TRUE)
hotspot_indices <- seq_len(800L)
cell_ids[hotspot_indices] <- sample(
  sprintf("cell%03d", seq_len(8L)),
  length(hotspot_indices),
  replace = TRUE
)
chromosomes[hotspot_indices] <- "chr1"
starts[hotspot_indices] <- sample(10000:10100, length(hotspot_indices), replace = TRUE)
fragment_widths <- sample.int(500L, n_fragments, replace = TRUE)

synthetic_GRanges <- GenomicRanges::GRanges(
  seqnames = chromosomes,
  ranges = IRanges::IRanges(start = starts, width = fragment_widths),
  name = cell_ids,
  cell_id = cell_ids
)
synthetic_fragments <- BPCells::convert_to_fragments(synthetic_GRanges)
synthetic_fragment_dir <- tempfile("amulet_bpcells_multichromosome_")
synthetic_fragments <- BPCells::write_fragments_dir(
  synthetic_fragments,
  synthetic_fragment_dir,
  overwrite = TRUE
)

expected_loci <- run_reference_fragment_overlaps(
  fragment_GRanges = synthetic_GRanges,
  arguments = list(
    regionsToExclude = NULL,
    minFrags = 0L,
    removeHighOverlapSites = FALSE,
    ret = "loci"
  )
)
observed_loci <- get_amulet_fragment_overlaps_BPCells(
  fragments = synthetic_fragments,
  regions_to_exclude = NULL,
  min_fragments = 0L,
  remove_high_overlap_sites = FALSE,
  cellranger_end_inclusive = FALSE,
  return_type = "loci",
  verbose = FALSE
)
expect_identical_result(
  observed = observed_loci,
  expected = expected_loci,
  label = "multi-chromosome loci and ordering"
)

expected_synthetic_metrics <- run_reference_amulet(
  fragment_input = synthetic_GRanges,
  arguments = list(
    regionsToExclude = NULL,
    minFrags = 0L,
    removeHighOverlapSites = TRUE
  )
)
observed_synthetic_metrics <- calculate_amulet_metrics_BPCells(
  fragments = synthetic_fragments,
  regions_to_exclude = NULL,
  min_fragments = 0L,
  remove_high_overlap_sites = TRUE,
  cellranger_end_inclusive = FALSE,
  verbose = FALSE
)
expect_identical_result(
  observed = observed_synthetic_metrics,
  expected = expected_synthetic_metrics,
  label = "multi-chromosome AMULET metrics"
)

unsupported_error <- tryCatch(
  {
    get_amulet_fragment_overlaps_BPCells(
      fragments = fragments,
      unique_fragments = FALSE,
      verbose = FALSE
    )
    NULL
  },
  error = function(error) conditionMessage(error)
)
if (
  is.null(unsupported_error) ||
    !grepl("PCR-duplicate counts", unsupported_error, fixed = TRUE)
) {
  fail("`unique_fragments = FALSE` did not fail with the expected explicit error")
}
message("PASS: unsupported PCR-duplicate expansion fails explicitly")

cat("BPCells-native AMULET parity validation passed.\n")
