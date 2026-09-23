convert_amulet_loci_to_GRanges <- function(loci, cell_names, chromosome_names, selected_cell_names) {
  if (length(loci$start) == 0L) {
    return(GenomicRanges::GRanges(
      seqnames = factor(character(), chromosome_names),
      ranges = IRanges::IRanges(start = integer(), end = integer()),
      name = factor(character(), selected_cell_names)
    ))
  }

  GenomicRanges::GRanges(
    seqnames = factor(chromosome_names[loci$chr], chromosome_names),
    ranges = IRanges::IRanges(start = loci$start, end = loci$end),
    name = factor(cell_names[loci$cell], selected_cell_names)
  )
}

remove_high_overlap_amulet_loci <- function(loci_GRanges, p_value_threshold = 0.01) {
  high_overlap_GRanges <- GenomicRanges::reduce(
    loci_GRanges,
    min.gapwidth = 0L,
    with.revmap = TRUE
  )
  overlap_counts <- lengths(high_overlap_GRanges$revmap)
  high_overlap_GRanges$p_value <- stats::ppois(
    overlap_counts,
    mean(overlap_counts),
    lower.tail = FALSE
  )
  indices_to_remove <- unlist(
    high_overlap_GRanges$revmap[high_overlap_GRanges$p_value < p_value_threshold],
    use.names = FALSE
  )
  if (length(indices_to_remove) > 0L) {
    loci_GRanges <- loci_GRanges[-indices_to_remove]
  }
  loci_GRanges
}

#' Calculate AMULET fragment overlaps from BPCells fragments
#'
#' Stream a BPCells `IterableFragments` object and reproduce the overlap counts
#' returned by `scDblFinder::getFragmentOverlaps()` without materializing all
#' fragments as `GRanges`.
#'
#' @param fragments A BPCells `IterableFragments` object with known chromosome
#'   and cell names.
#' @param barcodes Optional character vector of cell names to retain.
#' @param min_fragments Minimum fragments required when `barcodes` is `NULL`.
#'   Values between zero and one are interpreted as a fraction of all fragments.
#' @param regions_to_exclude `GRanges` of regions to exclude.
#' @param max_fragment_size Maximum fragment size to retain.
#' @param remove_high_overlap_sites Whether to remove loci unexpectedly covered
#'   in many cells, matching scDblFinder's AMULET implementation.
#' @param cellranger_end_inclusive Whether `fragments` was opened with
#'   `BPCells::open_fragments_10x()`'s Cell Ranger default. `TRUE` removes the
#'   extra end-coordinate base before comparison with scDblFinder's BED import.
#' @param return_type Return barcode statistics or the loci covered by more than
#'   two fragments.
#' @param verbose Whether to report progress.
#' @param native_source_file Path to the tracked native C++ implementation.
#'
#' @return A data frame of fragment-overlap statistics or a `GRanges` object.
#' @keywords internal

get_amulet_fragment_overlaps_BPCells <- function(
  fragments,
  barcodes = NULL,
  min_fragments = 500L,
  regions_to_exclude = GenomicRanges::GRanges(
    c("M", "chrM", "MT", "X", "Y", "chrX", "chrY"),
    IRanges::IRanges(1L, width = 10^8)
  ),
  max_fragment_size = 1000L,
  remove_high_overlap_sites = TRUE,
  cellranger_end_inclusive = TRUE,
  return_type = c("stats", "loci"),
  verbose = TRUE,
  native_source_file = file.path(get_project_root(), "src", "amulet_bpcells.cpp")
) {
  return_type <- match.arg(return_type)
  dll_name <- load_native_library(native_source_file, "multiomeR_amulet_bpcells")
  if (isTRUE(cellranger_end_inclusive)) {
    fragments <- BPCells::shift_fragments(fragments, shift_end = -1L)
  }
  fragments <- BPCells::subset_lengths(
    fragments,
    max_len = as.integer(max_fragment_size)
  )
  if (!is.null(regions_to_exclude) && length(regions_to_exclude) > 0L) {
    fragments <- BPCells::select_regions(
      fragments,
      regions_to_exclude,
      invert_selection = TRUE
    )
  }

  cell_names <- BPCells::cellNames(fragments)
  chromosome_names <- BPCells::chrNames(fragments)

  if (isTRUE(verbose)) {
    message(format(Sys.time(), "%X"), " - Computing BPCells-native AMULET overlaps")
  }
  fragment_iterator <- get("iterate_fragments", envir = asNamespace("BPCells"))(fragments)
  fragment_counts <- .Call(
    "multiomeR_bpcells_fragment_counts",
    fragment_iterator,
    as.integer(length(cell_names)),
    PACKAGE = dll_name
  )
  names(fragment_counts) <- cell_names

  nonzero_cells <- fragment_counts > 0L
  if (min_fragments > 0 && min_fragments < 1) {
    min_fragments <- round(min_fragments * sum(fragment_counts))
  }
  selected_cells <- if (is.null(barcodes)) {
    nonzero_cells & fragment_counts >= min_fragments
  } else {
    missing_barcodes <- setdiff(barcodes, cell_names[nonzero_cells])
    if (length(missing_barcodes) > 0L && isTRUE(verbose)) {
      warning(
        length(missing_barcodes),
        " requested barcode(s) are absent from the BPCells fragments.",
        call. = FALSE
      )
    }
    nonzero_cells & cell_names %in% barcodes
  }

  selected_cell_names <- levels(as.factor(cell_names[selected_cells]))
  empty_stats_df <- data.frame(
    nFrags = integer(),
    uniqFrags = integer(),
    nAbove2 = integer(),
    total.nAbove2 = integer(),
    row.names = character()
  )
  if (!any(selected_cells)) {
    if (identical(return_type, "stats")) {
      return(empty_stats_df)
    }
    return(convert_amulet_loci_to_GRanges(
      loci = list(start = integer()),
      cell_names = cell_names,
      chromosome_names = chromosome_names,
      selected_cell_names = selected_cell_names
    ))
  }

  loci <- .Call(
    "multiomeR_bpcells_amulet_loci",
    fragment_iterator,
    as.logical(selected_cells),
    PACKAGE = dll_name
  )
  loci_GRanges <- convert_amulet_loci_to_GRanges(
    loci = loci,
    cell_names = cell_names,
    chromosome_names = chromosome_names,
    selected_cell_names = selected_cell_names
  )
  loci_GRanges <- loci_GRanges[order(
    loci_GRanges$name,
    as.integer(GenomicRanges::seqnames(loci_GRanges)),
    IRanges::start(loci_GRanges),
    IRanges::end(loci_GRanges)
  )]
  if (identical(return_type, "loci")) {
    return(loci_GRanges)
  }

  stats_df <- data.frame(
    nFrags = as.integer(fragment_counts[selected_cell_names]),
    uniqFrags = as.integer(fragment_counts[selected_cell_names]),
    nAbove2 = 0L,
    total.nAbove2 = 0L,
    row.names = selected_cell_names
  )
  total_overlap_counts <- table(loci_GRanges$name)
  stats_df[names(total_overlap_counts), "total.nAbove2"] <- as.integer(total_overlap_counts)

  if (isTRUE(remove_high_overlap_sites) && length(loci_GRanges) > 0L) {
    loci_GRanges <- remove_high_overlap_amulet_loci(loci_GRanges)
  }
  overlap_counts <- table(loci_GRanges$name)
  stats_df[names(overlap_counts), "nAbove2"] <- as.integer(overlap_counts)
  stats_df
}

#' Calculate BPCells-native AMULET metrics
#'
#' @inheritParams get_amulet_fragment_overlaps_BPCells
#' @param ... Further arguments passed to `get_amulet_fragment_overlaps_BPCells()`.
#'
#' @return A data frame with fragment counts, loci covered by more than two
#'   fragments, and AMULET p- and q-values for each retained cell.
#' @keywords internal

calculate_amulet_metrics_BPCells <- function(fragments, ...) {
  metrics_df <- get_amulet_fragment_overlaps_BPCells(fragments, ..., return_type = "stats")
  metrics_df$p.value <- stats::ppois(
    metrics_df$nAbove2,
    mean(metrics_df$nAbove2),
    lower.tail = FALSE
  )
  metrics_df$q.value <- stats::p.adjust(metrics_df$p.value, method = "BH")
  metrics_df
}
