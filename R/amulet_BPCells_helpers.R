if (!exists("amulet_BPCells_native_state_env", inherits = FALSE)) {
  amulet_BPCells_native_state_env <- new.env(parent = emptyenv())
  amulet_BPCells_native_state_env$dll_name <- NULL
}

amulet_BPCells_supported_ref <- "28759cdd512578b6cbe549e226e1cd52a2d2308c"

validate_amulet_BPCells_native_abi <- function() {
  BPCells_description <- utils::packageDescription("BPCells")
  installed_ref <- BPCells_description[["RemoteSha"]]
  if (is.null(installed_ref)) {
    installed_ref <- BPCells_description[["GithubSHA1"]]
  }

  if (!identical(installed_ref, amulet_BPCells_supported_ref)) {
    installed_ref_label <- if (is.null(installed_ref)) "unknown" else installed_ref
    stop(
      "The BPCells-native AMULET helper uses a private fragment iterator ABI and only supports BPCells commit ",
      amulet_BPCells_supported_ref,
      ". The installed BPCells commit is ",
      installed_ref_label,
      ". Run `pixi run install-r-github-packages` or revalidate the native interface before updating the pin.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

load_amulet_BPCells_native_library <- function(native_source_file) {
  if (!is.null(amulet_BPCells_native_state_env$dll_name)) {
    return(amulet_BPCells_native_state_env$dll_name)
  }

  validate_amulet_BPCells_native_abi()

  build_dir <- tempfile("multiomeR_amulet_bpcells_")
  dir.create(build_dir)
  build_source_file <- file.path(build_dir, basename(native_source_file))
  if (!file.copy(native_source_file, build_source_file)) {
    stop("Could not copy the BPCells-native AMULET source into the temporary build directory.", call. = FALSE)
  }

  shared_library_file <- file.path(
    build_dir,
    paste0("multiomeR_amulet_bpcells", .Platform$dynlib.ext)
  )
  build_result <- processx::run(
    command = file.path(R.home("bin"), "R"),
    args = c("CMD", "SHLIB", "-o", shared_library_file, build_source_file),
    wd = build_dir,
    echo = FALSE,
    error_on_status = FALSE
  )
  if (build_result$status != 0L) {
    stop(
      "Could not compile the BPCells-native AMULET helper:\n",
      paste(c(build_result$stdout, build_result$stderr), collapse = "\n"),
      call. = FALSE
    )
  }

  loaded_library <- dyn.load(shared_library_file)
  amulet_BPCells_native_state_env$dll_name <- loaded_library[["name"]]
  amulet_BPCells_native_state_env$dll_name
}

iterate_BPCells_fragments <- function(fragments) {
  get("iterate_fragments", envir = asNamespace("BPCells"))(fragments)
}

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
#' @param unique_fragments Must be `TRUE`. BPCells fragment objects do not retain
#'   the Cell Ranger PCR-duplicate count needed to expand non-unique fragments.
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
  unique_fragments = TRUE,
  cellranger_end_inclusive = TRUE,
  return_type = c("stats", "loci"),
  verbose = TRUE,
  native_source_file = file.path(get_project_root(), "src", "amulet_bpcells.cpp")
) {
  return_type <- match.arg(return_type)
  if (!isTRUE(unique_fragments)) {
    stop(
      "`unique_fragments = FALSE` is unsupported because BPCells does not retain Cell Ranger PCR-duplicate counts.",
      call. = FALSE
    )
  }
  if (!methods::is(fragments, "IterableFragments")) {
    stop("`fragments` must be a BPCells IterableFragments object.", call. = FALSE)
  }
  if (!is.null(barcodes) && !is.character(barcodes)) {
    stop("`barcodes` must be NULL or a character vector.", call. = FALSE)
  }

  dll_name <- load_amulet_BPCells_native_library(native_source_file)
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
  if (is.null(cell_names) || is.null(chromosome_names)) {
    stop("BPCells fragments must have known cell and chromosome names.", call. = FALSE)
  }
  if (anyNA(cell_names) || anyNA(chromosome_names)) {
    stop("BPCells fragments cannot contain missing cell or chromosome names.", call. = FALSE)
  }

  if (isTRUE(verbose)) {
    message(format(Sys.time(), "%X"), " - Computing BPCells-native AMULET overlaps")
  }
  fragment_iterator <- iterate_BPCells_fragments(fragments)
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
#'
#' @return A data frame with fragment counts, loci covered by more than two
#'   fragments, and AMULET p- and q-values for each retained cell.
#' @keywords internal

calculate_amulet_metrics_BPCells <- function(
  fragments,
  barcodes = NULL,
  min_fragments = 500L,
  regions_to_exclude = GenomicRanges::GRanges(
    c("M", "chrM", "MT", "X", "Y", "chrX", "chrY"),
    IRanges::IRanges(1L, width = 10^8)
  ),
  max_fragment_size = 1000L,
  remove_high_overlap_sites = TRUE,
  unique_fragments = TRUE,
  cellranger_end_inclusive = TRUE,
  verbose = TRUE,
  native_source_file = file.path(get_project_root(), "src", "amulet_bpcells.cpp")
) {
  metrics_df <- get_amulet_fragment_overlaps_BPCells(
    fragments = fragments,
    barcodes = barcodes,
    min_fragments = min_fragments,
    regions_to_exclude = regions_to_exclude,
    max_fragment_size = max_fragment_size,
    remove_high_overlap_sites = remove_high_overlap_sites,
    unique_fragments = unique_fragments,
    cellranger_end_inclusive = cellranger_end_inclusive,
    return_type = "stats",
    verbose = verbose,
    native_source_file = native_source_file
  )
  metrics_df$p.value <- stats::ppois(
    metrics_df$nAbove2,
    mean(metrics_df$nAbove2),
    lower.tail = FALSE
  )
  metrics_df$q.value <- stats::p.adjust(metrics_df$p.value, method = "BH")
  metrics_df
}
