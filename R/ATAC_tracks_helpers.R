#' Visualize peak calling
#'
#' Import BED peaks and BigWig cut-site coverage for selected chromosomes and draw Gviz tracks over the requested window.
#'
#' @param chroms Chromosome names to plot, including the `chr` prefix expected by Gviz and imported genomic tracks.
#' @param peaks_bed_file BED/narrowPeak file imported as a peak annotation track.
#' @param cuts_BW_file BigWig cut-site coverage file imported as a quantitative
#'   coverage track.
#' @param from 1-based genomic start coordinate for each chromosome window shown in the diagnostic plot.
#' @param to Genomic end coordinate for each chromosome window shown in the diagnostic plot.
#' @param do_plot Logical; when TRUE, open a PDF device and write the Gviz tracks as a side effect.
#' @param name Suffix appended to the optional diagnostic PDF filename.
#' @return A list of Gviz track plot return values, one per chromosome; optionally also writes a PDF.
#' @keywords internal

visualize_peak_calling <- function(
  chroms = c(1:22, "X") %>% paste0("chr", .),
  peaks_bed_file,
  cuts_BW_file,
  from = 100E3,
  to = 400E3,
  do_plot = TRUE,
  name = ""
) {
  options(ucscChromosomeNames = FALSE)
  if (do_plot) {
    grDevices::pdf(paste0("peak_calling", name, ".pdf"))
  }

  plots <- purrr::map(
    chroms,
    ~ {
      bed_data <- rtracklayer::import(peaks_bed_file, format = "bed", which = GenomicRanges::GRanges(.x, IRanges::IRanges(from, to)))
      bed_track <- Gviz::AnnotationTrack(bed_data, name = "Peaks", genome = "hg38")

      bw_data <- rtracklayer::import(cuts_BW_file, format = "bigWig", which = GenomicRanges::GRanges(.x, IRanges::IRanges(from, to)))
      bw_track <- Gviz::DataTrack(bw_data, type = "l", name = "Cutsites", genome = "hg38")

      tracks_list <- list(bed_track, bw_track)
      track_plot <- Gviz::plotTracks(
        trackList = tracks_list,
        from = from,
        to = to,
        chromosome = .x,
        main = paste0(.x, ":", as.integer(from), "-", as.integer(to))
      )
      return(track_plot)
    }
  )
  if (do_plot) {
    (grDevices::dev.off())
  }
  return(plots)
}

#' Get coverage regions tibble BPCells
#'
#' Prepare ATAC fragment coverage regions as one row per plot.
#'
#' @param formatted_peak_tibble Peak tibble with chromosome, start, end, and
#'   display fields already normalized for track plotting.
#' @param gene_GRanges Gene annotation GRanges with gene-name metadata used to resolve gene-centered regions.
#' @param marker_genes Character vector of validated marker-gene names.
#' @param n_top Number of top significant peak regions to include.
#' @param n_bottom Number of bottom significant peak regions to include.
#' @return A tibble with `name` and `region` columns for downstream dynamic branching.
#' @keywords internal

get_coverage_regions_tibble_BPCells <- function(formatted_peak_tibble, gene_GRanges, marker_genes, n_top = 3, n_bottom = 3) {
  sorted_peak_tibble <- dplyr::arrange(formatted_peak_tibble, dplyr::desc(neg_log10qvalue_summit))

  marker_genes_requested <- marker_genes
  gene_names <- if ("gene_name" %in% colnames(S4Vectors::mcols(gene_GRanges))) {
    as.character(gene_GRanges$gene_name)
  } else {
    as.character(names(gene_GRanges))
  }
  missing_marker_genes <- setdiff(marker_genes_requested, gene_names)
  if (length(marker_genes_requested) == 0) {
    stop("No coverage marker genes were configured.", call. = FALSE)
  }
  if (length(missing_marker_genes) > 0) {
    stop(
      "Coverage marker gene(s) were not found in the gene annotation names: ",
      paste(missing_marker_genes, collapse = ", "),
      call. = FALSE
    )
  }

  dplyr::bind_rows(
    sorted_peak_tibble |>
      dplyr::slice_head(n = n_top) |>
      dplyr::transmute(name = paste0("top_", dplyr::row_number(), "_", .data$region_vec), region = .data$region_vec),
    sorted_peak_tibble |>
      dplyr::slice_tail(n = n_bottom) |>
      dplyr::transmute(name = paste0("bottom_", dplyr::row_number(), "_", .data$region_vec), region = .data$region_vec),
    tibble::tibble(name = paste0("marker_", marker_genes_requested), region = marker_genes_requested)
  ) |>
    dplyr::distinct(.data$name, .keep_all = TRUE)
}

#' Get region GRanges for BPCells track
#'
#' Construct or transform genomic ranges with coordinates compatible with downstream ATAC helpers.
#'
#' @param region_id Region identifier, either a peak string like `chr1:1-100` or a gene symbol resolvable in `gene_GRanges`.
#' @param gene_GRanges Gene annotation GRanges with gene-name metadata used to resolve gene-centered regions.
#' @param extend_upstream Number of bases to extend upstream when converting a gene name into a plotting interval.
#' @param extend_downstream Number of bases to extend downstream when converting a gene name into a plotting interval.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for saving or composition.
#' @keywords internal

get_region_GRanges_for_BPCells_track <- function(region_id, gene_GRanges = NULL, extend_upstream = 5e4, extend_downstream = 5e4) {
  if (!is.null(gene_GRanges)) {
    gene_names <- if ("gene_name" %in% colnames(S4Vectors::mcols(gene_GRanges))) {
      as.character(gene_GRanges$gene_name)
    } else {
      as.character(names(gene_GRanges))
    }
    gene_idx <- which(gene_names == region_id)
  }
  if (!is.null(gene_GRanges) && length(gene_idx) == 1) {
    region <- gene_GRanges[gene_idx]
    GenomicRanges::start(region) <- pmax(1L, GenomicRanges::start(region) - extend_upstream)
    GenomicRanges::end(region) <- GenomicRanges::end(region) + extend_downstream
    return(region)
  }
  if (!is.null(gene_GRanges) && length(gene_idx) > 1) {
    stop("Coverage region gene identifier matches multiple genes: ", region_id, call. = FALSE)
  }

  region_parts <- stringr::str_match(region_id, "^([^:-]+)[:-]([0-9,]+)-([0-9,]+)$")
  if (anyNA(region_parts)) {
    stop("Could not parse coverage region: ", region_id)
  }

  GenomicRanges::GRanges(
    seqnames = region_parts[, 2],
    ranges = IRanges::IRanges(
      start = as.integer(stringr::str_remove_all(region_parts[, 3], ",")),
      end = as.integer(stringr::str_remove_all(region_parts[, 4], ","))
    )
  )
}

#' Plot coverage at region BPCells
#'
#' Plot BPCells fragment coverage tracks across one requested genomic region.
#'
#' @param fragments BPCells fragment object, already subsettable by cell barcode and genomic interval.
#' @param metadata_tibble Tibble with one row per cell or pseudobulk sample; must contain the barcode/grouping columns referenced by the helper arguments.
#' @param region_id Region ID resolved before coverage is summarized.
#' @param collapsed_peak_tibble Peak tibble with collapsed region coordinates used for track overlays.
#' @param gene_GRanges Gene annotation GRanges with gene-name metadata used to resolve gene-centered regions.
#' @param group_cells_by_col Metadata column used to split cells into separate
#'   coverage tracks.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for saving or composition.
#' @keywords internal

plot_coverage_at_region_BPCells <- function(
  fragments,
  metadata_tibble,
  region_id,
  collapsed_peak_tibble,
  gene_GRanges,
  group_cells_by_col = "LSI_harmony_SNN_cluster_cell_type"
) {
  fragments <- BPCells::select_cells(fragments, metadata_tibble$barcode_w_prefix)
  fragment_cell_names <- BPCells::cellNames(fragments)
  metadata <- metadata_tibble |>
    dplyr::filter(.data$barcode_w_prefix %in% fragment_cell_names) |>
    dplyr::arrange(match(.data$barcode_w_prefix, fragment_cell_names))

  groups <- metadata[[group_cells_by_col]]
  cell_read_counts <- if ("atac_fragments" %in% colnames(metadata)) metadata$atac_fragments else metadata$nCount_ATAC

  collapsed_peaks <- collapsed_peak_tibble |>
    dplyr::transmute(
      chr = as.character(seqnames),
      start = start,
      end = end
    )
  peak_color_col <- NULL
  if ("cluster" %in% colnames(collapsed_peak_tibble)) {
    collapsed_peaks <- collapsed_peaks |>
      dplyr::mutate(cluster = as.character(collapsed_peak_tibble$cluster))
    peak_color_col <- "cluster"
  }

  region <- get_region_GRanges_for_BPCells_track(region_id, gene_GRanges = gene_GRanges)
  coverage_track <- BPCells::trackplot_coverage(
    fragments = fragments,
    region = region,
    groups = groups,
    cell_read_counts = cell_read_counts,
    group_order = gtools::mixedsort(unique(as.character(groups))),
    bins = 500
  )
  peak_track <- BPCells::trackplot_genome_annotation(
    loci = collapsed_peaks,
    region = region,
    color_by = peak_color_col,
    track_label = "Collapsed peaks"
  )

  BPCells::trackplot_combine(
    tracks = list(coverage_track, peak_track),
    title = region_id
  )
}

#' Make BPCells ATAC coverage track from tibble
#'
#' Prepare or plot ATAC fragment coverage over named genomic regions using BPCells trackplot data structures.
#'
#' @param coverage_tibble BPCells trackplot coverage tibble with genomic position, group, and coverage columns.
#' @param region Genomic region accepted by BPCells trackplot helpers; normalized internally to a single plotting interval.
#' @param clip_quantile Upper coverage quantile used to cap extreme track values; set near 1 to keep almost all signal.
#' @param colors Named or positional colors passed to BPCells/ggplot track layers for groups or links.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for saving or composition.
#' @keywords internal

make_BPCells_ATAC_coverage_track_from_tibble <- function(
  coverage_tibble,
  region,
  clip_quantile = 0.999,
  colors = BPCells:::discrete_palette("stallion")
) {
  region <- BPCells:::normalize_ranges(region)
  ymax <- stats::quantile(coverage_tibble$normalized_insertions, clip_quantile)
  ymax_accuracy <- 10^floor(log10(0.01 * ymax))
  range_label <- sprintf(
    "[0-%s]",
    scales::label_comma(accuracy = ymax_accuracy, big.mark = " ")(ymax)
  )
  coverage_tibble <- coverage_tibble |>
    dplyr::mutate(normalized_insertions = pmin(.data$normalized_insertions, ymax))

  if (is.null(names(colors))) {
    names(colors) <- levels(coverage_tibble$group)
  }
  colors <- colors[seq_len(length(levels(coverage_tibble$group)))]

  BPCells:::wrap_trackplot(
    ggplot2::ggplot(coverage_tibble) +
      ggplot2::geom_area(
        ggplot2::aes(
          x = .data$pos,
          y = .data$normalized_insertions,
          fill = .data$group
        )
      ) +
      ggplot2::scale_fill_manual(values = colors, drop = FALSE) +
      ggplot2::scale_x_continuous(
        limits = c(region$start, region$end),
        expand = c(0, 0),
        labels = scales::label_number()
      ) +
      ggplot2::scale_y_continuous(
        limits = c(0, ymax),
        expand = c(0, 0)
      ) +
      ggplot2::annotate(
        "text",
        x = region$start,
        y = ymax,
        label = range_label,
        vjust = 1.5,
        hjust = -0.1,
        size = 11 * 0.8 / ggplot2::.pt
      ) +
      ggplot2::labs(
        x = "Genomic Position (bp)",
        y = "Insertions (RPKM)"
      ) +
      ggplot2::guides(y = "none", fill = "none") +
      ggplot2::facet_wrap("group", ncol = 1, strip.position = "left") +
      BPCells:::trackplot_theme(),
    ggplot2::unit(length(levels(coverage_tibble$group)), "null"),
    takes_sideplot = TRUE,
    region = region
  )
}

#' Make BPCells peak gene loop track from tibble
#'
#' Prepare peak-gene loop coordinates for BPCells track plotting.
#'
#' @param loop_data_tibble Tibble of peak-gene loop records with source/target coordinates and scores.
#' @param region Genomic region accepted by BPCells trackplot helpers; normalized internally to a single plotting interval.
#' @param track_label Text label shown on the rendered BPCells track.
#' @param color_label Metadata column or legend label used to map loop/link colors.
#' @param colors Named or positional colors passed to BPCells/ggplot track layers for groups or links.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for saving or composition.
#' @keywords internal

make_BPCells_peak_gene_loop_track_from_tibble <- function(
  loop_data_tibble,
  region,
  track_label = "Peak-gene links",
  color_label = "cell_group",
  colors = BPCells:::discrete_palette("tableau", length(levels(loop_data_tibble$color)))
) {
  region <- BPCells:::normalize_ranges(region)
  if (nrow(loop_data_tibble) == 0L) {
    return(BPCells:::trackplot_empty(region, track_label))
  }

  ymin <- loop_data_tibble |>
    dplyr::filter(.data$end <= region$end, .data$start >= region$start) |>
    dplyr::summarise(ymin = min(.data$y), .groups = "drop") |>
    dplyr::pull(.data$ymin)
  loop_data_tibble <- loop_data_tibble |>
    dplyr::mutate(
      y = pmax(.data$y, 1.05 * ymin),
      x = pmax(region$start, pmin(region$end, .data$x))
    )

  BPCells:::wrap_trackplot(
    ggplot2::ggplot(
      loop_data_tibble,
      ggplot2::aes(
        x = .data$x,
        y = .data$y,
        group = .data$loop_id,
        color = .data$color,
        linewidth = .data$abs_correlation
      )
    ) +
      ggplot2::geom_line() +
      ggplot2::scale_color_manual(values = colors) +
      ggplot2::scale_linewidth(range = c(0.2, 1.5), guide = "none") +
      ggplot2::scale_x_continuous(
        limits = c(region$start, region$end),
        expand = c(0, 0),
        labels = scales::label_number()
      ) +
      ggplot2::scale_y_continuous(
        labels = NULL,
        breaks = NULL,
        expand = c(0.05, 0, 0, 0)
      ) +
      ggplot2::guides(size = "none") +
      ggplot2::labs(
        x = "Genomic Position (bp)",
        y = NULL,
        color = color_label
      ) +
      ggplot2::facet_wrap("facet_label", strip.position = "left") +
      BPCells:::trackplot_theme(),
    ggplot2::unit(1, "null"),
    region = region
  )
}
