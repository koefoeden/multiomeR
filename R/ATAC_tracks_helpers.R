#' Get coverage regions tibble BPCells
#'
#' Prepare ATAC fragment coverage regions as one row per plot.
#'
#' @param formatted_peak_tibble Peak tibble with chromosome, start, end, and
#'   display fields already normalized for track plotting.
#' @param marker_genes Character vector of validated marker-gene names.
#' @param n_top Number of top significant peak regions to include.
#' @param n_bottom Number of bottom significant peak regions to include.
#' @return A tibble with `name` and `region` columns for downstream dynamic branching.
#' @keywords internal

# Marker genes are validated against the gene annotation upstream.
get_coverage_regions_tibble_BPCells <- function(formatted_peak_tibble, marker_genes, n_top = 3, n_bottom = 3) {
  sorted_peak_tibble <- dplyr::arrange(formatted_peak_tibble, dplyr::desc(neg_log10qvalue_summit))

  dplyr::bind_rows(
    sorted_peak_tibble |>
      dplyr::slice_head(n = n_top) |>
      dplyr::transmute(name = paste0("top_", dplyr::row_number(), "_", .data$region_vec), region = .data$region_vec),
    sorted_peak_tibble |>
      dplyr::slice_tail(n = n_bottom) |>
      dplyr::transmute(name = paste0("bottom_", dplyr::row_number(), "_", .data$region_vec), region = .data$region_vec),
    tibble::tibble(name = paste0("marker_", marker_genes), region = marker_genes)
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

# Gene names are unique in the annotation; other region IDs are coordinates.
get_region_GRanges_for_BPCells_track <- function(region_id, gene_GRanges, extend_upstream = 5e4, extend_downstream = 5e4) {
  gene_idx <- which(gene_GRanges$gene_name == region_id)
  if (length(gene_idx) == 1) {
    region <- gene_GRanges[gene_idx]
    GenomicRanges::start(region) <- pmax(1L, GenomicRanges::start(region) - extend_upstream)
    GenomicRanges::end(region) <- GenomicRanges::end(region) + extend_downstream
    return(region)
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
#' @return A combined BPCells trackplot whose reachable data are limited to
#'   summarized coverage bins and peaks overlapping the displayed region. The
#'   plot must not retain the fragment object or full peak table.
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
  group_order <- gtools::mixedsort(unique(as.character(groups)))
  coverage_colors <- stats::setNames(
    grDevices::hcl.colors(length(group_order), palette = "Dark 3"),
    group_order
  )
  collapsed_peaks <- collapsed_peak_tibble |>
    dplyr::transmute(
      chr = as.character(seqnames),
      start = start,
      end = end
    )

  region <- get_region_GRanges_for_BPCells_track(region_id, gene_GRanges = gene_GRanges)
  coverage_tibble <- BPCells::trackplot_coverage(
    fragments = fragments,
    region = region,
    groups = groups,
    cell_read_counts = metadata$atac_fragments,
    group_order = group_order,
    colors = coverage_colors,
    bins = 500,
    return_data = TRUE
  )
  coverage_track <- make_BPCells_ATAC_coverage_track_from_tibble(
    coverage_tibble = coverage_tibble,
    region = region,
    colors = coverage_colors
  )
  collapsed_peaks <- collapsed_peaks |>
    dplyr::filter(
      .data$chr == as.character(GenomeInfoDb::seqnames(region)),
      .data$end > GenomicRanges::start(region),
      .data$start < GenomicRanges::end(region)
    )
  peak_track <- BPCells::trackplot_genome_annotation(
    loci = collapsed_peaks,
    region = region,
    track_label = "Collapsed peaks"
  )

  BPCells::trackplot_combine(
    tracks = list(coverage_track, peak_track),
    title = paste("ATAC coverage:", region_id)
  ) + patchwork::plot_annotation(
    subtitle = stringr::str_wrap("Compare local accessibility across groups; overlapping peaks do not prove regulation of a nearby gene.", width = 100),
    caption = stringr::str_wrap(paste("Coverage uses 500 bins, normalized by bin width and summed group depth; extremes are clipped at the 99.9th percentile.",
      "Depth denominator: total ATAC fragments.",
      "Groups:", label_plot_variable(group_cells_by_col), ". Tracks share the displayed coverage scale. Collapsed peaks are shown below; gene-centred windows extend 50 kb on each side."), width = 110))
}

#' Make BPCells ATAC coverage track from tibble
#'
#' Prepare or plot ATAC fragment coverage over named genomic regions using BPCells trackplot data structures.
#'
#' @param coverage_tibble BPCells trackplot coverage tibble with genomic position, group, and coverage columns.
#' @param region Genomic region accepted by BPCells trackplot helpers; normalized internally to a single plotting interval.
#' @param clip_quantile Upper coverage quantile used to cap extreme track values; set near 1 to keep almost all signal.
#' @param background_layers Optional ggplot layers drawn behind the coverage signal.
#' @param colors Named or positional colors passed to BPCells/ggplot track layers for groups or links.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for saving or composition.
#' @keywords internal

make_BPCells_ATAC_coverage_track_from_tibble <- function(
  coverage_tibble,
  region,
  clip_quantile = 0.999,
  colors = grDevices::hcl.colors(nlevels(coverage_tibble$group), palette = "Dark 3"),
  background_layers = NULL
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
    ggplot2::ggplot(coverage_tibble) + background_layers +
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

get_peak_gene_track_colors <- function(cell_groups, focal_cell_group) {
  colors <- stats::setNames(
    grDevices::hcl.colors(length(cell_groups) + 1L, "Dark 3")[-1], cell_groups
  )
  colors[focal_cell_group] <- grDevices::hcl.colors(1, "Dark 3")
  colors
}

#' Prepare clipped, strand-aware gene bodies with nonoverlapping lanes
prepare_genomic_gene_bodies <- function(gene_ranges, region) {
  selected <- IRanges::subsetByOverlaps(gene_ranges, region)
  tibble::tibble(start = pmax(GenomicRanges::start(selected), GenomicRanges::start(region)),
    end = pmin(GenomicRanges::end(selected), GenomicRanges::end(region)),
    gene = as.character(selected$gene_name), strand = as.character(GenomicRanges::strand(selected))) |>
    dplyr::arrange(.data$start, .data$end) |>
    dplyr::mutate(lane = IRanges::disjointBins(IRanges::IRanges(.data$start, .data$end)))
}

#' Shared strand arrows and labels for genomic gene tracks
genomic_gene_body_layers <- function(genes) {
  list(
    ggplot2::geom_segment(data = genes, ggplot2::aes(
      x = ifelse(strand == "-", end, start), xend = ifelse(strand == "-", start, end),
      y = -lane, yend = -lane), arrow = grid::arrow(length = grid::unit(1.3, "mm")), linewidth = .5),
    ggplot2::geom_text(data = genes, ggplot2::aes((start + end) / 2, -lane, label = gene),
      vjust = -.7, size = 2.5, check_overlap = TRUE)
  )
}
