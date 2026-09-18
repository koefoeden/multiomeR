# Genomic context and rendering for hierarchical top-link diagnostics.

prepare_peak_gene_top_link_plot_data <- function(link, context, gene_ranges, consensus_peaks) {
  if (!isTRUE(link$is_analyzable_link[[1]])) return(list(link = link))
  context <- context |> dplyr::filter(.data$TargetGeneID == link$TargetGeneID[[1]],
    .data$chr == link$chr[[1]])
  region <- GenomicRanges::GRanges(link$chr[[1]], IRanges::IRanges(
    max(1L, min(context$start, link$start, link$TargetGeneTSS) - 2000L),
    max(context$end, link$end, link$TargetGeneTSS) + 2000L))
  peaks <- IRanges::subsetByOverlaps(consensus_peaks, region)
  list(link = link, region = region,
    genes = prepare_genomic_gene_bodies(gene_ranges, region),
    peaks = tibble::tibble(start = pmax(GenomicRanges::start(peaks), GenomicRanges::start(region)),
      end = pmin(GenomicRanges::end(peaks), GenomicRanges::end(region))),
    evidence = context |> dplyr::filter(is.finite(.data$hierarchical_pvalue),
      .data$hierarchical_pvalue >= 0, .data$hierarchical_pvalue <= 1) |>
      dplyr::mutate(evidence = -log10(pmax(.data$hierarchical_pvalue, .Machine$double.xmin))))
}

compute_peak_gene_focal_coverage <- function(plot_data, metadata, fragments) {
  if (!isTRUE(plot_data$link$is_analyzable_link[[1]])) return(tibble::tibble())
  focal_group <- plot_data$link$cell_group[[1]]
  metadata <- metadata |>
    dplyr::filter(.data$WNN_harmony_SNN_cluster_cell_type == focal_group) |>
    dplyr::distinct(.data$barcode_w_prefix, .keep_all = TRUE)
  fragment_names <- BPCells::cellNames(fragments)
  metadata <- metadata |> dplyr::filter(.data$barcode_w_prefix %in% fragment_names) |>
    dplyr::arrange(match(.data$barcode_w_prefix, fragment_names))
  fragments <- BPCells::select_cells(fragments, metadata$barcode_w_prefix)
  BPCells::trackplot_coverage(fragments = fragments, region = plot_data$region,
    groups = rep(focal_group, nrow(metadata)), cell_read_counts = metadata$atac_fragments,
    group_order = focal_group, bins = 500, return_data = TRUE)
}

plot_peak_gene_correlation_top_link <- function(plot_data, coverage, aggregate_scatter, cell_types) {
  top <- plot_data$link
  if (!isTRUE(top$is_analyzable_link[[1]])) {
    return(make_empty_peak_gene_correlation_plot("No candidate enhancer-gene links"))
  }
  region <- plot_data$region
  lo <- GenomicRanges::start(region)
  hi <- GenomicRanges::end(region)
  fits <- plot_data$evidence
  cell_colors <- get_peak_gene_track_colors(cell_types, top$cell_group[[1]])
  peak_shading <- ggplot2::geom_rect(data = dplyr::filter(plot_data$peaks,
      !(.data$start == top$start[[1]] & .data$end == top$end[[1]])),
    ggplot2::aes(xmin = start, xmax = end), ymin = -Inf, ymax = Inf,
    inherit.aes = FALSE, fill = "grey70", alpha = .45, colour = NA)
  focal_shading <- ggplot2::annotate("rect", xmin = top$start[[1]], xmax = top$end[[1]],
    ymin = -Inf, ymax = Inf, fill = "red", alpha = .4, colour = NA)
  background_layers <- list(peak_shading, focal_shading)
  tss_line <- ggplot2::geom_vline(xintercept = top$TargetGeneTSS[[1]],
    linetype = "dashed", colour = "black", linewidth = .45)
  gene_plot <- ggplot2::ggplot() + background_layers +
    genomic_gene_body_layers(plot_data$genes) + tss_line +
    ggplot2::scale_x_continuous(limits = c(lo, hi), expand = c(0, 0),
      position = "top", labels = scales::label_number(scale = 1e-6, accuracy = .01)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = c(.3, .7))) +
    ggplot2::labs(x = "Genomic position (Mb)", y = "Gene bodies") +
    ggplot2::guides(y = "none") + ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  coverage_plot <- make_BPCells_ATAC_coverage_track_from_tibble(coverage, region,
    colors = cell_colors[top$cell_group[[1]]], background_layers = background_layers) +
    tss_line + ggplot2::labs(x = NULL) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(5.5, 5.5, 0, 5.5))
  bar_plot <- ggplot2::ggplot(dplyr::arrange(fits, dplyr::desc(evidence))) +
    ggplot2::geom_rect(ggplot2::aes(xmin = start, xmax = end, ymin = 0,
      ymax = evidence, fill = cell_group, colour = cell_group),
      alpha = .45, linewidth = .2) +
    ggplot2::scale_fill_manual(values = cell_colors, name = "Cell type") +
    ggplot2::scale_colour_manual(values = cell_colors, guide = "none") +
    ggplot2::scale_x_continuous(limits = c(lo, hi), expand = c(0, 0),
      labels = scales::label_number()) +
    ggplot2::scale_y_reverse(limits = c(NA, 0),
      expand = ggplot2::expansion(mult = c(.03, 0))) +
    ggplot2::labs(x = NULL, y = expression(-log[10](p["hierarchical"]))) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(0, 5.5, 5.5, 5.5))
  scatter_plot <- plot_peak_gene_correlation_aggregate_scatter(aggregate_scatter)
  patchwork::wrap_plots(gene_plot, coverage_plot, bar_plot, scatter_plot,
    ncol = 1, heights = c(1, .8, 4, 3), guides = "keep") +
    patchwork::plot_annotation(
      title = paste(top$TargetGene[[1]], "— peak–gene associations across the candidate window"),
      subtitle = stringr::str_wrap(sprintf(
        "%s links (%s in %s) with valid hierarchical p-values. Bars extend downward from zero; longer bars indicate stronger evidence. Overlaid colours compare cell types at each peak.",
        nrow(fits), sum(fits$cell_group == top$cell_group[[1]]), top$cell_group[[1]]), 145),
      caption = stringr::str_wrap(paste(
        "Bars: -log10 of hierarchical nominal p-values; peak intervals are overlaid, not stacked. Unavailable p-values are omitted; no p-value, slope or promoter cutoff is applied to contextual bars.",
        "Grey backgrounds: consensus peaks; red: focal peak; dashed line: focal gene TSS. Shared candidate window with 2 kb padding. Coverage: focal WNN cell type, 500 bins normalized by bin width and group depth, clipped at the 99.9th percentile.",
        "Pooled coverage and statistical associations do not establish donor replication, physical contacts or causality."), 145)
    ) & ggplot2::theme(legend.position = "bottom", legend.box = "vertical")
}
