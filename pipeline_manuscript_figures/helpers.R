# Figure 1 uses data targets, never pipeline ggplot RDS files.
figure_1_theme <- function() {
  ggplot2::theme_minimal(base_size = 8, base_family = "NimbusSan") +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 10, hjust = .5),
      plot.title.position = "plot",
      legend.title = ggplot2::element_text(size = 7),
      legend.text = ggplot2::element_text(size = 6))
}

# All mappings use data columns; detach their evaluation environments and the
# ggplot creation environment from transient target inputs before serialization.
compact_figure_1_plot <- function(plot) {
  clean_mapping <- function(mapping) {
    for (name in names(mapping)) {
      if (rlang::is_quosure(mapping[[name]])) {
        mapping[[name]] <- rlang::quo_set_env(mapping[[name]], globalenv())
      }
    }
    mapping
  }
  plot$plot_env <- globalenv()
  plot$mapping <- clean_mapping(plot$mapping)
  for (i in seq_along(plot$layers)) {
    plot$layers[[i]]$mapping <- clean_mapping(plot$layers[[i]]$mapping)
  }
  if (inherits(plot, "patchwork")) {
    plot$patches$plots <- lapply(plot$patches$plots, compact_figure_1_plot)
  }
  plot
}

select_figure_1_link <- function(links) {
  selected <- links |>
    dplyr::filter(.data$is_analyzable_link, is.finite(.data$hierarchical_pvalue)) |>
    dplyr::arrange(.data$hierarchical_pvalue, .data$cell_group, .data$TargetGeneID, .data$peak) |>
    dplyr::slice_head(n = 1)
  if (nrow(selected) != 1L) stop("No estimable peak-gene example available.")
  selected
}

# The figure shows assigned cell types; unassigned clusters stay in the module outputs.
is_figure_1_cell_type <- function(x) !startsWith(as.character(x), "Unassigned")

# Assigned cell type-study pairs enriched at q < 0.05, strongest first.
figure_1_enrichments <- function(scores) {
  scores |>
    dplyr::filter(is_figure_1_cell_type(.data$cluster), .data$z_q < .05) |>
    dplyr::arrange(dplyr::desc(.data$z), .data$GWAS_ID)
}

# Figure 1D shows broadly supported enrichments: q < 0.05 and a top locus below half
# of the summed contribution, which leaves out pairs driven by one or a few loci; the
# configured studies already have at least 20 effective peaks. In decreasing z, a pair
# is shown unless its study is shown already or its cell type shows a study of the
# same category.
select_figure_1_traits <- function(scores, loci) {
  loci <- dplyr::summarise(loci, n_loci = dplyr::n(),
    top_locus_share = max(.data$relative_deviation_contribution) / sum(.data$relative_deviation_contribution),
    .by = c("GWAS_ID", "cluster"))
  candidates <- figure_1_enrichments(scores) |>
    dplyr::left_join(loci, by = c("GWAS_ID", "cluster")) |>
    dplyr::filter(.data$top_locus_share < .5)
  if (!nrow(candidates)) stop("No broadly supported cell-type/GWAS enrichments.")
  selected <- candidates[0, ]
  for (i in seq_len(nrow(candidates))) {
    pair <- candidates[i, ]
    if (!pair$GWAS_ID %in% selected$GWAS_ID &&
      !any(selected$cluster == pair$cluster & selected$Category == pair$Category)) {
      selected <- dplyr::bind_rows(selected, pair)
    }
  }
  selected
}

# D and E show two of them: the strongest, LDL cholesterol in hepatocytes, and the
# strongest outside the liver in an epithelium, colorectal cancer in colonocytes.
select_figure_1_showcase <- function(traits) {
  showcase <- traits[traits$GWAS_ID == "ColorectalCancer_Law2019", ]
  stopifnot(nrow(showcase) == 1L)
  showcase
}

# Sentence case without author, year and ancestry, keeping acronyms such as LDL
# and HbA1c: LungAdenocarcinoma_McKay2017 -> Lung adenocarcinoma.
figure_1_trait_labels <- function(x) {
  x <- gsub("_", " ", sub("_[^_]+[0-9]{4}(_[A-Z]+)?$", "", x))
  x <- gsub("([a-z])(?=[A-Z][A-Za-z])", "\\1 ", x, perl = TRUE)
  gsub(" ([A-Z])(?=[a-z])", " \\L\\1", x, perl = TRUE)
}

# Deviation clusters use hyphens and metadata cell types underscores.
figure_1_cell_type_labels <- function(x) gsub("[-_]", " ", x)

# Cell types as the plural nouns of the manuscript's prose.
figure_1_cell_type_nouns <- function(x) {
  dplyr::coalesce(c(T = "T cells", B = "B cells", NK = "NK cells", Mono = "monocytes")[x],
    paste0(tolower(figure_1_cell_type_labels(x)), "s"))
}

# Lanes from gene bodies widened to their labels' approximate extent, so labels in
# one lane do not overlap; char_fraction is one character's share of the facet width.
figure_1_gene_lanes <- function(genes, window_width, char_fraction = .018) {
  half_label <- (nchar(genes$gene) + 2) * char_fraction * window_width / 2
  middle <- (genes$start + genes$end) / 2
  genes <- dplyr::mutate(genes, label_start = floor(pmin(.data$start, middle - half_label)),
    label_end = ceiling(pmax(.data$end, middle + half_label))) |>
    dplyr::arrange(.data$label_start, .data$label_end)
  genes$lane <- IRanges::disjointBins(IRanges::IRanges(genes$label_start, genes$label_end))
  dplyr::select(genes, -"label_start", -"label_end")
}

# A: the configured WNN marker gene of each cell type whose expression and promoter
# accessibility are most specific to it, from cell-type pseudobulks summed over
# donors, as log1p counts per million scaled across cell types.
figure_1_marker_profiles <- function(markers, pseudobulk_GEX, pseudobulk_ATAC, peaks) {
  markers <- tibble::enframe(markers, name = "cell_type", value = "gene") |>
    tidyr::unnest("gene") |>
    dplyr::filter(!endsWith(.data$gene, "-"))
  promoters <- tibble::as_tibble(as.data.frame(S4Vectors::mcols(peaks))) |>
    dplyr::filter(.data$geneName %in% markers$gene, .data$txType %in% c("promoter", "TSS"),
      .data$region_vec %in% rownames(pseudobulk_ATAC))
  markers <- dplyr::filter(markers, .data$gene %in% intersect(rownames(pseudobulk_GEX), promoters$geneName))
  # Pseudobulk columns are "<cell type, - for _>_<donor>".
  cell_type_CPM <- function(pseudobulk, features) {
    cell_type <- gsub("-", "_", sub("_[^_]+$", "", colnames(pseudobulk)))
    keep <- cell_type %in% markers$cell_type
    totals <- tapply(colSums(pseudobulk[, keep]), cell_type[keep], sum)
    counts <- rowsum(t(as.matrix(pseudobulk[features, keep])), cell_type[keep])
    counts / as.numeric(totals[rownames(counts)]) * 1e6
  }
  RNA <- cell_type_CPM(pseudobulk_GEX, unique(markers$gene))
  ATAC <- t(rowsum(t(cell_type_CPM(pseudobulk_ATAC, promoters$region_vec)), promoters$geneName))
  profiles <- dplyr::bind_rows(
    tibble::as_tibble(RNA, rownames = "cell_type") |> dplyr::mutate(modality = "Gene expression"),
    tibble::as_tibble(ATAC, rownames = "cell_type") |> dplyr::mutate(modality = "Promoter accessibility")) |>
    tidyr::pivot_longer(-c("cell_type", "modality"), names_to = "gene", values_to = "CPM") |>
    dplyr::mutate(z = as.numeric(scale(log1p(.data$CPM))), .by = c("modality", "gene"))
  selected <- markers |>
    dplyr::inner_join(profiles, by = c("cell_type", "gene")) |>
    dplyr::summarise(specificity = sum(.data$z), .by = c("cell_type", "gene")) |>
    dplyr::slice_max(.data$specificity, n = 1, with_ties = FALSE, by = "cell_type")
  cell_types <- unique(markers$cell_type)
  profiles |>
    dplyr::filter(.data$gene %in% selected$gene, .data$cell_type %in% selected$cell_type) |>
    dplyr::mutate(cell_type = factor(.data$cell_type, rev(cell_types)),
      gene = factor(.data$gene, selected$gene[order(match(selected$cell_type, cell_types))]))
}

plot_figure_1_markers <- function(profiles) {
  ggplot2::ggplot(profiles, ggplot2::aes(.data$gene, .data$cell_type, fill = .data$z)) +
    ggplot2::geom_tile(colour = "white", linewidth = .3) +
    ggplot2::facet_wrap(ggplot2::vars(.data$modality)) +
    ggplot2::scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", name = "Scaled\nlog CPM") +
    ggplot2::scale_y_discrete(labels = figure_1_cell_type_labels) + figure_1_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(face = "italic", angle = 90, hjust = 1, vjust = .5, size = 6),
      axis.text.y = ggplot2::element_text(size = 6),
      panel.grid = ggplot2::element_blank(), strip.text = ggplot2::element_text(size = 8)) +
    ggplot2::labs(title = "Cell-type markers in the same nuclei", x = NULL, y = NULL)
}

# B: inspired by top_link_aggregate_scatter_plots. Separate compact genomic
# panel, focal cell type only, selected peak/TSS and link arc; omit the large
# cross-cell-type evidence stack. Zoom to the peak and TSS with 25 kb flanks; retain cached coverage bins.
plot_figure_1_link_tracks <- function(record, coverage) {
  link <- record$link
  bounds <- c(max(1, min(link$start, link$TargetGeneTSS) - 25000),
    max(link$end, link$TargetGeneTSS) + 25000)
  # Other genes whose centred labels would run past the window keep their arrows but
  # lose their labels; half_label estimates a label's half-width as in figure_1_gene_lanes().
  genes <- record$genes |>
    dplyr::filter(.data$end >= bounds[[1]], .data$start <= bounds[[2]]) |>
    dplyr::mutate(start = pmax(.data$start, bounds[[1]]),
      end = pmin(.data$end, bounds[[2]]),
      half_label = (nchar(.data$gene) + 2) * .018 * diff(bounds) / 2,
      gene = dplyr::if_else(.data$gene != link$TargetGene[[1]] &
        ((.data$start + .data$end) / 2 - .data$half_label < bounds[[1]] |
          (.data$start + .data$end) / 2 + .data$half_label > bounds[[2]]), "", .data$gene),
      lane = IRanges::disjointBins(IRanges::IRanges(.data$start, .data$end))) |>
    dplyr::arrange(.data$gene != link$TargetGene[[1]])
  coverage <- dplyr::filter(coverage, .data$pos >= bounds[[1]], .data$pos <= bounds[[2]])
  base <- ggplot2::ggplot() + figure_1_theme() +
    ggplot2::annotate("rect", xmin = link$start, xmax = link$end,
      ymin = -Inf, ymax = Inf, fill = "#B40426", alpha = 0.25) +
    ggplot2::geom_vline(xintercept = link$TargetGeneTSS, linetype = 2, linewidth = 0.3) +
    ggplot2::scale_x_continuous(limits = bounds, expand = c(0, 0),
      labels = scales::label_number(scale = 1e-6, accuracy = .01))
  genes <- base + genomic_gene_body_layers(genes) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = c(.3, .8))) +
    ggplot2::labs(x = NULL, y = "Genes") + ggplot2::guides(y = "none") +
    ggplot2::theme(axis.text.x = ggplot2::element_blank())
  coverage <- dplyr::mutate(coverage,
    normalized_insertions = pmin(.data$normalized_insertions,
      stats::quantile(.data$normalized_insertions, .999)))
  accessibility <- base + ggplot2::geom_area(data = coverage,
    ggplot2::aes(.data$pos, .data$normalized_insertions), fill = "#446B8C") +
    ggplot2::scale_y_continuous(breaks = \(limits) c(0, floor(limits[[2]] / 10) * 10)) +
    ggplot2::labs(x = NULL, y = "ATAC") +
    ggplot2::theme(axis.text.x = ggplot2::element_blank())
  arc <- tibble::tibble(position = seq(link$peak_center, link$TargetGeneTSS, length.out = 101),
    height = sin(seq(0, pi, length.out = 101)))
  connection <- base + ggplot2::geom_line(data = arc,
    ggplot2::aes(.data$position, .data$height), colour = "#B40426", linewidth = .5) +
    ggplot2::labs(x = paste0(link$chr, " (Mb)"), y = "Link") + ggplot2::guides(y = "none")
  patchwork::wrap_plots(genes, accessibility, connection, ncol = 1,
    heights = c(1.4, 1, .5)) + patchwork::plot_annotation(
      title = paste(link$TargetGene, "in", figure_1_cell_type_nouns(link$cell_group)),
      theme = figure_1_theme())
}

# The nominal P-value as a plotmath label, "Nominal P = 1.4 \u00d7 10^-4".
figure_1_p_label <- function(p) {
  exponent <- floor(log10(p))
  if (exponent >= -2) return(sprintf("Nominal~italic(P)==%.2g", p))
  sprintf("Nominal~italic(P)==%.1f%%*%%10^%d", p / 10^exponent, exponent)
}

# C: same residuals as the parent scatter, separate panel; donor colours without
# a large legend. The line is descriptive, not the hierarchical model fit.
plot_figure_1_scatter <- function(data) {
  ggplot2::ggplot(data, ggplot2::aes(.data$peak_accessibility_residual,
    .data$gene_expression_residual)) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$donor_id), size = .9, alpha = .7) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
      colour = "black", linewidth = .4) + figure_1_theme() +
    ggplot2::annotate("text", x = -Inf, y = Inf, hjust = -.08, vjust = 1.3,
      label = figure_1_p_label(data$hierarchical_pvalue[[1]]), parse = TRUE,
      family = "NimbusSan", size = 2.5) +
    ggplot2::guides(colour = "none") +
    ggplot2::labs(title = "Donor\u2013state pseudobulks",
      x = "ATAC residual (log1p CPM)", y = "GEX residual\n(log1p CPM)")
}

# D: raw_deviation_unscaled inspiration; the selected traits, no metadata tracks or
# clustering. All assigned cell types retained, common raw-deviation scale. Stars
# mark enrichments at q < 0.05.
# Columns follow the configured cell types, as in A, and rows the column of each
# trait's selected cell type, so the selected pairs form a staircase.
plot_figure_1_heatmap <- function(data, traits, cell_types) {
  cell_types <- intersect(cell_types, gsub("-", "_", data$cluster))
  traits <- traits[order(match(gsub("-", "_", traits$cluster), cell_types)), ]
  data <- data |> dplyr::mutate(cluster = gsub("-", "_", .data$cluster)) |>
    dplyr::filter(.data$GWAS_ID %in% traits$GWAS_ID, .data$cluster %in% cell_types) |>
    dplyr::mutate(GWAS_ID = factor(.data$GWAS_ID, levels = rev(traits$GWAS_ID)),
      cluster = factor(.data$cluster, levels = cell_types),
      label = dplyr::if_else(.data$z_q < .05, "*", ""))
  limit <- max(abs(data$deviation), na.rm = TRUE)
  ggplot2::ggplot(data, ggplot2::aes(.data$cluster, .data$GWAS_ID, fill = .data$deviation)) +
    ggplot2::geom_tile(colour = "white", linewidth = .3) +
    ggplot2::geom_text(ggplot2::aes(label = .data$label), size = 3) +
    ggplot2::scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426",
      limits = c(-limit, limit), name = "Deviation") +
    ggplot2::scale_x_discrete(labels = figure_1_cell_type_labels) +
    ggplot2::scale_y_discrete(labels = figure_1_trait_labels) + figure_1_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1, size = 6),
      panel.grid = ggplot2::element_blank(), legend.position = "right") +
    ggplot2::labs(title = "Trait-associated accessibility", x = NULL, y = NULL)
}

# Small tracks show 0 and the largest pretty break within the expanded range.
figure_1_track_breaks <- function(limits) {
  breaks <- scales::breaks_pretty(n = 2)(limits)
  c(0, max(breaks[breaks <= limits[[2]]]))
}

# E: aligned locus facets with independent y scales and shared axis titles
# and independent genomic x ranges. PIP retains the same scale in both columns.
plot_figure_1_loci <- function(records) {
  labels <- vapply(records, function(x) paste0(
    figure_1_trait_labels(x$selection$GWAS_ID[[1]]), " — ",
    figure_1_cell_type_labels(x$selection$cluster[[1]])), character(1))
  bounds <- purrr::map2_dfr(records, labels, function(x, label) tibble::tibble(
    locus = label, position = c(GenomicRanges::start(x$record$region), GenomicRanges::end(x$record$region))))
  variants <- purrr::map2_dfr(records, labels, function(x, label) {
    dplyr::mutate(x$record$variant_tibble, locus = label)
  })
  genes <- purrr::map2_dfr(records, labels, function(x, label) {
    figure_1_gene_lanes(x$genes, GenomicRanges::width(x$record$region)) |> dplyr::mutate(locus = label)
  })
  coverage <- purrr::map2_dfr(records, labels, function(x, label) {
    x$record$coverage_tibble |>
      dplyr::filter(.data$group == x$selection$cluster[[1]]) |>
      dplyr::mutate(locus = label)
  })
  peaks <- purrr::map2_dfr(records, labels, function(x, label) {
    tibble::as_tibble(as.data.frame(x$record$consensus_peak_GRanges)) |>
      dplyr::mutate(locus = label)
  })
  bounds$locus <- factor(bounds$locus, levels = labels)
  variants$locus <- factor(variants$locus, levels = labels)
  genes$locus <- factor(genes$locus, levels = labels)
  coverage$locus <- factor(coverage$locus, levels = labels)
  peaks$locus <- factor(peaks$locus, levels = labels)
  chromosomes <- stats::setNames(vapply(records, function(x)
    as.character(GenomicRanges::seqnames(x$record$region)), character(1)), labels)
  base <- ggplot2::ggplot() + figure_1_theme() +
    ggplot2::geom_blank(data = bounds, ggplot2::aes(.data$position, 0)) +
    ggplot2::facet_wrap(ggplot2::vars(locus), nrow = 1, scales = "free") +
    ggplot2::scale_x_continuous(expand = c(0, 0),
      labels = scales::label_number(scale = 1e-6, accuracy = .01)) +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 8), strip.clip = "off",
      axis.title.y = ggplot2::element_text(size = 6.5, angle = 0, hjust = 1, vjust = .5),
      axis.text.y = ggplot2::element_text(family = "mono"),
      panel.spacing.x = grid::unit(6, "mm"), legend.position = "right")
  contribution <- base +
    ggplot2::geom_hline(yintercept = 0, colour = "grey70") +
    ggplot2::geom_segment(data = variants, ggplot2::aes(x = .data$position,
      xend = .data$position, y = 0, yend = .data$relative_deviation_contribution), colour = "#B40426") +
    ggplot2::geom_point(data = variants, ggplot2::aes(.data$position,
      .data$relative_deviation_contribution, size = .data$posteriorProbability), colour = "#B40426") +
    ggplot2::scale_size_area(max_size = 3, name = "PIP", limits = c(0, 1), breaks = c(.1, .5, 1)) +
    ggplot2::scale_y_continuous(breaks = figure_1_track_breaks,
      labels = function(x) formatC(x, format = "f", digits = 1, width = 5),
      expand = ggplot2::expansion(mult = c(.05, .35))) +
    ggplot2::labs(x = NULL, y = "Genetic\ncontribution") +
    ggplot2::theme(axis.text.x = ggplot2::element_blank())
  # A window without protein-coding genes names the locus's Open Targets locus-to-gene prediction.
  empty <- dplyr::summarise(bounds, position = mean(.data$position), .by = "locus") |>
    dplyr::mutate(label = paste0("No protein-coding genes in this window\nLocus-to-gene prediction: ",
      vapply(records, \(x) x$locus$likely_genes_label[[1]], character(1)))) |>
    dplyr::filter(!.data$locus %in% genes$locus)
  gene_track <- base + genomic_gene_body_layers(genes) +
    ggplot2::geom_text(data = empty, ggplot2::aes(.data$position, -1, label = .data$label),
      size = 2.5, lineheight = .9, vjust = 0) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = c(.3, .8))) +
    ggplot2::labs(x = NULL, y = "Genes") + ggplot2::guides(y = "none") +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), strip.text = ggplot2::element_blank())
  accessibility <- base +
    ggplot2::geom_area(data = coverage, ggplot2::aes(.data$pos, .data$normalized_insertions), fill = "#446B8C") +
    ggplot2::geom_segment(data = peaks, ggplot2::aes(x = .data$start, xend = .data$end,
      y = 0, yend = 0), linewidth = 1, colour = "#B40426") +
    ggplot2::facet_wrap(ggplot2::vars(locus), nrow = 1, scales = "free", strip.position = "bottom",
      labeller = ggplot2::labeller(locus = ggplot2::as_labeller(chromosomes))) +
    ggplot2::scale_y_continuous(breaks = figure_1_track_breaks,
      labels = function(x) formatC(x, format = "f", digits = 1, width = 5)) +
    ggplot2::labs(x = "Genomic position (Mb)", y = "ATAC") +
    ggplot2::theme(strip.placement = "outside")
  # Keep each axis label beside its track while aligning the faceted panels.
  patchwork::wrap_plots(lapply(list(contribution, gene_track, accessibility),
    function(plot) patchwork::free(plot, type = "label", side = "l")),
    ncol = 1, heights = c(1, 1.2, 1), guides = "collect") +
    patchwork::plot_annotation(title = "Leading GWAS locus", theme = figure_1_theme())
}

# F: TRS_UMAPs inspiration; rebuild WNN UMAP from coordinates and scores, with
# WNN cell-type centroids, as in D, no score clipping or serialized upstream plot mutation.
# Only the cell types enriched for the trait are labelled, each placed outward from
# the embedding's centre with a line to its nuclei's median position.
plot_figure_1_TRS <- function(data, trait, cell_types) {
  labels <- data |> dplyr::summarise(x = stats::median(.data$WNN_UMAP_1),
    y = stats::median(.data$WNN_UMAP_2), .by = "WNN_harmony_SNN_cluster_cell_type") |>
    dplyr::filter(.data$WNN_harmony_SNN_cluster_cell_type %in% cell_types) |>
    dplyr::mutate(label = figure_1_cell_type_labels(.data$WNN_harmony_SNN_cluster_cell_type),
      label_x = .data$x + .15 * diff(range(data$WNN_UMAP_1)) * sign(.data$x - stats::median(data$WNN_UMAP_1)),
      label_y = .data$y + .12 * diff(range(data$WNN_UMAP_2)) * sign(.data$y - stats::median(data$WNN_UMAP_2)))
  ggplot2::ggplot(data, ggplot2::aes(.data$WNN_UMAP_1, .data$WNN_UMAP_2)) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$score), size = .15) +
    ggplot2::geom_segment(data = labels, ggplot2::aes(.data$x, .data$y, xend = .data$label_x, yend = .data$label_y),
      linewidth = .3, colour = "grey20") +
    ggplot2::geom_label(data = labels, ggplot2::aes(.data$label_x, .data$label_y, label = .data$label),
      size = 2.8, label.size = .2, label.padding = grid::unit(1.5, "pt"), fill = "white") +
    ggplot2::scale_colour_viridis_c(name = "TRS", na.value = "grey85") +
    figure_1_theme() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text = ggplot2::element_blank(),
      legend.position = "right") +
    ggplot2::guides(colour = ggplot2::guide_colourbar(barwidth = grid::unit(2, "mm"),
      barheight = grid::unit(25, "mm"))) +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::labs(title = "SCAVENGE trait relevance", subtitle = figure_1_trait_labels(trait),
      x = "UMAP 1", y = "UMAP 2")
}

# Supplementary Figures S5 and S6: a module heatmap of every scored study and assigned cell
# type, with Figure 1's cell-type labels and stars at q < 0.05 on the cell type's `value`.
# Values beyond the limit take the end colours. The metadata tracks must hold exactly the
# scored studies. Of them, only trait category, whose rows also name the study, and the
# listed `bar_tracks` are kept; Supplementary Table S3 gives each study's sample size,
# method and ancestry.
plot_trait_cell_type_heatmap <- function(scores, GWAS_metadata_tracks_plot, value_label, limit, bar_tracks = character()) {
  bars <- which(vapply(seq_along(GWAS_metadata_tracks_plot), \(i) "track" %in% names(GWAS_metadata_tracks_plot[[i]]$data), logical(1)))
  stopifnot(length(bars) == 1L)
  category <- GWAS_metadata_tracks_plot[[1]] +
    ggplot2::scale_y_discrete(drop = FALSE, expand = c(0, 0),
      labels = \(x) paste0(figure_1_trait_labels(x), " (", sub("^.*_([A-Za-z]+)([0-9]{4}).*$", "\\1 \\2", x), ")"))
  tracks <- list(category)
  if (length(bar_tracks)) {
    tracks[[2]] <- GWAS_metadata_tracks_plot[[bars]]
    tracks[[2]]$data <- dplyr::filter(tracks[[2]]$data, .data$track %in% bar_tracks)
  }
  ordered <- scores |>
    dplyr::mutate(cluster = figure_1_cell_type_labels(.data$cluster),
      value = pmax(pmin(.data$value, limit), -limit), label = dplyr::if_else(.data$q < .05, "*", "")) |>
    get_ordered_GWAS_score_plot_data(NULL, score_col = "value")
  heatmap <- plot_GWAS_feature_heatmap(ordered$scores, ordered$clusters, "cluster", "value", value_label,
    fill_limits = c(-limit, limit), support_label_col = "label")
  n <- length(tracks)
  patchwork::wrap_plots(c(tracks, list(heatmap, plot_GWAS_feature_support_tracks(ordered$clusters))),
    design = paste0(strrep("#", n), LETTERS[n + 2], "\n", paste(LETTERS[seq_len(n + 1)], collapse = "")),
    widths = c(.35, 1.6, 9)[c(seq_len(n), 3)], heights = grid::unit(c(24, 1), c("mm", "null")),
    guides = "collect") &
    ggplot2::theme(legend.position = "bottom", legend.box = "horizontal", legend.justification = "left")
}

# Supplementary Table S3: each configured study's Open Targets record and credible
# sets, with the loci and variants overlapping consensus peaks and the peaks they weight.
GWAS_study_tibble <- function(inputs, variant_weights, study_dataset_path) {
  studies <- arrow::open_dataset(study_dataset_path) |>
    dplyr::filter(studyId %in% inputs$studyId) |>
    dplyr::select("studyId", "pubmedId", "traitFromSource", "initialSampleSize") |>
    dplyr::collect()
  ancestry <- inputs |>
    tidyr::pivot_longer(tidyselect::starts_with("ancestry_"), names_to = "ancestry", names_prefix = "ancestry_") |>
    dplyr::filter(.data$value > 0) |>
    dplyr::arrange(dplyr::desc(.data$value)) |>
    dplyr::summarise(ancestry = paste0(.data$ancestry, " ", round(100 * .data$value), "%", collapse = ", "),
      .by = "GWAS_ID")
  overlaps <- dplyr::bind_rows(variant_weights) |>
    dplyr::summarise(overlapping_loci = dplyr::n_distinct(.data$studyLocusId),
      overlapping_variants = dplyr::n_distinct(.data$variantId),
      weighted_peaks = dplyr::n_distinct(.data$peak_name[.data$peak_weight > 0]), .by = "GWAS_ID")
  inputs |>
    dplyr::left_join(studies, by = "studyId", relationship = "one-to-one") |>
    dplyr::left_join(ancestry, by = "GWAS_ID") |>
    dplyr::left_join(overlaps, by = "GWAS_ID") |>
    dplyr::transmute(.data$GWAS_ID, study_id = .data$studyId, pubmed_id = .data$pubmedId,
      trait = .data$traitFromSource, category = .data$Category, source = .data$sourceType,
      release = .data$sourceRelease, finemapping_method = .data$finemappingMethod, .data$ancestry,
      discovery_sample = .data$initialSampleSize, .data$sample_size, credible_set_loci = .data$n_credible_set_loci,
      .data$overlapping_loci, .data$overlapping_variants, .data$weighted_peaks)
}

# Wrapping every panel prevents one panel's long axis labels from shifting the
# headings in other rows. Tags belong to the assembly, not the plot titles.
compose_figure_1 <- function(A, B, C, D, E, F) {
  plots <- lapply(list(A, B, C, D, E, F), function(plot) {
    if (inherits(plot, "patchwork")) {
      title <- plot$patches$annotation$title
      plot <- plot + patchwork::plot_annotation(title = NULL)
    } else {
      title <- plot$labels$title
      plot <- plot + ggplot2::labs(title = NULL)
    }
    patchwork::wrap_elements(panel = plot) + ggplot2::labs(title = title) +
      ggplot2::theme(plot.title = ggplot2::element_text(family = "NimbusSan", size = 10,
          hjust = .5, margin = ggplot2::margin(b = 5.5)),
        plot.title.position = "plot")
  })
  (patchwork::wrap_plots(plots, design = "AA\nBC\nDD\nEF",
    heights = c(2.3, 1.6, 1.5, 2.7)) +
    patchwork::plot_annotation(tag_levels = "A")) &
    ggplot2::theme(plot.tag = ggplot2::element_text(family = "NimbusSan", face = "bold",
        size = 12, hjust = 0, vjust = 1), plot.tag.position = c(0, 1))
}
