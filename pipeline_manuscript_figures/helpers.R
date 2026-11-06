# Figure 1 uses data targets, never pipeline ggplot RDS files. Its font is DejaVu Sans with
# its bold face, which targets that save a figure register first (register_benchmark_font()).
# One type scale: panel titles 10 pt, axis titles and subtitles 7 pt, heatmap row and
# column names 5.5 pt, so that they fit their rows, and all other text 6 pt, which text
# layers get as figure_1_text_size (in mm).
figure_1_text_size <- 6 / ggplot2::.pt
figure_1_theme <- function() {
  ggplot2::theme_minimal(base_size = 7, base_family = BENCHMARK_FONT) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 10, hjust = .5),
      plot.title.position = "plot",
      plot.subtitle = ggplot2::element_text(size = 7),
      axis.text = ggplot2::element_text(size = 6),
      strip.text = ggplot2::element_text(size = 7),
      legend.title = ggplot2::element_text(size = 6),
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

# G and H show colorectal cancer, which lacks genome-wide summary statistics, in its
# most strongly enriched cell type, with its loci that have variants in peaks and the
# largest locus's share of their summed contribution.
select_figure_1_showcase <- function(scores, loci) {
  showcase <- figure_1_enrichments(scores) |>
    dplyr::filter(.data$GWAS_ID == "ColorectalCancer_Law2019") |>
    dplyr::slice_head(n = 1)
  stopifnot(nrow(showcase) == 1L)
  contribution <- dplyr::semi_join(loci, showcase, by = c("GWAS_ID", "cluster"))$relative_deviation_contribution
  dplyr::mutate(showcase, n_loci = length(contribution), top_locus_share = max(contribution) / sum(contribution))
}

# Sentence case without author, year and ancestry, keeping acronyms such as LDL
# and HbA1c: LungAdenocarcinoma_McKay2017 -> Lung adenocarcinoma.
figure_1_trait_labels <- function(x) {
  x <- gsub("_", " ", sub("_[^_]+[0-9]{4}(_[A-Z]+)?$", "", x))
  x <- gsub("([a-z])(?=[A-Z][A-Za-z])", "\\1 ", x, perl = TRUE)
  x <- gsub(" ([A-Z])(?=[a-z])", " \\L\\1", x, perl = TRUE)
  # The spirometry ratio reads FEV1/FVC.
  sub("FEV1FVC", "FEV1/FVC", x, fixed = TRUE)
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

# A: multiomeR read from left to right, without a title, on a 0-100 by 0-20 canvas:
# inputs, the primary module, the optional modules (dashed) and outputs, with the key
# ideas underneath. Grey letters name the panels that show each module's outputs.
plot_figure_1_overview <- function() {
  boxes <- tibble::tribble(
    ~xmin, ~xmax, ~ymin, ~ymax, ~label, ~kind, ~panels,
    0, 15, 10, 18.5, "Cell Ranger ARC\noutputs, many samples", "data", "",
    0, 15, 0, 8.5, "Plain-text\nconfiguration", "data", "",
    21, 49, 0, 18.5, "", "primary", "B, C",
    56, 79, 13, 18.5, "Differential analyses", "optional", "",
    56, 79, 6.5, 12, "Peak\u2013gene correlation", "optional", "D, E",
    56, 79, 0, 5.5, "Genetic enrichment", "optional", "F\u2013H",
    85, 100, 4, 14.5, "Results, plots,\nSeurat/Signac objects", "data", "")
  steps <- tibble::tibble(x = c(25.5, 35, 44.5), label = c("Quality\ncontrol", "WNN\nintegration", "Cell\ntypes"))
  arrows <- tibble::tribble(
    ~x, ~xend, ~y, ~yend,
    15, 21, 14.25, 14.25,
    15, 21, 4.25, 4.25,
    49, 56, 9.25, 15.75,
    49, 56, 9.25, 9.25,
    49, 56, 9.25, 2.75,
    79, 85, 9.25, 9.25,
    steps$x[1] + 4, steps$x[2] - 4.5, 7.6, 7.6,
    steps$x[2] + 4.5, steps$x[3] - 3.5, 7.6, 7.6)
  text <- \(...) ggplot2::annotate("text", family = BENCHMARK_FONT, lineheight = .9, ...)
  ggplot2::ggplot() +
    ggplot2::geom_rect(data = boxes, ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin, ymax = .data$ymax,
      fill = .data$kind, linetype = .data$kind), colour = "#446B8C", linewidth = .4) +
    ggplot2::geom_text(data = boxes, ggplot2::aes((.data$xmin + .data$xmax) / 2, (.data$ymin + .data$ymax) / 2, label = .data$label),
      family = BENCHMARK_FONT, size = 2.3, lineheight = .9) +
    ggplot2::geom_text(data = boxes, ggplot2::aes(.data$xmax - .5, .data$ymax - .7, label = .data$panels),
      family = BENCHMARK_FONT, fontface = "bold", size = 2.1, colour = "grey50", hjust = 1, vjust = 1) +
    text(x = 35, y = 16.4, label = "Primary module", fontface = "bold", size = 2.5) +
    text(x = 35, y = 13.4, label = "expression and chromatin accessibility", size = 2.1, colour = "grey30") +
    ggplot2::geom_label(data = steps, ggplot2::aes(.data$x, 7.6, label = .data$label), family = BENCHMARK_FONT, size = 2.1,
      lineheight = .85, fill = "white", linewidth = .25, label.padding = grid::unit(1.5, "pt")) +
    text(x = 35, y = 2.6, label = "review checkpoint after each stage", size = 2.1, colour = "grey30") +
    text(x = 67.5, y = 20.6, label = "Optional modules, donors as replicates", size = 2.1, colour = "grey30") +
    ggplot2::geom_segment(data = arrows, ggplot2::aes(.data$x, .data$y, xend = .data$xend, yend = .data$yend),
      arrow = grid::arrow(length = grid::unit(1.2, "mm"), type = "closed"), linewidth = .35, colour = "grey30") +
    text(x = 50, y = -3.2, size = 2.3, fontface = "italic",
      label = "Every choice recorded  \u00b7  a change reruns only the steps it affects  \u00b7  data stay on disk") +
    ggplot2::scale_fill_manual(values = c(data = "grey95", primary = "#E3EAF1", optional = "white"), guide = "none") +
    ggplot2::scale_linetype_manual(values = c(data = "solid", primary = "solid", optional = "22"), guide = "none") +
    ggplot2::coord_cartesian(xlim = c(0, 100), ylim = c(-4.8, 22), expand = FALSE, clip = "off") +
    ggplot2::theme_void(base_family = BENCHMARK_FONT) +
    ggplot2::theme(plot.margin = ggplot2::margin(t = 2, l = 14))
}

# The bins of a cell type's main island: the connected group of grid bins in which it
# is the majority, touching at edges or corners, that holds most of its nuclei.
figure_1_main_island <- function(bins) {
  key <- paste(bins$bin_x, bins$bin_y)
  edges <- tidyr::expand_grid(from = seq_along(key), dx = -1:1, dy = -1:1) |>
    dplyr::mutate(to = match(paste(bins$bin_x[.data$from] + .data$dx, bins$bin_y[.data$from] + .data$dy), key)) |>
    dplyr::filter(!is.na(.data$to))
  island <- igraph::components(igraph::graph_from_data_frame(edges[c("from", "to")], vertices = seq_along(key)))$membership
  bins[island == which.max(tapply(bins$n, island, sum)), ]
}

# Places each label box, in bin units, on its island if the island is at least as wide
# as half the label and as tall as it, or beside it: right, left, above or below, and
# further out only if those are taken. A box scores 10 per bin of other islands it
# covers, 3 per bin of other islands within two bins of it, 3 per bin of its own
# island when beside it, 4 per bin of distance, 500 if another label's island is
# nearer than its own and it has no leader line, and 1,000 for another label or the
# panel edge. The most constrained label is placed first. Returns each label's centre
# and leader line.
figure_1_place_labels <- function(islands, occupied, limits, width, height) {
  grid <- matrix(0L, diff(limits[1:2]) + 1, diff(limits[3:4]) + 1)
  index <- \(bins) cbind(bins$bin_x - limits[[1]] + 1, bins$bin_y - limits[[3]] + 1)
  grid[index(occupied)] <- 1L
  candidates <- purrr::imap_dfr(islands, \(island, label) {
    centre_x <- stats::weighted.mean(island$bin_x, island$n)
    centre_y <- stats::weighted.mean(island$bin_y, island$n)
    w <- width[[label]]
    h <- height
    rows <- island$bin_y[abs(island$bin_y - centre_y) <= h / 2]
    columns <- island$bin_x[abs(island$bin_x - centre_x) <= w / 2]
    right <- max(island$bin_x[island$bin_y %in% rows], island$bin_x[!length(rows)])
    left <- min(island$bin_x[island$bin_y %in% rows], island$bin_x[!length(rows)])
    top <- max(island$bin_y[island$bin_x %in% columns], island$bin_y[!length(columns)])
    bottom <- min(island$bin_y[island$bin_x %in% columns], island$bin_y[!length(columns)])
    fits <- right - left + 1 >= w / 2 && top - bottom + 1 >= h
    tidyr::expand_grid(side = c("on", "right", "left", "top", "bottom"), shift = c(0, -.5, .5),
      distance = c(0, 1, 2, 4, 7, 10)) |>
      dplyr::filter(.data$side %in% c("top", "bottom") | .data$shift == 0, .data$side != "on" | (fits & .data$distance == 0)) |>
      dplyr::mutate(label = label,
        x = dplyr::case_when(.data$side == "right" ~ right + 1 + .data$distance + w / 2,
          .data$side == "left" ~ left - 1 - .data$distance - w / 2, .default = centre_x + .data$shift * w),
        y = dplyr::case_when(.data$side == "top" ~ top + 1 + .data$distance + h / 2,
          .data$side == "bottom" ~ bottom - 1 - .data$distance - h / 2, .default = centre_y),
        anchor_x = dplyr::case_when(.data$side == "right" ~ right, .data$side == "left" ~ left, .default = centre_x),
        anchor_y = dplyr::case_when(.data$side == "top" ~ top, .data$side == "bottom" ~ bottom, .default = centre_y),
        w = w, h = h)
  })
  covered <- \(x, y, w, h, label, side, distance) {
    xs <- floor(x - w / 2):ceiling(x + w / 2) - limits[[1]] + 1
    ys <- floor(y - h / 2):ceiling(y + h / 2) - limits[[3]] + 1
    if (min(xs) < 1 || max(xs) > nrow(grid) || min(ys) < 1 || max(ys) > ncol(grid)) return(1000)
    own <- index(islands[[label]])
    others <- \(xs, ys) sum(grid[xs, ys]) - sum(own[, 1] %in% xs & own[, 2] %in% ys & grid[own] == 1L)
    near <- others(max(1, min(xs) - 2):min(nrow(grid), max(xs) + 2), max(1, min(ys) - 2):min(ncol(grid), max(ys) + 2))
    inside <- others(xs, ys)
    own_inside <- sum(grid[xs, ys]) - inside
    gap <- \(bins) min(pmax(abs(bins$bin_x - x) - w / 2, abs(bins$bin_y - y) - h / 2, 0))
    misread <- side != "on" && distance < 7 && gap(dplyr::bind_rows(islands[names(islands) != label])) < gap(islands[[label]])
    10 * inside + 3 * (near - inside) + (side != "on") * 3 * own_inside + 500 * misread
  }
  candidates$base <- purrr::pmap_dbl(candidates[c("x", "y", "w", "h", "label", "side", "distance")], covered) +
    4 * candidates$distance
  # Each label box and, for labels seven or more bins out, its leader line from the box
  # to the island's edge, as rectangles that later labels and leaders must not overlap.
  candidates <- dplyr::mutate(candidates, leader = .data$distance >= 7,
    from_x = dplyr::case_when(.data$side == "right" ~ .data$x - .data$w / 2, .data$side == "left" ~ .data$x + .data$w / 2,
      .default = .data$anchor_x),
    from_y = dplyr::case_when(.data$side == "top" ~ .data$y - .data$h / 2, .data$side == "bottom" ~ .data$y + .data$h / 2,
      .default = .data$anchor_y))
  rectangles <- \(rows) dplyr::bind_rows(
    dplyr::transmute(rows, x0 = .data$x - .data$w / 2, x1 = .data$x + .data$w / 2, y0 = .data$y - .data$h / 2, y1 = .data$y + .data$h / 2),
    dplyr::filter(rows, .data$leader) |>
      dplyr::transmute(x0 = pmin(.data$from_x, .data$anchor_x) - .2, x1 = pmax(.data$from_x, .data$anchor_x) + .2,
        y0 = pmin(.data$from_y, .data$anchor_y) - .2, y1 = pmax(.data$from_y, .data$anchor_y) + .2))
  placed <- candidates[0, ]
  while (nrow(candidates)) {
    taken <- rectangles(placed)
    clash <- vapply(seq_len(nrow(candidates)), \(i) {
      own <- rectangles(candidates[i, ])
      any(outer(own$x0, taken$x1, `<`) & outer(own$x1, taken$x0, `>`) & outer(own$y0, taken$y1, `<`) & outer(own$y1, taken$y0, `>`))
    }, logical(1))
    candidates$score <- candidates$base + 1000 * clash
    best <- dplyr::slice_min(candidates, .data$score, n = 1, with_ties = FALSE, by = "label")
    free <- dplyr::count(dplyr::filter(candidates, .data$score < 10), .data$label)
    best <- dplyr::left_join(best, free, by = "label") |>
      dplyr::arrange(dplyr::coalesce(.data$n, 0L), .data$score) |>
      dplyr::slice(1) |>
      dplyr::select(-"n")
    placed <- dplyr::bind_rows(placed, best)
    candidates <- dplyr::filter(candidates, .data$label != best$label)
  }
  placed
}

# B: the WNN UMAP coloured by cell type in shades of its compartment's colour, with
# unassigned nuclei grey. A cell type from several tissues forms several islands, so
# each label sits beside its main island, as figure_1_place_labels() describes. Label
# sizes assume the panel's printed width in Figure 1, about 48 mm, and labels may use
# the 15 mm of free space to its right.
plot_figure_1_UMAP <- function(data, cell_types) {
  cell_types <- intersect(cell_types, data$WNN_harmony_SNN_cluster_cell_type)
  compartments <- get_cluster_compartments(cell_types, figure_1_compartments)
  hues <- get_compartment_colours(names(figure_1_compartments))
  colours <- unlist(unname(lapply(split(cell_types, compartments, drop = TRUE), \(members) {
    ramp <- grDevices::colorRampPalette(c("black", hues[[as.character(compartments[match(members[[1]], cell_types)])]], "white"))(11)
    stats::setNames(ramp[if (length(members) == 1) 6 else round(seq(4, 8, length.out = length(members)))], members)
  })))
  data <- dplyr::mutate(data, cell_type = factor(.data$WNN_harmony_SNN_cluster_cell_type, cell_types)) |>
    dplyr::arrange(!is.na(.data$cell_type))
  bin <- diff(range(data$WNN_UMAP_1)) / 80
  bin_mm <- 48 / (80 * 1.04)
  bins <- dplyr::mutate(data, bin_x = round(.data$WNN_UMAP_1 / bin), bin_y = round(.data$WNN_UMAP_2 / bin)) |>
    dplyr::summarise(n = dplyr::n(), .by = c("cell_type", "bin_x", "bin_y"))
  islands <- dplyr::filter(bins, !is.na(.data$cell_type), .data$n > sum(.data$n) / 2, .by = c("bin_x", "bin_y")) |>
    split(~ cell_type, drop = TRUE) |>
    lapply(figure_1_main_island)
  names(islands) <- figure_1_cell_type_labels(names(islands))
  register_benchmark_font()
  width_mm <- systemfonts::string_width(names(islands), family = BENCHMARK_FONT, size = 6) * 25.4 / 72
  labels <- figure_1_place_labels(lapply(islands, \(island) dplyr::filter(island, .data$n > 2)),
    occupied = dplyr::filter(dplyr::summarise(bins, n = sum(.data$n), .by = c("bin_x", "bin_y")), .data$n > 2),
    limits = c(range(bins$bin_x), range(bins$bin_y)) + c(-1, 1 + round(15 / bin_mm), -1, 1),
    width = stats::setNames(width_mm / bin_mm + 1, names(islands)), height = 2.6 / bin_mm) |>
    dplyr::mutate(dplyr::across(c("x", "y", "from_x", "from_y", "anchor_x", "anchor_y"), \(value) value * bin))
  ggplot2::ggplot(data, ggplot2::aes(.data$WNN_UMAP_1, .data$WNN_UMAP_2, colour = .data$cell_type)) +
    ggplot2::geom_point(size = .1, shape = 16) +
    ggplot2::geom_segment(data = dplyr::filter(labels, .data$leader), ggplot2::aes(x = .data$from_x, y = .data$from_y,
      xend = .data$anchor_x, yend = .data$anchor_y), inherit.aes = FALSE, linewidth = .25) +
    ggrepel::geom_text_repel(data = labels, ggplot2::aes(.data$x, .data$y, label = .data$label), inherit.aes = FALSE,
      family = BENCHMARK_FONT, size = figure_1_text_size, bg.colour = "white", bg.r = .12, force = 0, force_pull = 0,
      max.iter = 0, min.segment.length = Inf, max.overlaps = Inf, xlim = c(-Inf, Inf), ylim = c(-Inf, Inf)) +
    ggplot2::scale_colour_manual(values = colours, na.value = "grey85", guide = "none") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = .02)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = .02)) +
    figure_1_theme() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(linewidth = .3, colour = "grey30"),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 2)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 2))) +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::labs(title = "Cell types from expression and accessibility", x = "UMAP 1", y = "UMAP 2")
}

# C: the configured WNN marker gene of each cell type whose expression and promoter
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

# Markers along x and cell types down y in the same order, so each cell type's own
# marker, boxed, lies on the diagonal. Each cell is split: expression in the upper-left
# triangle and promoter accessibility in the lower-right one, so a cell where both
# agree reads as one square. Values below the cell-type mean are white. A key sits
# in the corner above the cell-type names.
plot_figure_1_markers <- function(profiles) {
  profiles <- droplevels(profiles)
  corners <- list(`Gene expression` = list(x = c(-.5, -.5, .5), y = c(-.5, .5, .5)),
    `Promoter accessibility` = list(x = c(-.5, .5, .5), y = c(-.5, -.5, .5)))
  triangles <- profiles |>
    dplyr::mutate(x = as.integer(.data$gene), y = as.integer(.data$cell_type), triangle = dplyr::row_number()) |>
    dplyr::reframe(x = .data$x + corners[[.data$modality]]$x, y = .data$y + corners[[.data$modality]]$y, z = .data$z,
      .by = c("triangle", "modality"))
  n <- nlevels(profiles$gene)
  key_y <- c(`Gene expression` = n + 2.9, `Promoter accessibility` = n + 1.5)
  key <- purrr::imap_dfr(corners, \(corner, modality) tibble::tibble(modality, y = key_y[[modality]] + corner$y, x = corner$x))
  ggplot2::ggplot(triangles, ggplot2::aes(.data$x, .data$y)) +
    ggplot2::geom_polygon(ggplot2::aes(group = .data$triangle, fill = .data$z)) +
    ggplot2::annotate("rect", xmin = seq_len(n) - .5, xmax = seq_len(n) + .5, ymin = rev(seq_len(n)) - .5,
      ymax = rev(seq_len(n)) + .5, fill = NA, colour = "grey10", linewidth = .35) +
    ggplot2::geom_polygon(data = key, ggplot2::aes(group = .data$modality), fill = "grey50") +
    ggplot2::annotate("rect", xmin = -.5, xmax = .5, ymin = key_y - .5, ymax = key_y + .5, fill = NA,
      colour = "grey50", linewidth = .25) +
    ggplot2::annotate("text", x = -.8, y = key_y, hjust = 1, family = BENCHMARK_FONT, size = figure_1_text_size,
      label = c("Expression", "Promoter accessibility")) +
    ggplot2::scale_x_continuous(breaks = seq_len(n), labels = levels(profiles$gene), position = "top") +
    ggplot2::scale_y_continuous(breaks = seq_len(n), labels = figure_1_cell_type_labels(levels(profiles$cell_type))) +
    ggplot2::scale_fill_gradient(low = "white", high = "#B40426", limits = c(0, max(triangles$z)), oob = scales::squish,
      name = "Scaled\nlog CPM") +
    ggplot2::coord_cartesian(xlim = c(.5, n + .5), ylim = c(.5, n + .5), expand = FALSE, clip = "off") +
    figure_1_theme() +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 5.5),
      axis.text.x = ggplot2::element_text(face = "italic", angle = 90, hjust = 0, vjust = .5),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(), legend.position = "right") +
    ggplot2::labs(title = "Cell-type markers in the same nuclei", x = NULL, y = NULL)
}

# D: inspired by top_link_aggregate_scatter_plots. Separate compact genomic
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
    ggplot2::theme(axis.title.y = ggplot2::element_text(angle = 0, hjust = 1, vjust = .5)) +
    ggplot2::annotate("rect", xmin = link$start, xmax = link$end,
      ymin = -Inf, ymax = Inf, fill = "#B40426", alpha = 0.25) +
    ggplot2::geom_vline(xintercept = link$TargetGeneTSS, linetype = 2, linewidth = 0.3) +
    ggplot2::scale_x_continuous(limits = bounds, expand = c(0, 0),
      labels = scales::label_number(scale = 1e-6, accuracy = .01))
  genes <- base + genomic_gene_body_layers(genes) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = c(.3, .8))) +
    ggplot2::coord_cartesian(clip = "off") +
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

# E: same residuals as the parent scatter, separate panel; donor colours without
# a large legend. The line is descriptive, not the hierarchical model fit.
plot_figure_1_scatter <- function(data) {
  ggplot2::ggplot(data, ggplot2::aes(.data$peak_accessibility_residual,
    .data$gene_expression_residual)) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$donor_id), size = .9, alpha = .7) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
      colour = "black", linewidth = .4) + figure_1_theme() +
    ggplot2::annotate("text", x = -Inf, y = Inf, hjust = -.08, vjust = 1.3,
      label = figure_1_p_label(data$hierarchical_pvalue[[1]]), parse = TRUE,
      family = BENCHMARK_FONT, size = figure_1_text_size) +
    ggplot2::guides(colour = "none") +
    ggplot2::labs(title = "Donor\u2013state pseudobulks",
      x = "ATAC residual (log1p CPM)", y = "GEX residual\n(log1p CPM)")
}

# F: S-LDSC tau* of each study with genome-wide summary statistics in each assigned cell
# type's share of accessibility, from the module's specificity model.
figure_1_heritability <- function(specificity, cell_types) {
  specificity |>
    dplyr::filter(.data$model == "peaks", .data$annotation == "peaks", !is.na(.data$tau_star)) |>
    dplyr::mutate(cluster = gsub("-", "_", .data$cluster)) |>
    dplyr::filter(.data$cluster %in% cell_types) |>
    dplyr::select("GWAS_ID", "Category", "cluster", "tau_star", "z", "p", "q")
}

# Small tracks show 0 and the largest pretty break within the expanded range.
figure_1_track_breaks <- function(limits) {
  breaks <- scales::breaks_pretty(n = 2)(limits)
  c(0, max(breaks[breaks <= limits[[2]]]))
}

# G: aligned locus facets with independent y scales and shared axis titles
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
    ggplot2::theme(strip.clip = "off", axis.title.y = ggplot2::element_text(angle = 0, hjust = 1, vjust = .5),
      panel.spacing.x = grid::unit(6, "mm"), legend.position = "right")
  contribution <- base +
    ggplot2::geom_hline(yintercept = 0, colour = "grey70") +
    ggplot2::geom_segment(data = variants, ggplot2::aes(x = .data$position,
      xend = .data$position, y = 0, yend = .data$relative_deviation_contribution), colour = "#B40426") +
    ggplot2::geom_point(data = variants, ggplot2::aes(.data$position,
      .data$relative_deviation_contribution, size = .data$posteriorProbability), colour = "#B40426") +
    ggplot2::scale_size_area(max_size = 3, name = "PIP", limits = c(0, 1), breaks = c(.1, .5, 1)) +
    ggplot2::scale_y_continuous(breaks = figure_1_track_breaks,
      labels = scales::label_number(drop0trailing = TRUE), expand = ggplot2::expansion(mult = c(.05, .35))) +
    ggplot2::labs(x = NULL, y = "Genetic\ncontribution") +
    ggplot2::theme(axis.text.x = ggplot2::element_blank())
  # A window without protein-coding genes names the locus's Open Targets locus-to-gene prediction.
  empty <- dplyr::summarise(bounds, position = mean(.data$position), .by = "locus") |>
    dplyr::mutate(label = paste0("No protein-coding genes in this window\nLocus-to-gene prediction: ",
      vapply(records, \(x) x$locus$likely_genes_label[[1]], character(1)))) |>
    dplyr::filter(!.data$locus %in% genes$locus)
  gene_track <- base + genomic_gene_body_layers(genes) +
    ggplot2::geom_text(data = empty, ggplot2::aes(.data$position, -1, label = .data$label),
      size = figure_1_text_size, lineheight = .9, vjust = 0) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = c(.3, .8))) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = NULL, y = "Genes") + ggplot2::guides(y = "none") +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), strip.text = ggplot2::element_blank())
  accessibility <- base +
    ggplot2::geom_area(data = coverage, ggplot2::aes(.data$pos, .data$normalized_insertions), fill = "#446B8C") +
    ggplot2::geom_segment(data = peaks, ggplot2::aes(x = .data$start, xend = .data$end,
      y = 0, yend = 0), linewidth = 1, colour = "#B40426") +
    ggplot2::facet_wrap(ggplot2::vars(locus), nrow = 1, scales = "free", strip.position = "bottom",
      labeller = ggplot2::labeller(locus = ggplot2::as_labeller(paste(chromosomes, "(Mb)") |> stats::setNames(labels)))) +
    ggplot2::scale_y_continuous(breaks = figure_1_track_breaks, labels = scales::label_number(drop0trailing = TRUE)) +
    ggplot2::labs(x = NULL, y = "ATAC") +
    ggplot2::theme(strip.placement = "outside")
  # Keep each axis label beside its track while aligning the faceted panels.
  patchwork::wrap_plots(lapply(list(contribution, gene_track, accessibility),
    function(plot) patchwork::free(plot, type = "label", side = "l")),
    ncol = 1, heights = c(1, 1.2, 1), guides = "collect") +
    patchwork::plot_annotation(title = "Leading GWAS locus", theme = figure_1_theme())
}

# H: TRS_UMAPs inspiration; rebuild WNN UMAP from coordinates and scores, with
# WNN cell-type centroids, as in B, no score clipping or serialized upstream plot mutation.
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
      family = BENCHMARK_FONT, size = figure_1_text_size, linewidth = .2, label.padding = grid::unit(1.2, "pt"), fill = "white") +
    ggplot2::scale_colour_viridis_c(name = "TRS", na.value = "grey85") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = .02)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = .02)) +
    figure_1_theme() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text = ggplot2::element_blank(),
      legend.position = "right") +
    ggplot2::guides(colour = ggplot2::guide_colourbar(barwidth = grid::unit(2, "mm"),
      barheight = grid::unit(25, "mm"))) +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::labs(title = paste(figure_1_trait_labels(trait), "relevance (SCAVENGE)"), x = NULL, y = NULL)
}

# The public data's cell types in compartments, in their default order.
figure_1_compartments <- list(Immune = "^(B|T|NK|Mono|pDC|Macrophage)$", Neural = "neuron|Oligodendrocyte",
  Hepatic = "Hepatocyte", Pancreatic = "Pancreatic|Ductal|Islet", Intestinal = "Colonocyte|Goblet",
  Pulmonary = "Alveolar", Cardiac = "Cardiomyocyte", Vascular = "Endothelial|Pericyte", Stromal = "Fibroblast")

# The configured trait categories in the figure's categories, in row order: the order of
# the compartments that each concerns, with the immune and the brain traits bundled.
figure_1_trait_categories <- c(blood_cell_counts = "Blood cell counts", immune_disease = "Immune disease",
  immune_gut = "Immune disease", lung_immune = "Immune disease", immune_brain = "Immune disease",
  neurodegeneration = "Brain traits", neuronal_disease = "Brain traits", neuronal_metabolic = "Brain traits",
  neuronal_complex = "Brain traits", hepatic_lipids = "Hepatic lipids", liver = "Liver",
  pancreatic_endocrine = "Pancreatic endocrine", colon_epithelial = "Colon epithelial", lung_function = "Lung function",
  cardiac = "Cardiac", vascular = "Vascular", fibroblast = "Fibroblast")

# Studies as Supplementary Table S3 names them: "Trait, Author Year".
figure_1_study_labels <- function(x) {
  paste0(figure_1_trait_labels(x), ", ", sub("^.*_([A-Za-z]+)([0-9]{4}).*$", "\\1 \\2", x))
}

# Figure 1F and Supplementary Figures S5 and S6: one value of each study (rows) in each
# assigned cell type (columns), in a diverging scale around `midpoint` clipped at `limits`.
# scores holds GWAS_ID, Category, cluster, value, p and q. Studies are grouped in the
# categories of figure_1_trait_categories, in its order, and by name within each. Each
# study's anchor, its largest value at q < 0.05, is boxed. Compartments, and cell types
# within them, follow the first significant tile they hold, taking anchors first and
# each in row order, so that the anchors step down from the top left; the rest keep
# their default order at the end. A header names each compartment of several cell
# types in bold above a line over its columns, above the cell-type names; a one-column
# compartment's cell type names it already. A side column names each category in bold
# beside a line along its rows, left of the study names.
# Stars: * P < 0.05, ** q < 0.05.
plot_trait_cell_type_grid <- function(scores, cell_types, value_label, limits, midpoint = 0,
                                      row_labels = figure_1_study_labels, title = NULL) {
  scores <- dplyr::filter(scores, .data$cluster %in% cell_types) |>
    dplyr::mutate(significant = .data$q < .05 & !is.na(.data$q),
      Category = factor(figure_1_trait_categories[.data$Category], unique(figure_1_trait_categories))) |>
    dplyr::arrange(!.data$significant, dplyr::desc(.data$value)) |>
    dplyr::mutate(anchor = .data$significant & dplyr::row_number() == 1, .by = "GWAS_ID")
  stopifnot(!anyNA(scores$Category))
  rows <- dplyr::distinct(scores, .data$GWAS_ID, .data$Category) |>
    dplyr::arrange(.data$Category, .data$GWAS_ID) |>
    dplyr::mutate(y = factor(.data$GWAS_ID, rev(.data$GWAS_ID)))
  steps <- dplyr::filter(scores, .data$significant) |> dplyr::arrange(!.data$anchor, match(.data$GWAS_ID, rows$GWAS_ID))
  columns <- tibble::tibble(cluster = intersect(cell_types, scores$cluster)) |>
    dplyr::mutate(compartment = get_cluster_compartments(.data$cluster, figure_1_compartments),
      step = match(.data$cluster, steps$cluster)) |>
    dplyr::mutate(first_step = min(.data$step, Inf, na.rm = TRUE), .by = "compartment") |>
    dplyr::arrange(.data$first_step, .data$compartment, .data$step) |>
    dplyr::mutate(compartment = factor(.data$compartment, unique(.data$compartment)), x = factor(.data$cluster, .data$cluster))
  tiles <- dplyr::mutate(scores, value = pmin(pmax(.data$value, limits[[1]]), limits[[2]]),
    x = factor(.data$cluster, levels(columns$x)), y = factor(.data$GWAS_ID, levels(rows$y)),
    compartment = columns$compartment[match(.data$cluster, columns$cluster)])
  named <- dplyr::count(columns, .data$compartment) |> dplyr::filter(.data$n > 1)
  # The three plots share their facets, scales, tile extents and spacing, so that their panels align.
  x_scale <- ggplot2::scale_x_discrete(position = "top", labels = figure_1_cell_type_labels, expand = ggplot2::expansion(add = .1))
  y_scale <- ggplot2::scale_y_discrete(labels = row_labels, expand = ggplot2::expansion())
  line <- \(..., data) ggplot2::geom_segment(ggplot2::aes(...), data = data, inherit.aes = FALSE, linewidth = .35, colour = "grey35")
  layout_theme <- figure_1_theme() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.title = ggplot2::element_blank(),
      axis.ticks.length = grid::unit(0, "pt"), strip.clip = "off", plot.margin = ggplot2::margin(0, 0, 0, 0),
      panel.spacing.x = grid::unit(2.4, "mm"), panel.spacing.y = grid::unit(.5, "mm"))
  header <- ggplot2::ggplot(columns, ggplot2::aes(.data$x, 0)) +
    ggplot2::geom_tile(width = .86, height = 0, fill = NA) +
    line(x = .62, xend = .data$n + .38, y = 0, yend = 0, data = named) +
    ggplot2::facet_grid(cols = ggplot2::vars(.data$compartment), scales = "free_x", space = "free_x",
      labeller = ggplot2::labeller(compartment = \(x) dplyr::if_else(x %in% named$compartment, x, ""))) +
    x_scale + ggplot2::scale_y_continuous(limits = c(0, 1), expand = ggplot2::expansion()) +
    ggplot2::coord_cartesian(clip = "off") + layout_theme +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
      strip.text.x.top = ggplot2::element_text(face = "bold", size = 6, margin = ggplot2::margin(b = 1.5)))
  side <- ggplot2::ggplot(rows, ggplot2::aes(0, .data$y)) +
    ggplot2::geom_tile(width = 0, height = .86, fill = NA) +
    line(x = 1, xend = 1, y = .62, yend = .data$n + .38, data = dplyr::count(rows, .data$Category)) +
    ggplot2::facet_grid(rows = ggplot2::vars(.data$Category), scales = "free_y", space = "free_y", switch = "y") +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = ggplot2::expansion()) + y_scale +
    ggplot2::coord_cartesian(clip = "off") + layout_theme +
    ggplot2::theme(axis.text = ggplot2::element_blank(), strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1, face = "bold", size = 5.5, margin = ggplot2::margin(r = 1.5)))
  grid <- ggplot2::ggplot(tiles, ggplot2::aes(.data$x, .data$y)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$value), width = .86, height = .86) +
    ggplot2::geom_tile(data = dplyr::filter(tiles, .data$anchor), fill = NA, colour = "black", linewidth = .45,
      width = .86, height = .86) +
    ggplot2::geom_text(ggplot2::aes(label = get_significance_stars(.data$p, .data$q)), size = 1.9, vjust = .78) +
    ggplot2::facet_grid(rows = ggplot2::vars(.data$Category), cols = ggplot2::vars(.data$compartment), scales = "free",
      space = "free") +
    x_scale + y_scale +
    ggplot2::scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = midpoint, limits = limits,
      name = value_label, guide = ggplot2::guide_colourbar(direction = "horizontal", position = "inside",
        theme = ggplot2::theme(legend.title.position = "top", legend.key.width = grid::unit(18, "mm"),
          legend.key.height = grid::unit(2.2, "mm")))) +
    ggplot2::coord_cartesian(clip = "off") + layout_theme +
    ggplot2::theme(strip.text = ggplot2::element_blank(), axis.text = ggplot2::element_text(size = 5.5),
      axis.text.x.top = ggplot2::element_text(angle = 90, hjust = 0, vjust = .5, margin = ggplot2::margin(t = 1.5, b = 2)),
      plot.margin = ggplot2::margin(0, 8, 0, 0),
      legend.position.inside = c(0, 1), legend.justification.inside = c(1, 0), legend.margin = ggplot2::margin(r = 6, b = 4))
  patchwork::wrap_plots(patchwork::plot_spacer(), header, side, grid, ncol = 2,
    widths = grid::unit(c(1.2, 1), c("mm", "null")), heights = grid::unit(c(1, 1), c("mm", "null"))) +
    patchwork::plot_annotation(title = title)
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
compose_figure_1 <- function(...) {
  plots <- lapply(list(...), function(plot) {
    if (inherits(plot, "patchwork")) {
      title <- plot$patches$annotation$title
      plot <- plot + patchwork::plot_annotation(title = NULL)
    } else {
      title <- plot$labels$title
      plot <- plot + ggplot2::labs(title = NULL)
    }
    patchwork::wrap_elements(panel = plot) + ggplot2::labs(title = title) +
      ggplot2::theme(plot.title = ggplot2::element_text(family = BENCHMARK_FONT, size = 10,
          hjust = 0, margin = ggplot2::margin(l = 14, b = 4)),
        plot.title.position = "plot")
  })
  (patchwork::wrap_plots(plots, design = "AA\nBC\nDE\nFF\nGH",
    heights = c(.85, 2.45, 1.1, 3.05, 1.85)) +
    patchwork::plot_annotation(tag_levels = "A")) &
    ggplot2::theme(plot.tag = ggplot2::element_text(family = BENCHMARK_FONT, face = "bold",
        size = 12, hjust = 0, vjust = 1), plot.tag.position = c(0, 1))
}
