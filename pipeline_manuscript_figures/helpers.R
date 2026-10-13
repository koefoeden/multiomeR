benchmark_aggregation_tibble <- function(aggregations) {
  GEM_well_tibble <- build_GEM_well_tibble()
  aggregation_tibble_all_from_yaml <- read_aggregation_config_tibble(config_file = configuration_path("cfg_aggregations.yaml"))
  aggregation_tibble <- build_aggregation_tibble(
    aggregation_tibble_all_from_yaml = aggregation_tibble_all_from_yaml,
    GEM_well_tibble = GEM_well_tibble
  )

  missing_aggregations <- setdiff(aggregations, aggregation_tibble$aggregation)
  if (length(missing_aggregations)) {
    stop(
      "Aggregation(s) not found in cfg_aggregations.yaml: ",
      paste(missing_aggregations, collapse = ", "),
      call. = FALSE
    )
  }

  aggregation_tibble |>
    dplyr::mutate(
      benchmark_cellranger_count_dirs = purrr::map(
        .data$aggregation_GEM_well_IDs,
        \(GEM_well_IDs) GEM_well_tibble$GEM_well_cellranger_arc_count_dir[match(GEM_well_IDs, GEM_well_tibble$GEM_well_ID)]
      )
    ) |>
    dplyr::slice(match(aggregations, .data$aggregation))
}


benchmark_cellranger_input_nuclei <- function(cellranger_count_dir) {
  cellranger_h5_file <- file.path(cellranger_count_dir, "outs", "filtered_feature_bc_matrix.h5")
  cellranger_h5_file_con <- hdf5r::H5File$new(cellranger_h5_file, mode = "r")
  on.exit(cellranger_h5_file_con$close_all())
  length(cellranger_h5_file_con[["matrix/barcodes"]][])
}


benchmark_aggregation_cellranger_input_nuclei <- function(aggregations) {
  aggregation_tibble <- benchmark_aggregation_tibble(aggregations)
  cellranger_input_nuclei <- purrr::map_int(
    aggregation_tibble$benchmark_cellranger_count_dirs,
    \(cellranger_count_dirs) {
      purrr::map_int(cellranger_count_dirs, benchmark_cellranger_input_nuclei) |>
        sum()
    }
  )

  stats::setNames(cellranger_input_nuclei, aggregation_tibble$aggregation)[aggregations]
}

#' Estimate wall time to a target from its recorded runtime critical path
#'
#' Pattern targets count their slowest dynamic branch, as if branches ran in
#' parallel. Targets matching `exclude_target_regex` keep their graph position
#' with zero runtime.
#' @param target_name Fully resolved static target name.
#' @param network Output from `targets::tar_network(targets_only = TRUE)`.
#' @param meta Output from `targets::tar_meta()`.
#' @param exclude_target_regex Regular expressions of excluded target names.
#' @return List with the critical-path hours and the ancestor targets.
#' @keywords internal
estimate_target_walltime <- function(target_name, network, meta, exclude_target_regex = character()) {
  graph <- igraph::graph_from_data_frame(network$edges, vertices = network$vertices["name"])
  ancestor_graph <- igraph::induced_subgraph(graph, igraph::subcomponent(graph, target_name, mode = "in"))
  ancestors <- network$vertices[match(igraph::as_ids(igraph::V(ancestor_graph)), network$vertices$name), c("name", "type")]
  seconds <- as.numeric(network$vertices$seconds[match(ancestors$name, network$vertices$name)])
  patterns <- ancestors$type == "pattern"
  seconds[patterns] <- purrr::map_dbl(ancestors$name[patterns], function(name) {
    child_seconds <- meta$seconds[match(unlist(meta$children[meta$name == name]), meta$name)]
    if (all(is.na(child_seconds))) NA_real_ else max(child_seconds, na.rm = TRUE)
  }) |>
    dplyr::coalesce(seconds[patterns])
  ancestors$excluded_from_benchmark <- Reduce(
    `|`, lapply(exclude_target_regex, grepl, x = ancestors$name, perl = TRUE), rep(FALSE, nrow(ancestors))
  )
  seconds[ancestors$excluded_from_benchmark] <- 0
  if (anyNA(seconds)) {
    stop("Missing recorded runtime for ancestor target(s): ", paste(ancestors$name[is.na(seconds)], collapse = ", "), call. = FALSE)
  }
  names(seconds) <- ancestors$name
  path_seconds <- seconds * 0
  for (node in igraph::as_ids(igraph::topo_sort(ancestor_graph))) {
    parents <- igraph::as_ids(igraph::neighbors(ancestor_graph, node, mode = "in"))
    path_seconds[[node]] <- seconds[[node]] + max(0, path_seconds[parents])
  }
  list(critical_path_hours = path_seconds[[target_name]] / 3600, ancestors = ancestors)
}

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

# Supplementary Figure S3: the module's unscaled deviation heatmap of every study
# and assigned cell type, with Figure 1D's colour scale, cell-type labels and stars;
# rows also name the study. Deviations beyond Figure 1D's range take its end colours.
plot_figure_S3 <- function(scores, GWAS_metadata_tracks_plot, traits) {
  limit <- max(abs(scores$deviation[scores$GWAS_ID %in% traits$GWAS_ID]))
  # Of the module's metadata tracks, only trait category and credible-set loci are
  # kept; Supplementary Table S3 gives each study's sample size, method and ancestry.
  bars <- which(vapply(seq_along(GWAS_metadata_tracks_plot), \(i) "track" %in% names(GWAS_metadata_tracks_plot[[i]]$data), logical(1)))
  stopifnot(length(bars) == 1L)
  loci <- GWAS_metadata_tracks_plot[[bars]]
  loci$data <- dplyr::filter(loci$data, .data$track == "Loci")
  category <- GWAS_metadata_tracks_plot[[1]] +
    ggplot2::scale_y_discrete(drop = FALSE, expand = c(0, 0),
      labels = \(x) paste0(figure_1_trait_labels(x), " (", sub("^.*_([A-Za-z]+)([0-9]{4}).*$", "\\1 \\2", x), ")"))
  ordered <- scores |>
    dplyr::mutate(cluster = figure_1_cell_type_labels(.data$cluster),
      deviation = pmax(pmin(.data$deviation, limit), -limit), label = dplyr::if_else(.data$z_q < .05, "*", "")) |>
    get_ordered_GWAS_score_plot_data(NULL, score_col = "deviation")
  heatmap <- plot_GWAS_feature_heatmap(ordered$scores, ordered$clusters, "cluster", "deviation", "Deviation",
    fill_limits = c(-limit, limit), support_label_col = "label")
  patchwork::wrap_plots(list(category, loci, heatmap, plot_GWAS_feature_support_tracks(ordered$clusters)),
    design = "##D\nABC", widths = c(.35, 1.6, 9), heights = grid::unit(c(24, 1), c("mm", "null")),
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

# Supplementary Figure S2: the workflows share one graph and store. Excluded
# targets keep their graph position with zero runtime and are left out of every
# resource sum.
SEURAT_SIGNAC_COMPARISON_SCRIPT <- "pipeline_manuscript_figures/targets.R"
SEURAT_SIGNAC_COMPARISON_ALLOCATIONS <- "pipeline_manuscript_figures/results/data/seurat_signac_comparison_allocations.tsv"
SEURAT_SIGNAC_COMPARISON_WORKFLOWS <- local({
  multiomeR_exclude <- c(
    # Cluster-fragment preparation is shared by both workflows.
    "^blacklist_GRanges[.]ATAC(\\.|$)", "^BCs_per_peak_cluster_list[.]ATAC(\\.|$)",
    "^peak_calling_cluster_(names|discovery_tibble)[.]ATAC(\\.|$)",
    "^fragments_per_(cluster|peak_calling_cluster_discovery)[.]fragments[.]ATAC(\\.|$)",
    # The conventional chain performs no cell-type annotation, doublet detection,
    # QC metrics, gene or peak annotation, chromatin-state projection or motif
    # analysis, and no ATAC-only or uncorrected RNA embedding.
    "^cluster_UCell_", "^(UCell_GEX_marker_genes_list|GEX_marker_genes_vec)(\\.|$)",
    "^metadata_w_cell_types(_unfiltered|_annotation)?_tibble[.]",
    "^(amulet|scDblFinder)_", "^(ATAC_qc_metrics|per_barcode_metrics|blacklist_counts)_tibble[.]",
    "^metadata_w_QC", "Ensembl_(gene_)?annotations?_GRanges_list(\\.|$)",
    "^(consensus_peak_annotated|signac_annotation)_GRanges(\\.|$)", "^chromHMMs_",
    "^(chromVAR|motif_family|peak_TF_motif|JASPAR)",
    "^UMAP_embeddings_tibble[.](ATAC|GEX_non_harmony)[.]", "^LSI_clusters[.]")
  list(
    # multiomeR's own result is its WNN cell metadata; the Seurat/Signac object
    # is an optional export built from it.
    multiomeR = list(endpoint = "metadata_w_clusters_tibble.WNN.", exclude = multiomeR_exclude),
    "multiomeR + Seurat exports" = list(endpoint = "multimodal_Seurat_object.8_multimodal_QC.",
      exclude = multiomeR_exclude),
    # Everything outside the conventional chain is imported from multiomeR. Its
    # peak calling uses multiomeR's clusters, so in the graph it can overlap the
    # RNA steps; a script runs the steps one after another.
    "Seurat/Signac" = list(endpoint = "comparison_seurat_signac_object.",
      exclude = "^(?!comparison_seurat_signac_)", sequential = TRUE))
})

# Measured use comes from the per-job history of slurm_monitor, which samples each
# crew worker's cgroup every 10-15 s together with its running target. RAM is the
# memory the kernel cannot reclaim, so page cache from reading and writing files is
# left out; shared pages count once and child processes are included.
read_slurm_monitor_samples <- function(since,
  history_dir = Sys.getenv("SLURM_MONITOR_HISTORY_DIR", "~/slurm_monitor_history")) {
  # Only schema 2 files (<job>.v2.tsv) hold unreclaimable RAM in ram_cur_gib; their
  # columns are found by header name. Rows from sstat hold RSS and are left out.
  awk_program <- paste(
    "FNR == 1 {delete col; for (i = 1; i <= NF; i++) col[$i] = i; next}",
    "$col[\"name\"] ~ /^crew-worker-/ && $col[\"ram_cur_gib\"] != \"NA\" && $col[\"scope\"] != \"sstat\"",
    "{print $col[\"timestamp\"] \"\\t\" $col[\"target\"] \"\\t\" $col[\"ram_cur_gib\"] \"\\t\" $col[\"cpu_cur\"]}")
  cmd <- sprintf("find %s -name '*.v2.tsv' -newermt '%s' -print0 | xargs -0 awk -F'\\t' %s",
    shQuote(path.expand(history_dir)), format(since, "%Y-%m-%d %H:%M:%S"), shQuote(awk_program))
  samples <- data.table::fread(cmd = cmd, sep = "\t", header = FALSE, na.strings = "NA",
    col.names = c("timestamp", "target", "RAM_GiB", "cpu_cores"), colClasses = list(character = 1))
  tibble::tibble(target = samples$target,
    time = as.POSIXct(samples$timestamp, tz = "Europe/Copenhagen"),
    RAM_GB = samples$RAM_GiB * 2^30 / 1e9, cpu_cores = samples$cpu_cores)
}

# Mean measured RAM and CPU of each execution within its recorded build window.
# Executions too short to be sampled get the median use of the sampled ones.
add_measured_resources <- function(executions, samples) {
  measured <- executions |>
    dplyr::select("name", "seconds", "time") |>
    dplyr::inner_join(samples, by = c(name = "target"), suffix = c("", "_sample"), relationship = "many-to-many") |>
    dplyr::filter(.data$time_sample >= .data$time - .data$seconds - 20, .data$time_sample <= .data$time + 20) |>
    dplyr::summarise(used_RAM_GB = mean(.data$RAM_GB), peak_RAM_GB = max(.data$RAM_GB),
      used_cores = mean(.data$cpu_cores, na.rm = TRUE), .by = "name")
  executions |>
    dplyr::left_join(measured, by = "name") |>
    dplyr::mutate(sampled = !is.na(.data$used_RAM_GB),
      used_RAM_GB = dplyr::coalesce(.data$used_RAM_GB, stats::median(.data$used_RAM_GB, na.rm = TRUE)),
      used_cores = dplyr::coalesce(.data$used_cores, stats::median(.data$used_cores, na.rm = TRUE)))
}

# Recorded runtimes and bytes from tar_meta(); allocated cores and RAM of each
# static target from `allocations`; measured CPU and RAM from slurm_monitor.
# Dynamic branches run concurrently in the critical path; a sequential workflow's
# critical path is the sum of its runtimes. Disk use counts target
# objects and file targets inside the store; the NA path of a disabled input
# file is outside it.
benchmark_resource_tibble <- function(workflows, aggregations, allocations, script = "_targets.R",
                                      store = targets::tar_config_get("store")) {
  network <- targets::tar_network(script = script, store = store,
    targets_only = TRUE, outdated = FALSE, reporter = "silent")
  meta <- targets::tar_meta(store = store)
  store_prefix <- paste0(fs::path_abs(store), "/")
  aggregation_tibble <- benchmark_aggregation_tibble(aggregations)
  nuclei <- benchmark_aggregation_cellranger_input_nuclei(aggregations)

  runs <- tidyr::expand_grid(aggregation = aggregations, workflow = names(workflows)) |>
    purrr::pmap(function(aggregation, workflow) {
      spec <- workflows[[workflow]]
      estimate <- estimate_target_walltime(paste0(spec$endpoint, aggregation),
        network = network, meta = meta, exclude_target_regex = spec$exclude)
      included <- dplyr::filter(estimate$ancestors, !.data$excluded_from_benchmark)
      patterns <- included$name[included$type == "pattern"]
      executions <- meta[match(c(setdiff(included$name, patterns),
        unlist(meta$children[match(patterns, meta$name)])), meta$name), ] |>
        dplyr::mutate(static_target = dplyr::if_else(.data$type == "branch", .data$parent, .data$name),
          in_store = .data$format != "file" | purrr::map_lgl(.data$path,
            \(path) isTRUE(all(startsWith(fs::path_abs(path), store_prefix))))) |>
        dplyr::left_join(allocations, by = "static_target")
      stopifnot(!anyNA(executions[c("seconds", "bytes", "cores", "RAM_GB")]),
        executions$repository == "local")
      list(aggregation = aggregation, workflow = workflow, estimate = estimate, executions = executions,
        sequential = isTRUE(spec$sequential))
    })
  first_start <- min(purrr::map_dbl(runs, \(run) min(as.numeric(run$executions$time) - run$executions$seconds)))
  samples <- read_slurm_monitor_samples(since = as.POSIXct(first_start - 3600, origin = "1970-01-01"))

  purrr::map_dfr(runs, \(run) {
    executions <- add_measured_resources(run$executions, samples)
    tibble::tibble(aggregation = run$aggregation,
      GEM_well_count = length(aggregation_tibble$aggregation_GEM_well_IDs[[match(run$aggregation, aggregation_tibble$aggregation)]]),
      cellranger_input_nuclei = nuclei[[run$aggregation]], workflow = run$workflow,
      critical_path_hours = if (run$sequential) sum(executions$seconds) / 3600 else run$estimate$critical_path_hours,
      allocated_CPU_hours = sum(executions$seconds * executions$cores) / 3600,
      allocated_RAM_GB_hours = sum(executions$seconds * executions$RAM_GB) / 3600,
      used_CPU_hours = sum(executions$seconds * executions$used_cores) / 3600,
      used_RAM_GB_hours = sum(executions$seconds * executions$used_RAM_GB) / 3600,
      peak_RAM_GB = max(executions$peak_RAM_GB, na.rm = TRUE),
      sampled_runtime_fraction = sum(executions$seconds[executions$sampled]) / sum(executions$seconds),
      disk_space_GB = sum(executions$bytes[executions$in_store]) / 1e9)
  })
}

# Every workflow also reads the same Cell Ranger outputs from outside the store,
# which count towards its disk space. Reserved resources use the worker sizes of
# the measured runs, which multiomeR 1.2.0 requested, not the current requests.
seurat_signac_comparison_resource_tibble <- function(aggregations, store = targets::tar_config_get("store")) {
  input_GB <- purrr::map_dbl(aggregations, \(aggregation)
    sum(file.size(comparison_input_files(comparison_build_config(aggregation)))) / 1e9)
  benchmark_resource_tibble(SEURAT_SIGNAC_COMPARISON_WORKFLOWS, aggregations,
    readr::read_tsv(SEURAT_SIGNAC_COMPARISON_ALLOCATIONS, show_col_types = FALSE),
    script = SEURAT_SIGNAC_COMPARISON_SCRIPT, store = store) |>
    dplyr::mutate(input_disk_space_GB = input_GB[match(.data$aggregation, aggregations)],
      disk_space_GB = .data$disk_space_GB + .data$input_disk_space_GB)
}

# One panel per measure for multiomeR with its Seurat exports and the Seurat/Signac
# baseline, with workflow names at the line ends. A dashed line and fold label mark each
# size whose difference spans at least 20% of the panel's height, so that the label fits in the gap.
plot_seurat_signac_comparison_resources <- function(data) {
  measures <- c(critical_path_hours = "Critical path (h)",
    used_CPU_hours = "CPU time (CPU-hours)",
    used_RAM_GB_hours = "Memory (GB-hours)",
    disk_space_GB = "Disk space (GB)")
  colours <- c(multiomeR = "#2a78d6", "Seurat/Signac" = "#eb6834")
  data <- data |>
    dplyr::filter(.data$workflow != "multiomeR") |>
    dplyr::mutate(workflow = dplyr::recode(.data$workflow, "multiomeR + Seurat exports" = "multiomeR"))
  largest_nuclei <- max(data$cellranger_input_nuclei)
  panels <- lapply(names(measures), function(measure) {
    values <- data |>
      dplyr::transmute(.data$cellranger_input_nuclei, .data$workflow, value = .data[[measure]])
    gaps <- values |>
      tidyr::pivot_wider(names_from = "workflow", values_from = "value") |>
      dplyr::filter(abs(.data$multiomeR - .data$`Seurat/Signac`) >= .2 * max(values$value)) |>
      dplyr::mutate(ratio = .data$`Seurat/Signac` / .data$multiomeR,
        label = sprintf("%.1f\u00d7\n%s", pmax(.data$ratio, 1 / .data$ratio),
          dplyr::if_else(.data$ratio > 1, "lower", "higher")))
    # Names closer than 8% of the panel's height move apart around their midpoint.
    space <- .08 * max(values$value)
    names <- dplyr::filter(values, .data$cellranger_input_nuclei == largest_nuclei) |>
      dplyr::arrange(.data$value) |>
      dplyr::mutate(x = 1.04 * largest_nuclei,
        y = if (diff(.data$value) < space) mean(.data$value) + c(-.5, .5) * space else .data$value)
    ggplot2::ggplot(values, ggplot2::aes(.data$cellranger_input_nuclei, .data$value, colour = .data$workflow)) +
      ggplot2::geom_segment(data = gaps, ggplot2::aes(x = .data$cellranger_input_nuclei,
        xend = .data$cellranger_input_nuclei, y = .data$multiomeR, yend = .data$`Seurat/Signac`),
        inherit.aes = FALSE, linetype = "dashed", linewidth = .3, colour = "grey45") +
      ggplot2::geom_line(linewidth = .5) + ggplot2::geom_point(size = 1.5) +
      ggplot2::geom_label(data = gaps, ggplot2::aes(.data$cellranger_input_nuclei,
        (.data$multiomeR + .data$`Seurat/Signac`) / 2, label = .data$label), inherit.aes = FALSE,
        family = "NimbusSan", size = 2, lineheight = .85, colour = "grey25",
        label.size = 0, label.padding = grid::unit(.6, "pt")) +
      ggplot2::geom_text(data = names, ggplot2::aes(.data$x, .data$y, label = .data$workflow),
        hjust = 0, family = "NimbusSan", size = 2.4) +
      ggplot2::scale_x_continuous(labels = scales::label_comma(), expand = ggplot2::expansion(mult = c(.02, .04))) +
      ggplot2::scale_y_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0, .12))) +
      ggplot2::scale_colour_manual(values = colours, guide = "none") +
      ggplot2::coord_cartesian(clip = "off") + figure_1_theme() +
      ggplot2::theme(plot.margin = ggplot2::margin(5.5, 48, 5.5, 5.5)) +
      ggplot2::labs(x = "Input nuclei", y = measures[[measure]])
  })
  patchwork::wrap_plots(panels, ncol = 2) + patchwork::plot_annotation(tag_levels = "A") &
    ggplot2::theme(plot.tag = ggplot2::element_text(family = "NimbusSan", face = "bold", size = 11))
}
