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

select_figure_1_combinations <- function(scores, n = 2L) {
  selected <- scores |>
    dplyr::filter(is.finite(.data$z), .data$z >= stats::qnorm(.95)) |>
    dplyr::arrange(dplyr::desc(.data$z), .data$GWAS_ID, .data$cluster) |>
    dplyr::slice_head(n = n)
  if (nrow(selected) < n) stop("Too few cell-type/GWAS combinations with raw Z >= qnorm(0.95).")
  selected
}

select_figure_1_traits <- function(scores, n = 6L) {
  selected <- scores |>
    dplyr::filter(is.finite(.data$z), .data$z >= stats::qnorm(.95)) |>
    dplyr::arrange(dplyr::desc(.data$z), .data$GWAS_ID, .data$cluster) |>
    dplyr::distinct(.data$GWAS_ID, .keep_all = TRUE) |>
    dplyr::slice_head(n = n)
  if (!nrow(selected)) stop("No GWAS traits meet raw Z >= qnorm(0.95).")
  selected
}

figure_1_trait_labels <- function(x) {
  gsub("([a-z])([A-Z])", "\\1 \\2", sub("_[^_]+[0-9]{4}$", "", x))
}

# Deviation clusters use hyphens and metadata cell types underscores.
figure_1_cell_type_labels <- function(x) gsub("[-_]", " ", x)

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

# A: inspired by top_link_aggregate_scatter_plots. Separate compact genomic
# panel, focal cell type only, selected peak/TSS and link arc; omit the large
# cross-cell-type evidence stack. Zoom to the peak and TSS with 25 kb flanks; retain cached coverage bins.
plot_figure_1_link_tracks <- function(record, coverage) {
  link <- record$link
  bounds <- c(max(1, min(link$start, link$TargetGeneTSS) - 25000),
    max(link$end, link$TargetGeneTSS) + 25000)
  genes <- record$genes |>
    dplyr::filter(.data$end >= bounds[[1]], .data$start <= bounds[[2]]) |>
    dplyr::mutate(start = pmax(.data$start, bounds[[1]]),
      end = pmin(.data$end, bounds[[2]]),
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
    ggplot2::labs(x = NULL, y = "ATAC") +
    ggplot2::theme(axis.text.x = ggplot2::element_blank())
  arc <- tibble::tibble(position = seq(link$peak_center, link$TargetGeneTSS, length.out = 101),
    height = sin(seq(0, pi, length.out = 101)))
  connection <- base + ggplot2::geom_line(data = arc,
    ggplot2::aes(.data$position, .data$height), colour = "#B40426", linewidth = .5) +
    ggplot2::labs(x = paste0(link$chr, " (Mb)"), y = "Link") + ggplot2::guides(y = "none")
  patchwork::wrap_plots(genes, accessibility, connection, ncol = 1,
    heights = c(1.4, 1, .5)) + patchwork::plot_annotation(
      title = paste(link$TargetGene, "in", link$cell_group),
      theme = figure_1_theme())
}

# B: same residuals as the parent scatter, separate panel; donor colours without
# a large legend. The line is descriptive, not the hierarchical model fit.
plot_figure_1_scatter <- function(data) {
  ggplot2::ggplot(data, ggplot2::aes(.data$peak_accessibility_residual,
    .data$gene_expression_residual)) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$donor_id), size = .9, alpha = .7) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
      colour = "black", linewidth = .4) + figure_1_theme() +
    ggplot2::annotate("text", x = -Inf, y = Inf, hjust = -.08, vjust = 1.3,
      label = sprintf("Hierarchical P = %.2g", data$hierarchical_pvalue[[1]]),
      family = "NimbusSan", size = 2.5) +
    ggplot2::guides(colour = "none") +
    ggplot2::labs(title = "Donor-state aggregates",
      x = "ATAC log1p CPM residual", y = "GEX log1p CPM residual")
}

# C: raw_deviation_unscaled inspiration; at most six traits passing the existing
# raw Z >= qnorm(0.95), no metadata tracks or clustering. All cell
# types retained, common raw-deviation scale. Stars use unadjusted upper-tail normal P <= 0.05 and <= 0.01.
plot_figure_1_heatmap <- function(data, traits) {
  data <- data |> dplyr::filter(.data$GWAS_ID %in% traits$GWAS_ID) |>
    dplyr::mutate(GWAS_ID = factor(.data$GWAS_ID, levels = rev(traits$GWAS_ID)),
      label = chromVAR_Z_support_labels(.data$z))
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

# D: aligned locus facets with independent y scales and shared axis titles
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
    ggplot2::theme(strip.text = ggplot2::element_text(size = 8),
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
    ggplot2::scale_y_continuous(breaks = scales::breaks_pretty(n = 2),
      labels = function(x) formatC(x, format = "f", digits = 1, width = 5),
      expand = ggplot2::expansion(mult = c(.05, .35))) +
    ggplot2::labs(x = NULL, y = "Genetic\ncontribution") +
    ggplot2::theme(axis.text.x = ggplot2::element_blank())
  gene_track <- base + genomic_gene_body_layers(genes) +
    ggplot2::geom_text(data = bounds |> dplyr::filter(!.data$locus %in% genes$locus) |>
      dplyr::summarise(position = mean(.data$position), .by = "locus"),
      ggplot2::aes(.data$position, -1), label = "No annotated genes in this window", size = 2.5) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = c(.3, .8))) +
    ggplot2::labs(x = NULL, y = "Genes") + ggplot2::guides(y = "none") +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), strip.text = ggplot2::element_blank())
  accessibility <- base +
    ggplot2::geom_area(data = coverage, ggplot2::aes(.data$pos, .data$normalized_insertions), fill = "#446B8C") +
    ggplot2::geom_segment(data = peaks, ggplot2::aes(x = .data$start, xend = .data$end,
      y = 0, yend = 0), linewidth = 1, colour = "#B40426") +
    ggplot2::facet_wrap(ggplot2::vars(locus), nrow = 1, scales = "free", strip.position = "bottom",
      labeller = ggplot2::labeller(locus = ggplot2::as_labeller(chromosomes))) +
    ggplot2::scale_y_continuous(breaks = scales::breaks_pretty(n = 2),
      labels = function(x) formatC(x, format = "f", digits = 1, width = 5)) +
    ggplot2::labs(x = "Genomic position (Mb)", y = "ATAC") +
    ggplot2::theme(strip.placement = "outside")
  # Keep each axis label beside its track while aligning the faceted panels.
  patchwork::wrap_plots(lapply(list(contribution, gene_track, accessibility),
    function(plot) patchwork::free(plot, type = "label", side = "l")),
    ncol = 1, heights = c(.9, 2.1, .7), guides = "collect") +
    patchwork::plot_annotation(title = "Leading GWAS loci", theme = figure_1_theme())
}

# E: TRS_UMAPs inspiration; rebuild WNN UMAP from coordinates and scores, with
# cell-type centroids, no score clipping or serialized upstream plot mutation.
plot_figure_1_TRS <- function(data, trait, score_limits) {
  labels <- data |> dplyr::summarise(x = stats::median(.data$WNN_UMAP_1),
    y = stats::median(.data$WNN_UMAP_2), .by = "PCA_harmony_SNN_cluster_cell_type") |>
    dplyr::mutate(PCA_harmony_SNN_cluster_cell_type =
      figure_1_cell_type_labels(.data$PCA_harmony_SNN_cluster_cell_type))
  # Unlabelled sampled nuclei repel text away from the point clouds, not just
  # from other centroids. A fixed seed keeps placement matched between traits.
  obstacles <- withr::with_seed(1L, dplyr::slice_sample(data, n = min(2000L, nrow(data)))) |>
    dplyr::transmute(x = .data$WNN_UMAP_1, y = .data$WNN_UMAP_2,
      PCA_harmony_SNN_cluster_cell_type = "")
  labels <- dplyr::bind_rows(labels, obstacles)
  ggplot2::ggplot(data, ggplot2::aes(.data$WNN_UMAP_1, .data$WNN_UMAP_2)) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$score), size = .15) +
    ggrepel::geom_text_repel(data = labels, ggplot2::aes(.data$x, .data$y,
      label = .data$PCA_harmony_SNN_cluster_cell_type), size = 2, seed = 1,
      max.overlaps = Inf, min.segment.length = 0,
      point.padding = .1, box.padding = .25, force = 1, force_pull = 2,
      segment.size = .2, segment.color = "grey40",
      max.iter = 3000, max.time = Inf) +
    ggplot2::scale_colour_viridis_c(name = "TRS", limits = score_limits, na.value = "grey85") +
    ggplot2::coord_equal() + figure_1_theme() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text = ggplot2::element_blank(),
      legend.position = "bottom") +
    ggplot2::guides(colour = ggplot2::guide_colourbar(direction = "horizontal",
      barwidth = grid::unit(35, "mm"), barheight = grid::unit(2, "mm"))) +
    ggplot2::labs(title = figure_1_trait_labels(trait), x = "UMAP 1", y = "UMAP 2")
}

# Wrapping every panel prevents one panel's long axis labels from shifting the
# headings in other rows. Tags belong to the assembly, not the plot titles.
compose_figure_1 <- function(A, B, C, D, E) {
  plots <- lapply(list(A, B, C, D, E), function(plot) {
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
  (patchwork::wrap_plots(plots, design = "AB\nCC\nDD\nEE",
    heights = c(1.15, .85, 1.4, 2.75)) +
    patchwork::plot_annotation(tag_levels = "A")) &
    ggplot2::theme(plot.tag = ggplot2::element_text(family = "NimbusSan", face = "bold",
        size = 12, hjust = 0, vjust = 1), plot.tag.position = c(0, 1))
}

# Supplementary Figure S2: both workflows share one graph and store. Excluded
# targets keep their graph position with zero runtime and are left out of every
# resource sum.
SEURAT_SIGNAC_COMPARISON_SCRIPT <- "pipeline_manuscript_figures/targets.R"
SEURAT_SIGNAC_COMPARISON_WORKFLOWS <- list(
  multiomeR = list(endpoint = "multimodal_Seurat_object.8_multimodal_QC.",
    exclude = c(
      # Cluster-fragment preparation is shared by both workflows.
      "^blacklist_GRanges[.]ATAC(\\.|$)", "^BCs_per_peak_cluster_list[.]ATAC(\\.|$)",
      "^peak_calling_cluster_(names|discovery_tibble)[.]ATAC(\\.|$)",
      "^fragments_per_(cluster|peak_calling_cluster_discovery)[.]fragments[.]ATAC(\\.|$)",
      # The conventional chain performs no cell-type annotation.
      "^cluster_UCell_", "^(UCell_GEX_marker_genes_list|GEX_marker_genes_vec)(\\.|$)",
      "^metadata_w_cell_types(_unfiltered|_annotation)?_tibble[.]")),
  # Everything outside the conventional chain is imported from multiomeR.
  "Seurat/Signac" = list(endpoint = "comparison_seurat_signac_object.",
    exclude = "^(?!comparison_seurat_signac_)"))

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

# Recorded runtimes and bytes from tar_meta(); allocated cores and RAM from each
# target's current controller tier; measured CPU and RAM from slurm_monitor.
# Dynamic branches run concurrently in the critical path. Disk use counts target
# objects and file targets inside the store; the NA path of a disabled input
# file is outside it.
benchmark_resource_tibble <- function(workflows, aggregations, script = "_targets.R", store = targets::tar_config_get("store")) {
  network <- targets::tar_network(script = script, store = store,
    targets_only = TRUE, outdated = FALSE, reporter = "silent")
  manifest <- targets::tar_manifest(script = script, fields = c("name", "resources"))
  allocations <- tibble::tibble(static_target = manifest$name,
    controller_name = purrr::map_chr(manifest$resources, \(x) x$crew$controller %||% NA_character_)) |>
    dplyr::left_join(multiomeR_core_runtime_state_env$controller_resources_tibble, by = "controller_name")
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
      list(aggregation = aggregation, workflow = workflow, estimate = estimate, executions = executions)
    })
  first_start <- min(purrr::map_dbl(runs, \(run) min(as.numeric(run$executions$time) - run$executions$seconds)))
  samples <- read_slurm_monitor_samples(since = as.POSIXct(first_start - 3600, origin = "1970-01-01"))

  purrr::map_dfr(runs, \(run) {
    executions <- add_measured_resources(run$executions, samples)
    tibble::tibble(aggregation = run$aggregation,
      GEM_well_count = length(aggregation_tibble$aggregation_GEM_well_IDs[[match(run$aggregation, aggregation_tibble$aggregation)]]),
      cellranger_input_nuclei = nuclei[[run$aggregation]], workflow = run$workflow,
      critical_path_hours = run$estimate$critical_path_hours,
      allocated_CPU_hours = sum(executions$seconds * executions$cores) / 3600,
      allocated_RAM_GB_hours = sum(executions$seconds * executions$RAM_GB) / 3600,
      used_CPU_hours = sum(executions$seconds * executions$used_cores) / 3600,
      used_RAM_GB_hours = sum(executions$seconds * executions$used_RAM_GB) / 3600,
      peak_RAM_GB = max(executions$peak_RAM_GB, na.rm = TRUE),
      sampled_runtime_fraction = sum(executions$seconds[executions$sampled]) / sum(executions$seconds),
      disk_space_GB = sum(executions$bytes[executions$in_store]) / 1e9)
  })
}

seurat_signac_comparison_resource_tibble <- function(aggregations, store = targets::tar_config_get("store")) {
  benchmark_resource_tibble(SEURAT_SIGNAC_COMPARISON_WORKFLOWS, aggregations,
    script = SEURAT_SIGNAC_COMPARISON_SCRIPT, store = store)
}

plot_seurat_signac_comparison_resources <- function(data) {
  measures <- c(critical_path_hours = "Critical path (h)",
    used_CPU_hours = "CPU time (h)",
    used_RAM_GB_hours = "RAM use (GB \u00d7 h)",
    disk_space_GB = "Disk space (GB)")
  colours <- c(multiomeR = "#2a78d6", "Seurat/Signac" = "#eb6834")
  data <- data |>
    tidyr::pivot_longer(tidyselect::all_of(names(measures)), names_to = "measure") |>
    dplyr::mutate(measure = factor(measures[.data$measure], measures),
      workflow = factor(.data$workflow, names(colours)))
  largest_nuclei <- max(data$cellranger_input_nuclei)
  labels <- data |>
    dplyr::filter(.data$cellranger_input_nuclei == largest_nuclei) |>
    dplyr::select("measure", "workflow", "value") |>
    tidyr::pivot_wider(names_from = "workflow", values_from = "value") |>
    dplyr::transmute(.data$measure,
      label = sprintf("  multiomeR %.1f\u00d7 lower at %s nuclei", .data$`Seurat/Signac` / .data$multiomeR, scales::comma(largest_nuclei)))
  ggplot2::ggplot(data, ggplot2::aes(.data$cellranger_input_nuclei, .data$value, colour = .data$workflow)) +
    ggplot2::geom_line(linewidth = .5) + ggplot2::geom_point(size = 1.5) +
    ggplot2::geom_text(data = labels, ggplot2::aes(x = -Inf, y = Inf, label = .data$label),
      inherit.aes = FALSE, hjust = 0, vjust = 1.5, size = 2.1) +
    ggplot2::facet_wrap(ggplot2::vars(.data$measure), scales = "free_y", strip.position = "left") +
    ggplot2::scale_x_continuous(labels = scales::label_comma()) +
    ggplot2::scale_y_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0, .12))) +
    ggplot2::scale_colour_manual(values = colours) +
    figure_1_theme() +
    ggplot2::theme(strip.placement = "outside", legend.position = "bottom",
      strip.text = ggplot2::element_text(size = 7)) +
    ggplot2::labs(x = "Input nuclei", y = NULL, colour = NULL)
}
