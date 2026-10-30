#' Count raw XIST and chromosome Y UMIs of Cell Ranger-called nuclei
#'
#' Read each GEM well's Cell Ranger filtered matrix, which holds every called
#' nucleus, including those that CellBender drops.
#'
#' @param cellranger_summary_files Cell Ranger summary CSV paths, one per GEM well.
#' @param GEM_well_IDs GEM well IDs in the same order.
#' @param gene_features_df Aggregation gene features in Cell Ranger matrix order.
#' @param cores Number of GEM wells read in parallel.
#' @return A tibble with `barcode_w_prefix`, `XIST`, `chrY` (all chromosome Y
#'   genes) and `total_UMIs`.
#' @keywords internal

count_cellranger_sex_gene_UMIs <- function(cellranger_summary_files, GEM_well_IDs, gene_features_df, cores = 1) {
  indicator <- cbind(
    XIST = as.numeric(toupper(gene_features_df$name) == "XIST"),
    chrY = as.numeric(gene_features_df$seqnames == "chrY"),
    total_UMIs = 1
  )
  parallel::mclapply(seq_along(GEM_well_IDs), \(well) {
    counts <- file.path(dirname(cellranger_summary_files[[well]]), "filtered_feature_bc_matrix.h5") |>
      BPCells::open_matrix_10x_hdf5(feature_type = "Gene Expression") |>
      BPCells::convert_matrix_type("double")
    tibble::as_tibble(t(t(indicator) %*% counts)) |>
      dplyr::mutate(barcode_w_prefix = paste0(GEM_well_IDs[[well]], "_", colnames(counts)), .before = 1)
  }, mc.cores = cores) |>
    dplyr::bind_rows()
}

#' Summarize QC metrics and sex-gene expression per sample
#'
#' @param metadata_tibble Nuclei that pass per-well QC, with `dataset`,
#'   `donor_id` (the sample) and the QC metric columns.
#' @param sex_gene_UMIs_tibble Output of `count_cellranger_sex_gene_UMIs()`.
#' @param donor_metadata_tibble Donor metadata; its optional `sex` column is the recorded sex.
#' @return One row per dataset and sample: nuclei, quartiles of TSS enrichment,
#'   fraction of fragments in peaks and mitochondrial UMI percentage, XIST and
#'   chromosome Y UMIs per 10,000, and the recorded sex.
#' @keywords internal

summarize_sample_QC <- function(metadata_tibble, sex_gene_UMIs_tibble, donor_metadata_tibble) {
  if (!"sex" %in% colnames(donor_metadata_tibble)) donor_metadata_tibble$sex <- NA_character_
  metadata_tibble |>
    dplyr::select("barcode_w_prefix", "dataset", "donor_id", "TSS.enrichment", "ATAC_peak_region_frac", "RNA_mito_percent") |>
    dplyr::inner_join(sex_gene_UMIs_tibble, by = "barcode_w_prefix") |>
    dplyr::summarise(
      n_nuclei = dplyr::n(),
      dplyr::across(c("TSS.enrichment", "ATAC_peak_region_frac", "RNA_mito_percent"), list(
        q25 = \(x) stats::quantile(x, 0.25, na.rm = TRUE, names = FALSE),
        median = \(x) stats::median(x, na.rm = TRUE),
        q75 = \(x) stats::quantile(x, 0.75, na.rm = TRUE, names = FALSE))),
      XIST_per_10k = 1e4 * sum(.data$XIST) / sum(.data$total_UMIs),
      chrY_per_10k = 1e4 * sum(.data$chrY) / sum(.data$total_UMIs),
      .by = c("dataset", "donor_id")
    ) |>
    dplyr::left_join(dplyr::select(donor_metadata_tibble, "donor_id", recorded_sex = "sex"), by = "donor_id")
}

#' Plot per-sample QC metrics and flag poor samples
#'
#' @param sample_QC_tibble Output of `summarize_sample_QC()`.
#' @return A ggplot.
#' @keywords internal

plot_sample_QC_metrics <- function(sample_QC_tibble) {
  metric_labels <- c(TSS.enrichment = "TSS enrichment", ATAC_peak_region_frac = "Fraction of fragments\nin peaks",
    RNA_mito_percent = "Mitochondrial UMIs (%)")
  poor_direction <- c(TSS.enrichment = -1, ATAC_peak_region_frac = -1, RNA_mito_percent = 1)
  plot_data <- sample_QC_tibble |>
    tidyr::pivot_longer(dplyr::starts_with(names(metric_labels)), names_to = c("metric", ".value"),
      names_pattern = "(.*)_(q25|median|q75)$") |>
    dplyr::mutate(
      MAD = stats::mad(.data$median),
      flagged = .data$MAD > 0 & poor_direction[.data$metric] * (.data$median - stats::median(.data$median)) > 3 * .data$MAD,
      .by = c("dataset", "metric")
    ) |>
    dplyr::mutate(metric = factor(metric_labels[.data$metric], levels = metric_labels))
  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$donor_id, y = .data$median, ymin = .data$q25, ymax = .data$q75,
    colour = .data$flagged)) +
    ggplot2::geom_pointrange(size = 0.2) +
    ggplot2::scale_colour_manual(values = c(`FALSE` = "grey30", `TRUE` = "#D55E00"), labels = c(`FALSE` = "No", `TRUE` = "Yes"),
      name = "Flagged") +
    ggplot2::facet_grid(rows = ggplot2::vars(.data$metric), cols = ggplot2::vars(.data$dataset), scales = "free",
      space = "free_x", switch = "y") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
      strip.text.x = ggplot2::element_text(angle = 90, hjust = 0), strip.placement = "outside") +
    ggplot2::labs(
      title = "Per-sample chromatin and RNA quality",
      subtitle = stringr::str_wrap(paste(
        "Look for flagged samples: their median lies more than 3 MADs on the poor side of the other samples in their dataset.",
        "Per-well rules judge nuclei one at a time, so the nuclei of a degraded sample can still pass them."), width = 100),
      x = "Sample", y = NULL,
      caption = stringr::str_wrap(paste(
        "Points show each sample's median over its nuclei that pass per-well QC, and bars the interquartile range.",
        "MAD is the median absolute deviation of the sample medians in the dataset; a metric whose medians do not vary",
        "within a dataset, as with one sample, flags nothing."), width = 120)
    )
}

#' Plot each sample's XIST and chromosome Y expression against its recorded sex
#'
#' A sample is labelled when it lies nearer the median of the other recorded
#' sex, on log scales, than the median of its own.
#'
#' @param sample_QC_tibble Output of `summarize_sample_QC()`.
#' @return A ggplot.
#' @keywords internal

plot_sample_sex_check <- function(sample_QC_tibble) {
  plot_data <- sample_QC_tibble |>
    dplyr::mutate(
      log_XIST = log10(.data$XIST_per_10k + 0.1), log_chrY = log10(.data$chrY_per_10k + 0.1),
      recorded_sex = tidyr::replace_na(as.character(.data$recorded_sex), "Not recorded")
    )
  centroids <- plot_data |>
    dplyr::filter(.data$recorded_sex != "Not recorded") |>
    dplyr::summarise(dplyr::across(c("log_XIST", "log_chrY"), stats::median), .by = "recorded_sex")
  plot_data$mismatch <- if (nrow(centroids) > 1) {
    distances <- outer(plot_data$log_XIST, centroids$log_XIST, "-")^2 + outer(plot_data$log_chrY, centroids$log_chrY, "-")^2
    plot_data$recorded_sex != "Not recorded" & centroids$recorded_sex[max.col(-distances, ties.method = "first")] != plot_data$recorded_sex
  } else {
    FALSE
  }
  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$XIST_per_10k + 0.1, y = .data$chrY_per_10k + 0.1)) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$recorded_sex), alpha = 0.8) +
    ggrepel::geom_text_repel(data = dplyr::filter(plot_data, .data$mismatch), ggplot2::aes(label = .data$donor_id),
      size = 3, min.segment.length = 0) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Sample sex from XIST and chromosome Y expression",
      subtitle = stringr::str_wrap(paste(
        "Look for samples whose colour differs from their neighbours'. Labelled samples lie nearer the other recorded sex,",
        "which points to a sample swap or a metadata error."), width = 100),
      x = "XIST UMIs per 10,000 (+0.1, log scale)", y = "Chromosome Y gene UMIs per 10,000 (+0.1, log scale)",
      colour = "Recorded sex",
      caption = stringr::str_wrap(paste(
        "Each point is a sample, summed over its nuclei that pass per-well QC, from raw Cell Ranger counts.",
        "Ambient RNA from other samples in a pooled well sets a floor on both axes.",
        "A sample is labelled when it is nearer the median of the other recorded sex than of its own."), width = 120)
    )
}

#' Count ATAC fragments of Cell Ranger-called nuclei in marker gene windows
#'
#' Chromatin marks a nucleus's identity even when its droplet carries RNA of
#' other cells, as many nuclei that per-well QC removes do. Each marker gene's
#' window is its body and 2 kb upstream.
#'
#' @param fragments_list BPCells fragments, one per GEM well.
#' @param GEM_well_IDs GEM well IDs in the same order.
#' @param cellranger_barcodes_tsvs Files of Cell Ranger-called barcodes, in the same order.
#' @param genes Gene annotation GRanges with a `gene_name` column.
#' @param marker_genes_list Named list of UCell-syntax marker sets; only positive markers count.
#' @param cores Number of GEM wells read in parallel.
#' @return A tibble with `GEM_well_ID`, `barcode_w_prefix`, one fragment count
#'   column per marker set and `total_fragments`, all fragments outside chrM.
#' @keywords internal

count_marker_window_fragments <- function(fragments_list, GEM_well_IDs, cellranger_barcodes_tsvs, genes, marker_genes_list,
  cores = 1) {
  marker_genes_list <- purrr::map(marker_genes_list, \(markers) stringr::str_remove(markers[!stringr::str_ends(markers, "-")], "[+]$"))
  windows <- genes[genes$gene_name %in% unlist(marker_genes_list)]
  windows <- GenomicRanges::resize(windows, GenomicRanges::width(windows) + 2000, fix = "end")
  parallel::mclapply(seq_along(GEM_well_IDs), \(well) {
    fragments <- BPCells::select_cells(fragments_list[[well]],
      paste0(GEM_well_IDs[[well]], "_", readLines(cellranger_barcodes_tsvs[[well]])))
    chromosomes <- setdiff(BPCells::chrNames(fragments), "chrM")
    gene_windows <- windows[as.character(GenomicRanges::seqnames(windows)) %in% chromosomes]
    # Whole-chromosome ranges count each nucleus's fragments.
    ranges <- tibble::tibble(
      chr = c(as.character(GenomicRanges::seqnames(gene_windows)), chromosomes),
      start = c(GenomicRanges::start(gene_windows) - 1L, rep(0L, length(chromosomes))),
      end = c(GenomicRanges::end(gene_windows), rep(1e9, length(chromosomes)))
    )
    range_genes <- c(gene_windows$gene_name, rep(NA, length(chromosomes)))
    indicator <- cbind(
      vapply(marker_genes_list, \(markers) as.numeric(range_genes %in% markers), numeric(nrow(ranges))),
      total_fragments = as.numeric(is.na(range_genes))
    )
    # peak_matrix() returns its ranges sorted, so sort them, and their indicator rows, first.
    ordering <- BPCells::order_ranges(ranges, BPCells::chrNames(fragments))
    counts <- BPCells::peak_matrix(fragments, ranges[ordering, ], mode = "fragments") |>
      BPCells::convert_matrix_type("double")
    tibble::as_tibble(t(t(indicator[ordering, ]) %*% counts)) |>
      dplyr::mutate(GEM_well_ID = GEM_well_IDs[[well]], barcode_w_prefix = BPCells::cellNames(fragments), .before = 1)
  }, mc.cores = cores) |>
    dplyr::bind_rows()
}

#' Classify nuclei by their share of fragments in each marker set's windows
#'
#' Each set's fragment share is scaled from the 90th percentile in nuclei with
#' other cluster labels (0) to the median in nuclei with its own label (1). A
#' nucleus takes the top-scoring set among the cell types labelled in its GEM
#' well, when that score reaches `min_score`. Sets that no cluster
#' carries, or whose labelled nuclei do not exceed that background, are not scored.
#'
#' @param marker_fragments_tibble Output of `count_marker_window_fragments()`.
#' @param labelled_metadata_tibble Nuclei with cluster labels, keyed by
#'   `barcode_w_prefix`, with `GEM_well_ID`.
#' @param set_names Marker set names, in display order.
#' @param label_col Column of `labelled_metadata_tibble` holding the cluster labels.
#' @param min_score Minimum scaled score for a nucleus to be classified.
#' @return A tibble with `barcode_w_prefix`, `scored_cell_type` (a set name or
#'   "Unclassified") and `cluster_label` (the nucleus's label when it is a scored set).
#' @keywords internal

classify_nuclei_by_marker_sets <- function(marker_fragments_tibble, labelled_metadata_tibble, set_names,
  label_col = "PCA_harmony_SNN_cluster_cell_type", min_score = 0.4) {
  all_labels <- as.character(labelled_metadata_tibble[[label_col]])
  labels <- all_labels[match(marker_fragments_tibble$barcode_w_prefix, labelled_metadata_tibble$barcode_w_prefix)]
  shares <- as.matrix(marker_fragments_tibble[set_names]) / pmax(marker_fragments_tibble$total_fragments, 1)
  labelled <- labels %in% set_names
  calibration <- tibble::tibble(set = intersect(set_names, labels)) |>
    dplyr::mutate(
      level = purrr::map_dbl(.data$set, \(set) stats::median(shares[labels %in% set, set])),
      background = purrr::map_dbl(.data$set, \(set) stats::quantile(shares[labelled & labels != set, set], 0.9, names = FALSE))
    ) |>
    dplyr::filter(.data$level > .data$background)
  scores <- sweep(sweep(shares[, calibration$set, drop = FALSE], 2, calibration$background), 2,
    calibration$level - calibration$background, "/")
  # A nucleus can only be of a cell type labelled in its well; wells without labelled nuclei get none.
  well_types <- unclass(table(labelled_metadata_tibble$GEM_well_ID, factor(all_labels, levels = calibration$set))) > 0
  allowed <- well_types[match(marker_fragments_tibble$GEM_well_ID, rownames(well_types)), , drop = FALSE]
  scores[is.na(allowed) | !allowed] <- -Inf
  best <- max.col(scores, ties.method = "first")
  tibble::tibble(
    barcode_w_prefix = marker_fragments_tibble$barcode_w_prefix,
    scored_cell_type = ifelse(scores[cbind(seq_along(best), best)] >= min_score, calibration$set[best], "Unclassified"),
    cluster_label = ifelse(labels %in% calibration$set, labels, NA_character_)
  )
}

#' Compare the scored cell types of nuclei that per-well QC removes and keeps
#'
#' @param scored_nuclei_tibble Output of `classify_nuclei_by_marker_sets()` for
#'   all Cell Ranger-called nuclei.
#' @param excluded_barcodes_lists Per GEM well, a named list of barcodes per exclusion rule.
#' @param set_names Marker set names, in display order.
#' @return A ggplot heatmap with one row per group of nuclei and one column per scored cell type.
#' @keywords internal

plot_excluded_nuclei_cell_types <- function(scored_nuclei_tibble, excluded_barcodes_lists, set_names) {
  removed <- purrr::map(excluded_barcodes_lists, \(rules) tibble::enframe(rules, name = "group", value = "barcode_w_prefix")) |>
    dplyr::bind_rows() |>
    tidyr::unnest("barcode_w_prefix")
  removed_any <- unique(removed$barcode_w_prefix)
  counts <- dplyr::bind_rows(
    tibble::tibble(group = "Kept by every rule",
      barcode_w_prefix = setdiff(scored_nuclei_tibble$barcode_w_prefix, removed_any)),
    tibble::tibble(group = "Removed by any rule", barcode_w_prefix = removed_any),
    removed
  ) |>
    dplyr::inner_join(scored_nuclei_tibble, by = "barcode_w_prefix") |>
    dplyr::summarise(n = dplyr::n(), .by = c("group", "scored_cell_type")) |>
    dplyr::mutate(group_n = sum(.data$n), percent = 100 * .data$n / .data$group_n, .by = "group")
  removed_percent <- counts |>
    dplyr::filter(.data$group %in% c("Kept by every rule", "Removed by any rule")) |>
    dplyr::summarise(percent = 100 * sum(.data$n[.data$group == "Removed by any rule"]) / sum(.data$n),
      .by = "scored_cell_type") |>
    tibble::deframe()
  type_levels <- intersect(c(set_names, "Unclassified"), counts$scored_cell_type)
  group_levels <- counts |>
    dplyr::distinct(.data$group, .data$group_n) |>
    dplyr::arrange(!.data$group %in% c("Kept by every rule", "Removed by any rule"), dplyr::desc(.data$group_n))
  plot_data <- counts |>
    dplyr::mutate(
      cell_type = factor(.data$scored_cell_type, levels = type_levels,
        labels = paste0(type_levels, " (", round(removed_percent[type_levels]), "%)")),
      group = factor(.data$group, levels = rev(group_levels$group),
        labels = rev(paste0(group_levels$group, " (", scales::comma(group_levels$group_n), ")")))
    )
  # Cluster labels exist only for nuclei that passed per-well QC.
  labelled <- dplyr::filter(scored_nuclei_tibble, !is.na(.data$cluster_label))
  classified <- dplyr::filter(labelled, .data$scored_cell_type != "Unclassified")
  agreement <- c(matching = mean(classified$scored_cell_type == classified$cluster_label),
    unclassified = 1 - nrow(classified) / nrow(labelled))
  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$cell_type, y = .data$group, fill = .data$percent)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(.data$percent >= 0.5, round(.data$percent), "")), size = 2.5) +
    ggplot2::scale_fill_distiller(palette = "Blues", direction = 1, name = "% of row") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), panel.grid = ggplot2::element_blank()) +
    ggplot2::labs(
      title = "Scored cell types of nuclei that per-well QC removes",
      subtitle = stringr::str_wrap(paste(
        "Compare each row with the kept nuclei: a cell type that makes up more of the removed nuclei than of the kept ones",
        "is lost selectively, which biases its composition estimates. Column labels give the share of each type's nuclei",
        "that any rule removes."), width = 110),
      x = "Scored cell type (share removed)", y = "Nuclei (number)",
      caption = stringr::str_wrap(paste0(
        "Cell Ranger-called nuclei, scored on chromatin, which marks a nucleus's identity even when its droplet carries ",
        "RNA of other cells, as many removed nuclei do. A nucleus's share of ATAC fragments in each marker set's gene ",
        "windows (gene body and 2 kb upstream) is scaled from the 90th percentile in nuclei with other cluster labels (0) ",
        "to the median in nuclei with that label (1). The nucleus takes the top set among the cell types labelled in its ",
        "GEM well when it reaches 0.4, and stays unclassified otherwise. Among kept nuclei with a cluster label, ",
        scales::percent(agreement[["matching"]], accuracy = 0.1), " of those classified match it and ",
        scales::percent(agreement[["unclassified"]], accuracy = 0.1), " stay unclassified, so compare rows rather ",
        "than trusting single nuclei. A nucleus that fails several rules counts in each of their rows."), width = 170)
    )
}
