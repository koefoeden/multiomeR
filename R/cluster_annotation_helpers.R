# Cluster annotation from marker lists and matched random signatures.
# Adjusted-score advantages and control tails are not identity probabilities.
# Controls are frozen per aggregation; cell types are never standardized against
# other clusters. Unassigned is an abstention, not a biological identity.
normalize_marker_panel <- function(markers, available_genes) {
  stopifnot(is.list(markers), length(markers) >= 2L, !is.null(names(markers)),
    !anyNA(names(markers)), all(nzchar(names(markers))), !anyDuplicated(names(markers)))
  markers <- lapply(markers, function(genes) unique(sub("\\+$", "", genes)))
  stopifnot(all(lengths(markers) > 0L), !anyNA(unlist(markers)))
  missing <- setdiff(sub("-$", "", unlist(markers)), available_genes)
  if (length(missing)) stop("Unavailable marker genes: ", paste(missing, collapse = ", "))
  markers
}

build_UCell_controls <- function(reference, markers, n_controls = 999L,
                                 neighbours = 50L, seed = 20260910L, exclude = character()) {
  withr::local_seed(seed)
  genes <- unique(sub("-$", "", unlist(markers, use.names = FALSE)))
  coordinates <- cbind(log1p(reference$abundance * 1e4), asin(sqrt(reference$detection)))
  spread <- apply(coordinates, 2L, stats::sd)
  spread[spread == 0] <- 1
  coordinates <- sweep(sweep(coordinates, 2L, colMeans(coordinates), `-`), 2L, spread, `/`)
  pool <- which(!reference$gene %in% c(genes, exclude) & reference$detection > 0)
  if (length(pool) < length(genes) + neighbours) stop("Insufficient eligible genes for matched controls.")
  candidate_count <- min(length(pool), max(4L * neighbours, 2L * length(genes)))
  candidates <- lapply(genes, function(gene) {
    distances <- rowSums(sweep(coordinates[pool, , drop = FALSE], 2L,
      coordinates[match(gene, reference$gene), ], `-`)^2)
    pool[head(order(distances), candidate_count)]
  })
  names(candidates) <- genes
  draws <- matrix(NA_integer_, length(genes), n_controls, dimnames = list(genes, NULL))
  for (replicate in seq_len(n_controls)) {
    used <- integer()
    for (gene in sample(genes)) {
      eligible <- head(setdiff(candidates[[gene]], used), neighbours)
      if (!length(eligible)) stop("Matched control pool exhausted for ", gene)
      chosen <- eligible[sample.int(length(eligible), 1L)]
      draws[gene, replicate] <- chosen
      used <- c(used, chosen)
    }
  }
  match_diagnostics <- do.call(rbind, lapply(genes, function(gene) {
    original <- match(gene, reference$gene)
    selected <- draws[gene, ]
    data.frame(gene, marker_abundance = reference$abundance[original],
      control_abundance = mean(reference$abundance[selected]),
      marker_detection = reference$detection[original],
      control_detection = mean(reference$detection[selected]),
      median_distance = stats::median(sqrt(rowSums(sweep(coordinates[selected, , drop = FALSE],
        2L, coordinates[original, ], `-`)^2))))
  }))
  list(markers = markers, draws = draws, reference_genes = reference$gene,
    seed = seed, n_controls = n_controls, neighbours = neighbours, excluded = exclude,
    match_diagnostics = match_diagnostics)
}

# Mean UCell = mean normalized capped-rank gap / its signature-length correction.
# The equality holds for positive signatures. Signed clipping would not commute with averaging.
score_UCell_rank_means <- function(rank_means, genes, max_rank = 1500L) {
  colMeans(rank_means[genes, , drop = FALSE]) / (1 - (length(genes) + 1) / (2 * max_rank))
}

score_UCell_controls <- function(rank_means, control, genes, max_rank = 1500L) {
  indices <- control$draws[genes, , drop = FALSE]
  indices <- match(control$reference_genes[indices], rownames(rank_means))
  stopifnot(!anyNA(indices))
  values <- rank_means[indices, , drop = FALSE]
  matrix(colMeans(matrix(values, nrow = length(genes))), nrow = control$n_controls,
    ncol = ncol(rank_means)) / (1 - (length(genes) + 1) / (2 * max_rank))
}

#' Score labels against their own matched background, without assignment gates.
score_UCell_group_evidence <- function(rank_means, control, detection, signed_means = list()) {
  markers <- control$markers
  dplyr::bind_rows(lapply(names(markers), function(label) {
    genes <- markers[[label]]
    if (any(grepl("-$", genes))) {
      scores <- signed_means[[label]]
      observed <- scores[1L, ]
      null <- scores[-1L, , drop = FALSE]
    } else {
      observed <- score_UCell_rank_means(rank_means, genes, control$max_rank)
      null <- score_UCell_controls(rank_means, control, genes, control$max_rank)
    }
    background <- apply(null, 2L, stats::quantile, probs = 0.95, names = FALSE)
    data.frame(cluster = colnames(rank_means), label,
      mean_score = unname(observed), control_mean = colMeans(null), control_q95 = background,
      excess = unname(observed) - background,
      tail_score = (1 + colSums(sweep(null, 2L, observed, `>=`))) / (nrow(null) + 1),
      broad_markers = colSums(detection[genes[!grepl("-$", genes)], , drop = FALSE] >= 0.1), marker_count = length(genes))
  }))
}

#' Assign the highest adjusted score only if it leads background and competitors.
#' Exact ties and scores at/below background remain unassigned even at delta zero.
assign_UCell_cluster_evidence <- function(evidence, min_advantage) {
  stopifnot(length(min_advantage) == 1L, is.finite(min_advantage), min_advantage >= 0)
  dplyr::bind_rows(lapply(unique(evidence$cluster), function(cluster) {
    current <- evidence[evidence$cluster == cluster, ]
    stopifnot(nrow(current) >= 2L, !anyDuplicated(current$label), all(is.finite(current$excess)))
    current <- current[order(-current$excess, current$label), ]
    best <- current$excess[1L]
    margin <- best - current$excess[2L]
    advantage <- min(best, margin)
    assigned <- best > 0 && margin > 0 && advantage >= min_advantage
    data.frame(cluster, label = if (assigned) current$label[1L] else NA_character_,
      status = if (assigned) "Assigned" else "Unassigned",
      reason = if (assigned) "sufficient_advantage" else if (best <= 0 || best < min_advantage)
        "insufficient_background_advantage" else "competing_labels",
      candidate = current$label[1L], alternative = current$label[2L],
      mean_score = current$mean_score[1L], excess = best, margin, advantage,
      competition_cutoff = best - min_advantage,
      broad_markers = current$broad_markers[1L], marker_count = current$marker_count[1L],
      single_marker = current$marker_count[1L] == 1L)
  }))
}

#' Freeze abundance/detection-matched controls using up to 50 cells per GEM well.
#' All modalities reuse this GEX reference and the same gene-to-control mappings.
prepare_cluster_UCell_controls <- function(counts_matrix, metadata_tibble, markers,
                                          reference_cells_per_well = 50L,
                                          seed = 20260910L) {
  withr::local_seed(seed)
  markers <- normalize_marker_panel(markers, rownames(counts_matrix))
  barcodes <- metadata_tibble$barcode_w_prefix
  stopifnot(length(barcodes) > 0L)
  stopifnot(!anyDuplicated(barcodes), all(barcodes %in% colnames(counts_matrix)),
    !anyNA(metadata_tibble$GEM_well_ID))
  rows <- unlist(lapply(split(seq_along(barcodes), metadata_tibble$GEM_well_ID), function(rows) {
    rows[sample.int(length(rows), min(reference_cells_per_well, length(rows)))]
  }), use.names = FALSE)
  abundance <- detection <- numeric(nrow(counts_matrix))
  chunks <- split(rows, ceiling(seq_along(rows) / 250L))
  for (chunk in chunks) {
    values <- as.matrix(counts_matrix[, barcodes[chunk], drop = FALSE])
    totals <- colSums(values)
    if (any(totals <= 0)) stop("Control reference includes cells with no GEX counts.")
    abundance <- abundance + rowSums(sweep(values, 2L, totals, `/`))
    detection <- detection + rowSums(values > 0)
  }
  reference <- data.frame(gene = rownames(counts_matrix),
    abundance = abundance / length(rows), detection = detection / length(rows))
  control <- build_UCell_controls(reference, markers, seed = seed)
  control$reference <- reference
  control$reference_barcodes <- barcodes[rows]
  control$max_rank <- min(1500L, nrow(counts_matrix))
  if (control$max_rank <= 1L || any(lengths(markers) > control$max_rank)) {
    stop("Marker signatures are too long for the UCell rank cutoff.")
  }
  control
}

#' Score one signed signature and its matched controls before averaging cells.
#' Negative genes retain their sign in every control and marker-deletion variant.
score_signed_UCell_cells <- function(gap, control, genes) {
  component <- function(selected) {
    if (!length(selected)) return(matrix(0, control$n_controls + 1L, ncol(gap)))
    rbind(score_UCell_rank_means(gap, selected, control$max_rank),
      score_UCell_controls(gap, control, selected, control$max_rank))
  }
  positive <- genes[!grepl("-$", genes)]
  negative <- sub("-$", "", genes[grepl("-$", genes)])
  scores <- component(positive) - component(negative)
  scores[] <- pmax(0, scores)
  scores
}

#' Aggregate cluster UCell statistics in bounded count chunks, clipping signed scores per cell.
#' A native kernel ranks each cell from its nonzero counts, because zero counts
#' tie below the rank cutoff, and keeps only marker/control genes. Chunks are
#' added in the same order and worker split as the reference R implementation,
#' so the sums are identical to it.
summarize_cluster_UCell_counts <- function(counts_matrix, metadata_tibble, control,
                                           cluster_column, chunk_size = 250L,
                                           workers = 2L, include_cell_scores = FALSE,
                                           native_source_file = file.path(get_project_root(), "src", "cluster_UCell_chunk.cpp")) {
  barcodes <- metadata_tibble$barcode_w_prefix
  clusters <- as.character(metadata_tibble[[cluster_column]])
  stopifnot(length(barcodes) > 0L, !anyNA(clusters), !anyDuplicated(barcodes),
    all(barcodes %in% colnames(counts_matrix)),
    identical(rownames(counts_matrix), control$reference_genes))
  library_name <- load_native_library(native_source_file, "multiomeR_UCell")
  columns <- match(barcodes, colnames(counts_matrix))
  cluster_names <- unique(clusters)
  cluster_index <- match(clusters, cluster_names) - 1L
  genes <- unique(c(sub("-$", "", unlist(control$markers, use.names = FALSE)),
    control$reference_genes[control$draws]))
  gene_position <- match(control$reference_genes, genes) - 1L
  gene_position[is.na(gene_position)] <- -1L
  replicates <- control$n_controls + 1L
  # 0-based gene indices of a signature (first column) and its matched controls.
  component <- function(selected) {
    if (!length(selected)) return(list(NULL, 1))
    index <- cbind(match(selected, genes),
      matrix(match(control$reference_genes[control$draws[selected, , drop = FALSE]], genes), length(selected)))
    stopifnot(!anyNA(index))
    list(index - 1L, 1 - (length(selected) + 1) / (2 * control$max_rank))
  }
  labels <- lapply(control$markers, function(markers) {
    c(component(markers[!grepl("-$", markers)]), component(sub("-$", "", markers[grepl("-$", markers)])))
  })
  signed_labels <- names(Filter(function(genes) any(grepl("-$", genes)), control$markers))
  chunks <- split(seq_along(barcodes), ceiling(seq_along(barcodes) / chunk_size))
  aggregate_chunks <- function(chunk_ids) {
    total <- list(
      rank_sum = matrix(0, length(genes), length(cluster_names), dimnames = list(genes, cluster_names)),
      detected_sum = matrix(0, length(genes), length(cluster_names), dimnames = list(genes, cluster_names)),
      signed_sum = lapply(stats::setNames(nm = signed_labels), function(label) {
        matrix(0, replicates, length(cluster_names), dimnames = list(NULL, cluster_names))
      }),
      cell_scores = list())
    # Read eight chunks per matrix access, but score and add each chunk separately.
    for (block in split(chunk_ids, ceiling(seq_along(chunk_ids) / 8L))) {
      values <- methods::as(counts_matrix[, columns[unlist(chunks[block])], drop = FALSE], "dgCMatrix")
      offsets <- cumsum(c(0L, lengths(chunks[block])))
      for (k in seq_along(block)) {
        rows <- chunks[[block[k]]]
        part <- .Call("multiomeR_UCell_chunk", values@p[(offsets[k] + 1L):(offsets[k + 1L] + 1L)],
          values@i, values@x, nrow(values), gene_position, length(genes), cluster_index[rows],
          length(cluster_names), labels, replicates, as.numeric(control$max_rank), PACKAGE = library_name)
        total$rank_sum <- total$rank_sum + part[[1L]]
        total$detected_sum <- total$detected_sum + part[[2L]]
        names(part[[3L]]) <- names(labels)
        for (label in signed_labels) total$signed_sum[[label]] <- total$signed_sum[[label]] + part[[3L]][[label]]
        if (include_cell_scores) {
          dimnames(part[[4L]]) <- list(barcodes[rows], names(control$markers))
          total$cell_scores[[length(total$cell_scores) + 1L]] <- part[[4L]]
        }
      }
    }
    total$cell_scores <- do.call(rbind, total$cell_scores)
    total
  }
  worker_chunks <- split(seq_along(chunks), rep(seq_len(min(workers, length(chunks))), length.out = length(chunks)))
  parts <- if (workers == 1L) lapply(worker_chunks, aggregate_chunks) else {
    parallel::mclapply(worker_chunks, aggregate_chunks, mc.cores = workers)
  }
  if (any(vapply(parts, inherits, logical(1L), "try-error"))) stop("UCell summary worker failed.")
  cells <- stats::setNames(tabulate(match(clusters, cluster_names)), cluster_names)
  cluster_means <- function(part_sums) sweep(Reduce(`+`, part_sums), 2L, cells, `/`)
  list(rank_means = cluster_means(lapply(parts, `[[`, "rank_sum")),
    detection = cluster_means(lapply(parts, `[[`, "detected_sum")),
    signed_means = lapply(stats::setNames(nm = signed_labels), function(label) {
      cluster_means(lapply(parts, function(part) part$signed_sum[[label]]))
    }),
    cells = cells,
    cell_scores = do.call(rbind, lapply(parts, `[[`, "cell_scores")))
}

#' Cache cluster score evidence independently of the assignment threshold.
score_cluster_UCell_summaries <- function(summaries, control) {
  evidence <- score_UCell_group_evidence(summaries$rank_means, control, summaries$detection, summaries$signed_means)
  evidence$cells <- unname(summaries$cells[evidence$cluster])
  list(evidence = evidence, cell_scores = summaries$cell_scores,
    settings = list(max_rank = control$max_rank, background_quantile = 0.95, diagnostic_marker_detection = 0.1))
}

prepare_cluster_UCell_evidence <- function(counts_matrix, metadata_tibble, control,
                                           cluster_column, include_cell_scores = FALSE, workers = 2L,
                                           native_source_file = file.path(get_project_root(), "src", "cluster_UCell_chunk.cpp")) {
  summaries <- summarize_cluster_UCell_counts(counts_matrix, metadata_tibble, control,
    cluster_column, workers = workers, include_cell_scores = include_cell_scores,
    native_source_file = native_source_file)
  score_cluster_UCell_summaries(summaries, control)
}

#' Apply the assignment threshold to cached cluster evidence.
evaluate_cluster_UCell_evidence <- function(scored, min_advantage) {
  decisions <- assign_UCell_cluster_evidence(scored$evidence, min_advantage)
  decisions$cells <- scored$evidence$cells[match(decisions$cluster, scored$evidence$cluster)]
  list(decisions = decisions, evidence = scored$evidence,
    settings = c(scored$settings, list(min_advantage = min_advantage)), cell_scores = scored$cell_scores)
}

#' Preserve existing metadata/score columns and attach explicit annotation status.
#' Abstaining clusters get distinct scDblFinder groups; they must not be collapsed
#' into one artificial unassigned population for doublet detection.
add_cluster_UCell_annotations <- function(metadata_tibble, annotation, cluster_column) {
  decisions <- annotation$decisions
  index <- match(as.character(metadata_tibble[[cluster_column]]), decisions$cluster)
  if (anyNA(index) || anyDuplicated(decisions$cluster)) stop("Cluster annotation does not cover metadata uniquely.")
  label <- ifelse(decisions$status == "Assigned", decisions$label, decisions$status)[index]
  metadata_tibble[[paste0(cluster_column, "_cell_type")]] <- get_mixsorted_factor(label)
  metadata_tibble[[paste0(cluster_column, "_named")]] <- get_mixsorted_factor(
    paste(metadata_tibble[[cluster_column]], label, sep = "-"))
  metadata_tibble[[paste0(cluster_column, "_annotation_status")]] <- decisions$status[index]
  metadata_tibble[[paste0(cluster_column, "_scDblFinder_group")]] <- ifelse(
    decisions$status[index] == "Assigned", label,
    paste0(label, "_cluster_", metadata_tibble[[cluster_column]]))
  if (!is.null(annotation$cell_scores)) {
    scores <- annotation$cell_scores[metadata_tibble$barcode_w_prefix, , drop = FALSE]
    for (name in colnames(scores)) metadata_tibble[[name]] <- scores[, name]
  }
  metadata_tibble
}

#' Summarize each marker set's support and competition across GEX clusters.
summarize_UCell_marker_sets <- function(annotation) {
  evidence <- annotation$evidence
  decisions <- annotation$decisions
  dplyr::bind_rows(lapply(unique(evidence$label), function(label) {
    current <- evidence[evidence$label == label, ]
    competitors <- evidence[evidence$label != label, ]
    competitors <- competitors[order(-competitors$excess, competitors$label), ]
    competitors <- competitors[match(current$cluster, competitors$cluster), ]
    advantage <- current$excess - pmax(0, competitors$excess)
    best <- order(-advantage, -current$excess, current$cluster)[1L]
    data.frame(marker_set = label, marker_count = current$marker_count[1L],
      clusters_assessed = nrow(current), clusters_above_background = sum(current$excess > 0),
      clusters_above_min_advantage = sum(current$excess > 0 & current$excess >= annotation$settings$min_advantage),
      clusters_leading = sum(decisions$candidate == label),
      clusters_assigned = sum(decisions$label == label, na.rm = TRUE),
      max_adjusted_score = max(current$excess),
      max_score_cluster = current$cluster[which.max(current$excess)],
      best_advantage = advantage[best], best_advantage_cluster = current$cluster[best],
      closest_competitor = competitors$label[best],
      competitor_adjusted_score = competitors$excess[best],
      min_advantage = annotation$settings$min_advantage)
  }))
}

#' Export review tables without duplicating per-cell score matrices.
save_cluster_UCell_diagnostics <- function(annotation, control,
                                            output_dir = get_structured_output_path(list_output = TRUE)) {
  fs::dir_create(output_dir)
  previous_files <- file.path(output_dir, c("clusters.tsv", "marker_evidence.tsv",
    "control_gene_matching.tsv", "settings.rds", "method.txt"))
  unlink(previous_files[file.exists(previous_files)])
  readr::write_tsv(annotation$decisions, file.path(output_dir, "clusters.tsv"))
  readr::write_tsv(annotation$evidence, file.path(output_dir, "marker_evidence.tsv"))
  readr::write_tsv(control$match_diagnostics, file.path(output_dir, "control_gene_matching.tsv"))
  saveRDS(list(settings = annotation$settings, markers = control$markers,
    seed = control$seed, n_controls = control$n_controls, neighbours = control$neighbours,
    reference_barcodes = control$reference_barcodes), file.path(output_dir, "settings.rds"))
  writeLines(c("Assignment uses the best adjusted UCell score's advantage over background and all other labels.",
    paste("Minimum advantage:", annotation$settings$min_advantage),
    "Adjusted score = observed mean UCell minus the label-specific control 95th percentile.",
    "Marker detection and control tails are diagnostic only.",
    "Scores are not calibrated identity probabilities. No mixture detection."), file.path(output_dir, "method.txt"))
  output_dir
}

#' Plot every label on shared axes; paginate clusters without changing scales.
plot_cluster_UCell_advantages <- function(annotation, clusters_per_page = 6L) {
  decisions <- annotation$decisions
  min_advantage <- annotation$settings$min_advantage
  cluster_order <- levels(get_mixsorted_factor(decisions$cluster))
  labels <- unique(annotation$evidence$label)
  data <- dplyr::left_join(annotation$evidence,
    decisions[, c("cluster", "candidate", "status", "label", "competition_cutoff")],
    by = "cluster", suffix = c("", "_assigned"))
  data$facet <- ifelse(data$status == "Assigned", paste(data$cluster, data$label_assigned, sep = "-"),
    paste0(data$cluster, "-Unassigned (candidate: ", data$candidate, ")"))
  data$label <- factor(data$label, levels = labels)
  data$role <- ifelse(as.character(data$label) == data$candidate, "Leading candidate",
    ifelse(data$excess > data$competition_cutoff |
      (min_advantage == 0 & data$excess == data$competition_cutoff), "Competing label", "Other label"))
  data$role <- factor(data$role, levels = c("Leading candidate", "Competing label", "Other label"))
  data$plot_score <- ifelse(data$role == "Leading candidate", data$excess, pmax(0, data$excess))
  limits <- range(c(0, data$plot_score, data$competition_cutoff))
  pages <- split(cluster_order, ceiling(seq_along(cluster_order) / clusters_per_page))
  plots <- lapply(pages, function(clusters) {
    current <- data[data$cluster %in% clusters, ]
    facet_order <- current$facet[match(clusters, current$cluster)]
    current$facet <- factor(current$facet, levels = facet_order)
    cutoffs <- unique(current[, c("facet", "competition_cutoff")])
    ggplot2::ggplot(current, ggplot2::aes(x = label, y = plot_score, fill = role)) +
      ggplot2::geom_col(width = 0.8) +
      ggplot2::geom_hline(yintercept = 0, color = "grey25", linewidth = 0.5) +
      ggplot2::geom_hline(data = cutoffs, ggplot2::aes(yintercept = competition_cutoff),
        linetype = "dashed", color = "grey25", linewidth = 0.5) +
      ggplot2::facet_wrap(ggplot2::vars(facet), ncol = 2, axes = "all_x") +
      ggplot2::scale_x_discrete(drop = FALSE, labels = function(x) gsub("_", " ", x)) +
      ggplot2::scale_y_continuous(limits = limits, expand = ggplot2::expansion(mult = 0.06)) +
      ggplot2::scale_fill_manual(values = c("Leading candidate" = "#0072B2", "Competing label" = "#D55E00", "Other label" = "grey75"), drop = FALSE) +
      ggplot2::labs(x = "Candidate cell type", y = "UCell score above matched background", fill = NULL,
        title = "Cluster assignment by adjusted-score advantage",
        subtitle = paste0(
          "Each panel compares cell-type marker evidence in one cluster: mean UCell score minus the 95th percentile of matched random-marker scores.\n",
          "Blue marks the leading candidate; the dashed line is its adjusted score minus the minimum advantage (", format(min_advantage), "). Zero marks matched background.\n",
          "Assignment requires the dashed line at or above zero and no competing bar above it. Orange bars flag close competitors; a line below zero flags weak evidence.\n",
          "Use these gaps to judge how clearly the markers distinguish a label and to calibrate the threshold; higher thresholds require stronger separation."),
        caption = paste0("Shared axes across all pages. Exact ties remain unassigned; adjusted scores are not assignment probabilities.\n",
          "Non-leading negative scores are displayed at zero; assignments use the original scores. Raising the threshold moves the dashed line downward.")) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 65, hjust = 1),
        legend.position = "bottom", panel.grid.major.x = ggplot2::element_blank())
  })
  stats::setNames(plots, sprintf("page_%02d", seq_along(plots)))
}
