# Cluster annotation from marker lists and matched random signatures.
# Adjusted-score advantages and control tails are not identity probabilities.
# Controls are frozen per aggregation; cell types are never standardized against
# other clusters. Unassigned is an abstention, not a biological identity.
annotation_group_seed <- function(key, seed = 20260910L) {
  value <- as.double(seed)
  for (code in utf8ToInt(enc2utf8(key))) value <- (131 * value + code) %% 2147483647
  as.integer(value)
}

make_annotation_blocks <- function(clusters, samples, n_blocks = 10L, seed = 20260910L) {
  withr::local_seed(seed)
  blocks <- integer(length(clusters))
  groups <- split(seq_along(clusters), paste(clusters, samples, sep = "/"))
  for (key in names(groups)) {
    rows <- groups[[key]]
    set.seed(annotation_group_seed(key, seed))
    assignments <- rep(sample(seq_len(n_blocks)), length.out = length(rows))
    blocks[rows] <- assignments[sample.int(length(rows))]
  }
  blocks
}

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
  markers <- normalize_marker_panel(markers, reference$gene)
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
  set.seed(seed)
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
  stopifnot(identical(dimnames(rank_means), dimnames(detection)),
    all(is.finite(rank_means)), all(rank_means >= 0), all(rank_means <= 1))
  markers <- normalize_marker_panel(control$markers, rownames(rank_means))
  if (control$max_rank <= 1L || any(lengths(markers) > control$max_rank)) {
    stop("Marker signatures are too long for the UCell rank cutoff.")
  }
  dplyr::bind_rows(lapply(names(markers), function(label) {
    genes <- markers[[label]]
    if (any(grepl("-$", genes))) {
      scores <- signed_means[[label]][[1L]]
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
assign_UCell_cluster_evidence <- function(evidence, min_advantage = 0.05) {
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

#' Aggregate UCell statistics in bounded count chunks, clipping signed scores per cell.
#' Cluster, stratified deletion-block, and GEM-well summaries share one rank pass.
#' Only marker/control genes are retained after ranking all genes in each cell.
summarize_cluster_UCell_counts <- function(counts_matrix, metadata_tibble, control,
                                           cluster_column, chunk_size = 250L,
                                           workers = 2L, include_cell_scores = FALSE,
                                           n_blocks = 10L) {
  barcodes <- metadata_tibble$barcode_w_prefix
  clusters <- as.character(metadata_tibble[[cluster_column]])
  wells <- as.character(metadata_tibble$GEM_well_ID)
  stopifnot(length(clusters) == length(barcodes), length(barcodes) > 0L,
    workers >= 1L, chunk_size >= 1L, n_blocks >= 2L,
    !anyNA(clusters), !anyNA(wells), !anyDuplicated(barcodes),
    all(barcodes %in% colnames(counts_matrix)),
    identical(rownames(counts_matrix), control$reference_genes))
  blocks <- make_annotation_blocks(clusters, wells, n_blocks, control$seed)
  group_rows <- rbind(
    data.frame(kind = "cluster", cluster = clusters, subgroup = ""),
    data.frame(kind = "block", cluster = clusters, subgroup = as.character(blocks)),
    data.frame(kind = "well", cluster = clusters, subgroup = wells)
  )
  groups <- unique(group_rows)
  groups$key <- as.character(seq_len(nrow(groups)))
  group_ids <- dplyr::left_join(group_rows, groups, by = c("kind", "cluster", "subgroup"))$key
  group_indices <- split(as.integer(group_ids), rep(seq_len(3L), each = length(barcodes)))
  groups$cells <- tabulate(as.integer(group_ids), nbins = nrow(groups))
  genes <- unique(c(sub("-$", "", unlist(control$markers, use.names = FALSE)),
    control$reference_genes[control$draws]))
  chunks <- split(seq_along(barcodes), ceiling(seq_along(barcodes) / chunk_size))
  signed_markers <- Filter(function(genes) any(grepl("-$", genes)), control$markers)
  signed_variants <- lapply(signed_markers, function(genes) {
    c(list(genes), if (length(genes) > 1L) lapply(genes, function(gene) setdiff(genes, gene)))
  })
  aggregate_chunks <- function(chunk_ids) {
    signed_sum <- lapply(signed_variants, function(variants) lapply(variants, function(genes) {
      matrix(0, control$n_controls + 1L, nrow(groups), dimnames = list(NULL, groups$key))
    }))
    rank_sum <- detected_sum <- matrix(0, length(genes), nrow(groups),
      dimnames = list(genes, groups$key))
    cell_scores <- list()
    for (chunk_id in chunk_ids) {
      rows <- chunks[[chunk_id]]
      values <- as.matrix(counts_matrix[, barcodes[rows], drop = FALSE])
      ranks <- rank_UCell_count_chunk(values)
      gap <- (control$max_rank - pmin(ranks[genes, , drop = FALSE], control$max_rank)) / control$max_rank
      detected <- values[genes, , drop = FALSE] > 0
      for (indices in group_indices) {
        sums <- rowsum(t(gap), indices[rows], reorder = FALSE)
        columns <- rownames(sums)
        rank_sum[, columns] <- rank_sum[, columns, drop = FALSE] + t(sums)
        detected_sum[, columns] <- detected_sum[, columns, drop = FALSE] +
          t(rowsum(t(detected * 1), indices[rows], reorder = FALSE))
      }
      for (label in names(signed_variants)) {
        for (variant in seq_along(signed_variants[[label]])) {
          scores <- score_signed_UCell_cells(gap, control, signed_variants[[label]][[variant]])
          for (indices in group_indices) {
            sums <- rowsum(t(scores), indices[rows], reorder = FALSE)
            columns <- rownames(sums)
            signed_sum[[label]][[variant]][, columns] <-
              signed_sum[[label]][[variant]][, columns, drop = FALSE] + t(sums)
          }
        }
      }
      if (include_cell_scores) {
        scores <- vapply(control$markers, function(marker_genes) {
          if (any(grepl("-$", marker_genes))) {
            score_signed_UCell_cells(gap, control, marker_genes)[1L, ]
          } else score_UCell_rank_means(gap, marker_genes, control$max_rank)
        }, numeric(length(rows)))
        dim(scores) <- c(length(rows), length(control$markers))
        dimnames(scores) <- list(barcodes[rows], names(control$markers))
        cell_scores[[as.character(chunk_id)]] <- scores
      }
    }
    list(rank_sum = rank_sum, detected_sum = detected_sum, signed_sum = signed_sum, cell_scores = do.call(rbind, cell_scores))
  }
  worker_chunks <- split(seq_along(chunks), rep(seq_len(min(workers, length(chunks))), length.out = length(chunks)))
  parts <- if (workers == 1L) lapply(worker_chunks, aggregate_chunks) else {
    parallel::mclapply(worker_chunks, aggregate_chunks, mc.cores = workers)
  }
  if (any(vapply(parts, inherits, logical(1L), "try-error"))) stop("UCell summary worker failed.")
  list(rank_means = sweep(Reduce(`+`, lapply(parts, `[[`, "rank_sum")), 2L, groups$cells, `/`),
    detection = sweep(Reduce(`+`, lapply(parts, `[[`, "detected_sum")), 2L, groups$cells, `/`),
    signed_means = stats::setNames(lapply(names(signed_variants), function(label) {
      lapply(seq_along(signed_variants[[label]]), function(variant) {
        sweep(Reduce(`+`, lapply(parts, function(part) part$signed_sum[[label]][[variant]])),
          2L, groups$cells, `/`)
      })
    }), names(signed_variants)),
    groups = groups, n_blocks = n_blocks,
    cell_scores = do.call(rbind, lapply(parts, `[[`, "cell_scores")))
}

#' Cache compact score evidence and perturbations independently of stringency.
score_cluster_UCell_summaries <- function(summaries, control) {
  groups <- summaries$groups
  main <- groups[groups$kind == "cluster", ]
  score_groups <- function(keys) score_UCell_group_evidence(
    summaries$rank_means[, keys, drop = FALSE], control, summaries$detection[, keys, drop = FALSE],
    lapply(summaries$signed_means, function(variants) lapply(variants, function(values) values[, keys, drop = FALSE])))
  evidence <- score_groups(main$key)
  evidence$cluster <- main$cluster[match(evidence$cluster, main$key)]
  evidence$cells <- main$cells[match(evidence$cluster, main$cluster)]
  marker_deletions <- dplyr::bind_rows(lapply(names(control$markers), function(label) {
    genes <- control$markers[[label]]
    if (length(genes) == 1L) return(NULL)
    dplyr::bind_rows(lapply(genes, function(omitted) {
      remaining <- setdiff(genes, omitted)
      if (any(grepl("-$", genes))) {
        scores <- summaries$signed_means[[label]][[match(omitted, genes) + 1L]][, main$key, drop = FALSE]
        observed <- scores[1L, ]
        null <- scores[-1L, , drop = FALSE]
      } else {
        observed <- score_UCell_rank_means(summaries$rank_means[, main$key, drop = FALSE], remaining, control$max_rank)
        null <- score_UCell_controls(summaries$rank_means[, main$key, drop = FALSE], control, remaining, control$max_rank)
      }
      data.frame(cluster = main$cluster, label, omitted,
        excess = unname(observed) - apply(null, 2L, stats::quantile, probs = 0.95, names = FALSE))
    }))
  }))
  cell_deletions <- dplyr::bind_rows(lapply(seq_len(summaries$n_blocks), function(block) {
    omitted_groups <- groups[groups$kind == "block" & groups$subgroup == as.character(block), ]
    omitted <- match(main$cluster, omitted_groups$cluster)
    sizes <- omitted_groups$cells[omitted]
    sizes[is.na(sizes)] <- 0L
    remaining <- main$cells - sizes
    assessed <- which(sizes > 0L & remaining > 0L)
    if (!length(assessed)) return(NULL)
    leave_out <- function(values) {
      means <- sweep(values[, main$key[assessed], drop = FALSE], 2L, main$cells[assessed], `*`) -
        sweep(values[, omitted_groups$key[omitted[assessed]], drop = FALSE], 2L, sizes[assessed], `*`)
      means <- sweep(means, 2L, remaining[assessed], `/`)
      means[] <- pmax(0, pmin(1, means))
      means
    }
    result <- score_UCell_group_evidence(leave_out(summaries$rank_means), control, leave_out(summaries$detection),
      lapply(summaries$signed_means, function(variants) lapply(variants, leave_out)))
    result$cluster <- main$cluster[match(result$cluster, main$key)]
    result$block <- block
    result
  }))
  if (!nrow(marker_deletions)) marker_deletions <- data.frame(cluster = character(), label = character(), excess = numeric())
  if (!nrow(cell_deletions)) cell_deletions <- data.frame(cluster = character(), block = integer())
  wells <- groups[groups$kind == "well" & groups$cells >= 25L, ]
  well_evidence <- if (nrow(wells)) score_groups(wells$key) else data.frame()
  list(evidence = evidence, marker_deletions = marker_deletions, cell_deletions = cell_deletions,
    well_evidence = well_evidence, wells = wells, cell_scores = summaries$cell_scores,
    settings = list(max_rank = control$max_rank, background_quantile = 0.95,
      diagnostic_marker_detection = 0.1, n_blocks = summaries$n_blocks, diagnostic_min_well_cells = 25L))
}

prepare_cluster_UCell_evidence <- function(counts_matrix, metadata_tibble, control,
                                           cluster_column, include_cell_scores = FALSE, workers = 2L) {
  summaries <- summarize_cluster_UCell_counts(counts_matrix, metadata_tibble, control,
    cluster_column, workers = workers, include_cell_scores = include_cell_scores)
  score_cluster_UCell_summaries(summaries, control)
}

#' Apply one threshold; perturbation and GEM-well agreement never veto assignments.
evaluate_cluster_UCell_evidence <- function(scored, min_advantage = 0.05) {
  decisions <- assign_UCell_cluster_evidence(scored$evidence, min_advantage)
  decisions$cells <- scored$evidence$cells[match(decisions$cluster, scored$evidence$cluster)]
  decisions$marker_stability <- decisions$cell_stability <- NA_real_
  decisions$cell_stability_replicates <- 0L
  decisions$sample_agreement <- NA_real_
  decisions$sample_groups_assessed <- 0L
  resampled <- dplyr::bind_rows(lapply(unique(scored$cell_deletions$block), function(block) {
    result <- assign_UCell_cluster_evidence(scored$cell_deletions[scored$cell_deletions$block == block, ], min_advantage)
    result$block <- block
    result
  }))
  if (!nrow(resampled)) resampled <- data.frame(cluster = character())
  sample_decisions <- data.frame()
  if (nrow(scored$well_evidence)) {
    sample_decisions <- assign_UCell_cluster_evidence(scored$well_evidence, min_advantage)
    index <- match(sample_decisions$cluster, scored$wells$key)
    sample_decisions$cluster <- scored$wells$cluster[index]
    sample_decisions$GEM_well_ID <- scored$wells$subgroup[index]
    sample_decisions$cells <- scored$wells$cells[index]
  }
  for (index in seq_len(nrow(decisions))) {
    cluster <- decisions$cluster[index]
    candidate <- decisions$candidate[index]
    omitted <- scored$marker_deletions[scored$marker_deletions$cluster == cluster &
      scored$marker_deletions$label == candidate, ]
    other <- scored$evidence$excess[scored$evidence$cluster == cluster & scored$evidence$label != candidate]
    if (nrow(omitted)) {
      advantage <- omitted$excess - max(0, other)
      decisions$marker_stability[index] <- mean(advantage >= min_advantage & advantage > 0)
    }
    cells <- resampled[resampled$cluster == cluster, , drop = FALSE]
    decisions$cell_stability_replicates[index] <- nrow(cells)
    if (nrow(cells) >= 2L) decisions$cell_stability[index] <- mean(
      cells$status == "Assigned" & cells$candidate == candidate)
    wells <- sample_decisions[sample_decisions$cluster == cluster, ]
    decisions$sample_groups_assessed[index] <- nrow(wells)
    if (nrow(wells) >= 2L) decisions$sample_agreement[index] <- mean(
      wells$status == "Assigned" & wells$candidate == candidate)
  }
  list(decisions = decisions, evidence = scored$evidence, sample_decisions = sample_decisions,
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
    "GEM_well_agreement.tsv", "control_gene_matching.tsv", "settings.rds", "method.txt"))
  unlink(previous_files[file.exists(previous_files)])
  readr::write_tsv(annotation$decisions, file.path(output_dir, "clusters.tsv"))
  readr::write_tsv(annotation$evidence, file.path(output_dir, "marker_evidence.tsv"))
  if (nrow(annotation$sample_decisions)) {
    readr::write_tsv(annotation$sample_decisions, file.path(output_dir, "GEM_well_agreement.tsv"))
  }
  readr::write_tsv(control$match_diagnostics, file.path(output_dir, "control_gene_matching.tsv"))
  saveRDS(list(settings = annotation$settings, markers = control$markers,
    seed = control$seed, n_controls = control$n_controls, neighbours = control$neighbours,
    reference_barcodes = control$reference_barcodes), file.path(output_dir, "settings.rds"))
  writeLines(c("Assignment uses the best adjusted UCell score's advantage over background and all other labels.",
    paste("Minimum advantage:", annotation$settings$min_advantage),
    "Adjusted score = observed mean UCell minus the label-specific control 95th percentile.",
    "Marker detection, marker deletion, cell stability, control tails and GEM-well agreement are diagnostic only.",
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
