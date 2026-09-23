# Parked cell-type annotation perturbation diagnostics.
#
# These functions reproduce the leave-one-marker-out, leave-one-block-out, and
# GEM-well agreement diagnostics that the primary module computed alongside its
# UCell cluster annotation until September 2026. They never changed an
# assignment and are not part of the target graph; run_diagnostics.R applies
# them to stored targets. They reuse the production helpers in
# R/cluster_annotation_helpers.R.

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

#' Aggregate UCell statistics in bounded count chunks, clipping signed scores per cell.
#' Cluster, stratified deletion-block, and GEM-well summaries share one rank pass.
#' Only marker/control genes are retained after ranking all genes in each cell.
summarize_cluster_UCell_diagnostic_counts <- function(counts_matrix, metadata_tibble, control,
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
score_cluster_UCell_diagnostic_summaries <- function(summaries, control) {
  groups <- summaries$groups
  main <- groups[groups$kind == "cluster", ]
  score_groups <- function(keys) score_UCell_group_evidence(
    summaries$rank_means[, keys, drop = FALSE], control, summaries$detection[, keys, drop = FALSE],
    lapply(summaries$signed_means, function(variants) variants[[1L]][, keys, drop = FALSE]))
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
      lapply(summaries$signed_means, function(variants) leave_out(variants[[1L]])))
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

#' Apply one threshold; perturbation and GEM-well agreement never veto assignments.
evaluate_cluster_UCell_diagnostics <- function(scored, min_advantage) {
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
