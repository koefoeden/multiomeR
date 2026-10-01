# Production-object parity and resource comparisons of the pipeline's algorithm
# reimplementations against their references (Supplementary Table S4). Fixtures
# are built from regular pipeline targets; each implementation runs in a fresh R
# subprocess while its target samples the resident memory of that process tree,
# since long-lived workers retain memory from earlier targets.

ALGORITHM_PARITY_PACKAGES <- c("BPCells", "UCell", "Seurat", "SeuratObject", "scDblFinder", "BiocParallel")

#' Installed versions of the packages the parity comparisons run
#'
#' Evaluated on every run, so comparisons rerun exactly when the lockfile changes
#' one of these packages.
get_algorithm_parity_package_versions <- function() {
  vapply(ALGORITHM_PARITY_PACKAGES, \(package) as.character(utils::packageVersion(package)), character(1))
}

get_UCell_parity_fixture <- function(counts, metadata, markers) {
  barcodes <- intersect(metadata$barcode_w_prefix, colnames(counts))
  stopifnot("No GEX metadata barcodes in the count matrix" = length(barcodes) > 0L)
  list(counts = counts[, barcodes, drop = FALSE], markers = markers,
    n_cells = length(barcodes), n_features = nrow(counts))
}

get_WNN_parity_fixture <- function(embeddings) {
  list(embeddings = embeddings, n_cells = nrow(embeddings[[1]]),
    n_features = sum(vapply(embeddings, ncol, integer(1))))
}

#' ATAC scDblFinder fixture: one GEM well's peak counts with the global LSI feature groups and clusters
get_scDblFinder_ATAC_parity_fixture <- function(branches, counts, feature_groups, loadings, GEM_well_ID) {
  branch <- branches[branches$GEM_well_ID == GEM_well_ID, , drop = FALSE]
  stopifnot("Expected one scDblFinder ATAC branch for the fixture well" = nrow(branch) == 1L)
  barcodes <- intersect(branch$barcode_vec[[1]], colnames(counts))
  stopifnot("No branch barcodes in the peak matrix" = length(barcodes) > 0L)
  list(counts = counts[, barcodes, drop = FALSE], loadings = loadings, feature_groups = feature_groups,
    branch = branch, n_cells = length(barcodes), n_features = nrow(counts),
    n_groups = length(unique(feature_groups)))
}

#' Write the pre-doublet ATAC fixture: insertion-counted consensus peaks for the
#' well's nuclei after ordinary QC, before doublet filters
write_scDblFinder_ATAC_predoublet_parity_matrix <- function(metadata, exclusions, consensus_peaks, fragments, dir) {
  doublet_filter <- grepl("amulet|vireo|doublet", names(exclusions), ignore.case = TRUE)
  QC_excluded <- unique(unlist(exclusions[!doublet_filter], use.names = FALSE))
  barcodes <- metadata$barcode_w_prefix[!metadata$barcode_w_prefix %in% QC_excluded]
  stopifnot("Missing or duplicated pre-doublet barcodes" = length(barcodes) > 0L && !anyDuplicated(barcodes))
  peaks <- GenomicRanges::as.data.frame(consensus_peaks) |>
    dplyr::transmute(chr = as.character(.data$seqnames), start = .data$start, end = .data$end)
  fragments <- BPCells::select_cells(fragments, barcodes)
  peak_order <- BPCells::order_ranges(peaks, BPCells::chrNames(fragments))
  counts <- BPCells::peak_matrix(fragments, peaks, mode = "insertions")
  rownames(counts) <- names(consensus_peaks)[peak_order]
  BPCells::write_matrix_dir(counts, dir, overwrite = TRUE)
  dir
}

get_scDblFinder_ATAC_predoublet_parity_fixture <- function(counts, feature_groups, GEM_well_ID) {
  feature_groups <- feature_groups[rownames(counts)]
  stopifnot("The global LSI feature groups do not cover the fixture peaks" = !anyNA(feature_groups))
  branch <- list(GEM_well_ID = GEM_well_ID, barcode_vec = list(colnames(counts)),
    cluster_vec = list(NULL), n_cells = ncol(counts))
  list(counts = counts, feature_groups = feature_groups, branch = branch,
    n_cells = ncol(counts), n_features = nrow(counts), n_groups = 50L)
}

get_amulet_parity_fixture <- function(fragments, barcode_file, summary_file, GEM_well_ID) {
  barcodes <- readLines(barcode_file)
  list(fragments = fragments, barcodes = barcodes, prefix = paste0(GEM_well_ID, "_"),
    fragment_file = file.path(dirname(summary_file), "atac_fragments.tsv.gz"), n_cells = length(barcodes))
}

#' SCAVENGE fixture: the final WNN SNN graph with one motif family's chromVAR
#' Z-scores as seed signal
get_SCAVENGE_parity_fixture <- function(WNN_results, chromVAR_results, motif_family) {
  graph <- get_SNN_matrix_from_WNN_results(WNN_results)
  z_score_matrix <- chromVAR_results$chromVAR_z_scores
  stopifnot("The chromVAR Z-scores lack the fixture motif family" = motif_family %in% rownames(z_score_matrix))
  z_score <- stats::setNames(as.numeric(z_score_matrix[motif_family, ]), colnames(z_score_matrix))
  cells <- intersect(rownames(graph), names(z_score))
  list(graph = graph[cells, cells, drop = FALSE], z_score = z_score[cells],
    motif_family = motif_family, n_cells = length(cells), n_features = length(graph@x))
}

# Runners execute in a fresh subprocess with the project runtime loaded; each
# returns the implementation's result and version.

run_UCell_parity_native <- function(fixture, cores) {
  list(result = calculate_BPCells_UCell_scores_from_matrix(counts_matrix = fixture$counts,
    features = fixture$markers, max_rank = 1500, chunk_size = 1000, workers = cores, missing_genes = "impute"),
    version = paste("BPCells", utils::packageVersion("BPCells")))
}

run_UCell_parity_reference <- function(fixture, cores) {
  scores <- UCell::ScoreSignatures_UCell(matrix = methods::as(fixture$counts, "dgCMatrix"),
    features = fixture$markers, maxRank = 1500, w_neg = 1, name = "", chunk.size = 1000,
    missing_genes = "impute", BPPARAM = BiocParallel::MulticoreParam(workers = cores), ncores = cores,
    ties.method = "average")
  list(result = as.data.frame(scores, check.names = FALSE), version = paste("UCell", utils::packageVersion("UCell")))
}

run_WNN_parity_native <- function(fixture, cores, native_source_file) {
  result <- weighted_nearest_neighbors_BPCells(embeddings_list = fixture$embeddings, k = 30, candidate_k = 200,
    threads = cores, ef = 500, native_source_file = native_source_file)
  list(result = list(weights = result$modality_weights, nn_idx = result$nn_idx),
    version = paste("BPCells", utils::packageVersion("BPCells")))
}

run_WNN_parity_reference <- function(fixture, cores) {
  cells <- rownames(fixture$embeddings[[1]])
  counts <- Matrix::sparseMatrix(i = rep(1L, length(cells)), j = seq_along(cells), x = 1,
    dims = c(1L, length(cells)), dimnames = list("dummy", cells))
  object <- SeuratObject::CreateSeuratObject(counts = counts)
  reductions <- paste0(tolower(names(fixture$embeddings)), "_benchmark")
  for (i in seq_along(fixture$embeddings)) {
    object[[reductions[[i]]]] <- SeuratObject::CreateDimReducObject(embeddings = fixture$embeddings[[i]],
      key = paste0(toupper(names(fixture$embeddings)[[i]]), "_"), assay = "RNA")
  }
  set.seed(847)
  object <- Seurat::FindMultiModalNeighbors(object = object, reduction.list = reductions,
    dims.list = lapply(fixture$embeddings, \(embedding) seq_len(ncol(embedding))), k.nn = 30, knn.range = 200,
    l2.norm = TRUE, modality.weight.name = paste0(names(fixture$embeddings), ".weight"), verbose = FALSE)
  list(result = list(weights = object[[]][, paste0(names(fixture$embeddings), ".weight"), drop = FALSE],
      nn_idx = object@neighbors$weighted.nn@nn.idx),
    version = paste("Seurat", utils::packageVersion("Seurat")))
}

run_scDblFinder_ATAC_parity <- function(feature_matrix, fixture, aggregate_features) {
  set.seed(713)
  run_scDblFinder_BPCells_GEM_well(feature_matrix = feature_matrix, scDblFinder_GEM_well_tibble = fixture$branch,
    output_suffix = "ATAC", dbr.sd = 1.0, aggregateFeatures = aggregate_features, nfeatures = fixture$n_groups,
    processing = "normFeatures")
}

run_scDblFinder_ATAC_parity_native <- function(fixture, cores) {
  if (!is.null(fixture$loadings)) {
    groups <- get_feature_groups_from_LSI_loadings(LSI_loadings_tibble = fixture$loadings, dims = 2:20,
      n_groups = fixture$n_groups, seed = 1)
    stopifnot("The global LSI grouping no longer reproduces the stored feature groups" =
      identical(groups, fixture$feature_groups))
  }
  feature_matrix <- aggregate_BPCells_rows_by_group(feature_matrix = fixture$counts,
    feature_groups = fixture$feature_groups, threads = cores)
  list(result = run_scDblFinder_ATAC_parity(feature_matrix, fixture, aggregate_features = FALSE),
    version = paste0("BPCells ", utils::packageVersion("BPCells"), " + scDblFinder ", utils::packageVersion("scDblFinder")))
}

run_scDblFinder_ATAC_parity_reference <- function(fixture, cores) {
  list(result = run_scDblFinder_ATAC_parity(fixture$counts, fixture, aggregate_features = TRUE),
    version = paste("scDblFinder", utils::packageVersion("scDblFinder")))
}

run_amulet_parity_native <- function(fixture, cores, native_source_file) {
  result <- calculate_amulet_metrics_BPCells(fragments = fixture$fragments,
    barcodes = paste0(fixture$prefix, fixture$barcodes), cellranger_end_inclusive = TRUE, verbose = FALSE,
    native_source_file = native_source_file)
  rownames(result) <- substring(rownames(result), nchar(fixture$prefix) + 1L)
  list(result = result, version = paste("BPCells", utils::packageVersion("BPCells")))
}

run_amulet_parity_reference <- function(fixture, cores) {
  result <- scDblFinder::amulet(x = fixture$fragment_file, barcodes = fixture$barcodes, uniqueFrags = TRUE,
    fullInMemory = TRUE, BPPARAM = BiocParallel::MulticoreParam(workers = cores), verbose = FALSE)
  list(result = result, version = paste0("scDblFinder ", utils::packageVersion("scDblFinder"), ", full-memory AMULET"))
}

run_SCAVENGE_parity_native <- function(fixture, cores) {
  record <- list(GWAS_ID = paste0("motif_family_", fixture$motif_family), z_score_vec = fixture$z_score)
  TRS_tibble <- get_SCAVENGE_TRS_tibble(chromVAR_z_score_record = record, NN_graph = fixture$graph,
    restart_prob = 0.05, seed_percent = 0.05)
  list(result = TRS_tibble[, c("barcode_w_prefix", "score")], version = "multiomeR native")
}

run_SCAVENGE_parity_reference <- function(fixture, cores, reference_file) {
  # The pinned SCAVENGE reference functions of the SCAVENGE parity tests.
  source(reference_file, local = TRUE)
  seed_idx <- reference_SCAVENGE_seed_index(fixture$z_score, seed_percent = 0.05)
  propagation <- reference_SCAVENGE_random_walk(graph = fixture$graph, seed_cells = names(fixture$z_score)[seed_idx])
  kept <- names(propagation)[propagation != 0]
  scores <- reference_SCAVENGE_scores(propagation_score = propagation[kept], z_score = fixture$z_score,
    scale_percent = 0.01)
  list(result = data.frame(barcode_w_prefix = kept, score = unname(scores), stringsAsFactors = FALSE),
    version = "SCAVENGE 1.0.2@8ee8b173d965")
}

read_process_tree_rss_kib <- function(root_pid) {
  output <- suppressWarnings(system2("ps", c("-e", "-o", "pid=,ppid=,rss="), stdout = TRUE, stderr = FALSE))
  fields <- strsplit(trimws(output), "[[:space:]]+")
  fields <- do.call(rbind, fields[lengths(fields) == 3L])
  if (is.null(fields)) return(0)
  pid <- as.integer(fields[, 1]); ppid <- as.integer(fields[, 2]); rss <- as.numeric(fields[, 3])
  tree <- as.integer(root_pid)
  repeat {
    expanded <- union(tree, pid[ppid %in% tree])
    if (length(expanded) == length(tree)) break
    tree <- expanded
  }
  sum(rss[pid %in% tree], na.rm = TRUE)
}

#' Run one implementation in a fresh R process and measure it
#'
#' Wall time covers data loading, package setup, representation conversion and
#' native compilation; peak RAM is the maximum summed resident set size of the
#' process and its descendants, sampled every 0.1 seconds.
#' @param package_versions Parity package versions, a dependency that reruns the
#'   case when the lockfile changes one of them.
measure_algorithm_parity_case <- function(algorithm, implementation, runner, fixture, cores, package_versions, ...) {
  # The project .Rprofile loads the pipeline runtime; the runners also use these helpers.
  process <- callr::r_bg(\(runner, fixture, cores, args) {
    grDevices::pdf(NULL)  # Reference packages may plot; keep Rplots.pdf out of the checkout.
    source("pipeline_manuscript_figures/algorithm_parity_helpers.R")
    do.call(runner, c(list(fixture, cores), args))
  }, args = list(runner = runner, fixture = fixture, cores = cores, args = list(...)), supervise = TRUE)
  started <- Sys.time()
  peak_rss_kib <- 0
  repeat {
    peak_rss_kib <- max(peak_rss_kib, read_process_tree_rss_kib(process$get_pid()))
    if (!process$is_alive()) break
    Sys.sleep(0.1)
  }
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  output <- process$get_result()
  list(result = output$result, measurement = tibble::tibble(algorithm = algorithm,
    implementation = implementation, version = output$version, n_cells = fixture$n_cells,
    n_features = fixture$n_features, cores = cores, elapsed_seconds = elapsed,
    peak_process_tree_rss_gb = peak_rss_kib / 1024^2, host = Sys.info()[["nodename"]]),
    package_versions = package_versions)
}

make_parity_row <- function(algorithm, metric, observed, operator, threshold) {
  passed <- switch(operator, "==" = isTRUE(observed == threshold), ">=" = isTRUE(observed >= threshold),
    "<=" = isTRUE(observed <= threshold))
  tibble::tibble(algorithm = algorithm, metric = metric, observed = as.numeric(observed),
    operator = operator, threshold = as.numeric(threshold), passed = passed)
}

compare_UCell_parity <- function(native, reference) {
  native <- as.matrix(native)
  reference <- as.matrix(reference[, colnames(native), drop = FALSE])
  max_delta <- if (identical(dim(native), dim(reference))) max(abs(native - reference)) else Inf
  dplyr::bind_rows(
    make_parity_row("ucell", "identical_values_dimensions_dimnames", identical(native, reference), "==", TRUE),
    make_parity_row("ucell", "maximum_absolute_delta", max_delta, "<=", 0))
}

compare_amulet_parity <- function(native, reference) {
  same_barcodes <- identical(rownames(native), rownames(reference))
  shared <- intersect(rownames(native), rownames(reference))
  native <- native[shared, , drop = FALSE]
  reference <- reference[shared, colnames(native), drop = FALSE]
  dplyr::bind_rows(
    make_parity_row("amulet", "same_retained_barcodes", same_barcodes, "==", TRUE),
    make_parity_row("amulet", "identical_metrics", identical(native, reference), "==", TRUE),
    make_parity_row("amulet", "maximum_absolute_delta", max(abs(as.matrix(native) - as.matrix(reference))), "<=", 0))
}

compare_WNN_parity <- function(native, reference) {
  modalities <- sub("\\.weight$", "", colnames(reference$weights))
  correlations <- vapply(modalities, \(modality) stats::cor(native$weights[[modality]],
    reference$weights[[paste0(modality, ".weight")]], method = "spearman"), numeric(1))
  k <- ncol(reference$nn_idx)
  overlap <- vapply(seq_len(nrow(reference$nn_idx)),
    \(cell) length(intersect(reference$nn_idx[cell, ], native$nn_idx[cell, ])) / k, numeric(1))
  dplyr::bind_rows(
    lapply(modalities, \(modality) make_parity_row("wnn", paste0(modality, "_weight_spearman"),
      correlations[[modality]], ">=", 0.85)),
    make_parity_row("wnn", "neighbor_overlap_mean", mean(overlap), ">=", 0.75),
    make_parity_row("wnn", "neighbor_overlap_q25", stats::quantile(overlap, 0.25, names = FALSE), ">=", 0.65))
}

summarize_scDblFinder_ATAC_parity <- function(native, reference, algorithm) {
  shared <- intersect(native$barcode_w_prefix, reference$barcode_w_prefix)
  native <- native[match(shared, native$barcode_w_prefix), , drop = FALSE]
  reference <- reference[match(shared, reference$barcode_w_prefix), , drop = FALSE]
  native_doublets <- native$barcode_w_prefix[native$scDblFinder.class_ATAC == "doublet"]
  reference_doublets <- reference$barcode_w_prefix[reference$scDblFinder.class_ATAC == "doublet"]
  union_doublets <- union(native_doublets, reference_doublets)
  tibble::tibble(algorithm = algorithm, n_native_cells = nrow(native), n_reference_cells = nrow(reference),
    n_shared_cells = length(shared), n_native_doublets = length(native_doublets),
    n_reference_doublets = length(reference_doublets),
    n_shared_doublets = length(intersect(native_doublets, reference_doublets)),
    score_rank_spearman = stats::cor(native$scDblFinder.score_ATAC, reference$scDblFinder.score_ATAC,
      method = "spearman"),
    doublet_jaccard = if (length(union_doublets)) length(intersect(native_doublets, reference_doublets)) /
      length(union_doublets) else 1,
    overall_class_agreement = mean(native$scDblFinder.class_ATAC == reference$scDblFinder.class_ATAC))
}

compare_scDblFinder_ATAC_parity <- function(native, reference, algorithm) {
  agreement <- summarize_scDblFinder_ATAC_parity(native, reference, algorithm)
  dplyr::bind_rows(
    make_parity_row(algorithm, "shared_cell_fraction", agreement$n_shared_cells / agreement$n_reference_cells, ">=", 1),
    make_parity_row(algorithm, "score_rank_spearman", agreement$score_rank_spearman, ">=", 0.6),
    make_parity_row(algorithm, "doublet_jaccard", agreement$doublet_jaccard, ">=", 0.45))
}

compare_SCAVENGE_parity <- function(native, reference) {
  shared <- intersect(native$barcode_w_prefix, reference$barcode_w_prefix)
  n_reference <- nrow(reference)
  native <- native[match(shared, native$barcode_w_prefix), , drop = FALSE]
  reference <- reference[match(shared, reference$barcode_w_prefix), , drop = FALSE]
  dplyr::bind_rows(
    make_parity_row("scavenge", "shared_cell_fraction", length(shared) / n_reference, ">=", 0.99),
    make_parity_row("scavenge", "maximum_absolute_score_delta", max(abs(native$score - reference$score)), "<=", 1e-12),
    make_parity_row("scavenge", "score_rank_spearman", stats::cor(native$score, reference$score, method = "spearman"),
      ">=", 0.999999))
}

#' Combine measurements and parity checks into one row per algorithm
summarize_algorithm_parity <- function(runs, parity) {
  measurements <- dplyr::bind_rows(lapply(runs, `[[`, "measurement"))
  native <- dplyr::filter(measurements, .data$implementation == "native")
  reference <- dplyr::filter(measurements, .data$implementation == "reference")
  dplyr::inner_join(native, reference, by = c("algorithm", "n_cells", "n_features", "cores"),
    suffix = c("_native", "_reference")) |>
    dplyr::mutate(time_ratio = .data$elapsed_seconds_native / .data$elapsed_seconds_reference,
      RAM_ratio = .data$peak_process_tree_rss_gb_native / .data$peak_process_tree_rss_gb_reference,
      parity_passed = vapply(.data$algorithm, \(a) all(parity$passed[parity$algorithm == a]), logical(1)))
}

format_parity_number <- function(x, digits = 3) formatC(x, format = "f", digits = digits, big.mark = ",")

#' Result-agreement text for one algorithm's parity rows
describe_algorithm_parity <- function(parity, agreement, algorithm) {
  p <- parity[parity$algorithm == algorithm, ]
  value <- \(metric) p$observed[p$metric == metric]
  switch(algorithm,
    ucell = , amulet = if (all(p$passed)) "Identical" else "Differs",
    scavenge = if (all(p$passed)) "Identical within floating-point precision" else "Differs",
    wnn = paste0(paste(sub("_weight_spearman$", "", p$metric[grepl("_weight_spearman$", p$metric)]), collapse = "/"),
      " weight Spearman ", paste(format_parity_number(p$observed[grepl("_weight_spearman$", p$metric)]),
        collapse = " / "), "#linebreak() Mean neighbour overlap ", format_parity_number(value("neighbor_overlap_mean")),
      "#linebreak() First-quartile neighbour overlap ", format_parity_number(value("neighbor_overlap_q25"))),
    with(agreement[agreement$algorithm == algorithm, ], paste0("Score-rank Spearman ",
      format_parity_number(score_rank_spearman), "#linebreak() Doublet-call Jaccard ",
      format_parity_number(doublet_jaccard), "#linebreak() Class agreement ",
      formatC(100 * overall_class_agreement, format = "f", digits = 2), "%")))
}

#' Typst table of Supplementary Table S4, one row per reported algorithm
format_algorithm_parity_table <- function(summary, parity, agreement) {
  rows <- tibble::tribble(
    ~algorithm, ~label, ~fixture,
    "ucell", "UCell-compatible scoring", "{cells} cells#linebreak() {features}-feature GEX matrix",
    "amulet", "AMULET", "{cells} cells#linebreak() Complete fragment input",
    "scdblfinder_atac_predoublet", "ATAC scDblFinder feature aggregation",
      "{cells} pre-doublet cells#linebreak() {features}-peak matrix",
    "wnn", "Weighted nearest neighbours", "{cells} cells#linebreak() {features}-dimensional RNA/ATAC embeddings",
    "scavenge", "SCAVENGE-style propagation", "{cells} cells#linebreak() WNN graph")
  count <- \(x) formatC(x, format = "d", big.mark = ",")
  cells <- vapply(rows$algorithm, \(a) {
    s <- summary[summary$algorithm == a, ]
    reference <- sub("SCAVENGE 1.0.2@8ee8b173d965", "SCAVENGE 1.0.2 (#raw(\"8ee8b173\"))", s$version_reference, fixed = TRUE)
    if (a == "scdblfinder_atac_predoublet") reference <- sub("$", " internal aggregation", reference)
    fixture <- rows$fixture[rows$algorithm == a]
    fixture <- sub("{cells}", count(s$n_cells), fixture, fixed = TRUE)
    fixture <- sub("{features}", count(s$n_features), fixture, fixed = TRUE)
    paste0("  [", rows$label[rows$algorithm == a], "],\n  [", fixture, "],\n  [", reference, "],\n  [",
      describe_algorithm_parity(parity, agreement, a), "],\n  [",
      sprintf("%s s / %s s (%.2f×)", format_parity_number(s$elapsed_seconds_native, 1),
        format_parity_number(s$elapsed_seconds_reference, 1), s$time_ratio), "],\n  [",
      sprintf("%.2f GB / %.2f GB (%.2f×)", s$peak_process_tree_rss_gb_native, s$peak_process_tree_rss_gb_reference,
        s$RAM_ratio), "],")
  }, character(1))
  c("```{=typst}", "#table(", "  columns: (1fr, 1.15fr, 1.05fr, 2.1fr, auto, auto),", "  stroke: none,",
    "  inset: (x: 4pt, y: 5pt),", "  align: (col, row) => if row == 0 {",
    "    if col >= 4 { right + horizon } else { left + horizon }", "  } else if col >= 4 {", "    right + top",
    "  } else {", "    left + top", "  },", "  table.header(", "    [#strong[Implementation]],",
    "    [#strong[Production fixture]],", "    [#strong[Reference]],", "    [#strong[Result agreement]],",
    "    [#strong[Time,#linebreak() reimplementation / reference]],",
    "    [#strong[Peak RAM,#linebreak() reimplementation / reference]],", "  ),", "  table.hline(stroke: 1pt),",
    cells, ")", "```")
}
