#' Make peak gene correlation gene TSS tibble
#'
#' Match Ensembl gene annotations to GEX matrix features and extract TSS positions.
#'
#' @param reference_Ensembl_annotations_GRanges_list List containing a `genes`
#'   GRanges with Ensembl gene metadata such as `gene_id`, `gene_name`, and
#'   `gene_biotype`.
#' @param GEX_counts_matrix Gene-by-cell count matrix; row names are gene IDs/names and column names are cell barcodes.
#' @return A tibble with one row per gene present in the GEX matrix, including
#'   Ensembl ID, display name, matched matrix feature, chromosome, strand, and TSS.
#' @keywords internal

make_peak_gene_correlation_gene_TSS_tibble <- function(
  reference_Ensembl_annotations_GRanges_list,
  GEX_counts_matrix
) {
  genes <- reference_Ensembl_annotations_GRanges_list$genes
  gex_feature_names <- rownames(GEX_counts_matrix)

  gene_tibble <- GenomicRanges::as.data.frame(genes) |>
    tibble::as_tibble() |>
    dplyr::mutate(
      gene_id = as.character(.data$gene_id),
      gene_name = as.character(.data$gene_name),
      gene_matrix_feature = dplyr::case_when(
        .data$gene_id %in% gex_feature_names ~ .data$gene_id,
        .data$gene_name %in% gex_feature_names ~ .data$gene_name,
        TRUE ~ NA_character_
      ),
      chr = as.character(.data$seqnames),
      TSS = dplyr::if_else(as.character(.data$strand) == "-", .data$end, .data$start),
      strand = as.character(.data$strand),
      gene_biotype = as.character(.data$gene_biotype)
    ) |>
    dplyr::filter(!is.na(.data$gene_matrix_feature)) |>
    dplyr::distinct(.data$gene_id, .keep_all = TRUE) |>
    dplyr::select(
      "gene_id",
      "gene_name",
      "gene_matrix_feature",
      "chr",
      "TSS",
      "strand",
      "gene_biotype"
    )

  gene_tibble
}

#' Make peak gene correlation candidate pairs
#'
#' Pair peak centers with genes whose TSS lies within a fixed genomic window.
#'
#' @param consensus_peak_GRanges Consensus peak GRanges whose names identify ATAC features in downstream matrices and track plots.
#' @param gene_TSS_tibble Tibble of gene TSS records with chromosome, TSS coordinate, gene ID/name, and promoter window columns.
#' @param max_distance Maximum absolute distance in bases between a peak center
#'   and target-gene TSS.
#' @return A candidate-pair tibble with peak coordinates, target-gene fields,
#'   signed TSS distance, and `isSelfPromoter` promoter-window flag.
#' @keywords internal

make_peak_gene_correlation_candidate_pairs <- function(
  consensus_peak_GRanges,
  gene_TSS_tibble,
  max_distance = 250000L
) {
  peak_tibble <- GenomicRanges::as.data.frame(consensus_peak_GRanges) |>
    tibble::as_tibble() |>
    dplyr::mutate(
      peak = names(consensus_peak_GRanges),
      chr = as.character(.data$seqnames),
      peak_center = as.integer(floor((.data$start + .data$end) / 2))
    ) |>
    dplyr::select("peak", "chr", "start", "end", "peak_center")

  gene_windows <- GenomicRanges::GRanges(
    seqnames = gene_TSS_tibble$chr,
    ranges = IRanges::IRanges(
      start = pmax(1L, gene_TSS_tibble$TSS - max_distance),
      end = gene_TSS_tibble$TSS + max_distance
    )
  )

  peak_centers <- GenomicRanges::GRanges(
    seqnames = peak_tibble$chr,
    ranges = IRanges::IRanges(start = peak_tibble$peak_center, width = 1L)
  )

  hits <- GenomicRanges::findOverlaps(peak_centers, gene_windows, ignore.strand = TRUE)
  if (length(hits) == 0L) {
    return(tibble::tibble(
      peak = character(),
      chr = character(),
      start = integer(),
      end = integer(),
      peak_center = integer(),
      TargetGeneID = character(),
      TargetGene = character(),
      gene_matrix_feature = character(),
      TargetGeneTSS = integer(),
      distance = integer(),
      isSelfPromoter = logical()
    ))
  }

  peak_idx <- S4Vectors::queryHits(hits)
  gene_idx <- S4Vectors::subjectHits(hits)

  candidate_pairs <- dplyr::bind_cols(
    peak_tibble[peak_idx, , drop = FALSE],
    gene_TSS_tibble[gene_idx, , drop = FALSE] |>
      dplyr::transmute(
        TargetGeneID = .data$gene_id,
        TargetGene = .data$gene_name,
        gene_matrix_feature = .data$gene_matrix_feature,
        TargetGeneTSS = .data$TSS,
        target_gene_strand = .data$strand
      )
  ) |>
    dplyr::mutate(
      distance = as.integer(.data$peak_center - .data$TargetGeneTSS),
      promoter_start = dplyr::if_else(
        .data$target_gene_strand == "-",
        .data$TargetGeneTSS - 500L,
        .data$TargetGeneTSS - 1500L
      ),
      promoter_end = dplyr::if_else(
        .data$target_gene_strand == "-",
        .data$TargetGeneTSS + 1500L,
        .data$TargetGeneTSS + 500L
      ),
      isSelfPromoter = .data$end >= .data$promoter_start & .data$start <= .data$promoter_end
    ) |>
    dplyr::select(
      "peak",
      "chr",
      "start",
      "end",
      "peak_center",
      "TargetGeneID",
      "TargetGene",
      "gene_matrix_feature",
      "TargetGeneTSS",
      "distance",
      "isSelfPromoter"
    ) |>
    dplyr::distinct(.data$peak, .data$TargetGeneID, .keep_all = TRUE)

  candidate_pairs
}

make_peak_gene_correlation_chromosome_tibble <- function(candidate_pairs_tibble) {
  candidate_pairs_tibble |>
    dplyr::distinct(.data$chr) |>
    dplyr::arrange(gtools::mixedorder(.data$chr))
}

#' Make peak gene correlation cell groups
#'
#' Select analyzable cell groups with barcodes present in GEX, ATAC, and embeddings.
#'
#' @param metadata_tibble Tibble with one row per cell or pseudobulk sample; must contain the barcode/grouping columns referenced by the helper arguments.
#' @param GEX_counts_matrix Gene-by-cell count matrix; row names are gene IDs/names and column names are cell barcodes.
#' @param ATAC_counts_matrix Peak-by-cell count matrix; row names are peak IDs and column names are cell barcodes.
#' @param embedding_matrix Numeric matrix with cells/barcodes in rows and embedding dimensions in columns; row names are carried into downstream coordinates.
#' @param cell_group_col Metadata column whose non-missing values define cell groups for branch construction or plotting.
#' @param min_cells Minimum cells required for a cell group to be retained in branch construction.
#' @return A tibble with one retained row per cell group, a list-column of
#'   member barcodes, and `n_cells`.
#' @keywords internal

make_peak_gene_correlation_cell_groups <- function(
  metadata_tibble,
  GEX_counts_matrix,
  ATAC_counts_matrix,
  embedding_matrix,
  cell_group_col = "PCA_harmony_SNN_cluster_cell_type",
  min_cells = 200L
) {
  available_barcodes <- Reduce(
    base::intersect,
    list(
      colnames(GEX_counts_matrix),
      colnames(ATAC_counts_matrix),
      rownames(embedding_matrix)
    )
  )

  metadata_tibble |>
    dplyr::distinct(.data$barcode_w_prefix, .keep_all = TRUE) |>
    dplyr::filter(.data$barcode_w_prefix %in% available_barcodes) |>
    dplyr::mutate(cell_group = as.character(.data[[cell_group_col]])) |>
    dplyr::filter(!is.na(.data$cell_group), .data$cell_group != "") |>
    dplyr::summarise(
      barcodes = list(.data$barcode_w_prefix),
      n_cells = dplyr::n(),
      .by = "cell_group"
    ) |>
    dplyr::filter(.data$n_cells >= min_cells) |>
    dplyr::arrange(.data$cell_group)
}

#' Make peak gene correlation cell group diagnostics
#'
#' Report cell groups that fail the minimum-cell threshold before correlation.
#'
#' @param metadata_tibble Tibble with one row per cell or pseudobulk sample; must contain the barcode/grouping columns referenced by the helper arguments.
#' @param GEX_counts_matrix Gene-by-cell count matrix; row names are gene IDs/names and column names are cell barcodes.
#' @param ATAC_counts_matrix Peak-by-cell count matrix; row names are peak IDs and column names are cell barcodes.
#' @param embedding_matrix Numeric matrix with cells/barcodes in rows and embedding dimensions in columns; row names are carried into downstream coordinates.
#' @param cell_group_col Metadata column whose non-missing values define cell groups for branch construction or plotting.
#' @param min_cells Minimum cells required for a cell group to be retained in branch construction.
#' @return A diagnostics tibble for too-small cell groups with `skipped_reason`
#'   set to `too_few_cells`.
#' @keywords internal

make_peak_gene_correlation_cell_group_diagnostics <- function(
  metadata_tibble,
  GEX_counts_matrix,
  ATAC_counts_matrix,
  embedding_matrix,
  cell_group_col = "PCA_harmony_SNN_cluster_cell_type",
  min_cells = 200L
) {
  available_barcodes <- Reduce(
    base::intersect,
    list(
      colnames(GEX_counts_matrix),
      colnames(ATAC_counts_matrix),
      rownames(embedding_matrix)
    )
  )

  metadata_tibble |>
    dplyr::distinct(.data$barcode_w_prefix, .keep_all = TRUE) |>
    dplyr::filter(.data$barcode_w_prefix %in% available_barcodes) |>
    dplyr::mutate(cell_group = as.character(.data[[cell_group_col]])) |>
    dplyr::filter(!is.na(.data$cell_group), .data$cell_group != "") |>
    dplyr::summarise(n_cells = dplyr::n(), .by = "cell_group") |>
    dplyr::filter(.data$n_cells < min_cells) |>
    dplyr::transmute(
      cell_group = .data$cell_group,
      chr = NA_character_,
      n_cells = .data$n_cells,
      n_aggregates = NA_integer_,
      n_candidate_pairs = NA_integer_,
      n_detected_genes = NA_integer_,
      n_accessible_peaks = NA_integer_,
      skipped_reason = "too_few_cells"
    )
}

#' Make peak gene correlation KNN aggregates
#'
#' Build overlapping KNN pseudobulk aggregates for one cell group.
#'
#' @param cell_group_tibble One-row tibble from `make_peak_gene_correlation_cell_groups()`
#'   containing `cell_group` and list-column `barcodes`.
#' @param embedding_matrix Numeric matrix with cells/barcodes in rows and embedding dimensions in columns; row names are carried into downstream coordinates.
#' @param dims Integer dimension indices to use; combined with `dim_prefix` to select columns such as `PCA_1` or `LSI_2`.
#' @param k Number of nearest neighbors to use for KNN/SNN construction.
#' @param seed_attempts Maximum number of candidate seed cells to sample for KNN
#'   aggregate construction.
#' @param overlap_cutoff Maximum allowed overlap fraction with already accepted
#'   aggregates; higher values retain more redundant aggregates.
#' @param seed Random seed passed to stochastic clustering, sampling, or embedding code for reproducibility.
#' @return A tibble with one accepted aggregate per row, including `cell_group`,
#'   `aggregate_id`, list-column `barcodes`, and `n_cells`.
#' @keywords internal

make_peak_gene_correlation_knn_aggregates <- function(
  cell_group_tibble,
  embedding_matrix,
  dims = 1:30,
  k = 100L,
  seed_attempts = 500L,
  overlap_cutoff = 0.8,
  seed = 1L
) {
  cell_group <- cell_group_tibble$cell_group[[1]]
  barcodes <- cell_group_tibble$barcodes[[1]]
  barcodes <- intersect(barcodes, rownames(embedding_matrix))

  dim_cols <- paste0("PCA_", dims)
  if (all(dim_cols %in% colnames(embedding_matrix))) {
    available_dims <- dim_cols
  } else {
    available_dims <- seq_len(min(max(dims), ncol(embedding_matrix)))
  }

  group_embedding <- embedding_matrix[barcodes, available_dims, drop = FALSE]
  aggregate_size <- min(as.integer(k), nrow(group_embedding))
  set.seed(seed)
  seed_indices <- sample(seq_len(nrow(group_embedding)), min(seed_attempts, nrow(group_embedding)))
  knn <- FNN::get.knnx(
    data = group_embedding,
    query = group_embedding[seed_indices, , drop = FALSE],
    k = aggregate_size,
    algorithm = "kd_tree"
  )
  accepted_barcodes <- list()

  for (seed_position in seq_along(seed_indices)) {
    seed_idx <- seed_indices[[seed_position]]
    neighbor_indices <- unique(c(seed_idx, knn$nn.index[seed_position, ]))
    neighbor_indices <- neighbor_indices[neighbor_indices >= 1L & neighbor_indices <= nrow(group_embedding)]
    candidate_barcodes <- rownames(group_embedding)[utils::head(neighbor_indices, aggregate_size)]

    overlaps <- if (length(accepted_barcodes) == 0L) {
      numeric()
    } else {
      vapply(
        accepted_barcodes,
        \(accepted) length(intersect(candidate_barcodes, accepted)) / aggregate_size,
        numeric(1)
      )
    }

    if (length(overlaps) == 0L || max(overlaps) <= overlap_cutoff) {
      accepted_barcodes[[length(accepted_barcodes) + 1L]] <- candidate_barcodes
    }
  }

  tibble::tibble(
    cell_group = cell_group,
    aggregate_id = paste0(make.names(cell_group), "_knn_", seq_along(accepted_barcodes)),
    barcodes = accepted_barcodes,
    n_cells = lengths(accepted_barcodes)
  )
}

make_peak_gene_correlation_group_chromosome_tibble <- function(
  knn_aggregates_tibble,
  chromosome_tibble,
  candidate_pairs_tibble
) {
  cell_groups <- knn_aggregates_tibble |>
    dplyr::distinct(.data$cell_group)

  tidyr::crossing(cell_groups, chromosome_tibble) |>
    dplyr::inner_join(
      candidate_pairs_tibble |> dplyr::count(.data$chr, name = "n_candidate_pairs"),
      by = "chr"
    ) |>
    dplyr::mutate(branch_id = paste(.data$cell_group, .data$chr, sep = "__")) |>
    dplyr::group_by(branch_id) |>
    targets::tar_group()
}

make_aggregate_membership_matrix <- function(knn_aggregates_tibble) {
  aggregate_ids <- knn_aggregates_tibble$aggregate_id
  aggregate_barcodes <- knn_aggregates_tibble$barcodes
  cell_barcodes <- unique(unlist(aggregate_barcodes, use.names = FALSE))

  Matrix::sparseMatrix(
    i = match(unlist(aggregate_barcodes, use.names = FALSE), cell_barcodes),
    j = rep(seq_along(aggregate_barcodes), lengths(aggregate_barcodes)),
    x = 1,
    dims = c(length(cell_barcodes), length(aggregate_ids)),
    dimnames = list(cell_barcodes, aggregate_ids)
  )
}

aggregate_BPCells_matrix_by_membership <- function(feature_matrix, features, membership_matrix) {
  features <- intersect(features, rownames(feature_matrix))
  cells <- intersect(rownames(membership_matrix), colnames(feature_matrix))
  membership_matrix <- membership_matrix[cells, , drop = FALSE]
  feature_matrix <- feature_matrix[features, cells, drop = FALSE]

  methods::as(feature_matrix %*% membership_matrix, "dgCMatrix")
}

#' Make peak gene correlation aggregate matrices
#'
#' Aggregate GEX and ATAC count matrices over accepted KNN pseudobulk groups.
#'
#' @param group_chromosome_tibble One-row branch tibble with `cell_group` and
#'   `chr` selecting the correlation branch.
#' @param knn_aggregates_tibble Aggregate membership tibble from
#'   `make_peak_gene_correlation_knn_aggregates()`.
#' @param candidate_pairs_tibble Candidate peak-gene pairs; only pairs on the
#'   branch chromosome are used to select features.
#' @param GEX_counts_matrix Gene-by-cell count matrix; row names are gene IDs/names and column names are cell barcodes.
#' @param ATAC_counts_matrix Peak-by-cell count matrix; row names are peak IDs and column names are cell barcodes.
#' @return A list with branch identifiers, aggregate-level GEX and ATAC count
#'   matrices, and an aggregate depth summary tibble.
#' @keywords internal

make_peak_gene_correlation_aggregate_matrices <- function(
  group_chromosome_tibble,
  knn_aggregates_tibble,
  candidate_pairs_tibble,
  GEX_counts_matrix,
  ATAC_counts_matrix
) {
  cell_group <- group_chromosome_tibble$cell_group[[1]]
  chr <- group_chromosome_tibble$chr[[1]]
  branch_pairs <- candidate_pairs_tibble |>
    dplyr::filter(.data$chr == !!chr)

  branch_aggregates <- knn_aggregates_tibble |>
    dplyr::filter(.data$cell_group == !!cell_group)
  membership_matrix <- make_aggregate_membership_matrix(branch_aggregates)

  GEX_counts <- aggregate_BPCells_matrix_by_membership(
    feature_matrix = GEX_counts_matrix,
    features = unique(branch_pairs$gene_matrix_feature),
    membership_matrix = membership_matrix
  )
  ATAC_counts <- aggregate_BPCells_matrix_by_membership(
    feature_matrix = ATAC_counts_matrix,
    features = unique(branch_pairs$peak),
    membership_matrix = membership_matrix
  )

  aggregate_depth_tibble <- tibble::tibble(
    cell_group = cell_group,
    chr = chr,
    aggregate_id = colnames(membership_matrix),
    n_cells = Matrix::colSums(membership_matrix),
    GEX_depth = Matrix::colSums(GEX_counts),
    ATAC_depth = Matrix::colSums(ATAC_counts)
  )

  list(
    cell_group = cell_group,
    chr = chr,
    GEX_counts = GEX_counts,
    ATAC_counts = ATAC_counts,
    aggregate_depth_tibble = aggregate_depth_tibble
  )
}

normalize_peak_gene_correlation_counts <- function(counts_matrix, scale_factor = 1e6) {
  counts_matrix <- methods::as(counts_matrix, "dgCMatrix")
  depth <- Matrix::colSums(counts_matrix)
  scale <- rep(0, length(depth))
  scale[depth > 0] <- scale_factor / depth[depth > 0]
  normalized <- counts_matrix %*% Matrix::Diagonal(x = scale)
  normalized@x <- log1p(normalized@x)
  Matrix::drop0(normalized)
}

normalize_peak_gene_correlation_aggregate_matrices <- function(
  aggregate_matrices,
  scale_factor = 1e6
) {
  GEX_norm <- normalize_peak_gene_correlation_counts(aggregate_matrices$GEX_counts, scale_factor)
  ATAC_norm <- normalize_peak_gene_correlation_counts(aggregate_matrices$ATAC_counts, scale_factor)

  list(
    cell_group = aggregate_matrices$cell_group,
    chr = aggregate_matrices$chr,
    GEX_norm = GEX_norm,
    ATAC_norm = ATAC_norm,
    gene_detected_frac = Matrix::rowSums(aggregate_matrices$GEX_counts > 0) / ncol(aggregate_matrices$GEX_counts),
    peak_accessible_frac = Matrix::rowSums(aggregate_matrices$ATAC_counts > 0) / ncol(aggregate_matrices$ATAC_counts),
    mean_gene_expression = Matrix::rowMeans(GEX_norm),
    mean_peak_accessibility = Matrix::rowMeans(ATAC_norm),
    aggregate_depth_tibble = aggregate_matrices$aggregate_depth_tibble
  )
}

#' Extract top-link aggregate values from one normalized branch
#'
#' Extract only the normalized GEX and ATAC rows needed for top-link scatter
#' plots from one cell-group/chromosome branch.
#'
#' @param normalized_aggregate_matrices List returned by
#'   `normalize_peak_gene_correlation_aggregate_matrices()`.
#' @param top_links_tibble Top peak-gene links, including scatter-plot names and
#'   matrix feature identifiers.
#' @return A compact tibble of aggregate-level values for links available in
#'   this branch.
#' @keywords internal

extract_peak_gene_correlation_top_link_aggregate_values <- function(
  normalized_aggregate_matrices,
  top_links_tibble
) {
  cell_group <- normalized_aggregate_matrices$cell_group
  chr <- normalized_aggregate_matrices$chr
  GEX_norm <- normalized_aggregate_matrices$GEX_norm
  ATAC_norm <- normalized_aggregate_matrices$ATAC_norm

  branch_links <- top_links_tibble |>
    dplyr::filter(
      .data$chr == !!chr,
      .data$gene_matrix_feature %in% rownames(GEX_norm),
      .data$peak %in% rownames(ATAC_norm)
    )

  if (nrow(branch_links) == 0L) {
    return(tibble::tibble(
      scatter_plot_name = character(),
      primary_cell_group = character(),
      cell_group = character(),
      chr = character(),
      peak = character(),
      TargetGeneID = character(),
      TargetGene = character(),
      correlation = numeric(),
      FDR = numeric(),
      rank_in_cell_group = integer(),
      rank_for_gene = integer(),
      aggregate_id = character(),
      n_cells = numeric(),
      GEX_depth = numeric(),
      ATAC_depth = numeric(),
      gene_expression_logCPM = numeric(),
      peak_accessibility_logCPM = numeric()
    ))
  }

  aggregate_depth_tibble <- normalized_aggregate_matrices$aggregate_depth_tibble |>
    dplyr::select("aggregate_id", "n_cells", "GEX_depth", "ATAC_depth")
  n_aggregates <- nrow(aggregate_depth_tibble)

  purrr::map_dfr(seq_len(nrow(branch_links)), \(index) {
    link_row <- branch_links[index, , drop = FALSE]

    dplyr::bind_cols(
      link_row[rep(1L, n_aggregates), , drop = FALSE] |>
        dplyr::transmute(
          scatter_plot_name = .data$scatter_plot_name,
          primary_cell_group = .data$cell_group,
          cell_group = !!cell_group,
          chr = .data$chr,
          peak = .data$peak,
          TargetGeneID = .data$TargetGeneID,
          TargetGene = .data$TargetGene,
          correlation = .data$correlation,
          FDR = .data$FDR,
          rank_in_cell_group = .data$rank_in_cell_group,
          rank_for_gene = .data$rank_for_gene
        ),
      aggregate_depth_tibble,
      tibble::tibble(
        gene_expression_logCPM = as.numeric(
          GEX_norm[link_row$gene_matrix_feature[[1]], ]
        ),
        peak_accessibility_logCPM = as.numeric(
          ATAC_norm[link_row$peak[[1]], ]
        )
      )
    )
  })
}

#' Empty peak gene correlation results tibble
#'
#' Return a correctly typed empty peak-gene correlation result table.
#'
#' @return A tibble with stable identifiers and derived columns consumed by downstream targets.
#' @keywords internal

empty_peak_gene_correlation_results_tibble <- function() {
  tibble::tibble(
    cell_group = character(),
    peak = character(),
    chr = character(),
    start = integer(),
    end = integer(),
    peak_center = integer(),
    TargetGeneID = character(),
    TargetGene = character(),
    TargetGeneTSS = integer(),
    distance = integer(),
    isSelfPromoter = logical(),
    n_aggregates = integer(),
    mean_gene_expression = numeric(),
    gene_detected_frac = numeric(),
    mean_peak_accessibility = numeric(),
    peak_accessible_frac = numeric(),
    correlation = numeric(),
    nominal_pvalue = numeric()
  )
}

#' Prepare peak gene correlation branch
#'
#' Apply branch-level detection filters and decide whether a branch can be scored.
#'
#' @param normalized_aggregate_matrices List returned by
#'   `normalize_peak_gene_correlation_aggregate_matrices()`.
#' @param candidate_pairs_tibble Candidate peak-gene pairs to filter for the
#'   branch chromosome and detected features.
#' @param min_gene_detection Minimum fraction of aggregates in which a gene must
#'   be detected to retain its candidate pairs.
#' @param min_peak_accessibility Minimum fraction of aggregates in which a peak
#'   must be accessible to retain its candidate pairs.
#' @param min_aggregates Minimum number of accepted aggregates needed to score a
#'   branch.
#' @return A list containing filtered candidate pairs, detected gene/peak sets,
#'   aggregate count, and `skipped_reason` when scoring should be skipped.
#' @keywords internal

prepare_peak_gene_correlation_branch <- function(
  normalized_aggregate_matrices,
  candidate_pairs_tibble,
  min_gene_detection = 0.05,
  min_peak_accessibility = 0.05,
  min_aggregates = 10L
) {
  cell_group <- normalized_aggregate_matrices$cell_group
  chr <- normalized_aggregate_matrices$chr
  n_aggregates <- ncol(normalized_aggregate_matrices$GEX_norm)

  branch_pairs <- candidate_pairs_tibble |>
    dplyr::filter(.data$chr == !!chr)

  detected_genes <- names(normalized_aggregate_matrices$gene_detected_frac)[
    normalized_aggregate_matrices$gene_detected_frac >= min_gene_detection
  ]
  accessible_peaks <- names(normalized_aggregate_matrices$peak_accessible_frac)[
    normalized_aggregate_matrices$peak_accessible_frac >= min_peak_accessibility
  ]

  filtered_pairs <- branch_pairs |>
    dplyr::filter(
      .data$gene_matrix_feature %in% detected_genes,
      .data$peak %in% accessible_peaks
    )

  skipped_reason <- dplyr::case_when(
    n_aggregates < min_aggregates ~ "too_few_accepted_knn_aggregates",
    nrow(branch_pairs) == 0L ~ "no_candidate_pairs",
    length(detected_genes) == 0L ~ "no_detected_genes",
    length(accessible_peaks) == 0L ~ "no_accessible_peaks",
    nrow(filtered_pairs) == 0L ~ "no_candidate_pairs_after_detection_filters",
    TRUE ~ NA_character_
  )

  list(
    cell_group = cell_group,
    chr = chr,
    n_aggregates = n_aggregates,
    candidate_pairs = filtered_pairs,
    detected_genes = detected_genes,
    accessible_peaks = accessible_peaks,
    skipped_reason = skipped_reason
  )
}

#' Diagnose peak gene correlation branch
#'
#' Summarize why a peak-gene correlation branch was scored or skipped.
#'
#' @param normalized_aggregate_matrices List returned by
#'   `normalize_peak_gene_correlation_aggregate_matrices()`.
#' @param candidate_pairs_tibble Candidate peak-gene pairs used to count retained
#'   branch pairs after detection filters.
#' @param min_gene_detection Minimum fraction of aggregates in which a gene must
#'   be detected.
#' @param min_peak_accessibility Minimum fraction of aggregates in which a peak
#'   must be accessible.
#' @param min_aggregates Minimum number of accepted aggregates required before a
#'   branch is considered scoreable.
#' @return One-row diagnostics tibble with branch size, retained feature counts,
#'   candidate-pair count, and optional skipped reason.
#' @keywords internal

diagnose_peak_gene_correlation_branch <- function(
  normalized_aggregate_matrices,
  candidate_pairs_tibble,
  min_gene_detection = 0.05,
  min_peak_accessibility = 0.05,
  min_aggregates = 10L
) {
  branch <- prepare_peak_gene_correlation_branch(
    normalized_aggregate_matrices = normalized_aggregate_matrices,
    candidate_pairs_tibble = candidate_pairs_tibble,
    min_gene_detection = min_gene_detection,
    min_peak_accessibility = min_peak_accessibility,
    min_aggregates = min_aggregates
  )

  tibble::tibble(
    cell_group = branch$cell_group,
    chr = branch$chr,
    n_cells = NA_integer_,
    n_aggregates = branch$n_aggregates,
    n_candidate_pairs = nrow(branch$candidate_pairs),
    n_detected_genes = length(branch$detected_genes),
    n_accessible_peaks = length(branch$accessible_peaks),
    skipped_reason = branch$skipped_reason
  )
}

#' Score peak gene correlations for cell group
#'
#' Compute aggregate-level Pearson correlations for candidate peak-gene pairs.
#'
#' @param normalized_aggregate_matrices List returned by
#'   `normalize_peak_gene_correlation_aggregate_matrices()`, including normalized
#'   GEX/ATAC matrices and feature detection summaries.
#' @param candidate_pairs_tibble Candidate peak-gene pairs to filter and score.
#' @param min_gene_detection Minimum fraction of aggregates in which a gene must
#'   be detected.
#' @param min_peak_accessibility Minimum fraction of aggregates in which a peak
#'   must be accessible.
#' @param min_aggregates Minimum number of accepted aggregates required before
#'   correlations are computed.
#' @return A result tibble with branch identifiers, feature means/detection
#'   fractions, Pearson correlation, and nominal p value for each retained pair.
#' @keywords internal

score_peak_gene_correlations_for_cell_group <- function(
  normalized_aggregate_matrices,
  candidate_pairs_tibble,
  min_gene_detection = 0.05,
  min_peak_accessibility = 0.05,
  min_aggregates = 10L
) {
  branch <- prepare_peak_gene_correlation_branch(
    normalized_aggregate_matrices = normalized_aggregate_matrices,
    candidate_pairs_tibble = candidate_pairs_tibble,
    min_gene_detection = min_gene_detection,
    min_peak_accessibility = min_peak_accessibility,
    min_aggregates = min_aggregates
  )

  if (!is.na(branch$skipped_reason)) {
    return(empty_peak_gene_correlation_results_tibble())
  }

  GEX_norm <- normalized_aggregate_matrices$GEX_norm
  ATAC_norm <- normalized_aggregate_matrices$ATAC_norm
  n_aggregates <- branch$n_aggregates
  df <- n_aggregates - 2L

  scored_pairs <- branch$candidate_pairs |>
    dplyr::group_by(gene_matrix_feature) |>
    dplyr::group_split() |>
    purrr::map_dfr(\(gene_pairs) {
      gene_feature <- gene_pairs$gene_matrix_feature[[1]]
      peak_features <- gene_pairs$peak

      gene_vec <- as.numeric(GEX_norm[gene_feature, ])
      gene_centered <- gene_vec - mean(gene_vec)
      gene_ss <- sum(gene_centered^2)

      peak_matrix <- as.matrix(ATAC_norm[peak_features, , drop = FALSE])
      peak_centered <- sweep(peak_matrix, 1, rowMeans(peak_matrix), "-")
      peak_ss <- rowSums(peak_centered^2)

      denominator <- sqrt(peak_ss * gene_ss)
      correlation <- as.numeric((peak_centered %*% gene_centered) / denominator)
      correlation[denominator == 0 | !is.finite(correlation)] <- NA_real_
      correlation <- pmax(pmin(correlation, 1), -1)

      t_statistic <- correlation * sqrt(df / pmax(1 - correlation^2, .Machine$double.eps))
      nominal_pvalue <- 2 * stats::pt(abs(t_statistic), df = df, lower.tail = FALSE)
      nominal_pvalue[is.na(correlation)] <- NA_real_

      gene_pairs |>
        dplyr::mutate(
          n_aggregates = n_aggregates,
          mean_gene_expression = normalized_aggregate_matrices$mean_gene_expression[.data$gene_matrix_feature],
          gene_detected_frac = normalized_aggregate_matrices$gene_detected_frac[.data$gene_matrix_feature],
          mean_peak_accessibility = normalized_aggregate_matrices$mean_peak_accessibility[.data$peak],
          peak_accessible_frac = normalized_aggregate_matrices$peak_accessible_frac[.data$peak],
          correlation = correlation,
          nominal_pvalue = nominal_pvalue
        )
    }) |>
    dplyr::mutate(cell_group = branch$cell_group, .before = 1) |>
    dplyr::select(-dplyr::all_of("gene_matrix_feature"))

  scored_pairs
}

#' Finalize peak gene correlation results
#'
#' Add FDR, ranks, and aggregation label to scored peak-gene correlations.
#'
#' @param results_tibble Raw scored-pair tibble from branch scoring helpers.
#' @param aggregation Aggregation label inserted into the first output column.
#' @return A standardized results tibble with BH FDR per cell group plus
#'   within-cell-group and within-gene ranks. Empty input returns the same schema.
#' @keywords internal

finalize_peak_gene_correlation_results <- function(results_tibble, aggregation) {
  aggregation_label <- aggregation

  if (nrow(results_tibble) == 0L) {
    return(empty_peak_gene_correlation_results_tibble() |>
      dplyr::mutate(
        aggregation = character(),
        FDR = numeric(),
        rank_in_cell_group = integer(),
        rank_for_gene = integer(),
        .before = 1
      ))
  }

  output_cols <- c(
    "aggregation",
    "cell_group",
    "peak",
    "chr",
    "start",
    "end",
    "peak_center",
    "TargetGeneID",
    "TargetGene",
    "TargetGeneTSS",
    "distance",
    "isSelfPromoter",
    "n_aggregates",
    "mean_gene_expression",
    "gene_detected_frac",
    "mean_peak_accessibility",
    "peak_accessible_frac",
    "correlation",
    "nominal_pvalue",
    "FDR",
    "rank_in_cell_group",
    "rank_for_gene"
  )

  data.table::setDT(results_tibble)
  results_tibble[, FDR := stats::p.adjust(nominal_pvalue, method = "BH"), by = cell_group]
  data.table::setorderv(
    results_tibble,
    c("cell_group", "correlation"),
    c(1L, -1L),
    na.last = TRUE
  )
  results_tibble[, rank_in_cell_group := seq_len(.N), by = cell_group]
  data.table::setorderv(
    results_tibble,
    c("cell_group", "TargetGeneID", "correlation"),
    c(1L, 1L, -1L),
    na.last = TRUE
  )
  results_tibble[, rank_for_gene := seq_len(.N), by = .(cell_group, TargetGeneID)]
  results_tibble[, aggregation := aggregation_label]
  extra_cols <- setdiff(names(results_tibble), output_cols)
  if (length(extra_cols) > 0L) {
    results_tibble[, (extra_cols) := NULL]
  }
  data.table::setcolorder(results_tibble, output_cols)
  tibble::as_tibble(results_tibble)
}

make_peak_gene_correlation_links <- function(results_tibble) {
  results_tibble |>
    dplyr::filter(
      .data$correlation > 0,
      .data$FDR < 0.05,
      !.data$isSelfPromoter
    ) |>
    dplyr::arrange(.data$cell_group, .data$rank_in_cell_group)
}

#' Summarize peak-gene correlations for a histogram
#'
#' @param results_tibble Peak-gene correlation results.
#' @param bin_width Correlation-bin width.
#' @return A compact tibble of pair counts by cell group and correlation bin.
#' @keywords internal

summarize_peak_gene_correlation_histogram <- function(
  results_tibble,
  bin_width = 0.025
) {
  results_tibble |>
    dplyr::select("cell_group", "correlation") |>
    dplyr::filter(!is.na(.data$correlation)) |>
    dplyr::mutate(
      correlation_bin = pmin(
        1 - bin_width,
        pmax(
          -1,
          floor((.data$correlation + 1) / bin_width) * bin_width - 1
        )
      ),
      correlation_mid = .data$correlation_bin + bin_width / 2
    ) |>
    dplyr::count(.data$cell_group, .data$correlation_mid, name = "n_pairs")
}

#' Plot a peak-gene correlation histogram
#'
#' @param plot_tibble Compact output from
#'   `summarize_peak_gene_correlation_histogram()`.
#' @param bin_width Correlation-bin width.
#' @return A ggplot ready for saving.
#' @keywords internal

plot_peak_gene_correlation_histogram <- function(
  plot_tibble,
  bin_width = 0.025
) {
  ggplot2::ggplot(
    plot_tibble,
    ggplot2::aes(x = .data$correlation_mid, y = .data$n_pairs)
  ) +
    ggplot2::geom_col(width = bin_width) +
    ggplot2::facet_wrap(~cell_group, scales = "free_y") +
    ggplot2::labs(
      x = "Pearson correlation",
      y = "Peak-gene pairs"
    )
}

#' Summarize peak-gene correlation support counts
#'
#' @param results_tibble Peak-gene correlation results.
#' @return A compact long-form tibble of support counts by cell group.
#' @keywords internal

summarize_peak_gene_correlation_support_counts <- function(results_tibble) {
  results_tibble |>
    dplyr::select(
      "cell_group",
      "FDR",
      "correlation",
      "isSelfPromoter"
    ) |>
    dplyr::summarise(
      tested_pairs = dplyr::n(),
      FDR_significant_pairs = sum(.data$FDR < 0.05, na.rm = TRUE),
      positive_non_promoter_links = sum(
        .data$correlation > 0 &
          .data$FDR < 0.05 &
          !.data$isSelfPromoter,
        na.rm = TRUE
      ),
      .by = "cell_group"
    ) |>
    tidyr::pivot_longer(
      cols = -"cell_group",
      names_to = "metric",
      values_to = "n"
    ) |>
    dplyr::mutate(n_for_plot = pmax(.data$n, 1))
}

#' Plot peak-gene correlation support counts
#'
#' @param plot_tibble Compact output from
#'   `summarize_peak_gene_correlation_support_counts()`.
#' @return A ggplot ready for saving.
#' @keywords internal

plot_peak_gene_correlation_support_counts <- function(plot_tibble) {
  ggplot2::ggplot(
    plot_tibble,
    ggplot2::aes(
      x = .data$n_for_plot,
      y = stats::reorder(.data$cell_group, .data$n_for_plot),
      fill = .data$metric
    )
  ) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::scale_x_log10() +
    ggplot2::labs(
      x = "Count, log10 scale",
      y = "Cell group",
      fill = NULL
    )
}

#' Summarize peak-gene correlations by TSS distance
#'
#' @param results_tibble Peak-gene correlation results.
#' @return A compact tibble of correlation summaries by cell group and distance
#'   bin.
#' @keywords internal

summarize_peak_gene_correlation_by_distance <- function(results_tibble) {
  results_tibble |>
    dplyr::select(
      "cell_group",
      "correlation",
      "isSelfPromoter",
      "distance",
      "FDR"
    ) |>
    dplyr::filter(!is.na(.data$correlation), !.data$isSelfPromoter) |>
    dplyr::mutate(
      abs_distance_bin = pmin(
        250000,
        floor(abs(.data$distance) / 5000) * 5000
      )
    ) |>
    dplyr::summarise(
      n_pairs = dplyr::n(),
      median_correlation = stats::median(.data$correlation, na.rm = TRUE),
      significant_fraction = mean(.data$FDR < 0.05, na.rm = TRUE),
      .by = c("cell_group", "abs_distance_bin")
    )
}

#' Plot peak-gene correlations by TSS distance
#'
#' @param plot_tibble Compact output from
#'   `summarize_peak_gene_correlation_by_distance()`.
#' @return A ggplot ready for saving.
#' @keywords internal

plot_peak_gene_correlation_by_distance <- function(plot_tibble) {
  ggplot2::ggplot(
    plot_tibble,
    ggplot2::aes(
      x = .data$abs_distance_bin / 1000,
      y = .data$median_correlation
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = 2,
      color = "grey70"
    ) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~cell_group) +
    ggplot2::labs(
      x = "Absolute TSS distance, kb",
      y = "Median correlation"
    )
}
