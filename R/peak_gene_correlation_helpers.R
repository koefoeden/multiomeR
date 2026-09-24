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
      gene_start = as.integer(.data$start),
      gene_end = as.integer(.data$end),
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
      "gene_start",
      "gene_end",
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
      TargetGeneStart = integer(),
      TargetGeneEnd = integer(),
      distance = integer(),
      isSelfPromoter = logical(),
      isTargetGeneBody = logical(),
      link_class = character()
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
        TargetGeneStart = .data$gene_start,
        TargetGeneEnd = .data$gene_end,
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
      isSelfPromoter = .data$end >= .data$promoter_start & .data$start <= .data$promoter_end,
      isTargetGeneBody = .data$end >= .data$TargetGeneStart & .data$start <= .data$TargetGeneEnd,
      link_class = dplyr::case_when(
        .data$isSelfPromoter ~ "self_promoter",
        .data$isTargetGeneBody ~ "target_gene_body",
        abs(.data$distance) <= 10000L ~ "proximal_nonpromoter",
        TRUE ~ "distal"
      )
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
      "TargetGeneStart",
      "TargetGeneEnd",
      "distance",
      "isSelfPromoter",
      "isTargetGeneBody",
      "link_class"
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
  cell_group_col = "WNN_harmony_SNN_cluster_cell_type",
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
  cell_group_col = "WNN_harmony_SNN_cluster_cell_type",
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
      n_donors = NA_integer_,
      n_state_bins = NA_integer_,
      n_aggregates = NA_integer_,
      n_candidate_pairs = NA_integer_,
      n_detected_genes = NA_integer_,
      n_accessible_peaks = NA_integer_,
      skipped_reason = "too_few_cells"
    )
}

#' Build the donor-aware peak gene correlation nuisance design
#'
#' @param aggregate_depth_tibble Donor-state aggregate metadata.
#' @return A full-rank numeric design matrix containing donor and
#'   library-depth terms.
#' @keywords internal

make_peak_gene_correlation_design_matrix <- function(aggregate_depth_tibble) {
  donor_id <- factor(aggregate_depth_tibble$donor_id)
  design <- if (nlevels(donor_id) > 1L) {
    stats::model.matrix(~ donor_id)
  } else {
    matrix(1, nrow(aggregate_depth_tibble), 1L, dimnames = list(NULL, "(Intercept)"))
  }
  for (depth_col in c("GEX_depth", "ATAC_depth")) {
    depth <- log1p(aggregate_depth_tibble[[depth_col]])
    if (all(is.finite(depth)) && isTRUE(stats::sd(depth) > 0)) {
      design <- cbind(design, scale(depth))
      colnames(design)[ncol(design)] <- depth_col
    }
  }
  design_qr <- qr(design)
  design[, design_qr$pivot[seq_len(design_qr$rank)], drop = FALSE]
}

residualize_peak_gene_correlation_matrix <- function(feature_matrix, design) {
  feature_matrix <- as.matrix(feature_matrix)
  t(qr.resid(qr(design), t(feature_matrix)))
}

#' Make donor by ATAC state pseudobulks for peak gene correlation
#'
#' Partition one broad cell group into mutually exclusive ATAC-defined states,
#' then pseudobulk cells within donor and state. The retained table requires
#' repeated state observations within donor. States need not be shared across
#' donors: candidate discovery is allowed even in a single donor.
#'
#' @param cell_group_tibble One-row tibble from
#'   `make_peak_gene_correlation_cell_groups()`.
#' @param metadata_tibble Cell metadata containing barcode, donor, and depth
#'   columns.
#' @param ATAC_embedding_matrix ATAC Harmony/LSI embedding with cells in rows.
#' @param donor_col Metadata column containing biological donor identifiers.
#' @param dims Embedding dimensions to use. Named `LSI_*` columns are preferred.
#' @param min_cells_per_donor_state Minimum cells in a retained donor-state
#'   pseudobulk.
#' @param min_donors Minimum retained donors required for the cell group.
#' @param min_donors_per_state Minimum donors supporting each retained state.
#' @param min_states_per_donor Minimum retained states required per donor.
#' @param max_state_bins Maximum number of ATAC state bins.
#' @param seed Random seed used for deterministic k-means initialization.
#' @return A list with `aggregates` and a one-row `diagnostics` tibble.
#' @keywords internal

make_peak_gene_correlation_donor_state_record <- function(
  cell_group_tibble,
  metadata_tibble,
  ATAC_embedding_matrix,
  donor_col = "donor_id",
  dims = 2:20,
  min_cells_per_donor_state = 20L,
  min_donors = 1L,
  min_donors_per_state = 1L,
  min_states_per_donor = 2L,
  max_state_bins = 20L,
  seed = 1L
) {
  cell_group <- cell_group_tibble$cell_group[[1]]
  empty_aggregates <- tibble::tibble(
    cell_group = character(),
    aggregate_id = character(),
    donor_id = character(),
    state_bin = character(),
    barcodes = list(),
    n_cells = integer(),
    GEX_depth = numeric(),
    ATAC_depth = numeric()
  )

  make_diagnostic <- function(
    skipped_reason = NA_character_,
    n_cells = 0L,
    n_donors = 0L,
    n_state_bins = 0L,
    n_aggregates = 0L
  ) {
    tibble::tibble(
      cell_group = cell_group,
      chr = NA_character_,
      n_cells = as.integer(n_cells),
      n_donors = as.integer(n_donors),
      n_state_bins = as.integer(n_state_bins),
      n_aggregates = as.integer(n_aggregates),
      n_candidate_pairs = NA_integer_,
      n_detected_genes = NA_integer_,
      n_accessible_peaks = NA_integer_,
      skipped_reason = skipped_reason
    )
  }

  if (!donor_col %in% colnames(metadata_tibble)) {
    return(list(
      aggregates = empty_aggregates,
      diagnostics = make_diagnostic("missing_donor_metadata")
    ))
  }

  requested_barcodes <- intersect(
    cell_group_tibble$barcodes[[1]],
    rownames(ATAC_embedding_matrix)
  )
  GEX_depth_cols <- intersect(c("nCount_RNA", "gex_umis_count"), colnames(metadata_tibble))
  ATAC_depth_cols <- intersect(c("nCount_ATAC", "atac_fragments"), colnames(metadata_tibble))
  if (length(GEX_depth_cols) == 0L || length(ATAC_depth_cols) == 0L) {
    return(list(
      aggregates = empty_aggregates,
      diagnostics = make_diagnostic("missing_library_depth_metadata")
    ))
  }
  GEX_depth_col <- GEX_depth_cols[[1]]
  ATAC_depth_col <- ATAC_depth_cols[[1]]
  group_metadata <- metadata_tibble |>
    dplyr::distinct(.data$barcode_w_prefix, .keep_all = TRUE) |>
    dplyr::filter(.data$barcode_w_prefix %in% requested_barcodes) |>
    dplyr::transmute(
      barcode_w_prefix = .data$barcode_w_prefix,
      donor_id = as.character(.data[[donor_col]]),
      GEX_cell_depth = dplyr::coalesce(as.numeric(.data[[GEX_depth_col]]), 0),
      ATAC_cell_depth = dplyr::coalesce(as.numeric(.data[[ATAC_depth_col]]), 0)
    ) |>
    dplyr::filter(!is.na(.data$donor_id), .data$donor_id != "")

  donor_counts <- group_metadata |>
    dplyr::count(.data$donor_id, name = "n_cells")
  eligible_donors <- donor_counts |>
    dplyr::filter(.data$n_cells >= min_cells_per_donor_state * min_states_per_donor) |>
    dplyr::pull(.data$donor_id)
  group_metadata <- group_metadata |>
    dplyr::filter(.data$donor_id %in% eligible_donors)

  if (length(eligible_donors) < min_donors) {
    return(list(
      aggregates = empty_aggregates,
      diagnostics = make_diagnostic(
        "too_few_donors_with_sufficient_cells",
        n_cells = nrow(group_metadata),
        n_donors = length(eligible_donors)
      )
    ))
  }

  median_cells_per_donor <- stats::median(
    donor_counts$n_cells[donor_counts$donor_id %in% eligible_donors]
  )
  n_state_bins <- min(
    as.integer(max_state_bins),
    as.integer(floor(median_cells_per_donor / min_cells_per_donor_state))
  )
  if (n_state_bins < min_states_per_donor) {
    return(list(
      aggregates = empty_aggregates,
      diagnostics = make_diagnostic(
        "too_few_adaptive_state_bins",
        n_cells = nrow(group_metadata),
        n_donors = length(eligible_donors),
        n_state_bins = n_state_bins
      )
    ))
  }

  available_dims <- intersect(paste0("LSI_", dims), colnames(ATAC_embedding_matrix))
  if (length(available_dims) == 0L) {
    available_dims <- seq_len(min(length(dims), ncol(ATAC_embedding_matrix)))
  }
  group_metadata <- group_metadata |>
    dplyr::arrange(match(.data$barcode_w_prefix, rownames(ATAC_embedding_matrix)))
  group_embedding <- ATAC_embedding_matrix[
    group_metadata$barcode_w_prefix,
    available_dims,
    drop = FALSE
  ]
  variable_dims <- apply(group_embedding, 2, stats::sd, na.rm = TRUE) > 0
  group_embedding <- scale(group_embedding[, variable_dims, drop = FALSE])
  if (ncol(group_embedding) == 0L || anyNA(group_embedding)) {
    return(list(
      aggregates = empty_aggregates,
      diagnostics = make_diagnostic(
        "invalid_ATAC_embedding",
        n_cells = nrow(group_metadata),
        n_donors = length(eligible_donors),
        n_state_bins = n_state_bins
      )
    ))
  }

  set.seed(seed)
  state_fit <- stats::kmeans(
    x = group_embedding,
    centers = n_state_bins,
    iter.max = 1000L,
    nstart = 1L,
    algorithm = "Lloyd"
  )
  if (isTRUE(state_fit$ifault == 2L)) {
    stop("ATAC state clustering did not converge for ", cell_group)
  }
  state_width <- nchar(as.character(n_state_bins))
  group_metadata$state_bin <- sprintf(
    paste0("ATAC_state_%0", state_width, "d"),
    state_fit$cluster
  )

  donor_state_cells <- group_metadata |>
    dplyr::summarise(
      barcodes = list(.data$barcode_w_prefix),
      n_cells = dplyr::n(),
      GEX_depth = sum(.data$GEX_cell_depth),
      ATAC_depth = sum(.data$ATAC_cell_depth),
      .by = c("donor_id", "state_bin")
    ) |>
    dplyr::filter(.data$n_cells >= min_cells_per_donor_state)

  repeat {
    previous_n <- nrow(donor_state_cells)
    retained_states <- donor_state_cells |>
      dplyr::summarise(n_donors = dplyr::n_distinct(.data$donor_id), .by = "state_bin") |>
      dplyr::filter(.data$n_donors >= min_donors_per_state) |>
      dplyr::pull(.data$state_bin)
    donor_state_cells <- donor_state_cells |>
      dplyr::filter(.data$state_bin %in% retained_states)
    retained_donors <- donor_state_cells |>
      dplyr::summarise(n_states = dplyr::n_distinct(.data$state_bin), .by = "donor_id") |>
      dplyr::filter(.data$n_states >= min_states_per_donor) |>
      dplyr::pull(.data$donor_id)
    donor_state_cells <- donor_state_cells |>
      dplyr::filter(.data$donor_id %in% retained_donors)
    if (nrow(donor_state_cells) == previous_n) {
      break
    }
  }

  retained_n_donors <- dplyr::n_distinct(donor_state_cells$donor_id)
  retained_n_states <- dplyr::n_distinct(donor_state_cells$state_bin)
  if (retained_n_donors < min_donors || retained_n_states < min_states_per_donor) {
    return(list(
      aggregates = empty_aggregates,
      diagnostics = make_diagnostic(
        "insufficient_repeated_donor_state_support",
        n_cells = sum(donor_state_cells$n_cells),
        n_donors = retained_n_donors,
        n_state_bins = retained_n_states,
        n_aggregates = nrow(donor_state_cells)
      )
    ))
  }

  aggregates <- donor_state_cells |>
    dplyr::arrange(.data$donor_id, .data$state_bin) |>
    dplyr::mutate(
      cell_group = cell_group,
      aggregate_id = paste(make.names(cell_group), make.names(.data$donor_id), .data$state_bin, sep = "__"),
      .before = 1
    )

  list(
    aggregates = aggregates,
    diagnostics = make_diagnostic(
      n_cells = sum(aggregates$n_cells),
      n_donors = retained_n_donors,
      n_state_bins = retained_n_states,
      n_aggregates = nrow(aggregates)
    )
  )
}

make_peak_gene_correlation_group_chromosome_tibble <- function(
  donor_state_aggregates_tibble,
  chromosome_tibble,
  candidate_pairs_tibble
) {
  cell_groups <- donor_state_aggregates_tibble |>
    dplyr::distinct(.data$cell_group)

  if (nrow(cell_groups) == 0L) {
    return(tibble::tibble(
      cell_group = "__no_analyzable_cell_group__",
      chr = NA_character_,
      n_candidate_pairs = 0L,
      branch_id = "__no_analyzable_branch__",
      analyzable = FALSE
    ) |>
      dplyr::group_by(.data$branch_id) |>
      targets::tar_group())
  }

  tidyr::crossing(cell_groups, chromosome_tibble) |>
    dplyr::inner_join(
      candidate_pairs_tibble |> dplyr::count(.data$chr, name = "n_candidate_pairs"),
      by = "chr"
    ) |>
    dplyr::mutate(
      branch_id = paste(.data$cell_group, .data$chr, sep = "__"),
      analyzable = TRUE
    ) |>
    dplyr::group_by(branch_id) |>
    targets::tar_group()
}

make_aggregate_membership_matrix <- function(aggregates_tibble) {
  aggregate_ids <- aggregates_tibble$aggregate_id
  aggregate_barcodes <- aggregates_tibble$barcodes
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
#' Aggregate GEX and ATAC count matrices over donor-state pseudobulk groups.
#'
#' @param group_chromosome_tibble One-row branch tibble with `cell_group` and
#'   `chr` selecting the correlation branch.
#' @param donor_state_aggregates_tibble Aggregate membership tibble from
#'   `make_peak_gene_correlation_donor_state_record()`.
#' @param candidate_pairs_tibble Candidate peak-gene pairs; only pairs on the
#'   branch chromosome are used to select features.
#' @param GEX_counts_matrix Gene-by-cell count matrix; row names are gene IDs/names and column names are cell barcodes.
#' @param ATAC_counts_matrix Peak-by-cell count matrix; row names are peak IDs and column names are cell barcodes.
#' @return A list with branch identifiers, aggregate-level GEX and ATAC count
#'   matrices, and an aggregate depth summary tibble.
#' @keywords internal

make_peak_gene_correlation_aggregate_matrices <- function(
  group_chromosome_tibble,
  donor_state_aggregates_tibble,
  candidate_pairs_tibble,
  GEX_counts_matrix,
  ATAC_counts_matrix
) {
  cell_group <- group_chromosome_tibble$cell_group[[1]]
  chr <- group_chromosome_tibble$chr[[1]]
  if (!isTRUE(group_chromosome_tibble$analyzable[[1]])) {
    empty_matrix <- Matrix::Matrix(
      matrix(numeric(), nrow = 0L, ncol = 0L),
      sparse = TRUE
    )
    return(list(
      cell_group = cell_group,
      chr = chr,
      GEX_counts = empty_matrix,
      ATAC_counts = empty_matrix,
      aggregate_depth_tibble = tibble::tibble(
        cell_group = character(),
        chr = character(),
        aggregate_id = character(),
        donor_id = character(),
        state_bin = character(),
        n_cells = integer(),
        GEX_depth = numeric(),
        ATAC_depth = numeric()
      )
    ))
  }
  branch_pairs <- candidate_pairs_tibble |>
    dplyr::filter(.data$chr == !!chr)

  branch_aggregates <- donor_state_aggregates_tibble |>
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

  aggregate_depth_tibble <- branch_aggregates |>
    dplyr::arrange(match(.data$aggregate_id, colnames(membership_matrix))) |>
    dplyr::select(
      "cell_group",
      "aggregate_id",
      "donor_id",
      "state_bin",
      "n_cells",
      "GEX_depth",
      "ATAC_depth"
    ) |>
    dplyr::mutate(chr = chr, .after = "cell_group")

  list(
    cell_group = cell_group,
    chr = chr,
    GEX_counts = GEX_counts,
    ATAC_counts = ATAC_counts,
    aggregate_depth_tibble = aggregate_depth_tibble
  )
}

normalize_peak_gene_correlation_counts <- function(
  counts_matrix,
  depth = Matrix::colSums(counts_matrix),
  scale_factor = 1e6
) {
  counts_matrix <- methods::as(counts_matrix, "dgCMatrix")
  stopifnot(length(depth) == ncol(counts_matrix))
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
  GEX_norm <- normalize_peak_gene_correlation_counts(
    aggregate_matrices$GEX_counts,
    depth = aggregate_matrices$aggregate_depth_tibble$GEX_depth,
    scale_factor = scale_factor
  )
  ATAC_norm <- normalize_peak_gene_correlation_counts(
    aggregate_matrices$ATAC_counts,
    depth = aggregate_matrices$aggregate_depth_tibble$ATAC_depth,
    scale_factor = scale_factor
  )

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
      .data$is_analyzable_link,
      .data$cell_group == !!cell_group,
      .data$chr == !!chr,
      .data$gene_matrix_feature %in% rownames(GEX_norm),
      .data$peak %in% rownames(ATAC_norm)
    )

  if (nrow(branch_links) == 0L) {
    return(tibble::tibble(
      scatter_plot_name = character(),
      is_analyzable_link = logical(),
      primary_cell_group = character(),
      cell_group = character(),
      chr = character(),
      peak = character(),
      TargetGeneID = character(),
      TargetGene = character(),
      hierarchical_pvalue = numeric(),
      hierarchical_df = numeric(),
      hierarchical_FDR = numeric(),
      rank_in_cell_group = integer(),
      aggregate_id = character(),
      donor_id = character(),
      state_bin = character(),
      n_cells = numeric(),
      GEX_depth = numeric(),
      ATAC_depth = numeric(),
      gene_expression_logCPM = numeric(),
      peak_accessibility_logCPM = numeric(),
      gene_expression_residual = numeric(),
      peak_accessibility_residual = numeric(),
      correlation = numeric()
    ))
  }

  aggregate_depth_tibble <- normalized_aggregate_matrices$aggregate_depth_tibble |>
    dplyr::select(
      "aggregate_id",
      "donor_id",
      "state_bin",
      "n_cells",
      "GEX_depth",
      "ATAC_depth"
    )
  n_aggregates <- nrow(aggregate_depth_tibble)
  design <- make_peak_gene_correlation_design_matrix(aggregate_depth_tibble)

  purrr::map_dfr(seq_len(nrow(branch_links)), \(index) {
    link_row <- branch_links[index, , drop = FALSE]

    dplyr::bind_cols(
      link_row[rep(1L, n_aggregates), , drop = FALSE] |>
        dplyr::transmute(
          scatter_plot_name = .data$scatter_plot_name,
          is_analyzable_link = .data$is_analyzable_link,
          primary_cell_group = .data$cell_group,
          cell_group = !!cell_group,
          chr = .data$chr,
          peak = .data$peak,
          TargetGeneID = .data$TargetGeneID,
          TargetGene = .data$TargetGene,
          hierarchical_pvalue = .data$hierarchical_pvalue,
          hierarchical_df = .data$hierarchical_df,
          hierarchical_FDR = .data$hierarchical_FDR,
          rank_in_cell_group = .data$rank_in_cell_group
        ),
      aggregate_depth_tibble,
      tibble::tibble(
        gene_expression_logCPM = as.numeric(
          GEX_norm[link_row$gene_matrix_feature[[1]], ]
        ),
        peak_accessibility_logCPM = as.numeric(
          ATAC_norm[link_row$peak[[1]], ]
        ),
        gene_expression_residual = as.numeric(qr.resid(qr(design), as.numeric(GEX_norm[link_row$gene_matrix_feature[[1]], ]))),
        peak_accessibility_residual = as.numeric(qr.resid(qr(design), as.numeric(ATAC_norm[link_row$peak[[1]], ])))
      )
    ) |>
      dplyr::mutate(correlation = if (stats::sd(.data$peak_accessibility_residual) > 0 &&
        stats::sd(.data$gene_expression_residual) > 0)
        stats::cor(.data$peak_accessibility_residual, .data$gene_expression_residual) else NA_real_)
  })
}

#' Plot primary-group donor/depth residuals corresponding to the reported score.
#' Raw log1p CPM values remain available in the compact aggregate-value table.
plot_peak_gene_correlation_aggregate_scatter <- function(plot_tibble) {
  if (!isTRUE(plot_tibble$is_analyzable_link[[1]])) {
    return(make_empty_peak_gene_correlation_plot("No exploratory candidate links"))
  }
  plot_tibble <- plot_tibble |>
    dplyr::filter(.data$cell_group == .data$primary_cell_group)
  ggplot2::ggplot(plot_tibble, ggplot2::aes(
    x = .data$peak_accessibility_residual,
    y = .data$gene_expression_residual
  )) +
    ggplot2::geom_point(ggplot2::aes(color = .data$donor_id), size = 1, alpha = 0.65) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.4) +
    ggplot2::labs(
      title = paste(plot_tibble$TargetGene[[1]], plot_tibble$peak[[1]], plot_tibble$primary_cell_group[[1]], sep = " - "),
      subtitle = stringr::str_wrap("Full-scan hierarchical model: check whether the association is supported across donors or driven by a few aggregates; ranking alone does not establish significance.", width = 100),
      x = "ATAC log1p CPM residual", y = "GEX log1p CPM residual", color = "Donor",
      caption = stringr::str_wrap(paste0(
        "Points: non-overlapping donor-state aggregates, residualized for donor and depth; line: descriptive fit. Adjusted r = ",
        round(plot_tibble$correlation[[1]], 3), ". Hierarchical Kenward-Roger p = ", signif(plot_tibble$hierarchical_pvalue[[1]], 3),
        " (df = ", round(plot_tibble$hierarchical_df[[1]], 2), "); model: donor fixed intercepts, depth covariates and Gaussian random donor slopes. Hierarchical BH FDR = ",
        signif(plot_tibble$hierarchical_FDR[[1]], 3), ", adjusted across all eligible pairs/cell type, counting unreliable fits in the family size. Ranked positive estimable nonpromoter slopes without a significance cutoff. ",
        "Aggregates are not independent donor replicates."), width = 110)
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
  aggregate_depth_tibble <- normalized_aggregate_matrices$aggregate_depth_tibble
  n_donors <- dplyr::n_distinct(aggregate_depth_tibble$donor_id)
  design <- if (n_aggregates == 0L) {
    matrix(numeric(), nrow = 0L, ncol = 0L)
  } else {
    make_peak_gene_correlation_design_matrix(aggregate_depth_tibble)
  }
  design_rank <- if (n_aggregates == 0L) 0L else qr(design)$rank
  residual_df <- n_aggregates - design_rank - 1L

  branch_pairs <- if (is.na(chr)) {
    candidate_pairs_tibble[0, , drop = FALSE]
  } else {
    candidate_pairs_tibble |>
      dplyr::filter(.data$chr == !!chr)
  }

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
    n_aggregates == 0L ~ "no_analyzable_donor_state_aggregates",
    n_aggregates < min_aggregates ~ "too_few_donor_state_aggregates",
    residual_df < 5L ~ "insufficient_residual_degrees_of_freedom",
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
    n_donors = n_donors,
    design = design,
    design_rank = design_rank,
    residual_df = residual_df,
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
#' @inheritParams prepare_peak_gene_correlation_branch
#' @return One-row diagnostics tibble with branch size, retained feature counts,
#'   candidate-pair count, and optional skipped reason.
#' @keywords internal

diagnose_peak_gene_correlation_branch <- function(normalized_aggregate_matrices, candidate_pairs_tibble) {
  branch <- prepare_peak_gene_correlation_branch(normalized_aggregate_matrices, candidate_pairs_tibble)

  tibble::tibble(
    cell_group = branch$cell_group,
    chr = branch$chr,
    n_cells = NA_integer_,
    n_donors = branch$n_donors,
    n_state_bins = dplyr::n_distinct(
      normalized_aggregate_matrices$aggregate_depth_tibble$state_bin
    ),
    n_aggregates = branch$n_aggregates,
    n_candidate_pairs = nrow(branch$candidate_pairs),
    n_detected_genes = length(branch$detected_genes),
    n_accessible_peaks = length(branch$accessible_peaks),
    skipped_reason = branch$skipped_reason
  )
}

make_empty_peak_gene_correlation_plot <- function(message = "No analyzable peak-gene associations") {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = message) +
    ggplot2::theme_void()
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
  if (nrow(plot_tibble) == 0L) {
    return(make_empty_peak_gene_correlation_plot())
  }
  ggplot2::ggplot(
    plot_tibble,
    ggplot2::aes(x = .data$correlation_mid, y = .data$n_pairs)
  ) +
    ggplot2::geom_col(width = bin_width) +
    ggplot2::facet_wrap(~cell_group, scales = "free_y") +
    ggplot2::labs(
      title = "Peak-gene correlation distributions by cell group",
      subtitle = stringr::str_wrap("Look for shifts or long tails before focusing on selected links; correlation alone does not establish regulation.", width = 100),
      caption = paste("All non-missing donor/depth-adjusted correlations in the result table; no FDR or enhancer filter is applied here. Bin width:", bin_width,
        ". Facets have independent count scales; pair counts are not independent biological replicates."),
      x = "Donor/depth-adjusted correlation (exploratory candidates)",
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
    dplyr::mutate(link = is_peak_gene_link(results_tibble)) |>
    dplyr::summarise(
      tested_pairs = dplyr::n(),
      FDR_significant_pairs = if (all(is.na(.data$hierarchical_FDR))) NA_integer_ else
        sum(.data$hierarchical_FDR < 0.05, na.rm = TRUE),
      candidate_enhancer_links = if (all(is.na(.data$hierarchical_FDR))) NA_integer_ else
        sum(.data$link, na.rm = TRUE),
      .by = "cell_group"
    ) |>
    tidyr::pivot_longer(
      cols = -"cell_group",
      names_to = "metric",
      values_to = "n"
    )
}


plot_peak_gene_correlation_support_counts <- function(plot_tibble) {
  if (nrow(plot_tibble) == 0L) {
    return(make_empty_peak_gene_correlation_plot())
  }
  max_significant <- max(c(0, plot_tibble$n[
    plot_tibble$metric == "FDR_significant_pairs"
  ]), na.rm = TRUE)
  display_limit <- if (max_significant > 0) 1.1 * max_significant else 1
  plot_tibble <- plot_tibble |>
    dplyr::mutate(
      metric = factor(.data$metric,
        levels = c("tested_pairs", "FDR_significant_pairs", "candidate_enhancer_links"),
        labels = c("Tested pairs", "FDR-significant pairs", "Candidate enhancer links")),
      truncated = !is.na(.data$n) & .data$n > display_limit,
      count_position = pmin(tidyr::replace_na(.data$n, 0), display_limit),
      count_label = dplyr::if_else(is.na(.data$n), "NA",
        paste0(scales::comma(.data$n, accuracy = 1),
          dplyr::if_else(.data$truncated, " \u00bb", "")))
    )
  ggplot2::ggplot(
    plot_tibble,
    ggplot2::aes(
      x = .data$count_position,
      y = .data$cell_group,
      fill = .data$metric
    )
  ) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::geom_text(ggplot2::aes(label = .data$count_label,
      hjust = dplyr::if_else(.data$truncated, 1.05, -0.15)),
      position = ggplot2::position_dodge(width = 0.9), size = 3) +
    ggplot2::scale_x_continuous(labels = scales::label_comma(),
      expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::coord_cartesian(xlim = c(0, display_limit), clip = "off") +
    ggplot2::labs(
      title = "Peak-gene association support by cell group",
      subtitle = stringr::str_wrap("Compare significant pairs and candidate links across cell types; these are exploratory associations, not causal links. Bars marked \u00bb exceed the displayed range; labels retain their full counts.", width = 95),
      caption = stringr::str_wrap("Significant: hierarchical Kenward-Roger BH FDR < 0.05 across tested pairs within each cell group. Candidates also require a reliable positive donor-slope coefficient and exclude self-promoter pairs. The linear axis ends 10% above the largest significant-pair count (at 1 if none are positive). Longer bars are truncated. NA marks groups with no available FDR values; categories overlap.", width = 110),
      x = "Pair count (linear scale; long bars truncated)",
      y = "Cell group",
      fill = NULL
    )
}

#' Summarize technical features and significant peak-gene pair counts
#'
#' @param results_tibble Peak-gene correlation results.
#' @param donor_state_aggregates_tibble Retained donor-state pseudobulk memberships.
#' @return A compact table of technical features and descriptive correlations.
#' @keywords internal

prepare_peak_gene_support_technical_features <- function(
  results_tibble,
  donor_state_aggregates_tibble
) {
  support_tibble <- summarize_peak_gene_correlation_support_counts(results_tibble) |>
    tidyr::pivot_wider(names_from = "metric", values_from = "n") |>
    dplyr::filter(!is.na(.data$FDR_significant_pairs))
  if (nrow(support_tibble) == 0L) {
    return(tibble::tibble())
  }
  aggregate_tibble <- donor_state_aggregates_tibble |>
    dplyr::summarise(
      n_nuclei = sum(.data$n_cells),
      n_donors = dplyr::n_distinct(.data$donor_id),
      n_aggregates = dplyr::n(),
      n_states = dplyr::n_distinct(.data$state_bin),
      median_nuclei = stats::median(.data$n_cells),
      GEX_depth = stats::median(.data$GEX_depth / .data$n_cells),
      ATAC_depth = stats::median(.data$ATAC_depth / .data$n_cells),
      .by = "cell_group"
    )
  donor_tibble <- donor_state_aggregates_tibble |>
    dplyr::summarise(n_cells = sum(.data$n_cells), .by = c("cell_group", "donor_id")) |>
    dplyr::summarise(largest_donor_fraction = max(.data$n_cells) / sum(.data$n_cells),
      .by = "cell_group")
  gene_tibble <- results_tibble |>
    dplyr::distinct(.data$cell_group, .data$TargetGeneID, .data$gene_detected_frac) |>
    dplyr::summarise(gene_detection = stats::median(.data$gene_detected_frac), .by = "cell_group")
  peak_tibble <- results_tibble |>
    dplyr::distinct(.data$cell_group, .data$peak, .data$peak_accessible_frac) |>
    dplyr::summarise(peak_detection = stats::median(.data$peak_accessible_frac), .by = "cell_group")
  design_tibble <- results_tibble |>
    dplyr::distinct(.data$cell_group, .data$chr, .data$residual_df) |>
    dplyr::summarise(residual_df = stats::median(.data$residual_df), .by = "cell_group")
  feature_labels <- c(
    n_nuclei = "Nuclei in retained pseudobulks",
    n_donors = "Retained donors",
    largest_donor_fraction = "Fraction of nuclei from largest donor",
    n_aggregates = "Retained donor-state pseudobulks",
    n_states = "Retained ATAC states",
    median_nuclei = "Median nuclei per pseudobulk",
    GEX_depth = "Median RNA depth per nucleus",
    ATAC_depth = "Median ATAC depth per nucleus",
    residual_df = "Median residual degrees of freedom",
    gene_detection = "Median gene detection fraction",
    peak_detection = "Median peak detection fraction",
    tested_pairs = "Tested peak-gene pairs"
  )
  plot_tibble <- support_tibble |>
    dplyr::inner_join(aggregate_tibble, by = "cell_group", relationship = "one-to-one") |>
    dplyr::left_join(donor_tibble, by = "cell_group", relationship = "one-to-one") |>
    dplyr::left_join(gene_tibble, by = "cell_group", relationship = "one-to-one") |>
    dplyr::left_join(peak_tibble, by = "cell_group", relationship = "one-to-one") |>
    dplyr::left_join(design_tibble, by = "cell_group", relationship = "one-to-one") |>
    tidyr::pivot_longer(dplyr::all_of(names(feature_labels)), names_to = "feature", values_to = "value") |>
    dplyr::filter(is.finite(.data$value))
  correlations <- plot_tibble |>
    dplyr::summarise(
      rho = if (dplyr::n() >= 3L && dplyr::n_distinct(.data$value) > 1L &&
          dplyr::n_distinct(.data$FDR_significant_pairs) > 1L) {
        stats::cor(.data$value, .data$FDR_significant_pairs, method = "spearman")
      } else NA_real_,
      n_cell_types = dplyr::n(), .by = "feature"
    ) |>
    dplyr::mutate(facet_label = paste0(feature_labels[.data$feature],
      "\nSpearman rho = ", ifelse(is.na(.data$rho), "NA", sprintf("%.2f", .data$rho)),
      "; cell types = ", .data$n_cell_types))
  plot_tibble <- plot_tibble |>
    dplyr::left_join(correlations, by = "feature", relationship = "many-to-one") |>
    dplyr::mutate(facet_label = factor(.data$facet_label,
      levels = correlations$facet_label[match(names(feature_labels), correlations$feature)]))
  plot_tibble
}

#' Plot significant peak-gene pair counts against technical features
#'
#' @param plot_tibble Output of `prepare_peak_gene_support_technical_features()`.
#' @return A faceted labelled scatter plot with descriptive Spearman correlations.
#' @keywords internal

plot_peak_gene_significant_pairs_vs_technical_features <- function(plot_tibble) {
  if (nrow(plot_tibble) == 0L) {
    return(make_empty_peak_gene_correlation_plot())
  }
  cell_groups <- sort(unique(plot_tibble$cell_group))
  ggplot2::ggplot(plot_tibble,
    ggplot2::aes(x = .data$value, y = .data$FDR_significant_pairs,
      colour = .data$cell_group)) +
    ggplot2::geom_point(size = 2) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = gsub("_", " ", .data$cell_group)),
      seed = 1, max.overlaps = Inf, box.padding = 0.4,
      min.segment.length = 0, size = 2.8
    ) +
    ggplot2::scale_colour_manual(values = stats::setNames(
      grDevices::hcl.colors(length(cell_groups), "Dark 3"), cell_groups
    ), guide = "none") +
    ggplot2::facet_wrap(~facet_label, scales = "free_x", ncol = 4) +
    ggplot2::scale_x_continuous(labels = scales::label_comma(), limits = c(0, NA),
      expand = ggplot2::expansion(mult = c(0, 0.2))) +
    ggplot2::scale_y_continuous(labels = scales::label_comma(), limits = c(0, NA),
      expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::labs(
      title = "Technical features associated with peak-gene discovery counts",
      subtitle = "Each point is a WNN-derived cell type. Compare detection yield with sampling and measurement quality; these associations do not separate technical effects from biology.",
      caption = stringr::str_wrap("Significant pairs have hierarchical BH FDR < 0.05 within each cell type, including self-promoter pairs. Sampling summaries use retained, mutually exclusive donor-state pseudobulks. Depth is the median across pseudobulks of total RNA UMIs or ATAC fragments divided by nuclei count. Detection fractions are medians across unique genes or peaks represented in tested pairs, measured as the fraction of pseudobulks with nonzero counts; they are conditional on feature filtering. Residual degrees of freedom are summarized across chromosomes. Facets use separate linear x-scales and a shared y-scale. Spearman correlations weight cell types equally and are descriptive; NA denotes insufficient variation or fewer than three points. Cell types with unavailable FDR are omitted; observed zeros are retained.", width = 190),
      x = "Technical feature value",
      y = "FDR-significant peak-gene pairs"
    )
}

#' Summarize peak-gene correlations and p-values by TSS distance
#'
#' @param results_tibble Peak-gene correlation results.
#' @return A long tibble of median correlation and median -log10 hierarchical
#'   p-value by cell group and 10 kb distance-bin midpoint.
#' @keywords internal

summarize_peak_gene_correlation_by_distance <- function(results_tibble) {
  results_tibble |>
    dplyr::select("cell_group", "correlation", "isSelfPromoter", "distance", "hierarchical_pvalue") |>
    dplyr::filter(!is.na(.data$correlation), !.data$isSelfPromoter) |>
    dplyr::mutate(distance_kb = pmin(24, abs(.data$distance) %/% 10000) * 10 + 5) |>
    dplyr::summarise(
      median_correlation = stats::median(.data$correlation),
      median_neg_log10_p = stats::median(-log10(.data$hierarchical_pvalue), na.rm = TRUE),
      .by = c("cell_group", "distance_kb")
    ) |>
    tidyr::pivot_longer(c("median_correlation", "median_neg_log10_p"),
      names_to = "metric", values_drop_na = TRUE)
}

#' Plot peak-gene correlations and p-values by TSS distance
#'
#' @param plot_tibble Compact output from
#'   `summarize_peak_gene_correlation_by_distance()`.
#' @return A ggplot ready for saving.
#' @keywords internal

plot_peak_gene_correlation_by_distance <- function(plot_tibble) {
  if (nrow(plot_tibble) == 0L) {
    return(make_empty_peak_gene_correlation_plot())
  }
  endpoint_tibble <- plot_tibble |>
    dplyr::slice_max(.data$distance_kb, n = 1, with_ties = FALSE, by = c("cell_group", "metric"))
  mean_tibble <- plot_tibble |>
    dplyr::summarise(value = mean(.data$value), .by = c("metric", "distance_kb"))
  # Uniform null p-values have median -log10(p) = log10(2).
  reference_tibble <- tibble::tibble(metric = c("median_correlation", "median_neg_log10_p"),
    value = c(0, log10(2)))
  cell_groups <- sort(unique(plot_tibble$cell_group))
  ggplot2::ggplot(
    plot_tibble,
    ggplot2::aes(
      x = .data$distance_kb,
      y = .data$value,
      colour = .data$cell_group
    )
  ) +
    ggplot2::geom_hline(
      data = reference_tibble,
      ggplot2::aes(yintercept = .data$value),
      linetype = 3,
      color = "grey70"
    ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_line(data = mean_tibble, colour = "black",
      linetype = "dashed", linewidth = 1) +
    ggrepel::geom_text_repel(
      data = endpoint_tibble,
      ggplot2::aes(label = gsub("_", " ", .data$cell_group)),
      nudge_x = 15, direction = "y", hjust = 0,
      box.padding = 0.4, min.segment.length = 0,
      max.overlaps = Inf, max.iter = 10000, seed = 1, size = 3.5
    ) +
    ggrepel::geom_text_repel(
      data = dplyr::slice_max(mean_tibble, .data$distance_kb, n = 1, by = "metric"),
      label = "Mean across cell types", colour = "black", fontface = "bold",
      nudge_x = -20, nudge_y = 0.015, min.segment.length = 0,
      seed = 1, size = 3.5
    ) +
    ggplot2::facet_wrap(~metric, ncol = 1, scales = "free_y", strip.position = "left",
      labeller = ggplot2::as_labeller(c(median_correlation = "Median correlation",
        median_neg_log10_p = "Median -log10(p)"))) +
    ggplot2::scale_colour_manual(values = stats::setNames(
      grDevices::hcl.colors(length(cell_groups), "Dark 3"), cell_groups
    ), guide = "none") +
    ggplot2::scale_x_continuous(breaks = seq(0, 250, 50), limits = c(0, NA),
      expand = ggplot2::expansion(mult = c(0.01, 0.3))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.15)) +
    ggplot2::theme(strip.placement = "outside", strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = ggplot2::rel(1))) +
    ggplot2::labs(
      title = "Peak-gene correlation and association p-values by distance from the gene TSS",
      subtitle = stringr::str_wrap("Compare distance-dependent trends across WNN-derived cell types; colours and endpoint labels identify each curve. -log10(p) also rises with each cell type's statistical power, so compare curve shapes rather than levels between cell types. Nearby peaks are not necessarily regulatory.", width = 110),
      caption = stringr::str_wrap("Medians among non-missing, non-self-promoter pairs, with no FDR filter. Correlation is the descriptive donor/depth-adjusted Pearson correlation; p is the nominal two-sided hierarchical Kenward-Roger p-value for the donor-varying peak slope, available only for reliable fits. Absolute peak-centre to TSS distances use 10 kb bins plotted at their midpoints; the final 240–250 kb bin includes the 250 kb boundary. Black dashed curves are means of the available cell-type medians in each bin, weighting cell types equally rather than pooling pairs. Pair counts vary between bins and cell types, and rows use separate y-scales. Grey dotted lines mark zero correlation and log10(2) ≈ 0.30, the median -log10(p) expected for null pairs.", width = 130),
      x = "Absolute TSS distance, kb",
      y = NULL
    )
}
