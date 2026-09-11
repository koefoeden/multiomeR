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

combine_peak_gene_correlation_donor_state_records <- function(records, component) {
  stopifnot(component %in% c("aggregates", "diagnostics"))
  purrr::map_dfr(records, component)
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
      correlation = numeric(),
      FDR = numeric(),
      rank_in_cell_group = integer(),
      rank_for_gene = integer(),
      aggregate_id = character(),
      donor_id = character(),
      state_bin = character(),
      n_cells = numeric(),
      GEX_depth = numeric(),
      ATAC_depth = numeric(),
      gene_expression_logCPM = numeric(),
      peak_accessibility_logCPM = numeric(),
      gene_expression_residual = numeric(),
      peak_accessibility_residual = numeric()
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
        ),
        gene_expression_residual = as.numeric(qr.resid(qr(design), as.numeric(GEX_norm[link_row$gene_matrix_feature[[1]], ]))),
        peak_accessibility_residual = as.numeric(qr.resid(qr(design), as.numeric(ATAC_norm[link_row$peak[[1]], ])))
      )
    )
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
      subtitle = paste0("Donor/depth-adjusted r = ", round(plot_tibble$correlation[[1]], 3), "; conditional BH FDR = ", signif(plot_tibble$FDR[[1]], 3)),
      x = "ATAC log1p CPM residual", y = "GEX log1p CPM residual", color = "Donor",
      caption = "Each point is a non-overlapping donor-by-state aggregate. No population-level significance claim."
    )
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
    TargetGeneStart = integer(),
    TargetGeneEnd = integer(),
    distance = integer(),
    isSelfPromoter = logical(),
    isTargetGeneBody = logical(),
    link_class = character(),
    n_aggregates = integer(),
    n_donors = integer(),
    design_rank = integer(),
    residual_df = integer(),
    mean_gene_expression = numeric(),
    gene_detected_frac = numeric(),
    mean_peak_accessibility = numeric(),
    peak_accessible_frac = numeric(),
    raw_correlation = numeric(),
    correlation = numeric(),
    coefficient = numeric(),
    n_informative_donors = integer(),
    n_positive_donors = integer(),
    donor_direction_agreement = numeric(),
    max_donor_covariance_fraction = numeric(),
    inference_status = character(),
    association_SE = numeric(),
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

#' Score peak gene correlations for cell group
#'
#' Compute donor-adjusted aggregate-level associations for candidate pairs.
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
#'   fractions, adjusted correlation, coefficient, donor-support diagnostics, and
#'   approximate HC3 t tests conditional on sampled donors. Non-overlapping
#'   aggregates are assumed to have independent errors after donor/depth
#'   adjustment. These tests do not establish population-level replication.
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
  donor_id <- normalized_aggregate_matrices$aggregate_depth_tibble$donor_id
  n_donors <- branch$n_donors
  design <- branch$design
  design_qr <- qr(design)
  residual_df <- branch$residual_df
  nuisance_leverage <- rowSums(qr.Q(design_qr)^2)

  scored_pairs <- branch$candidate_pairs |>
    dplyr::group_by(gene_matrix_feature) |>
    dplyr::group_split() |>
    purrr::map_dfr(\(gene_pairs) {
      gene_feature <- gene_pairs$gene_matrix_feature[[1]]
      peak_features <- gene_pairs$peak

      gene_vec <- as.numeric(GEX_norm[gene_feature, ])
      gene_centered <- gene_vec - mean(gene_vec)
      raw_gene_ss <- sum(gene_centered^2)
      gene_residual <- as.numeric(qr.resid(design_qr, gene_vec))
      residual_gene_ss <- sum(gene_residual^2)

      peak_matrix <- as.matrix(ATAC_norm[peak_features, , drop = FALSE])
      peak_centered <- sweep(peak_matrix, 1, rowMeans(peak_matrix), "-")
      raw_peak_ss <- rowSums(peak_centered^2)
      peak_residual <- residualize_peak_gene_correlation_matrix(peak_matrix, design)
      residual_peak_ss <- rowSums(peak_residual^2)
      # Avoid treating floating-point residuals of fully explained features as
      # biological variation.
      variable_pairs <- residual_peak_ss > 1e-12 * pmax(raw_peak_ss, 1) &
        residual_gene_ss > 1e-12 * max(raw_gene_ss, 1)

      raw_denominator <- sqrt(raw_peak_ss * raw_gene_ss)
      raw_correlation <- as.numeric((peak_centered %*% gene_centered) / raw_denominator)
      raw_correlation[raw_denominator == 0 | !is.finite(raw_correlation)] <- NA_real_
      raw_correlation <- pmax(pmin(raw_correlation, 1), -1)

      cross_product <- as.numeric(peak_residual %*% gene_residual)
      denominator <- sqrt(residual_peak_ss * residual_gene_ss)
      correlation <- cross_product / denominator
      correlation[!variable_pairs | !is.finite(correlation)] <- NA_real_
      correlation <- pmax(pmin(correlation, 1), -1)

      coefficient <- cross_product / residual_peak_ss
      coefficient[!variable_pairs | !is.finite(coefficient)] <- NA_real_
      association_residual <- matrix(gene_residual, nrow(peak_residual), ncol(peak_residual), byrow = TRUE) -
        sweep(peak_residual, 1, coefficient, "*")
      leverage <- sweep(sweep(peak_residual^2, 1, residual_peak_ss, "/"), 2, nuisance_leverage, "+")
      HC3_scores <- peak_residual * association_residual / (1 - leverage)
      association_SE <- sqrt(rowSums(HC3_scores^2)) / residual_peak_ss
      association_SE[!variable_pairs | apply(leverage >= 1 - 1e-8, 1, any) | !is.finite(association_SE)] <- NA_real_
      nominal_pvalue <- 2 * stats::pt(abs(coefficient / association_SE), df = residual_df, lower.tail = FALSE)
      # Describe donor support independently of inferential eligibility. The
      # score uses within-donor variation; donors need >=3 aggregates and
      # nonzero residual variation in both features to inform direction.
      donor_cross_products <- rowsum(
        t(sweep(peak_residual, 2, gene_residual, "*")),
        group = donor_id, reorder = FALSE
      )
      donor_peak_ss <- rowsum(t(peak_residual^2), group = donor_id, reorder = FALSE)
      donor_gene_ss <- rowsum(matrix(gene_residual^2, ncol = 1), group = donor_id, reorder = FALSE)[, 1]
      donor_n <- rowsum(matrix(1, n_aggregates, 1), group = donor_id, reorder = FALSE)[, 1]
      informative <- donor_peak_ss > .Machine$double.eps &
        donor_gene_ss > .Machine$double.eps & donor_n >= 3L
      informative[, !variable_pairs] <- FALSE
      n_informative_donors <- colSums(informative)
      n_positive_donors <- colSums(informative & donor_cross_products > 0)
      direction_agreement <- colSums(informative &
        sweep(sign(donor_cross_products), 2, sign(cross_product), "==")) / n_informative_donors
      direction_agreement[n_informative_donors == 0L] <- NA_real_
      absolute_contribution <- abs(donor_cross_products)
      max_donor_covariance_fraction <- apply(absolute_contribution, 2, max) / colSums(absolute_contribution)
      max_donor_covariance_fraction[!variable_pairs | !is.finite(max_donor_covariance_fraction)] <- NA_real_

      gene_pairs |>
        dplyr::mutate(
          n_aggregates = n_aggregates,
          n_donors = n_donors,
          design_rank = branch$design_rank,
          residual_df = residual_df,
          mean_gene_expression = normalized_aggregate_matrices$mean_gene_expression[.data$gene_matrix_feature],
          gene_detected_frac = normalized_aggregate_matrices$gene_detected_frac[.data$gene_matrix_feature],
          mean_peak_accessibility = normalized_aggregate_matrices$mean_peak_accessibility[.data$peak],
          peak_accessible_frac = normalized_aggregate_matrices$peak_accessible_frac[.data$peak],
          raw_correlation = raw_correlation,
          correlation = correlation,
          coefficient = as.numeric(coefficient),
          n_informative_donors = as.integer(n_informative_donors),
          n_positive_donors = as.integer(n_positive_donors),
          donor_direction_agreement = direction_agreement,
          max_donor_covariance_fraction = max_donor_covariance_fraction,
          inference_status = "conditional_aggregate_HC3",
          association_SE = as.numeric(association_SE),
          nominal_pvalue = as.numeric(nominal_pvalue)
        )
    }) |>
    dplyr::mutate(cell_group = branch$cell_group, .before = 1) |>
    dplyr::select(-dplyr::all_of("gene_matrix_feature"))

  scored_pairs
}

#' Finalize peak gene correlation results
#'
#' Add exploratory status, ranks, and aggregation label to peak-gene correlations.
#'
#' @param results_tibbles List of chromosome-level scored-pair tibbles.
#' @param aggregation Aggregation label inserted into the first output column.
#' @return A standardized results tibble with within-cell-group BH FDR plus
#'   within-cell-group and within-gene ranks. Empty input returns the same schema.
#' @keywords internal

finalize_peak_gene_correlation_results <- function(results_tibbles, aggregation) {
  # Bypass targets/vctrs concatenation and bind once, without recycling the
  # named numeric columns carried by the chromosome-level tibbles.
  results_tibble <- data.table::rbindlist(results_tibbles, use.names = TRUE)
  aggregation_label <- aggregation

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
    "TargetGeneStart",
    "TargetGeneEnd",
    "distance",
    "isSelfPromoter",
    "isTargetGeneBody",
    "link_class",
    "n_aggregates",
    "n_donors",
    "design_rank",
    "residual_df",
    "mean_gene_expression",
    "gene_detected_frac",
    "mean_peak_accessibility",
    "peak_accessible_frac",
    "raw_correlation",
    "correlation",
    "coefficient",
    "n_informative_donors",
    "n_positive_donors",
    "donor_direction_agreement",
    "max_donor_covariance_fraction",
    "inference_status",
    "association_SE",
    "nominal_pvalue",
    "FDR",
    "rank_in_cell_group",
    "rank_for_gene"
  )

  if (nrow(results_tibble) == 0L) {
    return(empty_peak_gene_correlation_results_tibble() |>
      dplyr::mutate(
        aggregation = character(),
        FDR = numeric(),
        rank_in_cell_group = integer(),
        rank_for_gene = integer(),
        .before = 1
      ) |> dplyr::select(dplyr::all_of(output_cols)))
  }

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
  # as_tibble.data.frame otherwise dispatches to as.data.frame.data.table,
  # which deep-copies every column before removing the data.table class.
  data.table::setDF(results_tibble)
  tibble::as_tibble(results_tibble)
}

make_peak_gene_correlation_links <- function(results_tibble) {
  results_tibble |>
    dplyr::filter(
      .data$correlation >= 0.15,
      .data$FDR < 0.05,
      !.data$isSelfPromoter
    ) |>
    dplyr::arrange(.data$cell_group, .data$rank_in_cell_group) |>
    dplyr::mutate(
      rank_in_cell_group = dplyr::row_number(),
      .by = "cell_group"
    )
}

make_peak_gene_correlation_top_links <- function(
  links_tibble,
  candidate_pairs_tibble,
  n_per_cell_group = 3L
) {
  top_links <- links_tibble |>
    dplyr::filter(.data$rank_in_cell_group <= n_per_cell_group) |>
    dplyr::arrange(.data$cell_group, .data$rank_in_cell_group) |>
    dplyr::mutate(
      scatter_plot_name = paste(
        make.names(.data$cell_group),
        sprintf("rank%03d", .data$rank_in_cell_group),
        make.names(.data$TargetGene),
        .data$chr,
        sep = "_"
      ),
      is_analyzable_link = TRUE
    ) |>
    dplyr::left_join(
      candidate_pairs_tibble |>
        dplyr::select("peak", "TargetGeneID", "gene_matrix_feature") |>
        dplyr::distinct(),
      by = c("peak", "TargetGeneID")
    )

  if (nrow(top_links) > 0L) {
    return(top_links)
  }

  links_tibble[NA_integer_, , drop = FALSE] |>
    dplyr::mutate(
      cell_group = "__no_analyzable_link__",
      chr = NA_character_,
      scatter_plot_name = "no_analyzable_peak_gene_link",
      is_analyzable_link = FALSE,
      gene_matrix_feature = NA_character_
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
    dplyr::select(
      "cell_group",
      "FDR",
      "correlation",
      "isSelfPromoter",
      "isTargetGeneBody"
    ) |>
    dplyr::summarise(
      tested_pairs = dplyr::n(),
      FDR_significant_pairs = if (all(is.na(.data$FDR))) NA_integer_ else sum(.data$FDR < 0.05, na.rm = TRUE),
      candidate_enhancer_links = sum(
        .data$correlation >= 0.15 &
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
  if (nrow(plot_tibble) == 0L) {
    return(make_empty_peak_gene_correlation_plot())
  }
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
      significant_fraction = if (all(is.na(.data$FDR))) NA_real_ else mean(.data$FDR < 0.05, na.rm = TRUE),
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
  if (nrow(plot_tibble) == 0L) {
    return(make_empty_peak_gene_correlation_plot())
  }
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
