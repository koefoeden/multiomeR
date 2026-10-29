#' Select the cell groups scored by activity by contact
#'
#' Keep labelled WNN cell types with at least `min_nuclei` nuclei, the floor the
#' genetic-enrichment module uses for a stable summed accessibility profile.
#' Every cell group the correlation can test has at least 200 nuclei, so all of
#' them are included.
#'
#' @param metadata_tibble Per-nucleus metadata with `barcode_w_prefix` and
#'   `WNN_harmony_SNN_cluster_cell_type`.
#' @param min_nuclei Minimum nuclei of a scored cell group.
#' @return A tibble with `cell_group` (the WNN label, as in the correlation
#'   results), `cluster` (the pseudobulk and genetic-enrichment name) and
#'   `n_nuclei`.
#' @keywords internal

make_peak_gene_ABC_cell_groups <- function(metadata_tibble, min_nuclei = 100L) {
  metadata_tibble |>
    dplyr::distinct(.data$barcode_w_prefix, .keep_all = TRUE) |>
    dplyr::count(cell_group = as.character(.data$WNN_harmony_SNN_cluster_cell_type), name = "n_nuclei") |>
    dplyr::filter(
      !is.na(.data$cell_group), .data$cell_group != "", !startsWith(.data$cell_group, "Unassigned"),
      .data$n_nuclei >= min_nuclei
    ) |>
    dplyr::mutate(cluster = stringr::str_replace_all(.data$cell_group, "_", "-"), .after = "cell_group")
}

#' Sum cell type by donor pseudobulks per cell type
#'
#' @param pseudobulk_matrix Feature-by-pseudobulk matrix with
#'   `<cluster>_<donor_id>` columns.
#' @param clusters Clusters to keep, in output column order.
#' @return A dense feature-by-cluster count matrix.
#' @keywords internal

sum_pseudobulks_by_cluster <- function(pseudobulk_matrix, clusters) {
  sample_clusters <- get_pseudobulk_sample_tibble(pseudobulk_matrix)$cluster
  is_kept <- sample_clusters %in% clusters
  BPCells::pseudobulk_matrix(pseudobulk_matrix[, is_kept], sample_clusters[is_kept])[, clusters, drop = FALSE]
}

#' Score candidate peak-gene pairs by activity by contact
#'
#' ABC(peak, gene) = A C / sum(A C), where the sum runs over all consensus peaks
#' whose centres lie within `window` of the gene's TSS, including its own
#' promoter peak (Fulco et al. 2019). Activity A is the cell group's ATAC counts
#' per million in the peak. Contact C is the ABC pipeline's power law
#' (max(d, 5 kb) + 1)^-gamma of the distance d from peak centre to TSS; its
#' scale cancels in the ratio. A gene is scored in a cell group when it is
#' expressed there by the S-LDSC gene rule: at least 10 counts and 1 count per
#' million.
#'
#' @param candidate_pairs_tibble The module's candidate peak-gene pairs.
#' @param gene_TSS_tibble The module's gene TSS table.
#' @param consensus_peak_GRanges Consensus peaks, named as the ATAC features.
#' @param ATAC_pseudobulk_matrix Consensus peak by cell type and donor counts.
#' @param GEX_pseudobulk_matrix Gene by cell type and donor counts.
#' @param cell_groups_tibble Output of `make_peak_gene_ABC_cell_groups()`.
#' @param window Distance from the TSS within which peaks enter the sum.
#' @param gamma Power-law exponent of contact with distance.
#' @param min_contact_distance Distance below which contact stays at its value
#'   at this distance.
#' @return A tibble with one row per candidate pair and cell group in which
#'   the gene is expressed: `cell_group`, `peak`, `TargetGeneID`, `TargetGene`,
#'   `link_class`, absolute `distance` and `ABC_score`.
#' @keywords internal

score_peak_gene_ABC <- function(candidate_pairs_tibble, gene_TSS_tibble, consensus_peak_GRanges,
                                ATAC_pseudobulk_matrix, GEX_pseudobulk_matrix, cell_groups_tibble,
                                window = 5000000L, gamma = 1.024238616787792, min_contact_distance = 5000L) {
  get_contact <- \(distance) (pmax(abs(distance), min_contact_distance) + 1)^-gamma
  ATAC_counts_matrix <- sum_pseudobulks_by_cluster(ATAC_pseudobulk_matrix, cell_groups_tibble$cluster)
  activity_matrix <- sweep(ATAC_counts_matrix, 2, colSums(ATAC_counts_matrix), "/") * 1e6
  GEX_counts_matrix <- sum_pseudobulks_by_cluster(GEX_pseudobulk_matrix, cell_groups_tibble$cluster)
  is_expressed_matrix <- GEX_counts_matrix >= 10 & sweep(GEX_counts_matrix, 2, colSums(GEX_counts_matrix), "/") >= 1e-6
  colnames(activity_matrix) <- colnames(is_expressed_matrix) <- cell_groups_tibble$cell_group

  peak_centres <- as.integer(floor((GenomicRanges::start(consensus_peak_GRanges) + GenomicRanges::end(consensus_peak_GRanges)) / 2))
  peak_chromosomes <- as.character(GenomicRanges::seqnames(consensus_peak_GRanges))
  gene_TSS_tibble <- gene_TSS_tibble[gene_TSS_tibble$gene_id %in% candidate_pairs_tibble$TargetGeneID, ]
  # Each gene's sum of A C over the peaks within the window, as a sparse
  # gene-by-peak contact matrix times the peak-by-cell-group activity.
  denominator_matrix <- split(seq_len(nrow(gene_TSS_tibble)), gene_TSS_tibble$chr) |>
    purrr::imap(\(gene_rows, chromosome) {
      peak_rows <- which(peak_chromosomes == chromosome)
      peak_rows <- peak_rows[order(peak_centres[peak_rows])]
      centres <- peak_centres[peak_rows]
      TSS <- gene_TSS_tibble$TSS[gene_rows]
      first_peak <- findInterval(TSS - window - 1L, centres) + 1L
      n_peaks <- findInterval(TSS + window, centres) - first_peak + 1L
      gene_index <- rep(seq_along(gene_rows), n_peaks)
      peak_index <- sequence(n_peaks, from = first_peak)
      contact_matrix <- Matrix::sparseMatrix(
        gene_index, peak_index, x = get_contact(centres[peak_index] - TSS[gene_index]),
        dims = c(length(gene_rows), length(peak_rows)), dimnames = list(gene_TSS_tibble$gene_id[gene_rows], NULL)
      )
      as.matrix(contact_matrix %*% activity_matrix[peak_rows, , drop = FALSE])
    }) |>
    do.call(what = rbind)

  pairs_tibble <- dplyr::transmute(
    candidate_pairs_tibble, .data$peak, .data$TargetGeneID, .data$TargetGene, .data$link_class,
    distance = abs(.data$distance)
  )
  pair_peak_index <- match(pairs_tibble$peak, rownames(activity_matrix))
  pair_gene_index <- match(pairs_tibble$TargetGeneID, rownames(denominator_matrix))
  pair_feature_index <- match(candidate_pairs_tibble$gene_matrix_feature, rownames(is_expressed_matrix))
  pair_contact <- get_contact(pairs_tibble$distance)
  purrr::map(cell_groups_tibble$cell_group, \(group) {
    is_scored <- is_expressed_matrix[pair_feature_index, group]
    tibble::tibble(
      cell_group = group, pairs_tibble[is_scored, ],
      ABC_score = activity_matrix[pair_peak_index[is_scored], group] * pair_contact[is_scored] /
        denominator_matrix[pair_gene_index[is_scored], group]
    )
  }) |>
    dplyr::bind_rows()
}
