# Locus attribution, an experimental route of the genetic-enrichment module, off by default: for each credible set of a
# GWAS, which cell types its posterior probability in accessible peaks points to, from the locus's own accessibility;
# and, for context, the cell-type weights that best explain all of the trait's loci.

#' Get a peak-by-cell-type accessibility likelihood
#'
#' A cell type takes a peak in proportion to the peak's share of its ATAC reads. Each profile is shrunk toward the
#' pooled profile of all cell types by the reads of `prior_nuclei` average nuclei, so that a peak without reads in a
#' small cell type counts as average rather than closed. The default of 50 nuclei adds about half a read per peak, at
#' most half of the reads of a cell type of 100 nuclei, the fewest the pseudobulk keeps. A uniform `background`
#' column takes peaks that no cell type explains better than reads spread evenly over all peaks: otherwise the cell
#' types with the most diffuse reads, such as neurons in a peripheral aggregation, collect every weak peak.
#'
#' @param pseudobulk_ATAC_data_matrix BPCells peak-by-cell-type ATAC counts.
#' @param peaks Row indices of the peaks to return.
#' @param support_tibble Nuclei per cell type, as `cell_type_pseudobulk_support_tibble.ATAC`.
#' @param prior_nuclei Average nuclei whose reads shrink each profile.
#' @return Dense matrix of likelihoods, one row per element of `peaks` and one column per cell type and `background`.
#' @keywords internal
get_locus_attribution_likelihood_matrix <- function(pseudobulk_ATAC_data_matrix, peaks, support_tibble, prior_nuclei = 50) {
  cell_type_reads <- BPCells::colSums(pseudobulk_ATAC_data_matrix)
  pooled_share <- BPCells::rowSums(pseudobulk_ATAC_data_matrix)[peaks] / sum(cell_type_reads)
  prior_reads <- prior_nuclei * sum(cell_type_reads) / sum(support_tibble$n_cells[match(names(cell_type_reads), support_tibble$cluster)])
  counts <- as.matrix(pseudobulk_ATAC_data_matrix[peaks, ])
  cbind(sweep(counts + prior_reads * pooled_share, 2, cell_type_reads + prior_reads, "/"),
    background = 1 / nrow(pseudobulk_ATAC_data_matrix))
}

#' Fit a trait's cell-type weights to its locus-peak pairs
#'
#' Expectation-maximisation of a mixture over cell types, with each pair weighted by its PIP and a Dirichlet prior of
#' `pseudo_loci` loci spread evenly over cell types.
#'
#' @param likelihood_matrix Pair-by-cell-type likelihoods.
#' @param PIP Posterior probability of each pair.
#' @param pseudo_loci Total prior weight, in loci.
#' @return Cell-type weights summing to 1.
#' @keywords internal
fit_locus_attribution_weights <- function(likelihood_matrix, PIP, pseudo_loci, tolerance = 1e-6, max_iterations = 1000L) {
  weights <- rep(1 / ncol(likelihood_matrix), ncol(likelihood_matrix))
  for (iteration in seq_len(max_iterations)) {
    # A cell type's PIP-weighted responsibility summed over pairs, in two matrix-vector products.
    responsibility <- weights * drop(crossprod(likelihood_matrix, PIP / drop(likelihood_matrix %*% weights)))
    updated <- (responsibility + pseudo_loci / length(weights)) / (sum(PIP) + pseudo_loci)
    converged <- max(abs(updated - weights)) < tolerance
    weights <- updated
    if (converged) break
  }
  rlang::set_names(weights, colnames(likelihood_matrix))
}

#' Attribute one GWAS's credible sets to cell types
#'
#' Each locus-peak pair carries the PIP of the locus's variants in the peak, as capped for chromVAR. Each pair's PIP is
#' split between the cell types and the uniform background in proportion to their likelihoods, so a locus's shares are
#' its own evidence. For context, `fit_locus_attribution_weights()` also fits the trait's cell-type weights to all of
#' its pairs, with the prior weight that best predicts the pairs of held-out chromosomes; weights fitted to random peaks
#' of the same chromVAR background bins, in GC content and accessibility, give the share each cell type takes by chance.
#' A GWAS whose weights predict held-out pairs no better than even weights keeps its loci but gets missing weights.
#'
#' @param peak_variant_weight_tibble Peak-variant weights of one GWAS.
#' @param GWAS_input_record The GWAS's input record, for its credible sets' total PIP, lead variants and strongest
#'   association.
#' @param pseudobulk_ATAC_data_matrix BPCells peak-by-cell-type ATAC counts.
#' @param chromVAR_background_record Cell-type pseudobulk chromVAR background with the peaks' bins.
#' @param support_tibble Nuclei per cell type.
#' @param n_null Random peak sets whose weights give the chance shares.
#' @param pseudo_loci_grid Candidate prior weights.
#' @param n_folds Chromosome folds of the prior weight's cross-validation.
#' @return List with one-row-per-cell-type `trait_tibble`, one-row-per-locus `locus_tibble` and long `share_tibble`,
#'   or `NULL` without credible variants in peaks.
#' @keywords internal
get_locus_attribution_record <- function(peak_variant_weight_tibble, GWAS_input_record, pseudobulk_ATAC_data_matrix,
                                         chromVAR_background_record, support_tibble, n_null = 10L,
                                         pseudo_loci_grid = c(0.1, 1, 10, 100), n_folds = 5L) {
  if (NROW(peak_variant_weight_tibble) == 0L) return(NULL)
  GWAS_ID <- peak_variant_weight_tibble$GWAS_ID[[1]]
  pair_tibble <- dplyr::summarise(peak_variant_weight_tibble, PIP = sum(peak_variant_weight),
    chromosome = peak_chromosome[[1]], .by = c(studyLocusId, peak_name))
  peaks <- match(pair_tibble$peak_name, chromVAR_background_record$peak_names)
  peak_bin <- chromVAR_background_record$background@peak2bin
  bin_peaks <- split(seq_along(peak_bin), peak_bin)
  # A seed of its own keeps the draws independent of the target's name.
  null_peaks <- withr::with_seed(1L, vapply(peaks, \(peak) {
    pool <- bin_peaks[[as.character(peak_bin[[peak]])]]
    pool[sample.int(length(pool), n_null, replace = TRUE)]
  }, integer(n_null))) |> matrix(nrow = n_null)
  used_peaks <- unique(c(peaks, null_peaks))
  likelihood_matrix <- get_locus_attribution_likelihood_matrix(pseudobulk_ATAC_data_matrix, used_peaks, support_tibble)
  pair_likelihood_matrix <- likelihood_matrix[match(peaks, used_peaks), , drop = FALSE]
  PIP <- pair_tibble$PIP

  fold <- match(pair_tibble$chromosome, unique(pair_tibble$chromosome)) %% n_folds
  held_out_log_likelihood <- vapply(pseudo_loci_grid, \(pseudo_loci) sum(vapply(unique(fold), \(held_out) {
    train <- fold != held_out
    weights <- fit_locus_attribution_weights(pair_likelihood_matrix[train, , drop = FALSE], PIP[train], pseudo_loci)
    sum(PIP[!train] * log(pair_likelihood_matrix[!train, , drop = FALSE] %*% weights))
  }, numeric(1))), numeric(1))
  pseudo_loci <- pseudo_loci_grid[[which.max(held_out_log_likelihood)]]
  # The weights count as estimable when they predict the pairs of held-out chromosomes at least one nat better than
  # even weights, which split each pair by its accessibility alone: positive evidence on the scale of Kass and Raftery
  # (1995). GWAS with few loci fall short. Even weights don't depend on the training chromosomes, so their held-out
  # likelihood is that of all pairs.
  weights_estimable <- max(held_out_log_likelihood) - sum(PIP * log(rowMeans(pair_likelihood_matrix))) >= 1
  weights <- fit_locus_attribution_weights(pair_likelihood_matrix, PIP, pseudo_loci)
  null_weight_matrix <- vapply(seq_len(n_null), \(draw) fit_locus_attribution_weights(
    likelihood_matrix[match(null_peaks[draw, ], used_peaks), , drop = FALSE], PIP, pseudo_loci
  ), numeric(length(weights)))

  # PIP mass of each locus per cell type and the background, in proportion to their likelihoods.
  share_matrix <- rowsum(pair_likelihood_matrix / rowSums(pair_likelihood_matrix) * PIP, pair_tibble$studyLocusId)
  PIP_in_peaks <- rowSums(share_matrix)

  credible_set_tibble <- tibble::as_tibble(as.data.frame(S4Vectors::mcols(GWAS_input_record$credible_set_GRanges))) |>
    dplyr::arrange(dplyr::desc(posteriorProbability), variantId) |>
    dplyr::summarise(credible_set_PIP = sum(posteriorProbability), lead_variantId = dplyr::first(variantId),
      lead_variant_PIP = dplyr::first(posteriorProbability), minus_log10_p = max(-(log10(pValueMantissa) + pValueExponent)),
      .by = studyLocusId)
  list(
    trait_tibble = tibble::tibble(GWAS_ID, cluster = names(weights),
      weight = if (weights_estimable) weights else NA_real_,
      matched_weight = if (weights_estimable) rowMeans(null_weight_matrix) else NA_real_,
      weights_estimable, pseudo_loci, n_loci = nrow(share_matrix), PIP_in_peaks = sum(PIP)),
    locus_tibble = tibble::tibble(GWAS_ID, studyLocusId = rownames(share_matrix), PIP_in_peaks,
      n_peaks = as.vector(table(pair_tibble$studyLocusId)[rownames(share_matrix)])) |>
      dplyr::left_join(credible_set_tibble, by = "studyLocusId", relationship = "one-to-one"),
    share_tibble = tibble::tibble(GWAS_ID, studyLocusId = rep(rownames(share_matrix), times = ncol(share_matrix)),
      cluster = rep(colnames(share_matrix), each = nrow(share_matrix)), share = as.vector(share_matrix / PIP_in_peaks))
  )
}

#' Get each locus's shares and folds per cell type
#'
#' A cell type's fold is its share among the cell types, leaving out the background, over an even share of 1 over the
#' number of cell types: 1 means the locus's accessibility does not favour the cell type.
#'
#' @param share_tibble Combined locus-by-cell-type shares, the background included.
#' @return The cell types' rows with their locus's background share and their fold.
#' @keywords internal
get_locus_attribution_cell_type_shares <- function(share_tibble) {
  n_cell_types <- dplyr::n_distinct(share_tibble$cluster) - 1L
  share_tibble |>
    dplyr::mutate(background_share = share[cluster == "background"], .by = c(GWAS_ID, studyLocusId)) |>
    dplyr::filter(cluster != "background") |>
    dplyr::mutate(fold = share / (1 - background_share) * n_cell_types)
}

#' Get each locus's shares and folds per cell-type compartment
#'
#' @inheritParams get_locus_attribution_cell_type_shares
#' @param compartments_patterns Optional named cell-type compartment patterns.
#' @return One row per locus and compartment with its number of cell types, its share of the locus's PIP in peaks and
#'   its fold, the mean of its cell types' folds.
#' @keywords internal
get_locus_attribution_compartment_shares <- function(share_tibble, compartments_patterns) {
  cell_type_tibble <- get_locus_attribution_cell_type_shares(share_tibble)
  clusters <- unique(cell_type_tibble$cluster)
  compartment_of <- rlang::set_names(as.character(get_cluster_compartments(clusters, compartments_patterns %||% list())), clusters)
  cell_type_tibble |>
    dplyr::mutate(compartment = compartment_of[cluster]) |>
    dplyr::summarise(n_cell_types = dplyr::n(), share = sum(share), fold = mean(fold), .by = c(GWAS_ID, studyLocusId, compartment))
}

#' Summarise each locus's attribution
#'
#' @param locus_tibble Combined per-locus records.
#' @param share_tibble Combined locus-by-cell-type shares.
#' @param trait_tibble Combined per-trait cell-type weights, for each trait's main cell type.
#' @param locus_to_gene_tibble Ranked Open Targets L2G predictions.
#' @param GWAS_inputs_tibble GWAS metadata, for categories.
#' @param compartments_patterns Optional named cell-type compartment patterns.
#' @return One row per locus with its association, likely genes, share in background-like peaks, top two cell types
#'   by share and top compartment by fold, with their shares and folds, and its trait's main cell type where the
#'   trait's weights are estimable.
#' @keywords internal
summarise_locus_attribution <- function(locus_tibble, share_tibble, trait_tibble, locus_to_gene_tibble, GWAS_inputs_tibble,
                                        compartments_patterns) {
  top_tibble <- get_locus_attribution_cell_type_shares(share_tibble) |>
    dplyr::slice_max(share, n = 2L, with_ties = FALSE, by = c(GWAS_ID, studyLocusId)) |>
    dplyr::mutate(rank = c("top", "second")[dplyr::row_number()], .by = c(GWAS_ID, studyLocusId)) |>
    tidyr::pivot_wider(id_cols = c(GWAS_ID, studyLocusId, background_share), names_from = rank,
      values_from = c(cluster, share, fold), names_glue = "{rank}_{.value}") |>
    dplyr::rename(top_cell_type = top_cluster, second_cell_type = second_cluster)
  # The top compartment is the one the locus favours most for its size, by fold, so that large compartments do not
  # win by their number of cell types alone.
  compartment_tibble <- get_locus_attribution_compartment_shares(share_tibble, compartments_patterns) |>
    dplyr::slice_max(fold, n = 1L, with_ties = FALSE, by = c(GWAS_ID, studyLocusId)) |>
    dplyr::select(GWAS_ID, studyLocusId, top_compartment = compartment, top_compartment_share = share,
      top_compartment_fold = fold)
  trait_cell_type_tibble <- dplyr::filter(trait_tibble, cluster != "background", weights_estimable) |>
    dplyr::slice_max(weight, n = 1L, with_ties = FALSE, by = GWAS_ID) |>
    dplyr::select(GWAS_ID, trait_cell_type = cluster)
  gene_tibble <- locus_to_gene_tibble |>
    dplyr::filter(L2G_rank <= 3L) |>
    dplyr::summarise(top_L2G_gene = dplyr::first(gene_symbol),
      likely_genes = paste0(gene_symbol, " (", sprintf("%.2f", L2G_score), ")", collapse = ", "), .by = c(GWAS_ID, studyLocusId))
  locus_tibble |>
    dplyr::left_join(dplyr::select(GWAS_inputs_tibble, GWAS_ID, Category), by = "GWAS_ID") |>
    dplyr::left_join(gene_tibble, by = c("GWAS_ID", "studyLocusId")) |>
    dplyr::left_join(top_tibble, by = c("GWAS_ID", "studyLocusId")) |>
    dplyr::left_join(compartment_tibble, by = c("GWAS_ID", "studyLocusId")) |>
    dplyr::left_join(trait_cell_type_tibble, by = "GWAS_ID") |>
    dplyr::relocate(GWAS_ID, Category, studyLocusId, lead_variantId, minus_log10_p, lead_variant_PIP, top_L2G_gene, likely_genes)
}

#' Summarise each locus's attribution by cell-type compartment
#'
#' @param locus_summary_tibble Per-locus summary from `summarise_locus_attribution()`.
#' @param share_tibble Combined locus-by-cell-type shares.
#' @param trait_tibble Combined per-trait cell-type weights, for each compartment's weight in the trait.
#' @param compartments_patterns Optional named cell-type compartment patterns.
#' @return One row per locus and compartment with the locus's identifiers and association, the compartment's number of
#'   cell types, share and fold, and its weight in the trait over all loci and by chance.
#' @keywords internal
summarise_locus_attribution_compartments <- function(locus_summary_tibble, share_tibble, trait_tibble, compartments_patterns) {
  trait_tibble <- dplyr::filter(trait_tibble, cluster != "background")
  clusters <- unique(trait_tibble$cluster)
  compartment_of <- rlang::set_names(as.character(get_cluster_compartments(clusters, compartments_patterns %||% list())), clusters)
  trait_weight_tibble <- trait_tibble |>
    dplyr::mutate(compartment = compartment_of[cluster]) |>
    dplyr::summarise(trait_weight = sum(weight), trait_chance_weight = sum(matched_weight), .by = c(GWAS_ID, compartment))
  get_locus_attribution_compartment_shares(share_tibble, compartments_patterns) |>
    dplyr::left_join(trait_weight_tibble, by = c("GWAS_ID", "compartment")) |>
    dplyr::left_join(dplyr::select(locus_summary_tibble, GWAS_ID, studyLocusId, lead_variantId, minus_log10_p, lead_variant_PIP,
      PIP_in_peaks, credible_set_PIP, top_L2G_gene), by = c("GWAS_ID", "studyLocusId")) |>
    dplyr::relocate(GWAS_ID, studyLocusId, lead_variantId, top_L2G_gene, minus_log10_p, lead_variant_PIP, PIP_in_peaks,
      credible_set_PIP, compartment, n_cell_types)
}

#' Plot each trait's cell-type weights against chance
#'
#' @param trait_tibble Combined per-trait cell-type weights with their chance shares.
#' @param GWAS_metadata_tracks_plot Shared GWAS metadata tracks.
#' @param support_tibble Nuclei per cell type, for the support track.
#' @param compartments_patterns Optional named cell-type compartment patterns.
#' @return Annotated GWAS-by-cell-type heatmap.
#' @keywords internal
plot_locus_attribution_trait_heatmap <- function(trait_tibble, GWAS_metadata_tracks_plot, support_tibble, compartments_patterns) {
  # GWAS without estimable weights keep the grey rows of the metadata tracks' GWAS without scores.
  heatmap_tibble <- trait_tibble |>
    dplyr::filter(cluster != "background", weights_estimable) |>
    dplyr::left_join(get_GWAS_track_tibble(GWAS_metadata_tracks_plot), by = "GWAS_ID") |>
    dplyr::left_join(dplyr::select(support_tibble, cluster, n_cells), by = "cluster") |>
    dplyr::mutate(excess = 100 * (weight - matched_weight), share_label = dplyr::if_else(weight >= 0.05, as.character(round(100 * weight)), ""))
  if (nrow(heatmap_tibble) == 0L) return(structure(list(), class = c("empty_plot_list", "list")))
  # A few large shares, such as microglia's in Alzheimer's disease, would otherwise wash out all others.
  colour_limit <- max(10, ceiling(stats::quantile(abs(heatmap_tibble$excess), 0.99, names = FALSE) / 10) * 10)
  heatmap_tibble$excess <- pmax(pmin(heatmap_tibble$excess, colour_limit), -colour_limit)
  plot_GWAS_by_cluster_heatmap(heatmap_tibble, GWAS_metadata_tracks_plot, compartments_patterns, fill_col = "excess",
    fill_label = "Share above chance\n(percentage points)", fill_limits = c(-colour_limit, colour_limit),
    support_label_col = "share_label") +
    patchwork::plot_annotation(
      title = "Cell types that each trait's fine-mapped signal in accessible peaks points to",
      subtitle = paste(
        "Red marks cell types that take more of the trait's credible-set PIP in peaks than they take of random peaks with",
        "the same GC content and accessibility; numbers give each cell type's share in percent where it is at least 5%.",
        "They summarise the trait; each locus's own shares, in the locus plots and tables, do not use them.",
        sep = "\n"
      ),
      caption = stringr::str_wrap(paste(
        "A mixture over cell types explains each trait's locus-peak pairs, weighted by the PIP of the credible variants in",
        "the peak: a cell type takes a peak in proportion to the peak's share of its pseudobulk ATAC reads, shrunk toward",
        "the pooled profile by the reads of 50 average nuclei. Expectation-maximisation fits the weights with a Dirichlet",
        "prior of a few pseudo-loci spread over all cell types, choosing among 0.1, 1, 10 and 100 by the likelihood of",
        "held-out chromosomes (five folds). Chance shares are the mean weights fitted to ten sets of random peaks drawn from",
        "the chromVAR background bins of the trait's peaks. A uniform background component, not shown, takes peaks that no",
        "cell type explains better than evenly spread reads. Grey rows mark GWAS without credible variants in peaks, or with",
        "too few loci to estimate weights: weights that predict held-out chromosomes less than one nat better than even",
        "weights. Shares describe where the fine-mapped signal in peaks lies, not significance; S-LDSC and chromVAR give that."
      ), width = 150),
      theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0))
    )
}

#' Plot the cell types whose accessibility each of a GWAS's best fine-mapped loci favours
#'
#' @param locus_summary_tibble Per-locus summary from `summarise_locus_attribution()`.
#' @param share_tibble Locus-by-cell-type shares.
#' @param trait_tibble Per-trait cell-type weights, for the comparison row.
#' @param compartments_patterns Optional named cell-type compartment patterns.
#' @param n_loci Loci per GWAS, those with the most PIP in peaks.
#' @return Named list of locus-by-cell-type heatmaps of folds, one per GWAS.
#' @keywords internal
plot_locus_attribution_heatmaps <- function(locus_summary_tibble, share_tibble, trait_tibble, compartments_patterns, n_loci = 30L) {
  cell_type_tibble <- get_locus_attribution_cell_type_shares(share_tibble)
  clusters <- unique(cell_type_tibble$cluster)
  compartment_of <- rlang::set_names(as.character(get_cluster_compartments(clusters, compartments_patterns %||% list())), clusters)
  compartment_levels <- intersect(c(names(compartments_patterns), "Other"), compartment_of)
  cluster_levels <- clusters[order(match(compartment_of[clusters], compartment_levels), clusters)]
  trait_row <- "Trait: weights over chance"
  GWAS_IDs <- unique(locus_summary_tibble$GWAS_ID)
  plots <- purrr::map(GWAS_IDs, \(GWAS_ID) {
    loci <- locus_summary_tibble |>
      dplyr::filter(.data$GWAS_ID == .env$GWAS_ID) |>
      dplyr::slice_max(dplyr::tibble(PIP_in_peaks, top_compartment_fold), n = n_loci, with_ties = FALSE) |>
      dplyr::arrange(match(top_compartment, compartment_levels), dplyr::desc(top_compartment_fold)) |>
      dplyr::mutate(row = paste0(dplyr::coalesce(top_L2G_gene, "no L2G gene"), "  ", lead_variantId),
        row = dplyr::if_else(duplicated(row) | duplicated(row, fromLast = TRUE), paste0(row, " (", substr(studyLocusId, 1, 6), ")"), row))
    tile_tibble <- dplyr::bind_rows(
      cell_type_tibble |>
        dplyr::inner_join(dplyr::select(loci, GWAS_ID, studyLocusId, row, group = top_compartment), by = c("GWAS_ID", "studyLocusId")) |>
        dplyr::select(row, group, cluster, fold),
      trait_tibble |>
        dplyr::filter(.data$GWAS_ID == .env$GWAS_ID, cluster != "background", weights_estimable) |>
        dplyr::transmute(row = trait_row, group = "Trait", cluster, fold = weight / matched_weight)
    ) |>
      dplyr::mutate(
        log2_fold = pmax(pmin(log2(fold), 3), -3),
        row = factor(row, levels = rev(c(trait_row, loci$row))),
        group = factor(group, levels = c("Trait", compartment_levels)),
        compartment = factor(compartment_of[cluster], levels = compartment_levels),
        cluster = factor(cluster, levels = cluster_levels)
      )
    n_total <- sum(locus_summary_tibble$GWAS_ID == GWAS_ID)
    # The subtitle's and the caption's notes on the trait row.
    trait_notes <- if (trait_row %in% tile_tibble$row) c(
      "The top row shows the trait's cell-type weights over all its loci against random matched peaks, for comparison.",
      paste("The trait row divides each cell type's weight in a mixture fitted to all the trait's loci by its weight for",
        "random peaks of the same GC content and accessibility.")
    ) else c(
      "The trait has too few loci to estimate its cell-type weights over all its loci.",
      paste("A mixture fitted to all the trait's loci predicts loci on held-out chromosomes less than one nat better than",
        "even weights, so its weights are not shown.")
    )
    ggplot2::ggplot(tile_tibble, ggplot2::aes(x = cluster, y = row, fill = log2_fold)) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
      ggplot2::facet_grid(group ~ compartment, scales = "free", space = "free") +
      ggplot2::scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0, limits = c(-3, 3),
        breaks = -3:3, labels = c("1/8", "1/4", "1/2", "1", "2", "4", "≥8"), name = "Fold") +
      ggplot2::labs(
        x = NULL, y = NULL,
        title = paste0("Cell types whose accessibility ", GWAS_ID, "'s loci favour"),
        subtitle = stringr::str_wrap(paste0(
          "Red marks cell types that take more of a locus's PIP in peaks than an even share of the ", length(clusters),
          " cell types, from the locus's own accessibility; loci are grouped, on the right, by the compartment they favour most. ",
          trait_notes[[1]]
        ), width = 130),
        caption = stringr::str_wrap(paste0(
          "Shown: the ", nrow(loci), " of ", n_total, " credible sets with credible variants (PIP > cutoff) in peaks that have the ",
          "most PIP in peaks; rows give the top Open Targets L2G gene and the lead variant. Within a credible set, each peak's PIP ",
          "is split between the cell types in proportion to the peak's share of each cell type's pseudobulk ATAC reads, and to a ",
          "uniform background for peaks no cell type explains better. A cell type's fold is its share among the cell types over ",
          "an even share; a compartment's is the mean of its cell types'. Folds assume that the causal variant acts through ",
          "accessible chromatin and do not test single loci, since a random peak is usually specific to some cell type too. ",
          trait_notes[[2]], " Colours stop at 8-fold either way."
        ), width = 150)
      ) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid = ggplot2::element_blank(),
        strip.text.x = ggplot2::element_text(angle = 90, hjust = 0), strip.text.y = ggplot2::element_text(angle = 0, hjust = 0),
        panel.spacing = grid::unit(2, "pt"), legend.position = "right",
        plot.title.position = "plot", plot.caption = ggplot2::element_text(hjust = 0)
      )
  })
  rlang::set_names(plots, GWAS_IDs)
}

#' Save a locus attribution table
#'
#' @param table Per-locus or per-locus-and-compartment summary.
#' @return Path of the written TSV.
#' @keywords internal
save_locus_attribution_table <- function(table) {
  path <- get_structured_file_path(filetype = "tsv")
  readr::write_tsv(table, path)
  path
}
