#' Binarise cell-type scores at matched shares of common SNPs
#'
#' Rank each cell type's autosomal territory by its score and keep the
#' highest-scoring ranges until they hold each share of the common (MAF > 5%)
#' reference SNPs. A SNP in overlapping ranges takes their largest score, as in
#' `get_SLDSC_LD_score_record()`, so the ranges are first split into disjoint
#' pieces. Pieces are kept whole, in genomic order among ties, so an annotation
#' can exceed its share by one piece; pieces without score are never kept. A
#' `random` annotation ranks the same ranges at random, as a reference for what
#' the territory alone concentrates.
#'
#' @param score_matrix Range-by-cell-type scores named by range.
#' @param ranges Ranges named like the score rows.
#' @param reference_dir Directory returned by `download_SLDSC_reference()`.
#' @param shares Shares of the common reference SNPs.
#' @return A list with the kept pieces, their binary piece-by-annotation matrix
#'   with annotations named `<cell type>:<share>`, and each annotation's sums of
#'   the baseline-LD annotations over its common SNPs, named like the baseline
#'   LD scores; its `baseL2` sum counts these SNPs.
#' @keywords internal

get_SLDSC_matched_annotation_record <- function(score_matrix, ranges, reference_dir, shares = c(0.005, 0.01, 0.02, 0.05)) {
  score_matrix <- cbind(score_matrix, random = withr::with_seed(1L, stats::runif(nrow(score_matrix))))
  pieces <- GenomicRanges::disjoin(ranges[rownames(score_matrix)], with.revmap = TRUE)
  pieces <- pieces[as.character(GenomicRanges::seqnames(pieces)) %in% paste0("chr", 1:22)]
  score_matrix <- apply(score_matrix, 2, \(values) max(IRanges::extractList(values, pieces$revmap)))
  chromosome_records <- purrr::map(1:22, \(chromosome) {
    baseline_tibble <- data.table::fread(
      file.path(reference_dir, "baselineLD_v2.2", paste0("baselineLD.", chromosome, ".annot.gz")),
      drop = c("CHR", "SNP", "CM")
    )
    common <- data.table::fread(file.path(reference_dir, "plink_files", paste0("1000G.EUR.hg38.", chromosome, ".frq")))$MAF > 0.05
    piece_idx <- GenomicRanges::findOverlaps(
      GenomicRanges::GRanges(paste0("chr", chromosome), IRanges::IRanges(baseline_tibble$BP[common], width = 1L)),
      pieces,
      select = "first"
    )
    baseline_matrix <- as.matrix(baseline_tibble[common, -"BP"])
    list(
      piece_sums = rowsum(baseline_matrix[!is.na(piece_idx), , drop = FALSE], piece_idx[!is.na(piece_idx)]),
      n_common = sum(common)
    )
  })
  piece_sums <- do.call(rbind, purrr::map(chromosome_records, "piece_sums"))
  piece_baseline_matrix <- matrix(0, length(pieces), ncol(piece_sums), dimnames = list(NULL, paste0(colnames(piece_sums), "L2")))
  piece_baseline_matrix[as.integer(rownames(piece_sums)), ] <- piece_sums
  n_common <- sum(purrr::map_int(chromosome_records, "n_common"))
  piece_n_common <- piece_baseline_matrix[, "baseL2"]
  annotation_matrix <- do.call(cbind, purrr::map(colnames(score_matrix), \(cluster) {
    ranked <- order(-score_matrix[, cluster])
    ranked <- ranked[score_matrix[ranked, cluster] > 0]
    SNPs_before <- cumsum(piece_n_common[ranked]) - piece_n_common[ranked]
    vapply(shares, \(share) seq_along(pieces) %in% ranked[SNPs_before < share * n_common], logical(length(pieces)))
  })) * 1
  dimnames(annotation_matrix) <- list(as.character(pieces), paste0(rep(colnames(score_matrix), each = length(shares)), ":", shares))
  kept <- rowSums(annotation_matrix) > 0
  list(
    ranges = stats::setNames(GenomicRanges::granges(pieces), as.character(pieces))[kept],
    annotation_matrix = annotation_matrix[kept, , drop = FALSE],
    baseline_overlap_matrix = crossprod(annotation_matrix, piece_baseline_matrix)
  )
}

#' Estimate S-LDSC heritability enrichment of candidate annotations
#'
#' Fit each candidate annotation with the shared annotations, which also set
#' the regression weights, as `fit_SLDSC_candidate_regressions()` does. A
#' candidate's share of heritability sums the fitted per-SNP heritability of
#' the whole model over its common SNPs, as ldsc's `--overlap-annot` Prop._h2
#' does, and its standard error comes from ldsc's ratio jackknife.
#'
#' @inheritParams fit_SLDSC_candidate_regressions
#' @param overlap_matrix Candidate-by-shared-annotation sums over each
#'   candidate's common SNPs.
#' @param candidates Binary candidate annotations, one per model.
#' @return One row per candidate with its coefficient, its share of the
#'   model's heritability and the model's total heritability, each with a
#'   jackknife standard error.
#' @keywords internal

fit_SLDSC_enrichment_regressions <- function(chi2, N, weight_L2, L2_matrix, M, overlap_matrix, shared_annotations, candidates, n_blocks = 200L) {
  annotations <- c(shared_annotations, candidates)
  blocks <- get_SLDSC_block_crossproducts(chi2, N, weight_L2, L2_matrix, M, shared_annotations, annotations, n_blocks)
  n_shared <- length(shared_annotations)
  purrr::map_dfr(seq_along(candidates), \(i) {
    model_cols <- c(seq_len(n_shared), n_shared + i, length(annotations) + 1L)
    fit <- jackknife_LDSC_blocks(
      blocks$block_crossproducts[model_cols, model_cols, , drop = FALSE],
      blocks$block_crossproducts[model_cols, length(annotations) + 2L, ]
    )
    annotation_idx <- seq_len(n_shared + 1L)
    model_M <- M[c(shared_annotations, candidates[[i]])]
    # Columns: the full data, then each deleted block.
    coefficients <- cbind(fit$estimate, fit$delete_values)[annotation_idx, ]
    h2 <- drop(model_M %*% coefficients)
    h2_share <- drop(c(overlap_matrix[candidates[[i]], shared_annotations], model_M[[n_shared + 1L]]) %*% coefficients) / h2
    pseudovalues <- n_blocks * h2_share[[1]] - (n_blocks - 1) * h2_share[-1]
    tibble::tibble(
      annotation = candidates[[i]],
      tau = fit$estimate[[n_shared + 1L]] / blocks$N_bar,
      tau_se = sqrt(fit$cov[n_shared + 1L, n_shared + 1L]) / blocks$N_bar,
      h2_share = h2_share[[1]],
      h2_share_se = stats::sd(pseudovalues) / sqrt(n_blocks),
      h2 = h2[[1]] / blocks$N_bar,
      h2_se = sqrt(drop(model_M %*% fit$cov[annotation_idx, annotation_idx] %*% model_M)) / blocks$N_bar
    )
  })
}

#' Measure the fine-mapped posterior probability inside range annotations
#'
#' Sum the posterior inclusion probabilities (PIPs) of autosomal credible-set
#' variants inside each annotation and divide by their total. Variants are
#' assigned by position, so indels and variants missing from the reference
#' count. Loci, the credible sets whose centres lie within 500 kb of each
#' other, are bootstrapped for uncertainty.
#'
#' @param credible_set_GRanges Credible-set variants with `studyLocusId` and
#'   `posteriorProbability`.
#' @param ranges Disjoint annotation ranges.
#' @param annotation_matrix Binary range-by-annotation matrix.
#' @param n_bootstraps Bootstrap samples of loci.
#' @return One row per annotation with its share of the PIP, the bootstrap
#'   standard error and 95% interval of that share, and the numbers of credible
#'   sets and loci.
#' @keywords internal

get_SLDSC_PIP_share_tibble <- function(credible_set_GRanges, ranges, annotation_matrix, n_bootstraps = 1000L) {
  variants <- credible_set_GRanges[as.character(GenomicRanges::seqnames(credible_set_GRanges)) %in% paste0("chr", 1:22)]
  set_tibble <- tibble::tibble(
    studyLocusId = variants$studyLocusId,
    chromosome = as.character(GenomicRanges::seqnames(variants)),
    position = GenomicRanges::start(variants)
  ) |>
    dplyr::summarise(chromosome = dplyr::first(chromosome), centre = (min(position) + max(position)) %/% 2L, .by = studyLocusId)
  set_centres <- GenomicRanges::GRanges(set_tibble$chromosome, IRanges::IRanges(set_tibble$centre, width = 1L))
  loci <- GenomicRanges::reduce(GenomicRanges::GRanges(
    set_tibble$chromosome,
    IRanges::IRanges(pmax(set_tibble$centre - 250000L, 1L), set_tibble$centre + 250000L)
  ))
  variant_locus <- GenomicRanges::findOverlaps(set_centres, loci, select = "first")[match(variants$studyLocusId, set_tibble$studyLocusId)]
  piece_idx <- GenomicRanges::findOverlaps(variants, ranges, select = "first")
  inside_matrix <- annotation_matrix[piece_idx, , drop = FALSE]
  inside_matrix[is.na(piece_idx), ] <- 0
  locus_PIP_matrix <- rowsum(variants$posteriorProbability * inside_matrix, variant_locus)
  locus_PIP <- rowsum(variants$posteriorProbability, variant_locus)[, 1]
  withr::local_seed(1L)
  bootstrap_weights <- t(stats::rmultinom(n_bootstraps, length(locus_PIP), rep(1, length(locus_PIP))))
  bootstrap_shares <- (bootstrap_weights %*% locus_PIP_matrix) / drop(bootstrap_weights %*% locus_PIP)
  tibble::tibble(
    annotation = colnames(annotation_matrix),
    PIP_share = colSums(locus_PIP_matrix) / sum(locus_PIP),
    PIP_share_se = apply(bootstrap_shares, 2, stats::sd),
    PIP_share_lower = apply(bootstrap_shares, 2, stats::quantile, 0.025, names = FALSE),
    PIP_share_upper = apply(bootstrap_shares, 2, stats::quantile, 0.975, names = FALSE),
    n_credible_sets = nrow(set_tibble),
    n_loci = length(loci)
  )
}

#' Evaluate one GWAS's enrichment in matched-size cell-type annotations
#'
#' Divide each annotation's share of S-LDSC heritability, fitted with the
#' baseline-LD model, and its share of fine-mapped PIP by its share of common
#' reference SNPs.
#'
#' @param GWAS_input_record One normalized GWAS input record.
#' @param sumstats_record The record of the same source from
#'   `get_SLDSC_sumstats_record()`.
#' @param baseline_record Record from `read_SLDSC_baseline_record()`.
#' @param annotation_record Record from `get_SLDSC_matched_annotation_record()`.
#' @param LD_score_record Record from `get_SLDSC_LD_score_record()` of its
#'   annotations.
#' @return One row per annotation with its cell type, target and actual share
#'   of common SNPs, both enrichments with their standard errors, the PIP
#'   enrichment's bootstrap interval, and the S-LDSC coefficient and total
#'   heritability; S-LDSC estimates are missing without summary statistics.
#' @keywords internal

get_SLDSC_matched_enrichment_tibble <- function(GWAS_input_record, sumstats_record, baseline_record, annotation_record, LD_score_record) {
  if (!identical(sumstats_record$sourceId, GWAS_input_record$sourceId)) {
    stop("S-LDSC summary statistics of ", sumstats_record$sourceId, " were paired with ", GWAS_input_record$GWAS_ID, ".")
  }
  M <- c(baseline_record$M_5_50, annotation_record$baseline_overlap_matrix[, "baseL2"])
  enrichment_tibble <- get_SLDSC_PIP_share_tibble(GWAS_input_record$credible_set_GRanges, annotation_record$ranges, annotation_record$annotation_matrix)
  sumstats_tibble <- sumstats_record$sumstats_tibble
  if (!is.null(sumstats_tibble)) {
    kept <- baseline_record$SNP_tibble$SNP %in% sumstats_tibble$SNP
    sumstats_idx <- match(baseline_record$SNP_tibble$SNP[kept], sumstats_tibble$SNP)
    SLDSC_tibble <- fit_SLDSC_enrichment_regressions(
      chi2 = sumstats_tibble$chi2[sumstats_idx],
      N = sumstats_tibble$N[sumstats_idx],
      weight_L2 = baseline_record$SNP_tibble$weight_L2[kept],
      L2_matrix = cbind(baseline_record$L2_matrix[kept, ], LD_score_record$L2_matrix[kept, ]),
      M = M,
      overlap_matrix = annotation_record$baseline_overlap_matrix,
      shared_annotations = colnames(baseline_record$L2_matrix),
      candidates = colnames(LD_score_record$L2_matrix)
    )
    enrichment_tibble <- dplyr::left_join(enrichment_tibble, SLDSC_tibble, by = "annotation", relationship = "one-to-one")
  }
  SNP_share <- M[enrichment_tibble$annotation] / baseline_record$M_5_50[["baseL2"]]
  enrichment_tibble |>
    dplyr::mutate(
      GWAS_ID = GWAS_input_record$GWAS_ID,
      finemappingMethod = GWAS_input_record$finemappingMethod,
      SNP_share = unname(SNP_share),
      dplyr::across(dplyr::any_of(c("h2_share", "h2_share_se")), \(x) x / SNP_share, .names = "SLDSC_{stringr::str_replace(.col, 'h2_share', 'enrichment')}"),
      dplyr::across(dplyr::starts_with("PIP_share"), \(x) x / SNP_share, .names = "{stringr::str_replace(.col, 'PIP_share', 'PIP_enrichment')}")
    ) |>
    tidyr::separate_wider_regex(annotation, c(cluster = ".*", ":", share = "[^:]+$")) |>
    dplyr::mutate(share = as.numeric(share)) |>
    dplyr::relocate(GWAS_ID, cluster, share, SNP_share, dplyr::starts_with("SLDSC_enrichment"), dplyr::starts_with("PIP_enrichment"))
}

#' Summarise matched-size enrichment across GWAS
#'
#' Take the median enrichment over GWAS-cell-type pairs of each strategy,
#' share and statistic, for all cell types together and for each one, among
#' the pairs that the single-modality S-LDSC models enrich and among all
#' pairs, and of the random reference over GWAS. Specificity is each GWAS's
#' median enrichment in its enriched cell types over that in its other cell
#' types, when both are positive, summarised by its median over GWAS.
#' Intervals are 2.5% and 97.5%
#' quantiles of the medians over bootstrap samples of GWAS.
#'
#' @param enrichment_tibble Combined matched-size enrichment of all GWAS.
#' @param n_bootstraps Bootstrap samples of GWAS.
#' @return One row per strategy, share, statistic, pair set and cell type
#'   (`All cell types` for all together) with the number of pairs and GWAS,
#'   the median and its interval.
#' @keywords internal

summarise_SLDSC_matched_enrichment <- function(enrichment_tibble, n_bootstraps = 1000L) {
  withr::local_seed(1L)
  bootstrap_median <- \(value, GWAS_ID) {
    GWAS_rows <- split(seq_along(value), GWAS_ID)
    medians <- replicate(n_bootstraps, stats::median(value[unlist(sample(GWAS_rows, replace = TRUE), use.names = FALSE)]))
    tibble::tibble(
      n_pairs = length(value), n_GWAS = length(GWAS_rows), median = stats::median(value),
      lower = stats::quantile(medians, 0.025, names = FALSE), upper = stats::quantile(medians, 0.975, names = FALSE)
    )
  }
  pair_tibble <- enrichment_tibble |>
    tidyr::pivot_longer(c(SLDSC_enrichment, PIP_enrichment), names_to = "statistic", values_to = "enrichment") |>
    dplyr::filter(!is.na(enrichment)) |>
    dplyr::mutate(statistic = factor(statistic, c("SLDSC_enrichment", "PIP_enrichment"), c("S-LDSC heritability", "Fine-mapped PIP")))
  pair_sets <- list(
    "Enriched pairs" = dplyr::filter(pair_tibble, enriched_in_models %in% TRUE),
    "All pairs" = dplyr::filter(pair_tibble, cluster != "random"),
    "Random reference" = dplyr::filter(pair_tibble, cluster == "random")
  )
  enrichment_summary <- purrr::imap_dfr(pair_sets, \(pairs, pair_set) {
    dplyr::bind_rows(dplyr::mutate(pairs, cluster = "All cell types"), if (pair_set != "Random reference") pairs) |>
      dplyr::reframe(bootstrap_median(enrichment, GWAS_ID), .by = c(strategy, share, statistic, cluster)) |>
      dplyr::mutate(pair_set = pair_set)
  })
  specificity_summary <- pair_tibble |>
    dplyr::filter(cluster != "random") |>
    dplyr::summarise(
      enriched = stats::median(enrichment[enriched_in_models %in% TRUE]),
      other = stats::median(enrichment[!enriched_in_models %in% TRUE]),
      .by = c(strategy, share, statistic, GWAS_ID)
    ) |>
    dplyr::filter(enriched > 0, other > 0) |>
    dplyr::mutate(enrichment = enriched / other) |>
    dplyr::reframe(bootstrap_median(enrichment, GWAS_ID), .by = c(strategy, share, statistic)) |>
    dplyr::mutate(cluster = "All cell types", pair_set = "Specificity")
  dplyr::bind_rows(enrichment_summary, specificity_summary) |>
    dplyr::relocate(strategy, share, statistic, pair_set, cluster)
}

# Shared strategy colour scale of the matched-size plots.
SLDSC_strategy_scale <- \(aesthetics = "colour") {
  ggplot2::scale_colour_manual(
    values = c(peaks = "#D95F02", genes = "#1B9E77"),
    labels = c(peaks = "Peaks by accessibility share", genes = "Gene windows by expression specificity"),
    name = "Strategy",
    aesthetics = aesthetics
  )
}

#' Plot matched-size enrichment against annotation size
#'
#' @param summary_tibble Summary from `summarise_SLDSC_matched_enrichment()`.
#' @return Patchwork of median enrichment and specificity curves.
#' @keywords internal

plot_SLDSC_matched_enrichment_curves <- function(summary_tibble) {
  summary_tibble <- dplyr::filter(summary_tibble, cluster == "All cell types")
  pair_sets <- c("Enriched pairs", "All pairs")
  curve_tibble <- dplyr::bind_rows(
    dplyr::filter(summary_tibble, pair_set %in% pair_sets) |> dplyr::mutate(annotation = "Cell-type annotation"),
    dplyr::filter(summary_tibble, pair_set == "Random reference") |>
      dplyr::select(-pair_set) |>
      tidyr::crossing(pair_set = pair_sets) |>
      dplyr::mutate(annotation = "Random reference")
  ) |>
    dplyr::mutate(pair_set = factor(pair_set, pair_sets))
  share_axis <- ggplot2::scale_x_log10(breaks = unique(summary_tibble$share), labels = scales::label_percent(accuracy = 0.1))
  curve_layers <- list(
    ggplot2::geom_hline(yintercept = 1, colour = "grey60", linewidth = 0.3),
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper, fill = strategy, group = interaction(strategy, annotation)), alpha = 0.15, colour = NA),
    ggplot2::geom_line(ggplot2::aes(colour = strategy, linetype = annotation)),
    ggplot2::geom_point(ggplot2::aes(colour = strategy), size = 1),
    share_axis,
    ggplot2::scale_y_log10(),
    SLDSC_strategy_scale(c("colour", "fill")),
    ggplot2::scale_linetype_manual(values = c("Cell-type annotation" = "solid", "Random reference" = "dashed"), name = NULL),
    ggplot2::theme_bw(base_size = 9)
  )
  enrichment_plot <- ggplot2::ggplot(curve_tibble, ggplot2::aes(share, median)) +
    curve_layers +
    ggplot2::facet_grid(ggplot2::vars(statistic), ggplot2::vars(pair_set), scales = "free_y") +
    ggplot2::labs(subtitle = "Median enrichment", x = "Annotation size (share of common SNPs)", y = "Enrichment (log scale)")
  specificity_plot <- dplyr::filter(summary_tibble, pair_set == "Specificity") |>
    dplyr::mutate(annotation = "Cell-type annotation") |>
    ggplot2::ggplot(ggplot2::aes(share, median)) +
    curve_layers +
    ggplot2::facet_grid(ggplot2::vars(statistic), scales = "free_y") +
    ggplot2::guides(linetype = "none") +
    ggplot2::labs(subtitle = "Specificity: enriched over other cell types", x = "Annotation size (share of common SNPs)", y = "Ratio of median enrichments (log scale)")
  n_pairs <- \(pair_set_name, statistic_name) max(summary_tibble$n_pairs[summary_tibble$pair_set == pair_set_name & summary_tibble$statistic == statistic_name])
  n_GWAS <- \(pair_set_name, statistic_name) max(summary_tibble$n_GWAS[summary_tibble$pair_set == pair_set_name & summary_tibble$statistic == statistic_name])
  patchwork::wrap_plots(enrichment_plot, specificity_plot, widths = c(2, 1), guides = "collect") +
    patchwork::plot_annotation(
      title = "Heritability and fine-mapped variant enrichment of cell-type annotations at matched sizes",
      subtitle = paste(
        "Compare strategies at equal size: a higher curve concentrates more of the signal in the same share of the genome.",
        "Dashed curves rank the same territory at random; the gap to them is what cell-type ranking adds.",
        sep = "\n"
      ),
      caption = stringr::str_wrap(paste(
        "Each annotation keeps a cell type's highest-scoring peaks (accessibility share) or gene windows (expression specificity, gene body +/- 100 kb) until they hold the given share of common (MAF > 5%) 1000 Genomes European SNPs.",
        "S-LDSC heritability: the annotation's share of heritability in a model with the 97 baseline-LD v2.2 annotations, over its share of common SNPs; SNPs with chi-square above 80 are excluded.",
        "Fine-mapped PIP: the annotation's share of the posterior inclusion probability of autosomal credible-set variants, over its share of common SNPs.",
        sprintf(
          "Enriched pairs: GWAS-cell-type pairs with q < 0.05 in the single-modality peak or gene S-LDSC model (%d for S-LDSC, %d for PIP); all pairs: every GWAS and cell type (%d and %d).",
          n_pairs("Enriched pairs", "S-LDSC heritability"), n_pairs("Enriched pairs", "Fine-mapped PIP"),
          n_pairs("All pairs", "S-LDSC heritability"), n_pairs("All pairs", "Fine-mapped PIP")
        ),
        "Most enriched pairs come from the peak model on the same summary statistics, which favours peaks in S-LDSC.",
        sprintf("Specificity: each GWAS's median enrichment in its enriched cell types over that in its other cell types, for the up to %d GWAS with an enriched pair and both medians positive.", n_GWAS("Specificity", "S-LDSC heritability")),
        "PIP enrichment rests on fine-mapping resolution: wide gene windows hold whole credible sets, while small peaks hold only the variants inside them.",
        "Lines: medians over pairs or GWAS; bands: 95% intervals from 1,000 bootstrap samples of GWAS."
      ), width = 180),
      theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0))
    ) &
    ggplot2::theme(legend.position = "bottom")
}

#' Compare annotation strategies per cell type at one matched size
#'
#' @param summary_tibble Summary from `summarise_SLDSC_matched_enrichment()`.
#' @param enrichment_tibble Combined matched-size enrichment of all GWAS.
#' @param share Annotation size to compare, as a share of common SNPs.
#' @return Dot plot of each cell type's median enrichment per strategy.
#' @keywords internal

plot_SLDSC_matched_size_comparison <- function(summary_tibble, enrichment_tibble, share = 0.01) {
  matched_share <- share
  summary_tibble <- dplyr::filter(summary_tibble, share == matched_share)
  pair_sets <- c("Enriched pairs", "All pairs")
  cell_type_tibble <- summary_tibble |>
    dplyr::filter(pair_set %in% pair_sets, cluster != "All cell types") |>
    dplyr::mutate(pair_set = factor(pair_set, pair_sets))
  cluster_order <- cell_type_tibble |>
    dplyr::filter(pair_set == "All pairs", statistic == "S-LDSC heritability") |>
    dplyr::summarise(median = max(median), .by = cluster) |>
    dplyr::arrange(median) |>
    dplyr::pull(cluster)
  reference_tibble <- summary_tibble |>
    dplyr::filter(pair_set == "Random reference") |>
    dplyr::select(-pair_set) |>
    tidyr::crossing(pair_set = factor(pair_sets, pair_sets))
  # The strategy with the higher enrichment in each pair, counted per statistic and pair set.
  best_tibble <- enrichment_tibble |>
    dplyr::filter(share == matched_share, cluster != "random") |>
    tidyr::pivot_longer(c(SLDSC_enrichment, PIP_enrichment), names_to = "statistic", values_to = "enrichment") |>
    dplyr::filter(!is.na(enrichment)) |>
    dplyr::summarise(best = strategy[which.max(enrichment)], enriched = dplyr::first(enriched_in_models) %in% TRUE, .by = c(statistic, GWAS_ID, cluster))
  count_best <- \(statistic_name, enriched_only) {
    best <- best_tibble$best[best_tibble$statistic == statistic_name & (best_tibble$enriched | !enriched_only)]
    paste(sprintf("%s %d", names(table(best)), table(best)), collapse = ", ")
  }
  cell_type_tibble |>
    dplyr::mutate(cluster = factor(cluster, cluster_order)) |>
    ggplot2::ggplot(ggplot2::aes(median, cluster, colour = strategy)) +
    ggplot2::geom_vline(xintercept = 1, colour = "grey60", linewidth = 0.3) +
    ggplot2::geom_vline(data = reference_tibble, ggplot2::aes(xintercept = median, colour = strategy), linetype = "dashed", linewidth = 0.4) +
    ggplot2::geom_line(ggplot2::aes(group = cluster), colour = "grey70", linewidth = 0.4) +
    ggplot2::geom_point(ggplot2::aes(size = n_pairs)) +
    ggplot2::facet_grid(cols = ggplot2::vars(statistic, pair_set), scales = "free_x") +
    ggplot2::scale_x_log10() +
    SLDSC_strategy_scale() +
    ggplot2::scale_size_area(max_size = 2.5, breaks = c(1, 3, 10, 40), name = "GWAS") +
    ggplot2::labs(
      title = sprintf("Cell-type enrichment of peak and gene annotations holding %s of common SNPs", scales::percent(matched_share, accuracy = 0.1)),
      subtitle = paste(
        "Compare strategies within each cell type: the farther right, the more heritability or fine-mapped probability the same share of SNPs holds.",
        "Enriched pairs show the cell types each GWAS points to; dashed lines are randomly ranked peaks or gene windows.",
        sep = "\n"
      ),
      caption = stringr::str_wrap(paste(
        "Points: median enrichment over GWAS of each cell type's annotation; point size counts the GWAS.",
        "Enriched pairs: GWAS-cell-type pairs with q < 0.05 in the single-modality peak or gene S-LDSC model; all pairs: every GWAS with the statistic.",
        sprintf("Strategy with the higher enrichment per pair, S-LDSC: enriched pairs %s; all pairs %s.", count_best("SLDSC_enrichment", TRUE), count_best("SLDSC_enrichment", FALSE)),
        sprintf("Fine-mapped PIP: enriched pairs %s; all pairs %s.", count_best("PIP_enrichment", TRUE), count_best("PIP_enrichment", FALSE)),
        "Most enriched pairs come from the peak model on the same summary statistics, which favours peaks in S-LDSC.",
        "Cell types are ordered by their highest median S-LDSC enrichment over all pairs; medians at or below zero cannot be shown on the log scale."
      ), width = 170),
      x = "Median enrichment (log scale)",
      y = NULL
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(legend.position = "bottom", plot.caption = ggplot2::element_text(hjust = 0))
}

#' Compare S-LDSC and fine-mapped PIP enrichment at one matched size
#'
#' @inheritParams plot_SLDSC_matched_size_comparison
#' @return Patchwork of the two enrichments per annotation and of their
#'   peak-over-gene ratios per GWAS-cell-type pair.
#' @keywords internal

plot_SLDSC_matched_enrichment_agreement <- function(enrichment_tibble, share = 0.01) {
  matched_share <- share
  annotation_tibble <- enrichment_tibble |>
    dplyr::filter(share == matched_share, cluster != "random", !is.na(SLDSC_enrichment)) |>
    dplyr::mutate(pair_set = dplyr::if_else(enriched_in_models %in% TRUE, "Enriched pair", "Other pair"))
  positive_tibble <- dplyr::filter(annotation_tibble, SLDSC_enrichment > 0, PIP_enrichment > 0)
  rho_label <- positive_tibble |>
    dplyr::summarise(rho = stats::cor(SLDSC_enrichment, PIP_enrichment, method = "spearman"), .by = strategy) |>
    dplyr::reframe(label = paste(sprintf("%s %.2f", strategy, rho), collapse = ", "))
  log_axis <- \(scale) scale(labels = scales::label_number(drop0trailing = TRUE))
  annotation_plot <- positive_tibble |>
    ggplot2::ggplot(ggplot2::aes(SLDSC_enrichment, PIP_enrichment, colour = strategy)) +
    ggplot2::geom_abline(colour = "grey60", linewidth = 0.3) +
    ggplot2::geom_point(ggplot2::aes(alpha = pair_set, size = pair_set)) +
    log_axis(ggplot2::scale_x_log10) +
    log_axis(ggplot2::scale_y_log10) +
    SLDSC_strategy_scale() +
    ggplot2::scale_alpha_manual(values = c("Enriched pair" = 0.9, "Other pair" = 0.25), name = NULL) +
    ggplot2::scale_size_manual(values = c("Enriched pair" = 1.6, "Other pair" = 0.6), name = NULL) +
    ggplot2::guides(alpha = ggplot2::guide_legend(override.aes = list(colour = "grey30"))) +
    ggplot2::labs(
      subtitle = paste0("Enrichment per GWAS-cell-type annotation\nSpearman ρ: ", rho_label$label),
      x = "S-LDSC heritability enrichment", y = "Fine-mapped PIP enrichment"
    ) +
    ggplot2::theme_bw(base_size = 9)
  ratio_tibble <- annotation_tibble |>
    dplyr::select(GWAS_ID, cluster, pair_set, strategy, SLDSC_enrichment, PIP_enrichment) |>
    tidyr::pivot_wider(names_from = strategy, values_from = c(SLDSC_enrichment, PIP_enrichment)) |>
    dplyr::filter(dplyr::if_all(dplyr::matches("_enrichment_"), \(x) x > 0)) |>
    dplyr::mutate(
      SLDSC_ratio = log2(SLDSC_enrichment_peaks / SLDSC_enrichment_genes),
      PIP_ratio = log2(PIP_enrichment_peaks / PIP_enrichment_genes)
    )
  quadrant_tibble <- ratio_tibble |>
    dplyr::count(pair_set, SLDSC = SLDSC_ratio > 0, PIP = PIP_ratio > 0) |>
    dplyr::mutate(
      x = dplyr::if_else(SLDSC, Inf, -Inf), y = dplyr::if_else(PIP, Inf, -Inf),
      hjust = dplyr::if_else(SLDSC, 1.2, -0.2), vjust = dplyr::if_else(PIP, 1.5, -0.5)
    )
  ratio_plot <- ratio_tibble |>
    ggplot2::ggplot(ggplot2::aes(SLDSC_ratio, PIP_ratio)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
    ggplot2::geom_point(colour = "grey30", size = 0.8, alpha = 0.5) +
    ggplot2::geom_text(data = quadrant_tibble, ggplot2::aes(x, y, label = n, hjust = hjust, vjust = vjust), size = 3) +
    ggplot2::facet_wrap(ggplot2::vars(pair_set)) +
    ggplot2::labs(
      subtitle = "Peaks over genes per GWAS-cell-type pair",
      x = "S-LDSC: log2 enrichment ratio, peaks over genes",
      y = "Fine-mapped PIP: log2 enrichment ratio, peaks over genes"
    ) +
    ggplot2::theme_bw(base_size = 9)
  patchwork::wrap_plots(A = annotation_plot, B = ratio_plot, widths = c(1, 1.6)) +
    patchwork::plot_annotation(
      title = sprintf("Agreement of S-LDSC and fine-mapped PIP enrichment at %s of common SNPs", scales::percent(matched_share, accuracy = 0.1)),
      subtitle = paste(
        "Where both statistics favour the same strategy, it targets the signal better; pairs in the off-diagonal quadrants of B",
        "favour peaks for one statistic and genes for the other, such as strong loci captured without the polygenic background.",
        sep = "\n"
      ),
      caption = stringr::str_wrap(paste(
        "S-LDSC heritability enrichment: the annotation's share of heritability in a model with the 97 baseline-LD v2.2 annotations, over its share of common SNPs; it excludes SNPs with chi-square above 80.",
        "Fine-mapped PIP enrichment: the annotation's share of the posterior inclusion probability of autosomal credible-set variants, over its share of common SNPs; it includes the strongest loci.",
        "Wide gene windows hold whole credible sets, while small peaks hold only the variants inside them, so imprecise fine-mapping lowers PIP enrichment of peaks more.",
        sprintf("Only GWAS with summary statistics and positive enrichments are shown (%d of %d annotations in A).", nrow(positive_tibble), nrow(annotation_tibble)),
        "Enriched pairs reach q < 0.05 in the single-modality peak or gene S-LDSC model. Numbers in B count the pairs in each quadrant."
      ), width = 170),
      tag_levels = "A",
      theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0))
    ) &
    ggplot2::theme(legend.position = "bottom", legend.box = "vertical")
}
