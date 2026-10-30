#' Score each cell type's peaks for comparing peak-gene strategies
#'
#' Scores belong to a peak of one cell type within one territory. The
#' candidate territory holds the peaks of activity-by-contact (ABC) candidate
#' pairs with a positive score, and the tested territory the peaks of pairs
#' whose correlation the hierarchical model estimated and calibrated.
#' Self-promoter pairs are left out of both. A peak takes the maximum of its
#' pairs' ABC scores. Its correlation evidence is the Šidák-adjusted smallest
#' one-sided calibrated P-value of its tested genes, so that peaks near many
#' genes gain nothing from multiplicity. A pair's one-sided P-value is half its
#' background P-value when it has a link's direction, a positive slope above
#' its gene's background, and one minus that half otherwise. Correlation with
#' ABC ranks a peak's pairs with a link's direction first, by ABC.
#'
#' Covariates are the peak's counts per million in the cell type's summed
#' donor pseudobulks, ABC's activity, its accessibility share as in the S-LDSC
#' peak annotation, and its distance to the nearest TSS among its pairs. Strata
#' cross quantile bins of the three covariates within each territory and cell
#' type, with as many bins per covariate (2 to 10) as leave about 100 peaks per
#' stratum.
#'
#' @param ABC_scores_tibble Pair scores from `score_peak_gene_ABC()`.
#' @param cell_groups_tibble Cell groups from `make_peak_gene_ABC_cell_groups()`.
#' @param hierarchical_results_tibble Results of
#'   `finalize_peak_gene_hierarchical_results()`.
#' @param ATAC_pseudobulk_matrix The ATAC pseudobulks ABC was scored from.
#' @param cell_type_pseudobulk_counts_matrix Peak-by-cell-type ATAC counts of
#'   the S-LDSC peak annotation.
#' @param peak_ranges Peak ranges of the S-LDSC annotations, named by peak.
#' @return A tibble with one row per territory (`tier`), cell type (`cluster`,
#'   as in genetic enrichment) and peak (`peak_index` into `peak_ranges`), with
#'   the covariates `activity`, `share` and `distance`, the scores `ABC`,
#'   `correlation` and `correlation_ABC` (missing outside the tested
#'   territory) and the covariate `stratum`.
#' @keywords internal

get_SLDSC_peak_gene_score_tibble <- function(ABC_scores_tibble, cell_groups_tibble, hierarchical_results_tibble,
                                             ATAC_pseudobulk_matrix, cell_type_pseudobulk_counts_matrix, peak_ranges) {
  ABC_pairs <- data.table::as.data.table(ABC_scores_tibble)[link_class != "self_promoter"]
  tested_pairs <- data.table::as.data.table(hierarchical_results_tibble)[
    hierarchical_status == "estimable" & !isSelfPromoter & is.finite(hierarchical_background_pvalue) &
      cell_group %in% cell_groups_tibble$cell_group,
    .(cell_group, peak, TargetGeneID, distance = abs(distance),
      is_positive = hierarchical_coefficient > 0 & hierarchical_gene_z > 0, p = hierarchical_background_pvalue / 2)
  ]
  tested_pairs[, p := data.table::fifelse(is_positive, p, 1 - p)]
  tested_pairs[ABC_pairs, ABC := i.ABC_score, on = c("cell_group", "peak", "TargetGeneID")]
  tested_pairs[is.na(ABC), ABC := 0]
  tested_pairs[, correlation_ABC := is_positive + data.table::frank(ABC) / .N, by = "cell_group"]
  score_tibble <- data.table::rbindlist(list(
    candidate = ABC_pairs[, .(ABC = max(ABC_score), distance = min(distance)), by = c("cell_group", "peak")][ABC > 0],
    tested = tested_pairs[, .(
      ABC = max(ABC), distance = min(distance), correlation_ABC = max(correlation_ABC),
      correlation = -log10(-expm1(.N * log1p(-min(p))))
    ), by = c("cell_group", "peak")]
  ), idcol = "tier", fill = TRUE)

  activity_matrix <- sum_pseudobulks_by_cluster(ATAC_pseudobulk_matrix, cell_groups_tibble$cluster)
  activity_matrix <- sweep(activity_matrix, 2, colSums(activity_matrix), "/") * 1e6
  share_matrix <- get_SLDSC_peak_annotation_matrix(cell_type_pseudobulk_counts_matrix)
  # Calling into GenomicRanges first lets names() dispatch on ranges just read from the store.
  peak_names <- names(GenomicRanges::granges(peak_ranges))
  score_tibble[, cluster := cell_groups_tibble$cluster[match(cell_group, cell_groups_tibble$cell_group)]]
  score_tibble[, `:=`(
    peak_index = match(peak, peak_names),
    activity = activity_matrix[cbind(match(peak, rownames(activity_matrix)), match(cluster, colnames(activity_matrix)))],
    share = data.table::fcoalesce(share_matrix[cbind(match(peak, rownames(share_matrix)), match(cluster, colnames(share_matrix)))], 0)
  )]
  score_tibble[, stratum := {
    n_bins <- as.integer(min(10, max(2, floor((.N / 100)^(1 / 3)))))
    ((dplyr::ntile(activity, n_bins) - 1L) * n_bins + dplyr::ntile(share, n_bins) - 1L) * n_bins + dplyr::ntile(-distance, n_bins)
  }, by = c("tier", "cluster")]
  score_tibble[!is.na(peak_index), .(
    tier = factor(tier), cluster = factor(cluster), peak_index, activity, share, distance, ABC, correlation, correlation_ABC, stratum
  )] |>
    tibble::as_tibble()
}

#' Rank each territory's peaks by one strategy
#'
#' @param score_tibble Output of `get_SLDSC_peak_gene_score_tibble()`.
#' @param ranking A score column, `distance` (nearest first), `territory`
#'   (every peak alike), or either of the first with `_conditional` to rank
#'   within covariate strata.
#' @param tiers Territories to rank.
#' @return `tier`, `cluster`, `peak_index` and `score`, the peak's percentile
#'   in its territory, or in its stratum for conditional rankings.
#' @keywords internal

rank_SLDSC_peak_gene_scores <- function(score_tibble, ranking, tiers) {
  strategy <- sub("_conditional$", "", ranking)
  score <- switch(strategy, territory = 1, distance = -score_tibble$distance, score_tibble[[strategy]])
  score_tibble |>
    dplyr::mutate(score = score) |>
    dplyr::filter(.data$tier %in% tiers, !is.na(.data$score)) |>
    dplyr::mutate(
      score = rank(.data$score) / dplyr::n(),
      .by = dplyr::all_of(c("tier", "cluster", if (strategy != ranking) "stratum"))
    ) |>
    dplyr::select("tier", "cluster", "peak_index", "score")
}

#' Arrange one strategy's peak ranks as an S-LDSC score matrix
#'
#' @inheritParams rank_SLDSC_peak_gene_scores
#' @param peak_ranges Peak ranges the score tibble indexes.
#' @return A peak-by-`<tier>/<cell type>` matrix of ranks, zero outside each
#'   column's territory, named for `get_SLDSC_matched_annotation_record()`.
#' @keywords internal

get_SLDSC_peak_gene_score_matrix <- function(score_tibble, ranking, tiers, peak_ranges) {
  ranked <- rank_SLDSC_peak_gene_scores(score_tibble, ranking, tiers)
  columns <- paste(ranked$tier, ranked$cluster, sep = "/")
  peaks <- sort(unique(ranked$peak_index))
  peak_names <- names(GenomicRanges::granges(peak_ranges))[peaks]
  score_matrix <- matrix(0, length(peaks), dplyr::n_distinct(columns), dimnames = list(peak_names, sort(unique(columns))))
  score_matrix[cbind(match(ranked$peak_index, peaks), match(columns, colnames(score_matrix)))] <- ranked$score
  score_matrix
}

#' Contrast peak-gene strategies across GWAS
#'
#' Each contrast divides a ranking's central enrichment over GWAS-cell-type
#' pairs by its reference's over the same pairs: nothing, the whole territory,
#' a baseline ranking of the same size, or ABC for correlation with ABC. The
#' centre is the median for S-LDSC, robust to its noisy small annotations, and
#' the mean for fine-mapped PIP, whose small annotations mostly hold none. A
#' ratio is missing when its reference's centre is not positive, as S-LDSC
#' heritability shares of small annotations can be, since its sign would
#' mislead. Contrasts are computed in enriched pairs, which a single-modality
#' S-LDSC model enriches at q < 0.05, and in other pairs; their quotient
#' measures specificity. Pairs without S-LDSC results enter neither. Intervals
#' are 2.5% and 97.5% quantiles over the bootstrap samples of study groups,
#' which resample near-duplicate GWAS together, shared by all contrasts, in
#' which the ratio is defined.
#'
#' @param enrichment_tibble Combined strategy enrichments with `tier`,
#'   `ranking`, `share`, `cluster`, `enriched_in_models` and `study_group`.
#' @param n_bootstraps Bootstrap samples of study groups.
#' @return One row per territory, ranking, reference, size, statistic, pair
#'   set and cell type (`All cell types` for all together; single cell types
#'   only against nothing or the territory) with the numbers of pairs and GWAS,
#'   the estimate, its interval and the share of bootstrap samples in which
#'   it is defined.
#' @keywords internal

summarise_SLDSC_peak_gene_enrichment <- function(enrichment_tibble, n_bootstraps = 1000L) {
  pair_tibble <- enrichment_tibble |>
    dplyr::filter(!is.na(.data$enriched_in_models)) |>
    tidyr::pivot_longer(dplyr::any_of(c("SLDSC_enrichment", "PIP_enrichment")), names_to = "statistic", values_to = "enrichment") |>
    dplyr::filter(!is.na(.data$enrichment)) |>
    dplyr::select("tier", "ranking", "share", "statistic", "GWAS_ID", "study_group", "cluster", "enriched_in_models", "enrichment")
  rankings <- unique(pair_tibble$ranking)
  strategies <- intersect(c("ABC", "ABC_conditional", "correlation", "correlation_conditional", "correlation_ABC"), rankings)
  contrast_tibble <- dplyr::bind_rows(
    tibble::tibble(ranking = rankings, reference = "none"),
    tibble::tibble(ranking = setdiff(rankings, "territory"), reference = "territory"),
    tidyr::expand_grid(ranking = strategies, reference = intersect(c("activity", "share", "distance"), rankings)),
    tibble::tibble(ranking = "correlation_ABC", reference = "ABC")
  )
  territory_tibble <- pair_tibble |>
    dplyr::filter(.data$ranking == "territory") |>
    dplyr::select("tier", "statistic", "GWAS_ID", "cluster", territory_enrichment = "enrichment")
  paired_tibble <- contrast_tibble |>
    dplyr::inner_join(pair_tibble, by = "ranking", relationship = "many-to-many") |>
    dplyr::left_join(
      dplyr::select(pair_tibble, "tier", reference = "ranking", "share", "statistic", "GWAS_ID", "cluster", reference_enrichment = "enrichment"),
      by = c("reference", "tier", "share", "statistic", "GWAS_ID", "cluster"), relationship = "many-to-one"
    ) |>
    dplyr::left_join(territory_tibble, by = c("tier", "statistic", "GWAS_ID", "cluster"), relationship = "many-to-one") |>
    dplyr::mutate(reference_enrichment = dplyr::case_when(
      .data$reference == "none" ~ 1,
      .data$reference == "territory" ~ .data$territory_enrichment,
      .default = .data$reference_enrichment
    )) |>
    dplyr::filter(!is.na(.data$reference_enrichment))

  study_groups <- sort(unique(paired_tibble$study_group))
  withr::local_seed(1L)
  bootstrap_groups <- matrix(sample.int(length(study_groups), length(study_groups) * n_bootstraps, replace = TRUE), ncol = n_bootstraps)
  summarise_contrast <- \(enrichment, reference_enrichment, enriched, GWAS_ID, study_group, statistic) {
    centre <- if (statistic[[1]] == "SLDSC_enrichment") stats::median else mean
    divide <- \(numerator, denominator) if (is.finite(denominator) && denominator > 0) numerator / denominator else NA_real_
    ratio <- \(rows) divide(centre(enrichment[rows]), centre(reference_enrichment[rows]))
    ratios <- \(rows) {
      in_enriched <- enriched[rows]
      enriched_ratio <- ratio(rows[in_enriched])
      other_ratio <- ratio(rows[!in_enriched])
      c(enriched_ratio, other_ratio, divide(enriched_ratio, other_ratio))
    }
    group_rows <- split(seq_along(enrichment), factor(study_group, study_groups))
    draws <- apply(bootstrap_groups, 2, \(sample) ratios(unlist(group_rows[sample], use.names = FALSE)))
    tibble::tibble(
      pair_set = c("Enriched pairs", "Other pairs", "Specificity"),
      n_pairs = c(sum(enriched), sum(!enriched), length(enriched)),
      n_GWAS = c(dplyr::n_distinct(GWAS_ID[enriched]), dplyr::n_distinct(GWAS_ID[!enriched]), dplyr::n_distinct(GWAS_ID)),
      estimate = ratios(seq_along(enrichment)),
      lower = apply(draws, 1, stats::quantile, 0.025, na.rm = TRUE, names = FALSE),
      upper = apply(draws, 1, stats::quantile, 0.975, na.rm = TRUE, names = FALSE),
      defined_share = rowMeans(!is.na(draws))
    )
  }
  dplyr::bind_rows(
    dplyr::mutate(paired_tibble, cluster = "All cell types"),
    dplyr::filter(paired_tibble, .data$reference %in% c("none", "territory"))
  ) |>
    dplyr::reframe(
      summarise_contrast(.data$enrichment, .data$reference_enrichment, .data$enriched_in_models, .data$GWAS_ID, .data$study_group, .data$statistic),
      .by = c("tier", "ranking", "reference", "share", "statistic", "cluster")
    ) |>
    dplyr::mutate(statistic = factor(.data$statistic, c("SLDSC_enrichment", "PIP_enrichment"), c("S-LDSC heritability", "Fine-mapped PIP")))
}

# Shared ranking colours and line types of the peak-gene strategy plots.
SLDSC_peak_gene_ranking_scales <- list(
  ggplot2::scale_colour_manual(
    values = c(
      ABC = "#D95F02", ABC_conditional = "#D95F02", correlation = "#7570B3", correlation_conditional = "#7570B3",
      correlation_ABC = "#1B9E77", activity = "#666666", share = "#A6761D", distance = "#E7298A", territory = "black"
    ),
    name = "Ranking", aesthetics = c("colour", "fill")
  ),
  ggplot2::scale_linetype_manual(values = c(Ranking = "solid", `Within strata` = "dashed"), name = NULL)
)

SLDSC_peak_gene_tier_labels <- c(candidate = "ABC candidate peaks", tested = "Correlation-tested peaks")

# Shared caption sentences on what the rankings and territories are.
SLDSC_peak_gene_method_caption <- paste(
  "Territories: in each cell type, the peaks of activity-by-contact (ABC) candidate pairs with a positive score, and the peaks of pairs that the donor-level hierarchical correlation tested; self-promoter pairs are excluded.",
  "Rankings order a territory's peaks by the maximum ABC score of their pairs; by correlation, the Šidák-adjusted best one-sided background-calibrated P-value of their tested genes; by correlation with ABC, ABC among pairs with a positive calibrated correlation first;",
  "and by the baselines activity (ABC's counts per million), share (the S-LDSC accessibility share) and distance (to the nearest TSS among their pairs).",
  "Within-strata (_conditional) rankings rank peaks within quantile bins of activity, share and distance, so their top peaks keep the territory's covariate profile.",
  "Each annotation keeps a ranking's top peaks until they hold the given share of the territory's common (MAF > 5%) 1000 Genomes European SNPs."
)

#' Count the ratios a log-scale plot cannot show
#'
#' @param ratio_tibble Plotted rows of `summarise_SLDSC_peak_gene_enrichment()`.
#' @return A caption sentence.
#' @keywords internal

describe_SLDSC_peak_gene_hidden_ratios <- function(ratio_tibble) {
  sprintf(
    "Not shown: %d of %d ratios undefined, because a pair set is empty or the reference's centre enrichment is not positive, and %d not positive.",
    sum(is.na(ratio_tibble$estimate)), nrow(ratio_tibble), sum(ratio_tibble$estimate <= 0, na.rm = TRUE)
  )
}

#' Plot strategy enrichment over the territory against annotation size
#'
#' @param summary_tibble Output of `summarise_SLDSC_peak_gene_enrichment()`.
#' @return A ggplot of enrichment ratios to the territory per size.
#' @keywords internal

plot_SLDSC_peak_gene_uplift_curves <- function(summary_tibble) {
  curve_tibble <- summary_tibble |>
    dplyr::filter(.data$cluster == "All cell types", .data$reference == "territory", .data$pair_set != "Specificity") |>
    dplyr::mutate(
      tier = SLDSC_peak_gene_tier_labels[as.character(.data$tier)],
      line = dplyr::if_else(endsWith(.data$ranking, "_conditional"), "Within strata", "Ranking")
    )
  n_pairs <- curve_tibble |>
    dplyr::summarise(n = max(.data$n_pairs), GWAS = max(.data$n_GWAS), .by = c("tier", "pair_set")) |>
    dplyr::mutate(label = sprintf("%s, %s: up to %d pairs of %d GWAS", .data$tier, .data$pair_set, .data$n, .data$GWAS))
  hidden <- describe_SLDSC_peak_gene_hidden_ratios(curve_tibble)
  ggplot2::ggplot(dplyr::filter(curve_tibble, .data$estimate > 0), ggplot2::aes(.data$share, .data$estimate, colour = .data$ranking, fill = .data$ranking)) +
    ggplot2::geom_hline(yintercept = 1, colour = "grey60", linewidth = 0.3) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = pmax(.data$lower, 1e-3), ymax = .data$upper, group = .data$ranking), alpha = 0.08, colour = NA) +
    ggplot2::geom_line(ggplot2::aes(linetype = .data$line)) +
    ggplot2::geom_point(size = 1) +
    ggplot2::facet_grid(ggplot2::vars(.data$statistic), ggplot2::vars(.data$tier, .data$pair_set), scales = "free_y") +
    ggplot2::scale_x_continuous(breaks = unique(curve_tibble$share), labels = scales::label_percent()) +
    ggplot2::scale_y_log10() +
    SLDSC_peak_gene_ranking_scales +
    ggplot2::labs(
      title = "GWAS signal in each peak-gene ranking's top peaks relative to its whole territory",
      subtitle = paste(
        "Above 1, a ranking concentrates more heritability or fine-mapped probability in its top peaks than the territory holds per SNP; smaller sizes should rise further.",
        "Dashed within-strata rankings above 1 show signal beyond accessibility and distance; enriched pairs should rise and other pairs stay near 1.",
        sep = "\n"
      ),
      caption = stringr::str_wrap(paste(
        SLDSC_peak_gene_method_caption,
        "Enrichment is an annotation's share of S-LDSC heritability (fitted alone with the 97 baseline-LD v2.2 annotations) or of autosomal credible-set PIP, over its share of common SNPs.",
        "Y: the ranking's median (S-LDSC) or mean (PIP) enrichment over GWAS-cell-type pairs, divided by the territory's over the same pairs.",
        "Enriched pairs reach q < 0.05 in the single-modality peak or gene S-LDSC model; other pairs do not; GWAS without summary statistics are excluded.",
        "Bands: 95% intervals from 1,000 bootstrap samples of study groups, which keep near-duplicate GWAS together, drawn from 0.001.",
        hidden,
        paste(n_pairs$label, collapse = "; ")
      ), width = 200),
      x = "Annotation size (share of the territory's common SNPs)",
      y = "Enrichment over the territory's (log scale)"
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(legend.position = "bottom", plot.caption = ggplot2::element_text(hjust = 0))
}

#' Plot the prespecified strategy contrasts at one size
#'
#' @inheritParams plot_SLDSC_peak_gene_uplift_curves
#' @param share Annotation size as a share of the territory.
#' @return A forest plot of enrichment ratios.
#' @keywords internal

plot_SLDSC_peak_gene_contrasts <- function(summary_tibble, share = 0.25) {
  matched_share <- share
  primary <- c("candidate ABC_conditional territory", "tested correlation_conditional territory", "tested correlation_ABC ABC")
  contrast_tibble <- summary_tibble |>
    dplyr::filter(.data$cluster == "All cell types", .data$share == matched_share, .data$reference != "none") |>
    dplyr::mutate(
      is_primary = paste(.data$tier, .data$ranking, .data$reference) %in% primary,
      label = paste0(.data$ranking, " / ", .data$reference, dplyr::if_else(.data$is_primary, " (primary)", "")),
      tier = SLDSC_peak_gene_tier_labels[as.character(.data$tier)]
    ) |>
    dplyr::arrange(.data$is_primary, .data$reference != "territory", .data$label) |>
    dplyr::mutate(label = factor(.data$label, unique(.data$label)))
  hidden <- describe_SLDSC_peak_gene_hidden_ratios(contrast_tibble)
  ggplot2::ggplot(dplyr::filter(contrast_tibble, .data$estimate > 0), ggplot2::aes(.data$estimate, .data$label, colour = .data$pair_set)) +
    ggplot2::geom_vline(xintercept = 1, colour = "grey60", linewidth = 0.3) +
    ggplot2::geom_pointrange(
      ggplot2::aes(xmin = pmax(.data$lower, 1e-3), xmax = .data$upper, shape = .data$is_primary),
      position = ggplot2::position_dodge(width = 0.6), size = 0.25
    ) +
    ggplot2::facet_grid(ggplot2::vars(.data$tier), ggplot2::vars(.data$statistic), scales = "free", space = "free_y") +
    ggplot2::scale_x_log10() +
    ggplot2::scale_colour_manual(values = c(`Enriched pairs` = "#D95F02", `Other pairs` = "grey50", Specificity = "#1B9E77"), name = NULL) +
    ggplot2::scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
    ggplot2::labs(
      title = sprintf("Peak-gene strategies against their references at %s of each territory", scales::percent(matched_share)),
      subtitle = paste(
        "Ratios above 1 favour the ranking over its reference. Filled points are the prespecified primary contrasts;",
        "specificity above 1 means the gain is larger in the cell types each GWAS points to than in the others.",
        sep = "\n"
      ),
      caption = stringr::str_wrap(paste(
        SLDSC_peak_gene_method_caption,
        "Ratio: the ranking's median (S-LDSC) or mean (fine-mapped PIP) enrichment over GWAS-cell-type pairs, divided by its reference's over the same pairs; the reference is the whole territory or a ranking of the same size.",
        "Enriched pairs reach q < 0.05 in the single-modality peak or gene S-LDSC model; specificity divides the enriched-pair ratio by the other-pair ratio.",
        "Lines: 95% intervals from 1,000 bootstrap samples of study groups, which keep near-duplicate GWAS together, the same samples for every contrast; lower ends are drawn from 0.001.",
        hidden
      ), width = 190),
      x = "Ratio of enrichments (log scale)",
      y = NULL
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(legend.position = "bottom", plot.caption = ggplot2::element_text(hjust = 0))
}

#' Plot each cell type's enrichment over its territory at one size
#'
#' @inheritParams plot_SLDSC_peak_gene_contrasts
#' @return A dot plot of each cell type's ratio per ranking.
#' @keywords internal

plot_SLDSC_peak_gene_cell_types <- function(summary_tibble, share = 0.25) {
  matched_share <- share
  cell_type_tibble <- summary_tibble |>
    dplyr::filter(
      .data$cluster != "All cell types", .data$share == matched_share, .data$reference == "territory",
      .data$pair_set != "Specificity"
    ) |>
    dplyr::mutate(tier = SLDSC_peak_gene_tier_labels[as.character(.data$tier)])
  hidden <- describe_SLDSC_peak_gene_hidden_ratios(cell_type_tibble)
  ggplot2::ggplot(dplyr::filter(cell_type_tibble, .data$estimate > 0), ggplot2::aes(.data$estimate, .data$cluster, colour = .data$ranking)) +
    ggplot2::geom_vline(xintercept = 1, colour = "grey60", linewidth = 0.3) +
    ggplot2::geom_point(ggplot2::aes(size = .data$n_GWAS), position = ggplot2::position_dodge(width = 0.7), alpha = 0.8) +
    ggplot2::facet_grid(ggplot2::vars(.data$tier), ggplot2::vars(.data$statistic, .data$pair_set), scales = "free", space = "free_y") +
    ggplot2::scale_x_log10() +
    ggplot2::scale_size_area(max_size = 2.5, breaks = c(1, 3, 10, 30), name = "GWAS") +
    SLDSC_peak_gene_ranking_scales[[1]] +
    ggplot2::labs(
      title = sprintf("Each cell type's GWAS signal in the top %s of its peak-gene rankings, relative to its territory", scales::percent(matched_share)),
      subtitle = paste(
        "Compare rankings within a cell type: the farther right, the more signal the top peaks hold per SNP than the territory.",
        "Cell types with few GWAS give noisy ratios; enriched pairs show the cell types each GWAS points to.",
        sep = "\n"
      ),
      caption = stringr::str_wrap(paste(
        SLDSC_peak_gene_method_caption,
        "Points: the ranking's median (S-LDSC) or mean (fine-mapped PIP) enrichment over the cell type's GWAS, divided by the territory's; point size counts the GWAS.",
        "Enriched pairs reach q < 0.05 in the single-modality peak or gene S-LDSC model.",
        hidden
      ), width = 190),
      x = "Enrichment over the territory's (log scale)",
      y = NULL
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(legend.position = "bottom", plot.caption = ggplot2::element_text(hjust = 0))
}

#' Plot the territories and the covariate balance of each ranking
#'
#' Balance is the standardised mean difference of a covariate between a
#' ranking's top quarter of a territory's peaks and the whole territory.
#'
#' @param score_tibble Output of `get_SLDSC_peak_gene_score_tibble()`.
#' @param enrichment_tibble Combined strategy enrichments, for territory sizes.
#' @param rankings_tibble Rankings with their `tiers`.
#' @return A patchwork of territory sizes and covariate balance.
#' @keywords internal

plot_SLDSC_peak_gene_design <- function(score_tibble, enrichment_tibble, rankings_tibble) {
  territory_tibble <- enrichment_tibble |>
    dplyr::filter(.data$ranking == "territory") |>
    dplyr::distinct(.data$tier, .data$cluster, .data$SNP_share) |>
    dplyr::left_join(dplyr::count(score_tibble, .data$tier, .data$cluster, name = "n_peaks") |>
      dplyr::mutate(dplyr::across(c("tier", "cluster"), as.character)), by = c("tier", "cluster")) |>
    dplyr::mutate(tier = SLDSC_peak_gene_tier_labels[.data$tier])
  covariate_tibble <- dplyr::transmute(score_tibble, .data$tier, .data$cluster, .data$peak_index,
    `log10 activity` = log10(.data$activity + 1), share = .data$share, `log10 distance` = log10(.data$distance + 1))
  balance_tibble <- purrr::pmap_dfr(dplyr::filter(rankings_tibble, .data$ranking != "territory"), \(ranking, tiers, ...) {
    rank_SLDSC_peak_gene_scores(score_tibble, ranking, tiers) |>
      dplyr::inner_join(covariate_tibble, by = c("tier", "cluster", "peak_index")) |>
      tidyr::pivot_longer(c("log10 activity", "share", "log10 distance"), names_to = "covariate") |>
      dplyr::summarise(
        difference = (mean(.data$value[.data$score > 0.75]) - mean(.data$value)) / stats::sd(.data$value),
        .by = c("tier", "cluster", "covariate")
      ) |>
      dplyr::mutate(ranking = ranking)
  }) |>
    dplyr::mutate(
      tier = SLDSC_peak_gene_tier_labels[as.character(.data$tier)],
      ranking = factor(.data$ranking, rev(rankings_tibble$ranking))
    )
  size_plot <- ggplot2::ggplot(territory_tibble, ggplot2::aes(.data$SNP_share, stats::reorder(.data$cluster, .data$SNP_share))) +
    ggplot2::geom_col(fill = "grey60") +
    ggplot2::geom_text(ggplot2::aes(label = scales::comma(.data$n_peaks)), hjust = -0.1, size = 2.5) +
    ggplot2::facet_grid(ggplot2::vars(.data$tier), scales = "free_y", space = "free_y") +
    ggplot2::scale_x_continuous(labels = scales::label_percent(accuracy = 0.01), expand = ggplot2::expansion(mult = c(0, 0.3))) +
    ggplot2::labs(subtitle = "Territory size (labels: peaks)", x = "Share of common SNPs", y = NULL) +
    ggplot2::theme_bw(base_size = 9)
  balance_plot <- ggplot2::ggplot(balance_tibble, ggplot2::aes(.data$difference, .data$ranking, colour = .data$covariate)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.6), size = 0.8, alpha = 0.6) +
    ggplot2::facet_grid(ggplot2::vars(.data$tier), scales = "free_y", space = "free_y") +
    ggplot2::scale_colour_brewer(palette = "Dark2", name = NULL) +
    ggplot2::labs(subtitle = "Covariates of each ranking's top quarter (points: cell types)", x = "Standardised mean difference from the territory", y = NULL) +
    ggplot2::theme_bw(base_size = 9)
  patchwork::wrap_plots(A = size_plot, B = balance_plot, widths = c(1, 1.4)) +
    patchwork::plot_annotation(
      title = "Territories and covariate balance of the peak-gene strategy comparison",
      subtitle = paste(
        "Small territories give noisy enrichments, so read the correlation-tested territory with A in mind.",
        "In B, within-strata rankings should sit near zero: their top peaks keep the territory's accessibility and distance.",
        sep = "\n"
      ),
      caption = stringr::str_wrap(paste(
        SLDSC_peak_gene_method_caption,
        "B compares the top quarter of each territory's peaks, counted by peaks rather than SNPs, with the whole territory, in units of the territory's standard deviation.",
        "Activity is in counts per million; distance is in base pairs."
      ), width = 190),
      tag_levels = "A",
      theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0))
    ) &
    ggplot2::theme(legend.position = "bottom")
}
