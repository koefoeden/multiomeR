#' Binarise cell-type scores at matched shares of common SNPs
#'
#' Rank each cell type's autosomal territory by its score and keep the
#' highest-scoring ranges until they hold each share of the common (MAF > 5%)
#' reference SNPs. A SNP in overlapping ranges takes their largest score, as in
#' `get_SLDSC_LD_score_record()`, so the ranges are first split into disjoint
#' pieces. Pieces are kept whole, in genomic order among ties, so an annotation
#' can exceed its share by one piece; pieces without score are never kept, and
#' a cell type whose scored pieces hold fewer common SNPs than a share gets no
#' annotation of that share. Ten `random_<seed>` annotations rank the same
#' ranges at random with seeds 1 to 10, as a reference for what the territory
#' alone concentrates. Within territories, each column's shares are of the
#' common SNPs in its own scored pieces instead, and a share of 1 marks the
#' whole territory, the reference for its rankings.
#'
#' @param score_matrix Range-by-cell-type scores named by range.
#' @param ranges Ranges named like the score rows.
#' @param reference_dir Directory returned by `download_SLDSC_reference()`.
#' @param shares Shares of the common reference SNPs.
#' @param within_territory Whether shares are of each column's scored common
#'   SNPs, without random annotations.
#' @return A list with the kept pieces, their binary piece-by-annotation matrix
#'   with annotations named `<cell type>:<share>`, each annotation's sums of
#'   the baseline-LD annotations over its common SNPs, named like the baseline
#'   LD scores (its `baseL2` sum counts these SNPs), and the `cluster` and
#'   `share` of the dropped annotations.
#' @keywords internal

get_SLDSC_matched_annotation_record <- function(score_matrix, ranges, reference_dir, shares = c(0.005, 0.01, 0.02, 0.05),
                                                within_territory = FALSE) {
  if (!within_territory) {
    random_matrix <- vapply(1:10, \(seed) withr::with_seed(seed, stats::runif(nrow(score_matrix))), numeric(nrow(score_matrix)))
    colnames(random_matrix) <- paste0("random_", 1:10)
    score_matrix <- cbind(score_matrix, random_matrix)
  }
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
  territory_SNPs <- colSums(piece_n_common * (score_matrix > 0))
  share_SNPs <- if (within_territory) territory_SNPs else rep(n_common, ncol(score_matrix))
  annotation_matrix <- do.call(cbind, purrr::map(seq_len(ncol(score_matrix)), \(column) {
    ranked <- order(-score_matrix[, column])
    ranked <- ranked[score_matrix[ranked, column] > 0]
    SNPs_before <- cumsum(piece_n_common[ranked]) - piece_n_common[ranked]
    vapply(shares, \(share) seq_along(pieces) %in% ranked[SNPs_before < share * share_SNPs[[column]]], logical(length(pieces)))
  })) * 1
  dimnames(annotation_matrix) <- list(as.character(pieces), paste0(rep(colnames(score_matrix), each = length(shares)), ":", shares))
  # A cell type's whole territory would otherwise stand in for a share it cannot fill.
  is_filled <- rep(territory_SNPs, each = length(shares)) >= as.vector(outer(shares, share_SNPs))
  annotation_matrix <- annotation_matrix[, is_filled, drop = FALSE]
  kept <- rowSums(annotation_matrix) > 0
  list(
    ranges = stats::setNames(GenomicRanges::granges(pieces), as.character(pieces))[kept],
    annotation_matrix = annotation_matrix[kept, , drop = FALSE],
    baseline_overlap_matrix = crossprod(annotation_matrix, piece_baseline_matrix),
    dropped_tibble = tibble::tibble(
      cluster = rep(colnames(score_matrix), each = length(shares)),
      share = rep(shares, times = ncol(score_matrix))
    )[!is_filled, ]
  )
}

#' Describe the matched-size annotations that cell types cannot fill
#'
#' @param dropped_tibble The `strategy`, `cluster` and `share` of the
#'   annotations `get_SLDSC_matched_annotation_record()` dropped.
#' @return A caption sentence listing them by strategy and share.
#' @keywords internal

describe_SLDSC_dropped_annotations <- function(dropped_tibble) {
  if (nrow(dropped_tibble) == 0L) {
    return("Every cell type's peaks and gene windows fill every size.")
  }
  groups <- dropped_tibble |>
    dplyr::summarise(clusters = toString(sort(.data$cluster)), n = dplyr::n(), .by = c("strategy", "share")) |>
    dplyr::arrange(.data$strategy, .data$share)
  paste0(
    nrow(dropped_tibble), " cell-type annotations are dropped because the cell type's scored peaks or gene windows hold fewer common SNPs than their size: ",
    paste(sprintf("%s at %s (%d: %s)", groups$strategy, scales::percent(groups$share, accuracy = 0.1), groups$n, groups$clusters), collapse = "; "), "."
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
#' @param candidates Binary candidate annotations, one per model. A candidate
#'   with LD scores in fewer than two blocks has no jackknife and no row.
#' @return One row per candidate with its coefficient, its share of the
#'   model's heritability with its delete-one-block values, and the model's
#'   total heritability, each with a jackknife standard error.
#' @keywords internal

fit_SLDSC_enrichment_regressions <- function(chi2, N, weight_L2, L2_matrix, M, overlap_matrix, shared_annotations, candidates, n_blocks = 200L) {
  annotations <- c(shared_annotations, candidates)
  blocks <- get_SLDSC_block_crossproducts(chi2, N, weight_L2, L2_matrix, M, shared_annotations, annotations, n_blocks)
  n_shared <- length(shared_annotations)
  blocks_with_scores <- rowSums(apply(blocks$block_crossproducts, 3L, diag) > 0)
  purrr::map_dfr(seq_along(candidates), \(i) {
    if (blocks_with_scores[[n_shared + i]] < 2L) return(NULL)
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
      h2_share_delete_values = list(h2_share[-1]),
      h2_share = h2_share[[1]],
      h2_share_se = stats::sd(pseudovalues) / sqrt(n_blocks),
      h2 = h2[[1]] / blocks$N_bar,
      h2_se = sqrt(drop(model_M %*% fit$cov[annotation_idx, annotation_idx] %*% model_M)) / blocks$N_bar
    )
  })
}

#' Place autosomal credible-set variants in range annotations
#'
#' Variants are assigned by position, so indels and variants missing from the
#' reference count.
#'
#' @param credible_set_GRanges Credible-set variants.
#' @param ranges Disjoint annotation ranges.
#' @param annotation_matrix Range-by-annotation matrix.
#' @return A list with the autosomal `variants` and their variant-by-annotation
#'   `inside_matrix`, zero outside every range.
#' @keywords internal

get_SLDSC_credible_set_annotation_matrix <- function(credible_set_GRanges, ranges, annotation_matrix) {
  variants <- credible_set_GRanges[as.character(GenomicRanges::seqnames(credible_set_GRanges)) %in% paste0("chr", 1:22)]
  piece_idx <- GenomicRanges::findOverlaps(variants, ranges, select = "first")
  inside_matrix <- annotation_matrix[piece_idx, , drop = FALSE]
  inside_matrix[is.na(piece_idx), ] <- 0
  list(variants = variants, inside_matrix = inside_matrix)
}

#' Measure the fine-mapped posterior probability inside range annotations
#'
#' Sum the posterior inclusion probabilities (PIPs) of autosomal credible-set
#' variants inside each annotation, placed by
#' `get_SLDSC_credible_set_annotation_matrix()`, and divide by their total.
#' Loci, the credible sets whose centres lie within 500 kb of each other, are
#' bootstrapped for uncertainty.
#'
#' @param credible_set_GRanges Credible-set variants with `studyLocusId` and
#'   `posteriorProbability`.
#' @param ranges Disjoint annotation ranges.
#' @param annotation_matrix Binary range-by-annotation matrix.
#' @param n_bootstraps Bootstrap samples of loci.
#' @return One row per annotation with its share of the PIP, the bootstrap
#'   standard error and 95% interval of that share, its bootstrap shares
#'   `PIP_bootstrap_shares`, paired across annotations, and the numbers of
#'   credible sets and loci.
#' @keywords internal

get_SLDSC_PIP_share_tibble <- function(credible_set_GRanges, ranges, annotation_matrix, n_bootstraps = 1000L) {
  placed <- get_SLDSC_credible_set_annotation_matrix(credible_set_GRanges, ranges, annotation_matrix)
  variants <- placed$variants
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
  locus_PIP_matrix <- rowsum(variants$posteriorProbability * placed$inside_matrix, variant_locus)
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
    PIP_bootstrap_shares = lapply(seq_len(ncol(bootstrap_shares)), \(i) bootstrap_shares[, i]),
    n_credible_sets = nrow(set_tibble),
    n_loci = length(loci)
  )
}

#' Evaluate one GWAS's enrichment in matched-size cell-type annotations
#'
#' Divide each annotation's share of S-LDSC heritability, fitted with the
#' baseline-LD model, and its share of fine-mapped PIP by its share of common
#' reference SNPs. The random draws of each share are averaged into one
#' `random` annotation: their PIP shares and common SNPs directly, and their
#' heritability shares within each jackknife block, which gives the mean a
#' jackknife standard error. Its coefficient and total heritability are means
#' over the draws, without standard errors.
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
#'   enrichment's bootstrap interval, the S-LDSC coefficient and total
#'   heritability, and the heritability share's delete-one-block values;
#'   S-LDSC estimates are missing without summary statistics.
#' @keywords internal

get_SLDSC_matched_enrichment_tibble <- function(GWAS_input_record, sumstats_record, baseline_record, annotation_record, LD_score_record) {
  if (!identical(sumstats_record$sourceId, GWAS_input_record$sourceId)) {
    stop("S-LDSC summary statistics of ", sumstats_record$sourceId, " were paired with ", GWAS_input_record$GWAS_ID, ".")
  }
  annotation_matrix <- annotation_record$annotation_matrix
  group_names <- sub("^random_[0-9]+:", "random:", colnames(annotation_matrix))
  annotation_groups <- split(colnames(annotation_matrix), factor(group_names, unique(group_names)))
  M <- c(baseline_record$M_5_50, annotation_record$baseline_overlap_matrix[, "baseL2"])
  enrichment_tibble <- get_SLDSC_PIP_share_tibble(
    GWAS_input_record$credible_set_GRanges,
    annotation_record$ranges,
    vapply(annotation_groups, \(group) rowMeans(annotation_matrix[, group, drop = FALSE]), numeric(nrow(annotation_matrix)))
  )
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
    reference_tibble <- SLDSC_tibble |>
      dplyr::filter(startsWith(annotation, "random_")) |>
      dplyr::mutate(annotation = sub("^random_[0-9]+:", "random:", annotation)) |>
      dplyr::summarise(
        dplyr::across(c(tau, h2_share, h2), mean),
        h2_share_delete_values = list(Reduce(`+`, h2_share_delete_values) / dplyr::n()),
        .by = annotation
      ) |>
      dplyr::mutate(h2_share_se = purrr::map2_dbl(h2_share, h2_share_delete_values, \(estimate, delete_values) {
        n_blocks <- length(delete_values)
        stats::sd(n_blocks * estimate - (n_blocks - 1) * delete_values) / sqrt(n_blocks)
      }))
    SLDSC_tibble <- dplyr::bind_rows(dplyr::filter(SLDSC_tibble, !startsWith(annotation, "random_")), reference_tibble)
    enrichment_tibble <- dplyr::left_join(enrichment_tibble, SLDSC_tibble, by = "annotation", relationship = "one-to-one")
  }
  SNP_share <- vapply(annotation_groups, \(group) mean(M[group]), numeric(1))[enrichment_tibble$annotation] / baseline_record$M_5_50[["baseL2"]]
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

#' Decompose one GWAS's fine-mapped PIP enrichment into credible-set contributions
#'
#' A credible set's contribution to an annotation's PIP enrichment is its
#' posterior inclusion probability inside the annotation, over the PIP of all
#' autosomal credible sets and the annotation's share of common SNPs, so the
#' contributions of all credible sets sum to the enrichment of
#' `get_SLDSC_matched_enrichment_tibble()`. The random annotations are left
#' out.
#'
#' @param GWAS_input_record One normalized GWAS input record.
#' @param annotation_record Record from `get_SLDSC_matched_annotation_record()`.
#' @param baseline_record Record from `read_SLDSC_baseline_record()`, which
#'   gives the number of common SNPs.
#' @return One row per cell type, size and credible set with a nonzero
#'   contribution, with the credible set's lead variant.
#' @keywords internal

get_SLDSC_PIP_locus_contribution_tibble <- function(GWAS_input_record, annotation_record, baseline_record) {
  cell_type_annotations <- grep("^random_", colnames(annotation_record$annotation_matrix), value = TRUE, invert = TRUE)
  placed <- get_SLDSC_credible_set_annotation_matrix(
    GWAS_input_record$credible_set_GRanges,
    annotation_record$ranges,
    annotation_record$annotation_matrix[, cell_type_annotations, drop = FALSE]
  )
  variants <- placed$variants
  SNP_share <- annotation_record$baseline_overlap_matrix[cell_type_annotations, "baseL2"] / baseline_record$M_5_50[["baseL2"]]
  contribution_matrix <- sweep(
    rowsum(variants$posteriorProbability * placed$inside_matrix, variants$studyLocusId) / sum(variants$posteriorProbability),
    2, SNP_share, "/"
  )
  nonzero <- which(contribution_matrix > 0, arr.ind = TRUE)
  lead_tibble <- tibble::tibble(studyLocusId = variants$studyLocusId, lead_variantId = variants$variantId, PIP = variants$posteriorProbability) |>
    dplyr::slice_max(PIP, n = 1, with_ties = FALSE, by = studyLocusId) |>
    dplyr::select(-PIP)
  tibble::tibble(
    GWAS_ID = GWAS_input_record$GWAS_ID,
    annotation = cell_type_annotations[nonzero[, "col"]],
    studyLocusId = rownames(contribution_matrix)[nonzero[, "row"]],
    PIP_enrichment_contribution = contribution_matrix[nonzero]
  ) |>
    tidyr::separate_wider_regex(annotation, c(cluster = ".*", ":", share = "[^:]+$")) |>
    dplyr::mutate(share = as.numeric(share)) |>
    dplyr::left_join(lead_tibble, by = "studyLocusId", relationship = "many-to-one")
}

#' Mark the GWAS-cell-type pairs that single-modality S-LDSC enriches
#'
#' @param SLDSC_cell_type_tibble Combined cell-type S-LDSC results.
#' @return `GWAS_ID`, `cluster` and `enriched_in_models`, whether the peak or
#'   gene model gives q < 0.05.
#' @keywords internal

get_SLDSC_enriched_pair_tibble <- function(SLDSC_cell_type_tibble) {
  SLDSC_cell_type_tibble |>
    dplyr::filter(.data$model != "joint") |>
    dplyr::summarise(enriched_in_models = any(.data$q < 0.05), .by = c("GWAS_ID", "cluster"))
}

#' Test matched-size enrichment and find each GWAS's top sets
#'
#' Test each annotation's S-LDSC heritability enrichment above 1, the
#' enrichment of a random variant, in a one-sided normal test with its
#' jackknife standard error. Benjamini-Hochberg q-values adjust the cell types'
#' P-values over GWAS and cell types within each strategy and size. Each GWAS's
#' top set of cell types at each strategy and size comes from
#' `get_top_set_tibble()` on enrichment, with paired delete-one-block values.
#' Fine-mapped PIP enrichment is tested in the same way, with its locus-bootstrap
#' standard error and paired bootstrap values, in GWAS with at least
#' `min_loci` loci, since fewer give no usable bootstrap. The random reference
#' is tested, but takes no part in q-values or top sets.
#'
#' @param enrichment_tibble Combined matched-size enrichment of all GWAS with
#'   the heritability shares' delete-one-block values and the PIP shares'
#'   bootstrap values.
#' @param min_loci Fewest loci whose bootstrap tests fine-mapped enrichment.
#' @return The tibble with `z`, `p`, `q`, `in_top_set` and `z_below_anchor`,
#'   and their fine-mapped counterparts prefixed `PIP_`, without the
#'   delete-one-block and bootstrap values.
#' @keywords internal

test_SLDSC_matched_enrichment <- function(enrichment_tibble, min_loci = 10L) {
  enrichment_tibble |>
    dplyr::mutate(z = (SLDSC_enrichment - 1) / SLDSC_enrichment_se, p = stats::pnorm(z, lower.tail = FALSE)) |>
    dplyr::mutate(q = stats::p.adjust(dplyr::if_else(cluster == "random", NA, p), method = "BH"), .by = c(strategy, share)) |>
    dplyr::mutate(
      get_top_set_tibble(
        dplyr::if_else(cluster == "random", NA, SLDSC_enrichment),
        do.call(rbind, h2_share_delete_values) / SNP_share,
        q < 0.05,
        jackknife = TRUE
      ),
      .by = c(strategy, GWAS_ID, share)
    ) |>
    dplyr::mutate(
      PIP_z = dplyr::if_else(n_loci >= min_loci, (PIP_enrichment - 1) / PIP_enrichment_se, NA),
      PIP_p = stats::pnorm(PIP_z, lower.tail = FALSE)
    ) |>
    dplyr::mutate(PIP_q = stats::p.adjust(dplyr::if_else(cluster == "random", NA, PIP_p), method = "BH"), .by = c(strategy, share)) |>
    dplyr::mutate(
      get_top_set_tibble(
        dplyr::if_else(cluster == "random", NA, PIP_enrichment),
        do.call(rbind, PIP_bootstrap_shares) / SNP_share,
        PIP_q < 0.05
      ) |>
        dplyr::rename(PIP_in_top_set = in_top_set, PIP_z_below_anchor = z_below_anchor),
      .by = c(strategy, GWAS_ID, share)
    ) |>
    dplyr::select(-h2_share_delete_values, -PIP_bootstrap_shares) |>
    dplyr::relocate(z, p, q, in_top_set, z_below_anchor, .after = SLDSC_enrichment_se) |>
    dplyr::relocate(PIP_z, PIP_p, PIP_q, PIP_in_top_set, PIP_z_below_anchor, .after = PIP_enrichment_upper)
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
    values = c(accessible_peaks = "#7570B3", specific_peaks = "#D95F02", specific_genes = "#1B9E77"),
    labels = c(
      accessible_peaks = "Most accessible peaks",
      specific_peaks = "Most specific peaks (share of accessibility)",
      specific_genes = "Most specific genes (share of expression, \u00b1100 kb)"
    ),
    name = "Strategy",
    aesthetics = aesthetics
  )
}

#' Get the shared limits of the enrichment curves and heatmaps
#'
#' The limits cap the outer `cap_quantile` of all cell-type estimates of the
#' `accessible_peaks` strategy at each end, with the lower end at most 0,
#' rounded outward to multiples of 5.
#'
#' @param enrichment_tibble Matched-size enrichment from
#'   `test_SLDSC_matched_enrichment()`.
#' @param cap_quantile Share of cell-type estimates beyond each limit.
#' @param estimate_col Enrichment column, S-LDSC or fine-mapped PIP.
#' @return The lower and upper limit.
#' @keywords internal

get_SLDSC_enrichment_limits <- function(enrichment_tibble, cap_quantile = 0.01, estimate_col = "SLDSC_enrichment") {
  estimates <- dplyr::filter(enrichment_tibble, strategy == "accessible_peaks", cluster != "random")[[estimate_col]]
  c(
    floor(min(0, stats::quantile(estimates, cap_quantile, na.rm = TRUE, names = FALSE)) / 5) * 5,
    ceiling(stats::quantile(estimates, 1 - cap_quantile, na.rm = TRUE, names = FALSE) / 5) * 5
  )
}

#' Plot each cell type's heritability enrichment in its most accessible peaks
#'
#' One plot per GWAS category, with one panel per GWAS with S-LDSC estimates,
#' shows each cell type's enrichment in the `accessible_peaks` strategy's
#' annotations against their size, coloured by compartment, on one y scale
#' shared by all plots, by default with the limits of
#' `get_SLDSC_enrichment_limits()`. Points beyond the limits are drawn at them as
#' triangles labelled with their value. Up to `n_labels` cell types per GWAS
#' are highlighted and labelled: those that anchor its top set at the most
#' sizes, then those in its top set at the most sizes, then the most enriched.
#'
#' @param enrichment_tibble Matched-size enrichment from
#'   `test_SLDSC_matched_enrichment()` with GWAS categories.
#' @param compartments_patterns Optional cell-type compartment patterns.
#' @inheritParams plot_SLDSC_matched_enrichment_curves
#' @param n_labels Cell types to highlight and label per GWAS.
#' @param n_columns Panels per row.
#' @param cap_quantile Share of cell-type estimates beyond each end of the
#'   shared y scale.
#' @param y_limits Optional fixed limits of the shared y scale instead.
#' @return A list of faceted line plots named by GWAS category.
#' @keywords internal

plot_SLDSC_enrichment_curves <- function(enrichment_tibble, compartments_patterns, dropped_tibble, n_labels = 4L, n_columns = 4L,
                                         cap_quantile = 0.01, y_limits = NULL) {
  without_sumstats <- enrichment_tibble |>
    dplyr::filter(all(is.na(SLDSC_enrichment)), .by = GWAS_ID) |>
    dplyr::distinct(Category, GWAS_ID)
  curve_tibble <- enrichment_tibble |>
    dplyr::filter(strategy == "accessible_peaks", !is.na(SLDSC_enrichment)) |>
    dplyr::arrange(Category, GWAS_ID)
  limits <- y_limits %||% get_SLDSC_enrichment_limits(enrichment_tibble, cap_quantile)
  capped <- \(values) pmin(pmax(values, limits[[1]]), limits[[2]])
  curve_tibble <- curve_tibble |>
    dplyr::mutate(
      drawn_enrichment = capped(SLDSC_enrichment),
      lower = capped(SLDSC_enrichment - stats::qnorm(0.975) * SLDSC_enrichment_se),
      upper = capped(SLDSC_enrichment + stats::qnorm(0.975) * SLDSC_enrichment_se),
      in_top_set = in_top_set %in% TRUE
    )
  reference_tibble <- dplyr::filter(curve_tibble, cluster == "random")
  cell_type_tibble <- curve_tibble |>
    dplyr::filter(cluster != "random") |>
    dplyr::mutate(compartment = get_cluster_compartments(cluster, compartments_patterns))
  label_tibble <- cell_type_tibble |>
    dplyr::summarise(
      anchor_sizes = sum(z_below_anchor %in% 0),
      top_set_sizes = sum(in_top_set),
      highest_enrichment = max(SLDSC_enrichment),
      share = min(share),
      .by = c(Category, GWAS_ID, cluster)
    ) |>
    dplyr::filter(top_set_sizes > 0) |>
    dplyr::arrange(dplyr::desc(anchor_sizes), dplyr::desc(top_set_sizes), dplyr::desc(highest_enrichment)) |>
    dplyr::slice_head(n = n_labels, by = GWAS_ID) |>
    dplyr::left_join(dplyr::select(cell_type_tibble, GWAS_ID, cluster, share, drawn_enrichment, compartment), by = c("GWAS_ID", "cluster", "share"))
  highlighted_tibble <- dplyr::semi_join(cell_type_tibble, label_tibble, by = c("GWAS_ID", "cluster"))
  beyond_cap_tibble <- cell_type_tibble |>
    dplyr::filter(SLDSC_enrichment != drawn_enrichment) |>
    dplyr::mutate(
      direction = dplyr::if_else(SLDSC_enrichment > limits[[2]], "above", "below"),
      value_label = dplyr::if_else(abs(SLDSC_enrichment) < 10, sprintf("%.1f", SLDSC_enrichment), sprintf("%.0f", SLDSC_enrichment))
    )
  median_top_set_size <- stats::median(dplyr::count(dplyr::filter(cell_type_tibble, in_top_set), GWAS_ID, share)$n)
  compartment_colours <- get_compartment_colours(setdiff(levels(cell_type_tibble$compartment), "Other"))
  categories <- unique(curve_tibble$Category)
  purrr::map(rlang::set_names(categories), \(category) {
    of_category <- \(tibble) dplyr::filter(tibble, Category %in% category)
    category_without_sumstats <- of_category(without_sumstats)$GWAS_ID
    ggplot2::ggplot(of_category(cell_type_tibble), ggplot2::aes(share, drawn_enrichment, colour = compartment, group = cluster)) +
      ggplot2::geom_hline(yintercept = 1, colour = "grey60", linewidth = 0.3) +
      ggplot2::geom_line(data = dplyr::anti_join(of_category(cell_type_tibble), label_tibble, by = c("GWAS_ID", "cluster")), linewidth = 0.3, alpha = 0.35) +
      ggplot2::geom_ribbon(data = of_category(reference_tibble), ggplot2::aes(share, ymin = lower, ymax = upper), inherit.aes = FALSE, fill = "grey20", alpha = 0.08) +
      ggplot2::geom_line(data = of_category(reference_tibble), ggplot2::aes(share, drawn_enrichment), inherit.aes = FALSE, colour = "grey20", linewidth = 0.5, linetype = "22") +
      ggplot2::geom_errorbar(data = of_category(highlighted_tibble), ggplot2::aes(ymin = lower, ymax = upper), width = 0, linewidth = 0.3, alpha = 0.6, show.legend = FALSE) +
      ggplot2::geom_line(data = of_category(highlighted_tibble), linewidth = 0.8) +
      ggplot2::geom_point(data = of_category(highlighted_tibble), shape = 21, fill = "white", size = 1.6, stroke = 0.7, show.legend = FALSE) +
      ggplot2::geom_point(data = dplyr::filter(of_category(highlighted_tibble), in_top_set), size = 1.9, show.legend = FALSE) +
      ggplot2::geom_point(data = of_category(beyond_cap_tibble), ggplot2::aes(shape = direction, fill = compartment), size = 1.8, stroke = 0.3, colour = "grey20", show.legend = FALSE) +
      ggrepel::geom_text_repel(
        data = of_category(beyond_cap_tibble), ggplot2::aes(label = value_label), colour = "grey30", size = 2, max.overlaps = Inf,
        segment.colour = "grey70", segment.size = 0.2, min.segment.length = 0.1, box.padding = 0.15, seed = 1
      ) +
      ggrepel::geom_text_repel(
        data = of_category(label_tibble), ggplot2::aes(label = cluster), colour = "grey15", size = 2.4, hjust = 1, direction = "y",
        nudge_x = -0.003, segment.colour = "grey60", segment.size = 0.25, min.segment.length = 0, seed = 1
      ) +
      ggplot2::scale_x_continuous(
        breaks = sort(unique(curve_tibble$share)), labels = scales::label_number(scale = 100, drop0trailing = TRUE),
        expand = ggplot2::expansion(mult = c(0.75, 0.04))
      ) +
      ggplot2::scale_colour_manual(values = compartment_colours, limits = names(compartment_colours), name = "Compartment") +
      ggplot2::scale_fill_manual(values = compartment_colours, limits = names(compartment_colours), guide = "none") +
      ggplot2::scale_shape_manual(values = c(above = 24, below = 25), guide = "none") +
      ggplot2::coord_cartesian(ylim = limits) +
      ggplot2::facet_wrap(ggplot2::vars(GWAS_ID), ncol = n_columns) +
      ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2, override.aes = list(linewidth = 0.8, alpha = 1))) +
      ggplot2::labs(
        title = paste("Heritability enrichment in each cell type's most accessible peaks:", category, "GWAS"),
        subtitle = paste(
          "Compare cell types at each size: a higher point holds more heritability per common variant than a random one (1, grey line).",
          "The dashed line and band are randomly chosen peaks. Bold, labelled lines lead each GWAS's top sets,",
          "with filled points at the sizes where they are in them.",
          sep = "\n"
        ),
        caption = stringr::str_wrap(paste(
          "Each cell type's annotation keeps its peaks with the highest depth-normalized pseudobulk accessibility, among peaks with at least 10 reads in it, until they hold the given share of common (MAF > 5%) 1000 Genomes European SNPs.",
          "Enrichment: the annotation's share of S-LDSC heritability, in a model with the 97 baseline-LD v2.2 annotations, over its share of common SNPs; SNPs with chi-square above 80 are excluded. Bars and band: 95% intervals from the 200-block jackknife.",
          "Random reference: the mean of 10 random rankings of the same peaks, averaged within jackknife blocks.",
          "Top set at each size: of the cell types with enrichment above 1 at Benjamini-Hochberg q < 0.05 (one-sided jackknife z-test, over all GWAS and cell types at that size), the most enriched one and every cell type whose shortfall from it is not significant at one-sided P < 0.05 divided by the number of other cell types, from paired jackknife differences.",
          sprintf("Highlighted and labelled: up to %d cell types per GWAS, first those that anchor its top set (the most enriched of those at q < 0.05) at the most sizes, then those in its top set at the most sizes, then the most enriched.", n_labels),
          sprintf("Top sets are wide: a median of %s cell types per GWAS and size, over all categories.", median_top_set_size),
          paste0(
            sprintf("All categories share the y scale from %s to %s", limits[[1]], limits[[2]]),
            if (is.null(y_limits)) sprintf(
              ": the %s%% and %s%% quantiles of all cell-type estimates, with the lower end at most 0, rounded outward to multiples of 5",
              100 * cap_quantile, 100 * (1 - cap_quantile)
            ),
            ". Estimates beyond it are drawn at its ends as triangles labelled with their value, and intervals are cut there."
          ),
          "Estimates below 0, which no true enrichment can take, come from noisy annotations.",
          if (length(category_without_sumstats)) paste0("Not shown, without usable genome-wide summary statistics: ", toString(sort(category_without_sumstats)), "."),
          describe_SLDSC_dropped_annotations(dplyr::filter(dropped_tibble, strategy == "accessible_peaks"))
        ), width = 165),
        x = "Annotation size (% of common SNPs)",
        y = "Heritability enrichment"
      ) +
      ggplot2::theme_bw(base_size = 9) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(size = 7),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_line(colour = "grey93", linewidth = 0.25),
        legend.position = "bottom",
        plot.caption = ggplot2::element_text(hjust = 0)
      )
  })
}

#' Plot enrichment in each cell type's most accessible peaks as heatmaps
#'
#' One heatmap per size of the `accessible_peaks` strategy, laid out as the
#' specificity heatmap: every GWAS and cell type of the peak model, with
#' cell types that cannot fill the size left missing. All sizes share the
#' statistic's limits from `get_SLDSC_enrichment_limits()`, which for S-LDSC
#' are those of the enrichment curves, and one cell-type order, from
#' clustering the statistic at all sizes within compartments. Stars mark
#' enrichment above 1 at q < 0.05 and outlines the top sets, from
#' `test_SLDSC_matched_enrichment()`.
#'
#' @param enrichment_tibble Matched-size enrichment from
#'   `test_SLDSC_matched_enrichment()`.
#' @param SLDSC_tibble Combined cell-type S-LDSC results of the peak model,
#'   which give the GWAS and cell types with their metadata.
#' @param statistic `SLDSC` for heritability enrichment or `PIP` for
#'   fine-mapped PIP enrichment.
#' @inheritParams plot_SLDSC_cell_type_heatmap
#' @inheritParams plot_SLDSC_matched_enrichment_curves
#' @return A list of annotated heatmaps named by size.
#' @keywords internal

plot_SLDSC_enrichment_heatmaps <- function(enrichment_tibble, SLDSC_tibble, GWAS_inputs_tibble, compartments_patterns, dropped_tibble,
                                           statistic = c("SLDSC", "PIP")) {
  statistic <- match.arg(statistic)
  is_SLDSC <- statistic == "SLDSC"
  estimate_col <- paste0(statistic, "_enrichment")
  limits <- get_SLDSC_enrichment_limits(enrichment_tibble, estimate_col = estimate_col)
  accessible_tibble <- enrichment_tibble |>
    dplyr::filter(strategy == "accessible_peaks", cluster != "random") |>
    dplyr::mutate(
      estimate = .data[[estimate_col]],
      # Stars mark only FDR-significant enrichment.
      significance_label = get_significance_stars(NA, if (is_SLDSC) q else PIP_q),
      in_top_set = if (is_SLDSC) in_top_set else PIP_in_top_set,
      z_below_anchor = if (is_SLDSC) z_below_anchor else PIP_z_below_anchor
    )
  sizes <- sort(unique(accessible_tibble$share))
  heatmap_tibbles <- purrr::map(sizes, \(size) {
    SLDSC_tibble |>
      # Without h2, PIP heatmaps draw every GWAS and add no note on summary statistics.
      dplyr::select(GWAS_ID, Category, cluster, dplyr::any_of(if (is_SLDSC) c("h2", "h2_se")), sample_size, dplyr::starts_with("ancestry_"), n_cells) |>
      dplyr::left_join(
        dplyr::select(dplyr::filter(accessible_tibble, share == size), GWAS_ID, cluster, estimate, significance_label, in_top_set, z_below_anchor),
        by = c("GWAS_ID", "cluster"),
        relationship = "one-to-one"
      ) |>
      dplyr::mutate(fill = pmin(pmax(estimate, limits[[1]]), limits[[2]]))
  })
  # One order for all sizes, in which cell types that cannot fill a size count as not enriched.
  cluster_levels <- dplyr::bind_rows(heatmap_tibbles, .id = "size") |>
    dplyr::mutate(GWAS_ID = paste(GWAS_ID, size), fill = dplyr::coalesce(fill, 1)) |>
    get_ordered_GWAS_score_plot_data(compartments_patterns, score_col = "fill") |>
    purrr::pluck("cluster_levels")
  untested <- sort(unique(dplyr::filter(accessible_tibble, is.na(PIP_z))$GWAS_ID))
  purrr::map2(heatmap_tibbles, sizes, \(heatmap_tibble, size) {
    size_label <- scales::percent(size, accuracy = 0.1, drop0trailing = TRUE)
    not_measurable <- dplyr::filter(dropped_tibble, strategy == "accessible_peaks", share == size)$cluster
    annotation_sentence <- sprintf("Each cell type's annotation keeps its peaks with the highest depth-normalized pseudobulk accessibility, among peaks with at least 10 reads in it, until they hold %s of common (MAF > 5%%) 1000 Genomes European SNPs.", size_label)
    limits_sentence <- sprintf(
      "Colour: linear around 1; values beyond %s and %s take the end colours. %s share these limits: the 1%% and 99%% quantiles of all cell-type estimates at all sizes, with the lower end at most 0, rounded outward to multiples of 5.",
      limits[[1]], limits[[2]],
      if (is_SLDSC) "The enrichment curves and the heatmaps of all sizes" else "The fine-mapped heatmaps of all sizes"
    )
    not_measurable_sentence <- if (length(not_measurable)) paste0("Grey columns cannot be measured at this size, since the cell type's peaks with at least 10 reads hold fewer common SNPs: ", toString(sort(not_measurable)), ".")
    text <- if (is_SLDSC) list(
      fill_label = "Heritability enrichment",
      title = paste("Heritability enrichment in each cell type's most accessible peaks, sized to", size_label, "of common SNPs (S-LDSC)"),
      subtitle = paste(
        "Red marks cell types whose most accessible peaks hold more heritability per common variant than a random variant (enrichment 1, white); all sizes share this scale.",
        "Stars mark FDR-significant enrichment above 1 (**); a black outline marks the most enriched of these cell types and grey outlines the cell types not significantly below it.",
        sep = "\n"
      ),
      caption = c(
        annotation_sentence,
        "Enrichment: the annotation's share of S-LDSC heritability over its share of common SNPs, in a model with the 97 baseline-LD v2.2 annotations fitted to HapMap3 SNP summary statistics from the GWAS Catalog or a configured release, using GRCh38 1000 Genomes European LD; SNPs with chi-square above 80 are excluded.",
        limits_sentence,
        "Estimates below 0, which no true enrichment can take, come from noisy annotations.",
        "Stars: ** enrichment above 1 at Benjamini-Hochberg q < 0.05, from a one-sided z-test with its 200-block jackknife standard error, over all GWAS and cell types at this size.",
        "Outlines: in GWAS with a cell type at q < 0.05, black marks the most enriched of these and grey every cell type whose shortfall from it is not significant at one-sided P < 0.05 divided by the number of other cell types, from paired jackknife differences.",
        not_measurable_sentence
      )
    ) else list(
      fill_label = "Fine-mapped PIP enrichment",
      title = paste("Fine-mapped variant enrichment in each cell type's most accessible peaks, sized to", size_label, "of common SNPs"),
      subtitle = paste(
        "Red marks cell types whose most accessible peaks hold more fine-mapped posterior probability per common variant than a random variant (enrichment 1, white); all sizes share this scale.",
        "Stars mark FDR-significant enrichment above 1 (**); a black outline marks the most enriched of these cell types and grey outlines the cell types not significantly below it.",
        sep = "\n"
      ),
      caption = c(
        annotation_sentence,
        "Enrichment: the annotation's share of the posterior inclusion probability (PIP) of autosomal credible-set variants, placed by position, over its share of common SNPs.",
        limits_sentence,
        "Stars: ** enrichment above 1 at Benjamini-Hochberg q < 0.05, from a one-sided z-test with the standard error of 1,000 bootstrap samples of loci (credible sets whose centres lie within 500 kb of each other), over all GWAS and cell types at this size.",
        "Outlines: in GWAS with a cell type at q < 0.05, black marks the most enriched of these and grey every cell type whose shortfall from it is not significant at one-sided P < 0.05 divided by the number of other cell types, from the paired bootstrap samples.",
        "With enough loci, most cell types' most accessible peaks hold more fine-mapped probability than a random variant, so the outlines say more than the stars.",
        "Only genome-wide significant loci are fine-mapped, and a small annotation holds only the credible-set variants inside it, so the enrichment reflects the strongest loci and the fine-mapping resolution.",
        if (length(untested)) paste0("Not tested, with fewer than 10 loci: ", toString(untested), "."),
        not_measurable_sentence
      )
    )
    plot_SLDSC_cell_type_heatmap(
      heatmap_tibble,
      GWAS_inputs_tibble,
      compartments_patterns,
      fill_label = text$fill_label,
      title = text$title,
      subtitle = text$subtitle,
      caption = text$caption,
      fill_limits = limits,
      fill_midpoint = 1,
      cluster_levels = cluster_levels
    ) +
      # save_plots_structured() titles a listed plot without a title by its name, and reads a
      # patchwork's title from its last panel, so that panel gets a hidden one.
      ggplot2::labs(title = size_label) + ggplot2::theme(plot.title = ggplot2::element_blank())
  }) |>
    rlang::set_names(paste0(100 * sizes, "_percent"))
}

#' Plot the credible sets behind each cell type's fine-mapped PIP enrichment
#'
#' One heatmap per GWAS of the contributions from
#' `get_SLDSC_PIP_locus_contribution_tibble()`: cell types, grouped by
#' compartment, against the credible sets that most set a cell type apart,
#' with the rest summed, one panel per size, and each cell type's total
#' enrichment as a number.
#'
#' @param locus_contribution_tibble Contributions from
#'   `get_SLDSC_PIP_locus_contribution_tibble()`.
#' @param enrichment_tibble Matched-size enrichment, which gives each cell
#'   type's total and the sizes it fills.
#' @param locus_to_gene_tibble Ranked Open Targets L2G predictions.
#' @param compartments_patterns Cell-type compartment patterns.
#' @param n_top_loci Credible sets shown individually per GWAS.
#' @return A list of heatmaps named by GWAS.
#' @keywords internal

plot_SLDSC_PIP_locus_contribution_heatmaps <- function(locus_contribution_tibble, enrichment_tibble, locus_to_gene_tibble,
                                                       compartments_patterns, n_top_loci = 12L) {
  total_tibble <- enrichment_tibble |>
    dplyr::filter(strategy == "accessible_peaks", cluster != "random") |>
    dplyr::select(GWAS_ID, cluster, share, PIP_enrichment)
  gene_tibble <- dplyr::filter(locus_to_gene_tibble, L2G_rank == 1L) |>
    dplyr::select(GWAS_ID, studyLocusId, top_L2G_gene = gene_symbol)
  size_label <- \(share) scales::percent(share, accuracy = 0.1, drop0trailing = TRUE)
  GWAS_IDs <- intersect(unique(total_tibble$GWAS_ID), unique(locus_contribution_tibble$GWAS_ID))
  purrr::map(rlang::set_names(GWAS_IDs), \(GWAS_ID) {
    locus_tibble <- locus_contribution_tibble |>
      dplyr::filter(.data$GWAS_ID == .env$GWAS_ID) |>
      dplyr::left_join(gene_tibble, by = c("GWAS_ID", "studyLocusId"), relationship = "many-to-one") |>
      dplyr::mutate(label = dplyr::if_else(is.na(top_L2G_gene), lead_variantId, paste0(lead_variantId, " (", top_L2G_gene, ")")))
    cells <- total_tibble |>
      dplyr::filter(.data$GWAS_ID == .env$GWAS_ID) |>
      dplyr::mutate(compartment = get_cluster_compartments(cluster, compartments_patterns), size = size_label(share))
    # Credible sets that set a cell type apart come first: the largest excess of a cell type's contribution over
    # the median across the cell types that fill the size, counting those without the credible set as zero.
    n_cells <- table(cells$share)
    top_tibble <- locus_tibble |>
      dplyr::mutate(
        excess = PIP_enrichment_contribution - stats::median(c(PIP_enrichment_contribution, numeric(n_cells[[as.character(share[[1]])]] - dplyr::n()))),
        .by = c(studyLocusId, share)
      ) |>
      dplyr::summarise(excess = max(excess), largest = max(PIP_enrichment_contribution), label = dplyr::first(label), .by = studyLocusId) |>
      dplyr::arrange(dplyr::desc(excess), dplyr::desc(largest)) |>
      dplyr::slice_head(n = n_top_loci)
    column_levels <- c(top_tibble$label, "Other loci", "Total")
    cluster_levels <- cells |>
      dplyr::summarise(largest = max(PIP_enrichment), .by = c(compartment, cluster)) |>
      dplyr::arrange(compartment, dplyr::desc(largest)) |>
      dplyr::pull(cluster)
    contribution_tibble <- locus_tibble |>
      dplyr::filter(studyLocusId %in% top_tibble$studyLocusId) |>
      dplyr::transmute(cluster, share, column = top_tibble$label[match(studyLocusId, top_tibble$studyLocusId)], contribution = PIP_enrichment_contribution) |>
      dplyr::right_join(
        tidyr::expand_grid(dplyr::select(cells, cluster, share, size, compartment), column = top_tibble$label),
        by = c("cluster", "share", "column")
      ) |>
      dplyr::mutate(contribution = dplyr::coalesce(contribution, 0))
    # The remaining loci and the total are numbers, so that the colours span the individual credible sets.
    top_sums <- dplyr::summarise(contribution_tibble, top = sum(contribution), .by = c(cluster, share))
    number_tibble <- cells |>
      dplyr::left_join(top_sums, by = c("cluster", "share")) |>
      dplyr::mutate(top = dplyr::coalesce(top, 0))
    number_tibble <- dplyr::bind_rows(
      dplyr::mutate(number_tibble, column = "Other loci", value = PIP_enrichment - top),
      dplyr::mutate(number_tibble, column = "Total", value = PIP_enrichment)
    ) |>
      dplyr::mutate(label = sprintf("%.1f", value))
    plot_factors <- \(data) dplyr::mutate(data,
      column = factor(column, column_levels),
      cluster = factor(cluster, rev(cluster_levels)),
      size = factor(size, size_label(sort(unique(cells$share))))
    )
    number_cells <- plot_factors(number_tibble)
    ggplot2::ggplot(plot_factors(contribution_tibble), ggplot2::aes(column, cluster)) +
      ggplot2::geom_tile(ggplot2::aes(fill = contribution), colour = "grey90", linewidth = 0.2) +
      ggplot2::geom_tile(data = number_cells, fill = "grey96", colour = "grey90", linewidth = 0.2) +
      ggplot2::geom_text(data = number_cells, ggplot2::aes(label = label), size = 2.2, colour = "grey15") +
      ggplot2::scale_fill_gradient(low = "white", high = "#B40426", name = "Contribution to\nPIP enrichment") +
      ggplot2::scale_x_discrete(drop = FALSE) +
      ggplot2::facet_grid(
        rows = ggplot2::vars(compartment), cols = ggplot2::vars(size), scales = "free_y", space = "free_y",
        labeller = ggplot2::labeller(compartment = \(x) stringr::str_replace_all(x, "_", " "))
      ) +
      ggplot2::labs(
        title = paste("Credible sets behind fine-mapped variant enrichment in each cell type's most accessible peaks:", GWAS_ID),
        subtitle = "Look for cell types whose enrichment rests on one or two credible sets rather than many; Total gives the enrichment, where 1 is a random variant.",
        x = "Credible set (lead variant and top Open Targets L2G gene)",
        y = NULL,
        caption = stringr::str_wrap(paste(
          "A credible set's contribution is its posterior inclusion probability (PIP) inside the cell type's most accessible peaks of that size,",
          "over the PIP of all autosomal credible sets and the annotation's share of common (MAF > 5%) SNPs, so a cell type's contributions sum to its enrichment (Total).",
          sprintf("Coloured: the %d credible sets that most set a cell type apart, by the largest excess of a cell type's contribution over the median across cell types at any size; Other loci gives the sum of the rest, including credible sets accessible in most cell types.", n_top_loci),
          "A cell type without a row at a size cannot fill it. Gene labels are Open Targets locus-to-gene prioritizations, not causal assignments."
        ), width = 170)
      ) +
      ggplot2::theme_minimal(base_size = 8) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 60, hjust = 1),
        panel.grid = ggplot2::element_blank(),
        strip.text.y = ggplot2::element_text(angle = 0, hjust = 0),
        legend.position = "bottom",
        plot.caption = ggplot2::element_text(hjust = 0)
      )
  })
}

#' Plot matched-size enrichment against annotation size
#'
#' @param summary_tibble Summary from `summarise_SLDSC_matched_enrichment()`.
#' @param dropped_tibble Annotations that cell types cannot fill, with their
#'   `strategy`, `cluster` and `share`.
#' @return Patchwork of median enrichment and specificity curves.
#' @keywords internal

plot_SLDSC_matched_enrichment_curves <- function(summary_tibble, dropped_tibble) {
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
  share_axis <- ggplot2::scale_x_continuous(breaks = unique(summary_tibble$share), labels = scales::label_percent(accuracy = 0.1))
  curve_layers <- list(
    ggplot2::geom_hline(yintercept = 1, colour = "grey60", linewidth = 0.3),
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper, fill = strategy, group = interaction(strategy, annotation)), alpha = 0.15, colour = NA),
    ggplot2::geom_line(ggplot2::aes(colour = strategy, linetype = annotation)),
    ggplot2::geom_point(ggplot2::aes(colour = strategy), size = 1),
    share_axis,
    SLDSC_strategy_scale(c("colour", "fill")),
    ggplot2::scale_linetype_manual(values = c("Cell-type annotation" = "solid", "Random reference" = "dashed"), name = NULL),
    ggplot2::theme_bw(base_size = 9)
  )
  enrichment_plot <- ggplot2::ggplot(curve_tibble, ggplot2::aes(share, median)) +
    curve_layers +
    ggplot2::facet_grid(ggplot2::vars(statistic), ggplot2::vars(pair_set)) +
    ggplot2::labs(subtitle = "Median enrichment", x = "Annotation size (share of common SNPs)", y = "Enrichment")
  specificity_plot <- dplyr::filter(summary_tibble, pair_set == "Specificity") |>
    dplyr::mutate(annotation = "Cell-type annotation") |>
    ggplot2::ggplot(ggplot2::aes(share, median)) +
    curve_layers +
    ggplot2::facet_grid(ggplot2::vars(statistic)) +
    ggplot2::guides(linetype = "none") +
    ggplot2::labs(subtitle = "Specificity: enriched over other cell types", x = "Annotation size (share of common SNPs)", y = "Ratio of median enrichments")
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
        "Each annotation keeps a cell type's highest-scoring peaks (depth-normalized accessibility, or its share over cell types) or gene windows (expression specificity, gene body +/- 100 kb) until they hold the given share of common (MAF > 5%) 1000 Genomes European SNPs.",
        "Random reference: the mean of 10 random rankings of the same peaks or gene windows, so both peak strategies share it.",
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
        "Lines: medians over pairs or GWAS; bands: 95% intervals from 1,000 bootstrap samples of GWAS.",
        describe_SLDSC_dropped_annotations(dropped_tibble)
      ), width = 180),
      theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0))
    ) &
    ggplot2::theme(legend.position = "bottom")
}

#' Compare annotation strategies per cell type at one matched size
#'
#' @param summary_tibble Summary from `summarise_SLDSC_matched_enrichment()`.
#' @param enrichment_tibble Combined matched-size enrichment of all GWAS.
#' @inheritParams plot_SLDSC_matched_enrichment_curves
#' @param share Annotation size to compare, as a share of common SNPs.
#' @return Dot plot of each cell type's median enrichment per strategy.
#' @keywords internal

plot_SLDSC_matched_size_comparison <- function(summary_tibble, enrichment_tibble, dropped_tibble, share = 0.01) {
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
        "Cell types are ordered by their highest median S-LDSC enrichment over all pairs; medians at or below zero cannot be shown on the log scale.",
        describe_SLDSC_dropped_annotations(dplyr::filter(dropped_tibble, share == matched_share))
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

plot_SLDSC_matched_enrichment_agreement <- function(enrichment_tibble, dropped_tibble, share = 0.01) {
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
    dplyr::filter(strategy %in% c("specific_peaks", "specific_genes")) |>
    dplyr::select(GWAS_ID, cluster, pair_set, strategy, SLDSC_enrichment, PIP_enrichment) |>
    tidyr::pivot_wider(names_from = strategy, values_from = c(SLDSC_enrichment, PIP_enrichment)) |>
    dplyr::filter(dplyr::if_all(dplyr::matches("_enrichment_"), \(x) x > 0)) |>
    dplyr::mutate(
      SLDSC_ratio = log2(SLDSC_enrichment_specific_peaks / SLDSC_enrichment_specific_genes),
      PIP_ratio = log2(PIP_enrichment_specific_peaks / PIP_enrichment_specific_genes)
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
      subtitle = "Most specific peaks over most specific genes per GWAS-cell-type pair",
      x = "S-LDSC: log2 ratio, specific peaks over genes",
      y = "Fine-mapped PIP: log2 ratio, specific peaks over genes"
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
        "Enriched pairs reach q < 0.05 in the single-modality peak or gene S-LDSC model. Numbers in B count the pairs in each quadrant.",
        describe_SLDSC_dropped_annotations(dplyr::filter(dropped_tibble, share == matched_share))
      ), width = 170),
      tag_levels = "A",
      theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0))
    ) &
    ggplot2::theme(legend.position = "bottom", legend.box = "vertical")
}
