#' Read the published human CollecTRI network
#'
#' Validate the checksum-pinned CollecTRI snapshot and convert it to the schema
#' consumed by the TF activity targets.
#'
#' @param network_csv Path to the checksum-verified CollecTRI CSV download.
#' @return A signed source-target network tibble with CollecTRI provenance.
#' @keywords internal

read_CollecTRI_human_network <- function(network_csv) {
  # The checksum-pinned download fixes the columns and content.
  readr::read_csv(
    network_csv,
    col_types = readr::cols(.default = readr::col_character(), weight = readr::col_double()),
    progress = FALSE
  ) |>
    dplyr::transmute(
      source = stringr::str_to_upper(source),
      target = stringr::str_to_upper(target),
      mor = weight,
      TF_category = TF.category,
      resources,
      PMIDs = PMID,
      sign_decision = sign.decision
    )
}

#' Get regulon target enrichment results
#'
#' Test each CollecTRI regulon as a signed gene set. For every contrast of a
#' gene-expression model, `limma::cameraPR()` compares the targets' statistics,
#' with the repressed targets' (`mor < 0`) flipped, against all other tested
#' genes. Its variance comes from the spread of the gene-level statistics and
#' the set's inter-gene correlation, which camera estimates for the signed set
#' and which is floored at 0. So, unlike a per-sample regulon score, an
#' arbitrary set is not called merely because many genes change. Each class fit
#' sets the gene universe. A significant regulon means that its targets moved
#' together in the expected directions, not that the TF itself changed, so each
#' result names its leading targets and the TF's own expression.
#'
#' @param pseudobulk_feature_matrix_fit Fitted gene-expression model.
#' @param pseudobulk_feature_dynamic_tibble Dynamic-branch metadata row describing the model.
#' @param results_tibble The model's gene results, whose moderated t is tested.
#' @param pseudobulk_counts_matrix The model's gene-by-pseudobulk counts, for the TFs' expression.
#' @param CollecTRI_network_tibble Signed CollecTRI network with `source`, `target` and `mor` columns.
#' @param organism_chr Pipeline organism identifier; the network is human.
#' @param min_targets Minimum number of tested targets per regulon.
#' @param n_leading_targets Number of leading targets reported per regulon.
#' @param min_TF_CPM TF expression below which a regulon is flagged.
#' @return One row per contrast and regulon, with cameraPR's `PValue`, BH `FDR`
#'   within the contrast, a signed normal quantile `z`, the leading targets, and
#'   `TF_CPM`, the TF's lowest pooled expression across the contrast's cell types.
#' @keywords internal

get_regulon_target_enrichment_results <- function(
  pseudobulk_feature_matrix_fit,
  pseudobulk_feature_dynamic_tibble,
  results_tibble,
  pseudobulk_counts_matrix,
  CollecTRI_network_tibble,
  organism_chr,
  min_targets = 10L,
  n_leading_targets = 5L,
  min_TF_CPM = 0.5
) {
  if (!identical(organism_chr, "Homo_sapiens")) {
    stop("Regulon target enrichment uses the human CollecTRI network, so it supports only Homo_sapiens aggregations, not ",
      organism_chr, ".", call. = FALSE)
  }
  model_list <- pseudobulk_feature_dynamic_tibble$model[[1]]
  pairing_variable <- model_list$pairing_variable %||% model_list$random_effect %||% "donor_id"
  regulons <- split(stats::setNames(sign(CollecTRI_network_tibble$mor), CollecTRI_network_tibble$target), CollecTRI_network_tibble$source)

  # Each TF's expression, pooled over each cell type's pseudobulks.
  cell_type_counts <- pseudobulk_counts_matrix %*% Matrix::t(Matrix::fac2sparse(factor(get_pseudobulk_sample_tibble(pseudobulk_counts_matrix)$cluster)))
  TF_rows <- match(names(regulons), toupper(rownames(cell_type_counts)))
  TF_CPM <- 1e6 * t(t(as.matrix(cell_type_counts[TF_rows[!is.na(TF_rows)], , drop = FALSE])) / Matrix::colSums(cell_type_counts))
  TF_CPM <- TF_CPM[match(names(regulons), names(regulons)[!is.na(TF_rows)]), , drop = FALSE]
  rownames(TF_CPM) <- names(regulons)

  contrast_vec_list <- get_contrast_vec_list(pseudobulk_feature_matrix_fit, model_list)
  purrr::imap_dfr(contrast_vec_list, function(contrast_vec, contrast_name) {
    test <- get_set_test_inputs(pseudobulk_feature_matrix_fit, contrast_vec, pairing_variable)
    statistic <- get_contrast_statistic(results_tibble, contrast_name, test$genes)
    tested <- purrr::map(regulons, function(signs) {
      index <- match(names(signs), toupper(test$genes))
      list(index = index[!is.na(index)], sign = unname(signs[!is.na(index)]))
    }) |>
      purrr::keep(\(regulon) length(regulon$index) >= min_targets)
    if (!length(tested)) {
      return(tibble::tibble())
    }
    # Negating a gene's expression negates its residuals, so camera estimates a
    # signed set's correlation as that of a plain set in the expression stacked
    # over its negation.
    elist <- test$elist[test$genes, ]
    n_genes <- length(test$genes)
    stacked <- new("EList", list(E = rbind(elist$E, -elist$E), weights = if (!is.null(elist$weights)) rbind(elist$weights, elist$weights)))
    correlation <- limma::camera(stacked, purrr::map(tested, \(regulon) regulon$index + n_genes * (regulon$sign < 0)),
      design = test$design, contrast = test$contrast, inter.gene.cor = NA, sort = FALSE)$Correlation |>
      pmax(0) |>
      stats::setNames(names(tested))

    purrr::imap_dfr(tested, function(regulon, TF) {
      signed <- statistic
      signed[regulon$index] <- regulon$sign * statistic[regulon$index]
      set_result <- limma::cameraPR(signed, list(regulon$index), inter.gene.cor = correlation[[TF]], sort = FALSE)
      leading <- regulon$index[order(if (set_result$Direction == "Up") -signed[regulon$index] else signed[regulon$index])]
      tibble::tibble(ID = TF, NGenes = length(regulon$index), Correlation = correlation[[TF]],
        Direction = set_result$Direction, PValue = set_result$PValue,
        leading_targets = paste(test$genes[utils::head(leading, n_leading_targets)], collapse = ";"),
        TF_CPM = min(TF_CPM[TF, test$cell_types]))
    }) |>
      format_set_test_results(contrast_name, "cameraPR", pseudobulk_feature_dynamic_tibble) |>
      dplyr::mutate(
        z = dplyr::if_else(Direction == "Up", 1, -1) * stats::qnorm(PValue / 2, lower.tail = FALSE),
        low_TF_expression = TF_CPM < min_TF_CPM,
        .after = PValue
      )
  })
}


#' Map CollecTRI regulators to JASPAR motif families
#'
#' Match individual TF symbols directly and associate the retained AP1 and NFKB
#' complex regulons with motif families represented by their canonical members.
#'
#' @param CollecTRI_network_tibble CollecTRI network target.
#' @param JASPAR_motif_family_members_tibble JASPAR motif-family membership target.
#' @return One row per CollecTRI source and mapped JASPAR motif family.
#' @keywords internal

get_CollecTRI_JASPAR_family_map <- function(
  CollecTRI_network_tibble,
  JASPAR_motif_family_members_tibble
) {
  sources <- CollecTRI_network_tibble |>
    dplyr::distinct(source, TF_category)
  motif_members <- JASPAR_motif_family_members_tibble |>
    dplyr::transmute(
      motif_family,
      motif_id,
      TF_name = stringr::str_to_upper(TF_name)
    )

  exact_matches <- sources |>
    dplyr::inner_join(motif_members, by = c("source" = "TF_name")) |>
    dplyr::mutate(mapping_type = "exact TF symbol")

  complex_members <- tibble::tribble(
    ~source, ~TF_name,
    "AP1", "FOS",
    "AP1", "FOSB",
    "AP1", "FOSL1",
    "AP1", "FOSL2",
    "AP1", "JUN",
    "AP1", "JUNB",
    "AP1", "JUND",
    "NFKB", "NFKB1",
    "NFKB", "NFKB2",
    "NFKB", "REL",
    "NFKB", "RELA",
    "NFKB", "RELB"
  )
  complex_matches <- complex_members |>
    dplyr::inner_join(sources, by = "source") |>
    dplyr::inner_join(motif_members, by = "TF_name") |>
    dplyr::mutate(mapping_type = "TF complex members")

  dplyr::bind_rows(exact_matches, complex_matches) |>
    dplyr::summarise(
      TF_category = dplyr::first(TF_category),
      mapping_type = dplyr::first(mapping_type),
      mapped_TF_names = paste(sort(unique(TF_name)), collapse = ";"),
      motif_ids = paste(sort(unique(motif_id)), collapse = ";"),
      n_motif_profiles = dplyr::n_distinct(motif_id),
      .by = c(source, motif_family)
    ) |>
    dplyr::arrange(source, motif_family)
}


#' Correlate CollecTRI regulators' own motifs with their family roots
#'
#' A family root can follow another sub-family than a member's own motif, as the
#' forkhead root of family 007 does for MEF2. Each own motif's and family root's
#' betterChromVAR deviations are computed on the pseudobulks summed per cell type
#' and correlated (Spearman) across cell types. A regulator is compared through
#' its own motif where that correlation is below `min_root_correlation`.
#'
#' @param pseudobulk_ATAC_counts_matrix Peak-by-pseudobulk ATAC counts, with samples named `<cell type>_<donor>`.
#' @param chromVAR_obj Single-nucleus chromVAR object providing peak ranges and GC bias.
#' @param TF_motif_matrix Peak-by-motif matches of the regulators' own motifs, named by motif feature.
#' @param motif_family_matrix Peak-by-family matches of the family roots.
#' @param motif_family_members_tibble JASPAR motif-family membership target.
#' @param min_root_correlation Root correlation below which a regulator's own motif is used.
#' @return One row per own motif, with its family, `root_correlation` and `use_own_motif`.
#' @keywords internal

get_TF_motif_root_correlation_tibble <- function(
  pseudobulk_ATAC_counts_matrix,
  chromVAR_obj,
  TF_motif_matrix,
  motif_family_matrix,
  motif_family_members_tibble,
  min_root_correlation = 0.5
) {
  cell_type <- factor(get_pseudobulk_sample_tibble(pseudobulk_ATAC_counts_matrix)$cluster)
  cell_type_counts <- methods::as(pseudobulk_ATAC_counts_matrix, "dgCMatrix") %*% Matrix::t(Matrix::fac2sparse(cell_type))
  colnames(cell_type_counts) <- levels(cell_type)
  deviations <- get_pseudobulk_motif_family_accessibility_matrix(
    pseudobulk_ATAC_data_matrix = cell_type_counts,
    chromVAR_obj = chromVAR_obj,
    chromVAR_motif_family_matrix = cbind(motif_family_matrix, TF_motif_matrix[rownames(motif_family_matrix), , drop = FALSE])
  )
  motif_family_members_tibble |>
    dplyr::filter(motif_feature %in% colnames(TF_motif_matrix)) |>
    dplyr::transmute(
      motif_feature,
      motif_id,
      motif_family,
      root_correlation = purrr::map2_dbl(motif_feature, motif_family, \(motif, family) {
        stats::cor(deviations[motif, ], deviations[family, ], method = "spearman")
      }),
      use_own_motif = !is.na(root_correlation) & root_correlation < min_root_correlation
    )
}


#' Match regulon target enrichment with motif accessibility and TF expression
#'
#' Each regulon is compared with its family root's accessibility, or with its own
#' motif's where `get_TF_motif_root_correlation_tibble()` chose it.
#'
#' @param regulon_results_tibble Regulon target enrichment results.
#' @param motif_family_accessibility_results_tibble Differential JASPAR motif-family accessibility results.
#' @param TF_motif_accessibility_results_tibble Differential accessibility results of the chosen own motifs.
#' @param gene_expression_results_tibble Differential gene-expression results.
#' @param CollecTRI_JASPAR_family_map CollecTRI-to-JASPAR crosswalk.
#' @param TF_motif_root_correlation_tibble Own-motif root correlations and choices.
#' @return A regulator-level cross-modality comparison tibble, with `motif_feature`
#'   naming the family or own motif compared.
#' @keywords internal

get_CollecTRI_JASPAR_comparison_tibble <- function(
  regulon_results_tibble,
  motif_family_accessibility_results_tibble,
  TF_motif_accessibility_results_tibble,
  gene_expression_results_tibble,
  CollecTRI_JASPAR_family_map,
  TF_motif_root_correlation_tibble
) {
  join_keys <- c("model", "contrast", "cell_type_subset")
  regulon_results <- regulon_results_tibble |>
    dplyr::bind_rows() |>
    dplyr::transmute(
      dplyr::across(dplyr::all_of(join_keys)),
      source = ID,
      z_regulon = z,
      PValue_regulon = PValue,
      FDR_regulon = FDR
    )
  motif_map <- CollecTRI_JASPAR_family_map |>
    tidyr::separate_longer_delim(motif_ids, ";") |>
    dplyr::left_join(
      TF_motif_root_correlation_tibble |>
        dplyr::filter(use_own_motif) |>
        dplyr::select(motif_ids = motif_id, own_motif = motif_feature, root_correlation),
      by = "motif_ids"
    ) |>
    dplyr::mutate(motif_feature = dplyr::coalesce(own_motif, motif_family)) |>
    dplyr::summarise(
      TF_category = dplyr::first(TF_category),
      mapping_type = dplyr::first(mapping_type),
      mapped_TF_names = dplyr::first(mapped_TF_names),
      motif_ids = paste(sort(motif_ids), collapse = ";"),
      root_correlation = dplyr::first(root_correlation),
      .by = c(source, motif_family, motif_feature)
    )
  motif_accessibility_results <- dplyr::bind_rows(c(motif_family_accessibility_results_tibble, TF_motif_accessibility_results_tibble)) |>
    dplyr::transmute(
      dplyr::across(dplyr::all_of(join_keys)),
      motif_feature = feature_id,
      logFC_motif_accessibility = logFC,
      t_motif_accessibility = t,
      PValue_motif_accessibility = PValue,
      FDR_motif_accessibility = FDR
    )
  gene_expression_results <- gene_expression_results_tibble |>
    dplyr::bind_rows()
  gene_expression_t_statistic <- if ("t" %in% names(gene_expression_results)) {
    gene_expression_results$t
  } else {
    rep(NA_real_, nrow(gene_expression_results))
  }
  gene_expression_results <- gene_expression_results |>
    dplyr::transmute(
      dplyr::across(dplyr::all_of(join_keys)),
      source = stringr::str_to_upper(feature_id),
      logFC_TF_expression = logFC,
      t_TF_expression = gene_expression_t_statistic,
      FDR_TF_expression = FDR
    )

  regulon_results |>
    dplyr::left_join(motif_map, by = "source", relationship = "many-to-many") |>
    dplyr::left_join(
      motif_accessibility_results,
      by = c(join_keys, "motif_feature"),
      relationship = "many-to-one"
    ) |>
    dplyr::left_join(
      gene_expression_results,
      by = c(join_keys, "source"),
      relationship = "many-to-one"
    ) |>
    dplyr::mutate(
      motif_mapping_status = dplyr::case_when(
        is.na(motif_feature) ~ "unmapped",
        motif_feature == motif_family ~ "family root",
        .default = "own motif"
      ),
      direction_concordant = dplyr::if_else(
        motif_mapping_status != "unmapped",
        sign(z_regulon) == sign(t_motif_accessibility),
        NA
      ),
      jointly_FDR_significant = FDR_regulon < 0.05 & FDR_motif_accessibility < 0.05
    ) |>
    dplyr::arrange(model, contrast, source, motif_feature)
}


#' Collapse regulon comparisons to motif features
#'
#' @param comparison_tibble Output from `get_CollecTRI_JASPAR_comparison_tibble()`.
#' @return One row per model, contrast, and compared motif family or own motif.
#' @keywords internal

get_CollecTRI_JASPAR_family_comparison_tibble <- function(comparison_tibble) {
  comparison_tibble |>
    dplyr::filter(motif_mapping_status != "unmapped", is.finite(t_motif_accessibility)) |>
    dplyr::summarise(
      z_regulon = stats::median(z_regulon, na.rm = TRUE),
      t_motif_accessibility = dplyr::first(t_motif_accessibility),
      t_TF_expression = if (all(is.na(t_TF_expression))) {
        NA_real_
      } else {
        stats::median(t_TF_expression, na.rm = TRUE)
      },
      min_source_FDR_regulon = min(FDR_regulon, na.rm = TRUE),
      any_source_FDR_significant = any(FDR_regulon < 0.05, na.rm = TRUE),
      FDR_motif_accessibility = dplyr::first(FDR_motif_accessibility),
      n_CollecTRI_sources = dplyr::n_distinct(source),
      CollecTRI_sources = paste(sort(unique(source)), collapse = ";"),
      .by = c(model, contrast, cell_type_subset, motif_feature)
    ) |>
    dplyr::mutate(
      direction_concordant = sign(z_regulon) == sign(t_motif_accessibility),
      jointly_FDR_significant = any_source_FDR_significant & FDR_motif_accessibility < 0.05
    )
}


#' Summarize regulon-motif concordance by contrast
#'
#' `null_rho_95` is the 95% range of a Spearman correlation between unrelated
#' rankings of as many motif features.
#'
#' @param family_comparison_tibble Motif-feature-level CollecTRI-JASPAR comparison.
#' @return Contrast-level mapping, direction, significance, and correlation summary.
#' @keywords internal

get_CollecTRI_JASPAR_concordance_tibble <- function(family_comparison_tibble) {
  family_comparison_tibble |>
    dplyr::summarise(
      n_motif_features = dplyr::n(),
      spearman_rho = if (dplyr::n() >= 3L) {
        stats::cor(z_regulon, t_motif_accessibility, method = "spearman", use = "complete.obs")
      } else {
        NA_real_
      },
      null_rho_95 = if (dplyr::n() >= 3L) stats::qnorm(0.975) / sqrt(dplyr::n() - 1) else NA_real_,
      direction_concordance = mean(direction_concordant, na.rm = TRUE),
      n_regulon_FDR = sum(any_source_FDR_significant, na.rm = TRUE),
      n_motif_accessibility_FDR = sum(FDR_motif_accessibility < 0.05, na.rm = TRUE),
      n_joint_FDR = sum(jointly_FDR_significant, na.rm = TRUE),
      .by = c(model, contrast, cell_type_subset)
    )
}


#' Plot regulon-motif concordance across contrasts
#'
#' @param concordance_tibble Contrast-level CollecTRI-JASPAR summary.
#' @return A list of ggplots split by pseudobulk model.
#' @keywords internal

plot_CollecTRI_JASPAR_concordance <- function(concordance_tibble) {
  concordance_tibble |>
    group_split_by("model") |>
    purrr::imap(\(model_tibble, model_name) {
      family_count_range <- range(model_tibble$n_motif_features)
      family_count_label <- if (family_count_range[[1]] == family_count_range[[2]]) {
        as.character(family_count_range[[1]])
      } else {
        paste(family_count_range, collapse = "-")
      }

      model_tibble |>
        dplyr::mutate(
          contrast_label = contrast |>
            stringr::str_replace_all("_", " ") |>
            stringr::str_wrap(width = 48) |>
            stats::reorder(spearman_rho),
          joint_label = dplyr::if_else(n_joint_FDR > 0L, paste0(n_joint_FDR, " joint"), "")
        ) |>
        ggplot2::ggplot(
          ggplot2::aes(
            x = spearman_rho,
            y = contrast_label,
            color = direction_concordance
          )
        ) +
        ggplot2::geom_linerange(ggplot2::aes(xmin = -null_rho_95, xmax = null_rho_95), color = "grey85", linewidth = 4) +
        ggplot2::geom_vline(xintercept = 0, color = "grey65", linetype = "dashed") +
        ggplot2::geom_point(size = 3) +
        ggplot2::geom_text(ggplot2::aes(label = joint_label), hjust = -0.25, size = 3) +
        ggplot2::scale_color_gradient2(
          low = "#3B4CC0",
          mid = "grey65",
          high = "#B40426",
          midpoint = 0.5,
          limits = c(0, 1),
          labels = scales::label_percent()
        ) +
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.2))) +
        ggplot2::labs(
          title = paste("Regulon-accessibility concordance:", stringr::str_replace_all(model_name, "_", " ")),
          subtitle = stringr::str_wrap("Look for correlations outside the grey band and consistent directions; motif-family agreement does not identify one causal TF.", width = 100),
          x = "Spearman correlation: regulon target enrichment versus motif accessibility",
          y = NULL,
          color = "Direction concordance",
          caption = stringr::str_wrap(paste0(
            "Regulon values are median signed z of regulon target enrichment over the CollecTRI regulons mapped to each of ",
            family_count_label,
            " motif features per contrast: a JASPAR family, or a regulator's own motif where its accessibility correlates below 0.5 with its family root's across cell types. ",
            "Grey bands: 95% range of the correlation between unrelated rankings of as many features. ",
            "Colour: fraction with matching effect directions; joint labels count features significant in both analyses. ",
            "Scores summarize feature-level concordance, not independent donor replication."
          ), width = 110)
        ) +
        ggplot2::theme(legend.position = "bottom", plot.title.position = "plot", plot.caption.position = "plot")
    })
}
