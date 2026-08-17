#' Download the published human CollecTRI network
#'
#' Download the published snapshot from the OmniPath rescue archive and reject
#' any content that does not match the expected SHA-256 checksum.
#'
#' @param network_url URL of the published CollecTRI CSV snapshot.
#' @param expected_sha256 Expected SHA-256 checksum of the downloaded file.
#' @return Normalized path to the checksum-verified CSV under the targets store.
#' @keywords internal

download_CollecTRI_human_network <- function(network_url, expected_sha256) {
  output_dir <- file.path(
    targets::tar_config_get("store"),
    "files",
    "CollecTRI"
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  destfile <- file.path(output_dir, "CollecTRI_human_2023.csv")

  if (
    file.exists(destfile) &&
      identical(
        digest::digest(file = destfile, algo = "sha256"),
        expected_sha256
      )
  ) {
    return(normalizePath(destfile, mustWork = TRUE))
  }

  tempfile <- paste0(destfile, ".tmp")
  if (file.exists(tempfile)) {
    unlink(tempfile)
  }
  on.exit(unlink(tempfile), add = TRUE)

  utils::download.file(network_url, tempfile, mode = "wb")
  if (!file.exists(tempfile) || file.info(tempfile)$size == 0) {
    stop("CollecTRI download failed or produced an empty file: ", network_url)
  }

  observed_sha256 <- digest::digest(file = tempfile, algo = "sha256")
  if (!identical(observed_sha256, expected_sha256)) {
    stop(
      "CollecTRI checksum mismatch: expected ",
      expected_sha256,
      ", observed ",
      observed_sha256,
      "."
    )
  }

  if (file.exists(destfile) && !file.remove(destfile)) {
    stop("Could not replace the existing CollecTRI download: ", destfile)
  }
  if (!file.rename(tempfile, destfile)) {
    stop("Could not move the verified CollecTRI download into place: ", destfile)
  }

  normalizePath(destfile, mustWork = TRUE)
}


#' Read the published human CollecTRI network
#'
#' Validate the checksum-pinned CollecTRI snapshot and convert it to the schema
#' consumed by the TF activity targets.
#'
#' @param network_csv Path returned by `download_CollecTRI_human_network()`.
#' @return A signed source-target network tibble with CollecTRI provenance.
#' @keywords internal

read_CollecTRI_human_network <- function(network_csv) {
  expected_columns <- c(
    "source",
    "target",
    "weight",
    "TF.category",
    "resources",
    "PMID",
    "sign.decision"
  )
  network_tibble <- readr::read_csv(
    network_csv,
    col_types = readr::cols(
      source = readr::col_character(),
      target = readr::col_character(),
      weight = readr::col_double(),
      TF.category = readr::col_character(),
      resources = readr::col_character(),
      PMID = readr::col_character(),
      sign.decision = readr::col_character()
    ),
    progress = FALSE
  )

  if (!identical(names(network_tibble), expected_columns)) {
    stop(
      "The CollecTRI network must contain exactly these columns: ",
      paste(expected_columns, collapse = ", "),
      "."
    )
  }
  if (nrow(readr::problems(network_tibble)) > 0L) {
    stop("The CollecTRI network contains CSV parsing problems.")
  }

  network_tibble <- network_tibble |>
    dplyr::transmute(
      source = stringr::str_to_upper(source),
      target = stringr::str_to_upper(target),
      mor = weight,
      TF_category = TF.category,
      resources,
      PMIDs = PMID,
      sign_decision = sign.decision
    )

  validation <- c(
    interaction_rows = nrow(network_tibble) == 43536L,
    regulator_sources = dplyr::n_distinct(network_tibble$source) == 1189L,
    complete_records = !anyNA(network_tibble),
    unique_source_targets = !anyDuplicated(network_tibble[c("source", "target")]),
    signed_interactions = all(network_tibble$mor %in% c(-1, 1))
  )
  if (!all(validation)) {
    stop(
      "The checksum-pinned CollecTRI network failed: ",
      paste(names(validation)[!validation], collapse = ", "),
      "."
    )
  }

  network_tibble
}


#' Infer CollecTRI TF activities from pseudobulk expression
#'
#' Convert cluster-donor pseudobulk counts to filtered log-CPM values and infer
#' one signed ULM activity score per CollecTRI regulator and sample.
#'
#' @param psbulk_GEX_counts_matrix Gene-by-pseudobulk-sample count matrix.
#' @param CollecTRI_network_tibble Signed CollecTRI network with `source`,
#'   `target`, and `mor` columns.
#' @param min_targets Minimum detected regulon targets required per source.
#' @return A regulator-by-pseudobulk-sample numeric activity matrix.
#' @keywords internal

get_pseudobulk_CollecTRI_TF_activity_matrix <- function(
  psbulk_GEX_counts_matrix,
  CollecTRI_network_tibble,
  min_targets = 5L
) {
  counts_matrix <- if (inherits(psbulk_GEX_counts_matrix, "IterableMatrix")) {
    methods::as(psbulk_GEX_counts_matrix, "dgCMatrix")
  } else {
    psbulk_GEX_counts_matrix
  }

  sample_tibble <- get_psbulk_sample_tibble(counts_matrix)
  DGE_list <- edgeR::DGEList(counts = counts_matrix)
  expressed_features <- edgeR::filterByExpr(
    DGE_list,
    group = sample_tibble$cluster
  )
  DGE_list <- edgeR::normLibSizes(
    DGE_list[expressed_features, , keep.lib.sizes = FALSE]
  )
  log_CPM_matrix <- edgeR::cpm(DGE_list, log = TRUE, prior.count = 2)

  activity_matrix <- decoupleR::run_ulm(
    mat = log_CPM_matrix,
    network = dplyr::select(CollecTRI_network_tibble, source, target, mor),
    minsize = min_targets
  ) |>
    dplyr::select(source, condition, score) |>
    tidyr::pivot_wider(names_from = condition, values_from = score) |>
    tibble::column_to_rownames("source") |>
    as.matrix()

  activity_matrix <- activity_matrix[, colnames(counts_matrix), drop = FALSE]
  assert_with_info(
    nrow(activity_matrix) > 0L && all(is.finite(activity_matrix)),
    glue_info = "CollecTRI ULM did not produce a finite TF activity matrix."
  )
  activity_matrix
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


#' Match expression- and accessibility-derived TF results
#'
#' @param CollecTRI_results_tibble Differential CollecTRI activity results.
#' @param DTFA_results_tibble Differential JASPAR motif-family accessibility results.
#' @param DGE_results_tibble Differential gene-expression results.
#' @param CollecTRI_JASPAR_family_map CollecTRI-to-JASPAR crosswalk.
#' @return A regulator-level cross-modality comparison tibble.
#' @keywords internal

get_CollecTRI_DTFA_comparison_tibble <- function(
  CollecTRI_results_tibble,
  DTFA_results_tibble,
  DGE_results_tibble,
  CollecTRI_JASPAR_family_map
) {
  join_keys <- c("model", "contrast", "cell_type_subset")
  CollecTRI_results <- CollecTRI_results_tibble |>
    dplyr::bind_rows() |>
    dplyr::transmute(
      dplyr::across(dplyr::all_of(join_keys)),
      source = feature_id,
      logFC_CollecTRI = logFC,
      t_CollecTRI = t,
      PValue_CollecTRI = PValue,
      FDR_CollecTRI = FDR
    )
  DTFA_results <- DTFA_results_tibble |>
    dplyr::bind_rows() |>
    dplyr::transmute(
      dplyr::across(dplyr::all_of(join_keys)),
      motif_family = feature_id,
      logFC_DTFA = logFC,
      t_DTFA = t,
      PValue_DTFA = PValue,
      FDR_DTFA = FDR
    )
  DGE_results <- DGE_results_tibble |>
    dplyr::bind_rows() |>
    dplyr::transmute(
      dplyr::across(dplyr::all_of(join_keys)),
      source = stringr::str_to_upper(feature_id),
      logFC_TF_expression = logFC,
      t_TF_expression = t,
      FDR_TF_expression = FDR
    )

  CollecTRI_results |>
    dplyr::left_join(
      CollecTRI_JASPAR_family_map,
      by = "source",
      relationship = "many-to-many"
    ) |>
    dplyr::left_join(
      DTFA_results,
      by = c(join_keys, "motif_family"),
      relationship = "many-to-one"
    ) |>
    dplyr::left_join(
      DGE_results,
      by = c(join_keys, "source"),
      relationship = "many-to-one"
    ) |>
    dplyr::mutate(
      motif_mapping_status = dplyr::if_else(is.na(motif_family), "unmapped", "mapped"),
      direction_concordant = dplyr::if_else(
        motif_mapping_status == "mapped",
        sign(t_CollecTRI) == sign(t_DTFA),
        NA
      ),
      jointly_FDR_significant = FDR_CollecTRI < 0.05 & FDR_DTFA < 0.05
    ) |>
    dplyr::arrange(model, contrast, source, motif_family)
}


#' Collapse CollecTRI comparisons to motif families
#'
#' @param comparison_tibble Output from `get_CollecTRI_DTFA_comparison_tibble()`.
#' @return One row per model, contrast, and mapped motif family.
#' @keywords internal

get_CollecTRI_DTFA_family_comparison_tibble <- function(comparison_tibble) {
  comparison_tibble |>
    dplyr::filter(motif_mapping_status == "mapped", is.finite(t_DTFA)) |>
    dplyr::summarise(
      t_CollecTRI = stats::median(t_CollecTRI, na.rm = TRUE),
      t_DTFA = dplyr::first(t_DTFA),
      t_TF_expression = if (all(is.na(t_TF_expression))) {
        NA_real_
      } else {
        stats::median(t_TF_expression, na.rm = TRUE)
      },
      min_source_FDR_CollecTRI = min(FDR_CollecTRI, na.rm = TRUE),
      any_source_FDR_significant = any(FDR_CollecTRI < 0.05, na.rm = TRUE),
      FDR_DTFA = dplyr::first(FDR_DTFA),
      n_CollecTRI_sources = dplyr::n_distinct(source),
      CollecTRI_sources = paste(sort(unique(source)), collapse = ";"),
      .by = c(model, contrast, cell_type_subset, motif_family)
    ) |>
    dplyr::mutate(
      direction_concordant = sign(t_CollecTRI) == sign(t_DTFA),
      jointly_FDR_significant = any_source_FDR_significant & FDR_DTFA < 0.05
    )
}


#' Summarize CollecTRI-DTFA concordance by contrast
#'
#' @param family_comparison_tibble Family-level CollecTRI-DTFA comparison.
#' @return Contrast-level mapping, direction, significance, and correlation summary.
#' @keywords internal

get_CollecTRI_DTFA_concordance_tibble <- function(family_comparison_tibble) {
  family_comparison_tibble |>
    dplyr::summarise(
      n_motif_families = dplyr::n(),
      spearman_rho = if (dplyr::n() >= 3L) {
        stats::cor(t_CollecTRI, t_DTFA, method = "spearman", use = "complete.obs")
      } else {
        NA_real_
      },
      direction_concordance = mean(direction_concordant, na.rm = TRUE),
      n_CollecTRI_FDR = sum(any_source_FDR_significant, na.rm = TRUE),
      n_DTFA_FDR = sum(FDR_DTFA < 0.05, na.rm = TRUE),
      n_joint_FDR = sum(jointly_FDR_significant, na.rm = TRUE),
      .by = c(model, contrast, cell_type_subset)
    )
}


#' Plot CollecTRI-DTFA concordance across contrasts
#'
#' @param concordance_tibble Contrast-level CollecTRI-DTFA summary.
#' @return A list of ggplots split by pseudobulk model.
#' @keywords internal

plot_CollecTRI_DTFA_concordance <- function(concordance_tibble) {
  concordance_tibble |>
    group_split_by("model") |>
    purrr::imap(\(model_tibble, model_name) {
      family_count_range <- range(model_tibble$n_motif_families)
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
          title = stringr::str_replace_all(model_name, "_", " "),
          x = "Spearman correlation: CollecTRI versus DTFA t-statistics",
          y = NULL,
          color = "Direction concordance",
          caption = paste0(
            "CollecTRI values are median regulator t-statistics for ",
            family_count_label,
            " mapped motif families per contrast"
          )
        ) +
        ggplot2::theme(legend.position = "bottom")
    })
}
