# Shared model/cohort contracts; response construction and fitting stay modality-specific.
normalize_differential_models <- function(models) {
  if (is.null(models)) return(list())
  if (!is.list(models) || is.null(names(models)) || anyNA(names(models)) ||
      any(!nzchar(names(models))) || anyDuplicated(names(models))) {
    stop("Differential models must be a uniquely named list.")
  }
  models <- purrr::discard(models, is.null)
  purrr::imap(models, function(model, name) {
    if (!is.list(model)) stop("Invalid model: ", name)
    contrasts <- model$contrast_specs_vec
    if (!is.null(contrasts) && (is.null(names(contrasts)) || anyNA(names(contrasts)) ||
        any(!nzchar(names(contrasts))) || anyDuplicated(names(contrasts)))) {
      stop("Contrasts must have unique, nonempty names in model: ", name)
    }
    if (is.null(contrasts) && is.null(model$contrast_functions)) {
      stop("No contrasts configured for model: ", name)
    }
    if (!is.null(contrasts)) model$contrast_specs_vec <- unlist(contrasts, use.names = TRUE)
    model
  })
}

get_differential_model_metadata_columns <- function(model) {
  formulas <- purrr::compact(c(model$formula, model$cell_type_formula))
  columns <- unlist(lapply(formulas, function(formula) {
    all.vars(stats::delete.response(stats::terms(stats::as.formula(formula))))
  }), use.names = FALSE)
  unique(setdiff(c(columns, model$pairing_variable, model$correlation_block,
    model$random_effect), "cluster"))
}

get_differential_model_cohort <- function(donor_metadata, model, available_donor_ids) {
  if (!is.list(model)) stop("A model specification is required for cohort selection.")
  required <- get_differential_model_metadata_columns(model)
  missing <- setdiff(required, names(donor_metadata))
  if (length(missing)) stop("Missing model metadata: ", paste(missing, collapse = ", "))
  if (anyNA(donor_metadata$donor_id) || anyDuplicated(donor_metadata$donor_id)) {
    stop("Model metadata requires unique, nonmissing donor IDs.")
  }
  if (length(setdiff(model$donor_ids, donor_metadata$donor_id))) stop("Unknown donor_ids in model.")
  cohort <- tibble::as_tibble(donor_metadata)
  canonical <- function(x) gsub("_", "-", x, fixed = TRUE)
  complete <- if (length(required)) stats::complete.cases(cohort[, required, drop = FALSE]) else rep(TRUE, nrow(cohort))
  cohort$exclusion_reason <- dplyr::case_when(
    !is.null(model$donor_ids) & !cohort$donor_id %in% model$donor_ids ~ "donor_not_selected",
    !canonical(cohort$donor_id) %in% canonical(available_donor_ids) ~ "no_selected_samples",
    !complete ~ "missing_model_metadata",
    .default = NA_character_
  )
  cohort$included <- is.na(cohort$exclusion_reason)
  cohort
}

validate_differential_design <- function(design) {
  if (anyNA(design) || any(!is.finite(design)) || qr(design)$rank < ncol(design)) {
    stop("Model design is nonfinite or rank deficient; check the selected cohort and formula.")
  }
  if (nrow(design) <= ncol(design)) stop("Model has no residual degrees of freedom.")
  invisible(design)
}

prepare_DCTC_model_data <- function(metadata, donor_metadata, model, cluster_col) {
  wells <- model$GEM_well_IDs
  if (!is.null(wells)) {
    missing <- setdiff(wells, unique(metadata$GEM_well_ID))
    if (length(missing)) stop("Selected GEM wells have no retained cells: ", paste(missing, collapse = ", "))
    metadata <- dplyr::filter(metadata, .data$GEM_well_ID %in% wells)
  }
  if (anyNA(metadata[[cluster_col]])) stop("Missing cell-type labels in composition denominator.")
  cohort <- get_differential_model_cohort(donor_metadata, model, unique(metadata$donor_id))
  totals <- metadata |>
    dplyr::group_by(donor_id) |>
    dplyr::summarise(n_total_nuclei = dplyr::n(),
      GEM_well_IDs = paste(sort(unique(GEM_well_ID)), collapse = ";"), .groups = "drop")
  cohort <- dplyr::left_join(cohort, totals, by = "donor_id")
  cohort$n_total_nuclei <- tidyr::replace_na(cohort$n_total_nuclei, 0L)
  selected <- metadata |> dplyr::filter(donor_id %in% cohort$donor_id[cohort$included])
  if (!nrow(selected)) stop("No eligible donors remain for composition model.")
  cell_types <- sort(unique(as.character(selected[[cluster_col]])))
  if (!is.null(model$cell_types_to_test)) {
    if (length(setdiff(model$cell_types_to_test, cell_types))) stop("Requested cell types are absent from the eligible cohort.")
    cell_types <- model$cell_types_to_test
  }
  counts <- selected |>
    dplyr::transmute(donor_id, cluster = as.character(.data[[cluster_col]])) |>
    dplyr::count(donor_id, cluster, name = "n_nuclei")
  # Complete zeros BEFORE selecting response cell types; totals include every label.
  counts <- tidyr::expand_grid(donor_id = sort(unique(selected$donor_id)), cluster = cell_types) |>
    dplyr::left_join(counts, by = c("donor_id", "cluster")) |>
    dplyr::mutate(n_nuclei = tidyr::replace_na(n_nuclei, 0L)) |>
    dplyr::left_join(dplyr::filter(cohort, included), by = "donor_id") |>
    dplyr::mutate(n_other_nuclei = n_total_nuclei - n_nuclei, prop = n_nuclei / n_total_nuclei)
  list(cohort = cohort, counts = counts)
}

fit_DCTC_model <- function(model_data, model, model_name) {
  formula <- stats::as.formula(model$formula)
  if (length(formula) != 2L) {
    stop("Composition models require a one-sided formula, e.g. ~ sexMale; the response is supplied automatically.")
  }
  formula <- stats::update.formula(formula, cbind(n_nuclei, n_other_nuclei) ~ .)
  if (grepl("|", model$formula, fixed = TRUE)) stop("Composition models currently support fixed effects only.")
  if (!is.null(model$design_matrix_func_name) || !is.null(model$contrast_functions) ||
      !is.null(model$random_effect) || !is.null(model$cell_type_formula) ||
      !is.null(model$pairing_variable) || !is.null(model$correlation_block)) {
    stop("Composition models support formula-based fixed effects and literal named contrasts only.")
  }
  results <- model_data$counts |>
    group_split_by("cluster") |>
    purrr::imap(function(data, cell_type) {
      design <- get_design_and_basis_matrices_from_model_list(data, model)$design_matrix
      validate_differential_design(design)
      contrasts <- get_contrast_vec_list(list(design = design), model)
      fit <- tryCatch(glmmTMB::glmmTMB(formula, data = data,
        family = glmmTMB::betabinomial(link = "logit")), error = identity)
      failure <- if (inherits(fit, "error")) conditionMessage(fit) else if (
        fit$fit$convergence != 0L || !isTRUE(fit$sdr$pdHess)) "Model did not converge with a positive-definite Hessian" else NA_character_
      purrr::imap_dfr(contrasts, function(contrast, contrast_name) {
        estimate <- std.error <- NA_real_
        if (is.na(failure)) {
          coefficients <- glmmTMB::fixef(fit)$cond
          covariance <- stats::vcov(fit)$cond
          stopifnot(length(contrast) == length(coefficients))
          estimate <- sum(contrast * coefficients)
          std.error <- sqrt(as.numeric(t(contrast) %*% covariance %*% contrast))
        }
        p <- if (is.finite(std.error) && std.error > 0) 2 * stats::pnorm(-abs(estimate / std.error)) else NA_real_
        tibble::tibble(model = model_name, contrast = contrast_name, cluster = cell_type,
          estimate, std.error, p.value = p, n_donors = nrow(data),
          baseline_frac = sum(data$n_nuclei) / sum(data$n_total_nuclei),
          status = if (is.na(failure) && is.finite(p)) "fitted" else "not_estimable",
          diagnostic = failure)
      })
    }) |> dplyr::bind_rows() |>
    dplyr::group_by(model, contrast) |>
    dplyr::mutate(FDR = stats::p.adjust(p.value, method = "BH")) |>
    dplyr::ungroup()
  results
}

plot_DCTC_model <- function(model_data, model, model_name, results) {
  caption <- paste0("Model: ", model_name, ". Formula: ", model$formula,
    ". Selected wells: ", paste(model$GEM_well_IDs %||% "all aggregation wells", collapse = ", "),
    ". Denominator: all retained nuclei in selected wells for each eligible donor, including unassigned labels. ",
    "Zero cell-type counts retained. Cell types tested: ", paste(model$cell_types_to_test %||% "all observed labels", collapse = ", "), ".")
  plots <- lapply(model$plot_phenotype_vars, function(phenotype) {
    counts <- model_data$counts
    counts$phenotype <- counts[[phenotype]]
    if (is.numeric(counts$phenotype) && length(unique(counts$phenotype)) <= 2L) counts$phenotype <- factor(counts$phenotype)
    plot <- ggplot2::ggplot(counts, ggplot2::aes(phenotype, prop)) +
      ggplot2::geom_point(ggplot2::aes(colour = if (is.null(model$color_by)) NULL else .data[[model$color_by]]), alpha = 0.7) +
      ggplot2::facet_wrap(~cluster, scales = "free_y") +
      ggplot2::scale_y_continuous(labels = scales::label_percent()) +
      ggplot2::labs(title = paste("Cell-type abundance and", phenotype),
        subtitle = "Each point is one donor; compare proportions within the selected population. Facets have separate y scales.",
        x = phenotype, y = "Within-donor proportion", colour = model$color_by,
        caption = stringr::str_wrap(caption, 120))
    plot
  })
  names(plots) <- model$plot_phenotype_vars
  plots$contrasts <- ggplot2::ggplot(dplyr::filter(results, status == "fitted"),
    ggplot2::aes(estimate, cluster, colour = FDR < 0.05)) +
    ggplot2::geom_point() +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = estimate - 1.96 * std.error,
      xmax = estimate + 1.96 * std.error), orientation = "y", width = 0.2) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::facet_wrap(~contrast, scales = "free_x") +
    ggplot2::labs(title = paste("Cell-type abundance contrasts:", model_name),
      subtitle = "Positive effects indicate increased relative abundance; intervals show 95% Wald uncertainty, not absolute cell counts.",
      x = "Log-odds contrast", y = "Cell type", colour = "BH FDR < 0.05",
      caption = stringr::str_wrap(paste(caption,
        "Separate beta-binomial fits per cell type. BH correction across tested cell types within each model/contrast.",
        sum(results$status != "fitted"), "non-estimable contrasts omitted; see the results table."), 120))
  if (!any(results$status == "fitted")) {
    plots$contrasts <- ggplot2::ggplot() + ggplot2::theme_void() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No estimable abundance contrasts; inspect the results table.") +
      ggplot2::labs(title = paste("Cell-type abundance contrasts:", model_name),
        caption = stringr::str_wrap(caption, 120))
  }
  plots
}

save_differential_model_table <- function(table, model_name) {
  table$model <- model_name
  path <- get_structured_file_path(filetype = "tsv", override_suffix = model_name)
  readr::write_tsv(table, path)
  path
}

get_feature_model_cohort <- function(matrix, filtered_matrix, fit, donor_metadata, model, metadata) {
  samples <- get_psbulk_sample_tibble(matrix)
  if (!is.null(model$cell_type_subset)) samples <- dplyr::filter(samples, cluster %in% model$cell_type_subset)
  cohort <- get_differential_model_cohort(donor_metadata, model, samples$donor_id)
  retained <- get_psbulk_sample_tibble(filtered_matrix) |>
    dplyr::count(donor_id, name = "n_samples")
  retained$donor_id <- cohort$donor_id[match(retained$donor_id, gsub("_", "-", cohort$donor_id, fixed = TRUE))]
  cohort <- dplyr::left_join(cohort, retained, by = "donor_id")
  cohort$n_samples <- tidyr::replace_na(cohort$n_samples, 0L)
  cohort$exclusion_reason[cohort$included & cohort$n_samples == 0L] <- "sample_depth_filter"
  fitted_ids <- rownames(fit$samples)
  fitted <- get_psbulk_sample_tibble(matrix[, fitted_ids, drop = FALSE]) |>
    dplyr::count(donor_id, name = "n_fitted_samples")
  fitted$donor_id <- cohort$donor_id[match(fitted$donor_id, gsub("_", "-", cohort$donor_id, fixed = TRUE))]
  cohort <- dplyr::left_join(cohort, fitted, by = "donor_id")
  cohort$n_fitted_samples <- tidyr::replace_na(cohort$n_fitted_samples, 0L)
  cohort$exclusion_reason[is.na(cohort$exclusion_reason) & cohort$n_fitted_samples == 0L] <- "model_sample_filter"
  cohort$included <- is.na(cohort$exclusion_reason)
  cells <- metadata |>
    dplyr::filter(paste(gsub("_", "-", as.character(PCA_harmony_SNN_cluster_cell_type), fixed = TRUE),
      gsub("_", "-", donor_id, fixed = TRUE), sep = "_") %in% fitted_ids) |>
    dplyr::group_by(donor_id) |>
    dplyr::summarise(n_nuclei = dplyr::n(), GEM_well_IDs = paste(sort(unique(GEM_well_ID)), collapse = ";"), .groups = "drop")
  dplyr::left_join(cohort, cells, by = "donor_id") |>
    dplyr::mutate(n_nuclei = tidyr::replace_na(n_nuclei, 0L))
}
