#' Get spline matrices w sex interaction
#'
#' Build model and spline-basis matrices with a sex-by-spline interaction.
#'
#' @param phenodata_tibble Sample metadata used for design-matrix construction.
#' @param var_col Numeric metadata column modeled with a natural spline basis.
#' @param binary_sex_col Numeric or binary-coded sex column multiplied by each
#'   spline basis column to form interaction terms.
#' @param additional_covar_cols Optional additional covariate columns included as
#'   main effects in the design matrix.
#' @param n_df Degrees of freedom passed to `splines::ns()`.
#' @return A list with the full `design_matrix` and the spline `basis_matrix`.
#' @keywords internal

get_spline_matrices_w_sex_interaction <- function(phenodata_tibble, var_col, binary_sex_col, additional_covar_cols = NULL, n_df) {
  additional_covar_cols_str <- if (is.null(additional_covar_cols)) "" else stringr::str_c("+ ", additional_covar_cols) %>% stringr::str_c(collapse = " ")
  # make basis matrix
  spline_colnames <- stringr::str_c(var_col, "_spline", 1:n_df)
  basis_matrix <- phenodata_tibble[[var_col]] %>% splines::ns(df = n_df) %>% magrittr::set_colnames(spline_colnames)

  # make sex interaction matrix
  sex_interaction_colnames <- stringr::str_c(binary_sex_col, ":", spline_colnames)
  sex_interaction_matrix <- sweep(basis_matrix, 1, phenodata_tibble[[binary_sex_col]], `*`) %>% magrittr::set_colnames(sex_interaction_colnames)

  covar_model_matrix <- stringr::str_glue("~ 1 + {binary_sex_col} {additional_covar_cols_str}") %>% stats::as.formula() %>% stats::model.matrix(data = phenodata_tibble)

  design_matrix <- cbind(covar_model_matrix, basis_matrix, sex_interaction_matrix)

  return(list(design_matrix = design_matrix, basis_matrix = basis_matrix))
}

#' Get sex weighted spline contrast vec
#'
#' Construct a contrast vector for a spline change with weighted sex interaction.
#'
#' @param psbulk_DX_DGE_list_fit Fitted pseudobulk model object containing
#'   `basis_matrix`, `design`, and sample metadata.
#' @param z_from_to_num Two numeric values defining the low/high values at which
#'   to evaluate the spline basis when `q_from_to_num` is `NULL`.
#' @param sex_binary_col_weight_num Weight applied to sex-interaction spline
#'   coefficients.
#' @param q_from_to_num Optional two quantiles of the modeled variable; when set,
#'   these replace `z_from_to_num`.
#' @param main_terms_weight Weight applied to main-effect spline coefficients.
#' @param sex_binary_col Prefix of the binary sex term in the design-matrix
#'   interaction column names.
#' @return Named numeric contrast vector aligned to columns of
#'   `psbulk_DX_DGE_list_fit$design`.
#' @keywords internal

get_sex_weighted_spline_contrast_vec <- function(
  psbulk_DX_DGE_list_fit,
  z_from_to_num = c(-0.5, 0.5),
  sex_binary_col_weight_num = 0.5,
  q_from_to_num = NULL,
  main_terms_weight = 1,
  sex_binary_col = "sexMale"
) {
  # overwrite default z if quantiles are specified
  var_col <- psbulk_DX_DGE_list_fit$basis_matrix %>% colnames() %>% stringr::str_remove("_spline.*") %>% unique()

  if (!is.null(q_from_to_num)) {
    z_from_to_num <- psbulk_DX_DGE_list_fit$samples[[var_col]] %>% stats::quantile(q_from_to_num, na.rm = TRUE)
  }

  Z_preds_low_high_matrix <- psbulk_DX_DGE_list_fit$basis_matrix %>% stats::predict(newx = z_from_to_num)
  Z_diff_per_spline <- (Z_preds_low_high_matrix[2, ] - Z_preds_low_high_matrix[1, ]) %>% as.numeric()

  main_terms_vec <- psbulk_DX_DGE_list_fit$basis_matrix %>% colnames()
  sex_interaction_terms_vec <- psbulk_DX_DGE_list_fit$design %>% colnames() %>% stringr::str_subset(stringr::str_glue("^{sex_binary_col}:{var_col}_spline"))

  # check for equal lengths
  stopifnot(
    length(main_terms_vec) == length(Z_diff_per_spline),
    length(sex_interaction_terms_vec) == length(Z_diff_per_spline)
  )

  # build contrast vec, incorporating specified sex-specific weights and range
  contrast_vec <- rep(0, ncol(psbulk_DX_DGE_list_fit$design)) %>% purrr::set_names(colnames(psbulk_DX_DGE_list_fit$design))
  contrast_vec[main_terms_vec] <- main_terms_weight * Z_diff_per_spline
  contrast_vec[sex_interaction_terms_vec] <- sex_binary_col_weight_num * Z_diff_per_spline

  contrast_vec
}
