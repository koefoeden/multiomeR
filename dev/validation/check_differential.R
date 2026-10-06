# Require actual results for the mixed_human_31x cardiomyocyte sex comparison.
for (analysis in c("gene_expression", "chromatin_accessibility", "motif_family_accessibility",
                   "transcription_factor_activity")) {
  result <- targets::tar_read_raw(paste0(
    "results_tibble.", analysis, ".differential_analyses.mixed_human_31x"
  )) |>
    dplyr::bind_rows()
  stopifnot(nrow(result) > 0L, all(c("PValue", "contrast") %in% names(result)))
  comparison <- result[result$contrast == "male_vs_female", ]
  stopifnot(nrow(comparison) > 0L, any(is.finite(comparison$PValue)))
}
composition <- targets::tar_read_raw(
  "model_results.cell_type_composition.differential_analyses.mixed_human_31x"
) |>
  dplyr::bind_rows()
stopifnot(all(c("contrast", "estimate", "p.value") %in% names(composition)))
sex_effect <- composition[composition$contrast == "male_vs_female", ]
stopifnot(nrow(sex_effect) > 0L, any(is.finite(sex_effect$estimate) & is.finite(sex_effect$p.value)))
