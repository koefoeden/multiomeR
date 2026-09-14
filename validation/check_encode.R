# Require actual results for the public cardiomyocyte sex comparison.
for (modality in c("DGE", "DCA", "DTFA", "DCTA")) {
  result <- targets::tar_read_raw(paste0(
    "results_tibble.", modality, ".differential_analyses.ENCODE_heart_LV_6x"
  )) |>
    dplyr::bind_rows()
  stopifnot(nrow(result) > 0L, all(c("PValue", "contrast") %in% names(result)))
  comparison <- result[result$contrast == "male_vs_female", ]
  stopifnot(nrow(comparison) > 0L, any(is.finite(comparison$PValue)))
}
composition <- targets::tar_read_raw(
  "model_results.DCTC.differential_analyses.ENCODE_heart_LV_6x"
)
stopifnot(all(c("term", "estimate", "p.value") %in% names(composition)))
sex_effect <- composition[composition$term == "sexMale", ]
stopifnot(nrow(sex_effect) > 0L, any(is.finite(sex_effect$estimate) & is.finite(sex_effect$p.value)))
