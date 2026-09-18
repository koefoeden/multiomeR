source_project_file("R/pseudobulk_helpers.R")

make_correlation_sampling_fixture <- function() {
  withr::with_seed(19L, {
    block <- rep(seq_len(8L), each = 3L)
    exposure <- rep(c(-1, 0, 1), 8L)
    design <- stats::model.matrix(~ factor(exposure))
    family_effects <- matrix(stats::rnorm(200L * 8L, sd = 0.5), 200L, 8L)
    means <- 40 * exp(family_effects[, block])
    counts <- matrix(stats::rpois(length(means), means), nrow(means)) + 1L
    dimnames(counts) <- list(paste0("feature_", seq_len(200L)), paste0("sample_", seq_len(24L)))
    rownames(design) <- colnames(counts)
    list(counts = counts, design = design, block = block)
  })
}

testthat::test_that("unlimited correlation fitting preserves edgeR including sparse-count DF", {
  input <- make_correlation_sampling_fixture()
  input$counts[1:5, seq.int(1L, 24L, by = 3L)] <- 0L
  reference <- do.call(edgeR::voomLmFit, input)
  unlimited <- do.call(fit_psbulk_voom, input)
  all_features <- do.call(fit_psbulk_voom, c(input, list(correlation_max_features = 200L)))
  sampled <- do.call(fit_psbulk_voom, c(input, list(correlation_max_features = 60L)))

  testthat::expect_equal(unlimited, reference, tolerance = 1e-12)
  testthat::expect_equal(all_features, reference, tolerance = 1e-12)
  testthat::expect_identical(sampled$df.residual, reference$df.residual)
  testthat::expect_gt(length(unique(reference$df.residual)), 1L)
})

testthat::test_that("sampling changes only correlation estimation and restores RNG and namespaces", {
  input <- make_correlation_sampling_fixture()
  original_function <- edgeR::voomLmFit
  original_estimator <- get("duplicateCorrelation", envir = environment(original_function))
  withr::local_seed(53L)
  seed_before <- .Random.seed
  fit <- do.call(fit_psbulk_voom, c(input, list(correlation_max_features = 60L)))
  testthat::expect_identical(.Random.seed, seed_before)
  rows <- withr::with_seed(732L, sort(sample.int(nrow(input$counts), 60L)))
  correlation <- limma::duplicateCorrelation(
    fit$EList$E[rows, , drop = FALSE], design = input$design, block = input$block,
    weights = fit$EList$weights[rows, , drop = FALSE]
  )$consensus.correlation

  testthat::expect_equal(fit$correlation, correlation, tolerance = 1e-12)
  testthat::expect_identical(rownames(fit$coefficients), rownames(input$counts))
  testthat::expect_identical(dim(fit$EList$weights), dim(input$counts))
  testthat::expect_identical(edgeR::voomLmFit, original_function)
  testthat::expect_identical(
    get("duplicateCorrelation", envir = environment(edgeR::voomLmFit)), original_estimator
  )
})

source_project_file("R/differential_analysis_helpers.R")

make_abundance_fixture <- function() {
  donors <- data.frame(donor_id = paste0("d", 1:12), sexMale = rep(0:1, each = 6))
  n <- c(0, 9, 17, 25, 12, 20, 30, 47, 32, 55, 42, 28)
  metadata <- do.call(rbind, lapply(seq_len(12), function(i) {
    data.frame(donor_id = donors$donor_id[i], GEM_well_ID = paste0("w", i),
      cell_type = rep(c("Cardiomyocyte", "Other"), c(n[i], 100 - n[i])))
  }))
  model <- list(formula = "cbind(n_nuclei, n_other_nuclei) ~ sexMale",
    contrast_specs_vec = c(male_vs_female = "sexMale", female_vs_male = "-sexMale"))
  list(donors = donors, metadata = metadata, model = model)
}

testthat::test_that("abundance subsets preserve denominators, zero donors and explicit exclusions", {
  input <- make_abundance_fixture()
  input$donors$sexMale[12] <- NA
  input$model$GEM_well_IDs <- paste0("w", c(1:8, 12))
  input$model$cell_types_to_test <- "Cardiomyocyte"
  data <- prepare_DCTC_model_data(input$metadata, input$donors, input$model, "cell_type")
  testthat::expect_equal(nrow(data$counts), 8L)
  testthat::expect_true(all(data$counts$n_total_nuclei == 100L))
  testthat::expect_equal(data$counts$n_nuclei[data$counts$donor_id == "d1"], 0L)
  testthat::expect_equal(data$counts$n_other_nuclei[data$counts$donor_id == "d1"], 100L)
  testthat::expect_equal(data$cohort$exclusion_reason[data$cohort$donor_id == "d12"], "missing_model_metadata")
  testthat::expect_equal(data$cohort$exclusion_reason[data$cohort$donor_id == "d9"], "no_selected_samples")
  input$model$cell_types_to_test <- NULL
  all_types <- prepare_DCTC_model_data(input$metadata, input$donors, input$model, "cell_type")
  testthat::expect_equal(nrow(all_types$counts), 16L)
  sums <- tapply(all_types$counts$prop, all_types$counts$donor_id, sum)
  testthat::expect_equal(as.numeric(sums), rep(1, 8))
})

testthat::test_that("abundance named contrasts agree with direct beta-binomial coefficients", {
  input <- make_abundance_fixture()
  data <- prepare_DCTC_model_data(input$metadata, input$donors, input$model, "cell_type")
  results <- fit_DCTC_model(data, input$model, "sex")
  reference <- glmmTMB::glmmTMB(cbind(n_nuclei, n_other_nuclei) ~ sexMale,
    data = subset(data$counts, cluster == "Cardiomyocyte"), family = glmmTMB::betabinomial())
  forward <- subset(results, cluster == "Cardiomyocyte" & contrast == "male_vs_female")
  reverse <- subset(results, cluster == "Cardiomyocyte" & contrast == "female_vs_male")
  testthat::expect_equal(forward$estimate, unname(glmmTMB::fixef(reference)$cond['sexMale']), tolerance = 1e-8)
  testthat::expect_equal(forward$std.error, sqrt(unname(stats::vcov(reference)$cond['sexMale','sexMale'])), tolerance = 1e-8)
  testthat::expect_equal(reverse$estimate, -forward$estimate)
  testthat::expect_equal(reverse$p.value, forward$p.value)
  testthat::expect_true(all(results$status == "fitted"))
  testthat::expect_equal(results$FDR[results$contrast == "male_vs_female"],
    stats::p.adjust(results$p.value[results$contrast == "male_vs_female"], "BH"))
})

testthat::test_that("shared cohort validation rejects missing variables and non-estimable designs", {
  input <- make_abundance_fixture()
  testthat::expect_error(get_differential_model_cohort(input$donors['donor_id'], input$model,
    input$donors$donor_id), "Missing model metadata")
  input$model$donor_ids <- input$donors$donor_id[1:6]
  data <- prepare_DCTC_model_data(input$metadata, input$donors, input$model, "cell_type")
  testthat::expect_error(fit_DCTC_model(data, input$model, "sex"), "rank deficient")
  testthat::expect_error(normalize_differential_models(stats::setNames(list(input$model, input$model), c("sex", "sex"))), "uniquely named")
})

testthat::test_that("feature cohort exports distinguish sample-depth and fitted-sample exclusions", {
  input <- make_abundance_fixture()
  names(input$metadata)[names(input$metadata) == "cell_type"] <- "PCA_harmony_SNN_cluster_cell_type"
  matrix <- matrix(1, 1, 12, dimnames = list("gene", paste0("Cardiomyocyte_", input$donors$donor_id)))
  model <- list(formula = "~ sexMale", cell_type_subset = "Cardiomyocyte")
  fit <- list(samples = data.frame(row.names = colnames(matrix)[1:6]))
  cohort <- get_feature_model_cohort(matrix, matrix[, 1:10, drop = FALSE], fit,
    input$donors, model, input$metadata)
  testthat::expect_true(all(cohort$included[1:6]))
  testthat::expect_true(all(cohort$exclusion_reason[7:10] == "model_sample_filter"))
  testthat::expect_true(all(cohort$exclusion_reason[11:12] == "sample_depth_filter"))
  testthat::expect_true(all(cohort$n_nuclei[7:12] == 0))
})
