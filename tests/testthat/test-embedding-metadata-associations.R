testthat::test_that("association summaries align cells and preserve variable roles", {
  load_project_test_runtime()
  embedding <- cbind(PCA_1 = 1:6, PCA_2 = rep(c(0, 1), each = 3))
  rownames(embedding) <- letters[1:6]
  metadata <- tibble::tibble(
    barcode_w_prefix = rev(letters[1:6]),
    linear = 6:1,
    group = rep(c("B", "A"), each = 3),
    constant = 1
  )
  args <- list(
    embedding_matrix = embedding, harmony_embedding_matrix = -embedding,
    metadata_tibble = metadata, dims = 1:3, dim_prefix = "PCA_",
    continuous_technical_cols = c("linear", "missing", "group"),
    categorical_technical_cols = "group",
    continuous_biological_cols = "constant"
  )
  summaries <- do.call(get_embedding_metadata_association_tibbles, args)
  testthat::expect_named(summaries, c("continuous_technical", "categorical_technical", "continuous_biological", "categorical_biological"))
  continuous <- summaries$continuous_technical
  testthat::expect_equal(unique(continuous$variable), "linear")
  testthat::expect_equal(unique(continuous$dim), 1:2)
  testthat::expect_equal(continuous$metric[continuous$dim == 1], c(1, 1))
  categorical <- summaries$categorical_technical
  testthat::expect_equal(categorical$metric[categorical$dim == 2], c(1, 1))
  testthat::expect_true(all(is.na(summaries$continuous_biological$metric)))
  testthat::expect_equal(nrow(summaries$categorical_biological), 0L)
  plots <- do.call(plot_embedding_metadata_association_barplots, args)
  testthat::expect_identical(plots$continuous_technical$data$metric, continuous$metric)
  testthat::expect_equal(plots$continuous_technical$scales$get_scales("x")$breaks, 1:2)
})
