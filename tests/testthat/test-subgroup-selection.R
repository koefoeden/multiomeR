testthat::test_that("subgroups count and select exact labels from the configured column", {
  load_project_test_runtime()
  manifest <- targets::tar_manifest(
    script = file.path(multiomeR_project_root, "extra_targets/subgroups.R"),
    callr_function = NULL
  )
  commands <- setNames(lapply(manifest$command, str2lang), manifest$name)
  metadata <- tibble::tibble(
    barcode = letters[1:6],
    PCA_harmony_SNN_cluster_cell_type = rep("other", 6),
    chosen_group = factor(c("T", "T", "Treg", "Treg", "B+", "B+"))
  )
  context <- list2env(list(
    metadata_w_cell_types_subgroup_tibble.WNN = metadata,
    aggregation_subgroups_col = "chosen_group",
    aggregation_subgroups_min_nuclei_filter = 2L
  ))
  groups <- eval(commands$subgroups_to_process_vec, context)
  testthat::expect_equal(groups, c("B+", "T", "Treg"))
  for (group in groups) {
    context$subgroups_to_process_vec <- group
    selected <- eval(commands$metadata_tibble.subgroups, context)
    testthat::expect_identical(selected, metadata[metadata$chosen_group == group, ])
  }
  context$aggregation_subgroups_min_nuclei_filter <- 3L
  testthat::expect_error(
    eval(commands$subgroups_to_process_vec, context),
    "Consider lowering aggregation_subgroups_min_nuclei_filter"
  )
})
