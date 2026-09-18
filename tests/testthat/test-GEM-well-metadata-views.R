source_project_file("packages/multiomeRCore/R/null_default.R")
source_project_file("R/general_helpers.R")
source_project_file("R/target_mapping_helpers.R")
source_project_file("R/configuration_helpers.R")

testthat::test_that("configuration selection is local, explicit and never falls back", {
  root <- normalizePath(withr::local_tempdir())
  dir.create(file.path(root, "configuration"))
  dir.create(file.path(root, "configuration_my_project"))
  filename <- "cfg_GEM_wells.tsv"
  default_file <- file.path(root, "configuration", filename)
  selected_file <- file.path(root, "configuration_my_project", filename)
  writeLines("default", default_file)
  writeLines("selected", selected_file)
  selector <- file.path(root, "configuration.local")
  testthat::expect_identical(configuration_path(filename, root), normalizePath(default_file))
  writeLines("configuration_my_project", selector)
  withr::local_dir(tempdir())
  testthat::expect_identical(configuration_path(filename, root), normalizePath(selected_file))
  writeLines(file.path(root, "configuration_my_project"), selector)
  testthat::expect_identical(configuration_path(filename, root), normalizePath(selected_file))
  unlink(selected_file)
  testthat::expect_error(configuration_path(filename, root), "Missing file in selected")
  testthat::expect_identical(configuration_path(filename, root, must_exist = FALSE), selected_file)
  writeLines("missing", selector)
  testthat::expect_error(configuration_path(filename, root), "directory does not exist")
  writeLines("", selector)
  testthat::expect_error(configuration_path(filename, root), "exactly one non-empty")
  writeLines(c("configuration", "configuration_my_project"), selector)
  testthat::expect_error(configuration_path(filename, root), "exactly one non-empty")
})

testthat::test_that("disabled modules do not read configuration files", {
  aggregations <- tibble::tibble(aggregation = "example")
  testthat::expect_identical(
    read_module_config_tibble("missing.yaml", "differential_analyses", aggregations[0, ], aggregations),
    aggregations[0, ]
  )
  testthat::expect_error(
    read_module_config_tibble("missing.yaml", "differential_analyses", aggregations, aggregations,
      manifest_file = file.path(multiomeR_project_root, "cfg_pipeline_parameters.tsv")),
    "missing.yaml"
  )
})

make_GEM_well_metadata_fixture <- function() {
  tibble::tibble(
    GEM_well_ID = c("well_B", "well_A", "well_unused"),
    GEM_well_dataset = c("study", "study", "other"),
    GEM_well_cellranger_arc_count_dir = c("B", "A", "unused"),
    GEM_well_QC_exclude_list = c("nCount_RNA < 10", "", ""),
    GEM_well_multiplex_batch = c("batch_2", "batch_1", "batch_unused"),
    GEM_well_run_harmony = c(TRUE, TRUE, FALSE),
    GEM_well_notes = c("note B", "note A", "unused note")
  )
}

make_aggregation_metadata_fixture <- function() {
  subset_keyed_metadata_tibble(
    make_GEM_well_metadata_fixture(),
    key_col = "GEM_well_ID",
    keys = c("well_A", "well_B"),
    source_label = "synthetic canonical metadata"
  )
}

testthat::test_that("keyed metadata subsets follow the requested key order", {
  aggregation_metadata <- make_aggregation_metadata_fixture()
  testthat::expect_identical(
    aggregation_metadata$GEM_well_ID,
    c("well_A", "well_B")
  )
})

testthat::test_that("legacy projections preserve input row and requested column order", {
  aggregation_metadata <- make_aggregation_metadata_fixture()
  reordered_metadata <- aggregation_metadata[
    rev(seq_len(nrow(aggregation_metadata))),
    rev(names(aggregation_metadata))
  ]
  legacy_projection <- project_keyed_metadata_tibble(
    reordered_metadata,
    key_col = "GEM_well_ID",
    requested_columns = c(
      "GEM_well_run_harmony",
      "GEM_well_multiplex_batch"
    )
  )

  testthat::expect_identical(
    legacy_projection$GEM_well_ID,
    c("well_B", "well_A")
  )
  testthat::expect_identical(
    names(legacy_projection),
    c(
      "GEM_well_ID",
      "GEM_well_run_harmony",
      "GEM_well_multiplex_batch"
    )
  )
})

testthat::test_that("canonical projections ignore input and request ordering", {
  aggregation_metadata <- make_aggregation_metadata_fixture()
  reordered_metadata <- aggregation_metadata[
    rev(seq_len(nrow(aggregation_metadata))),
    rev(names(aggregation_metadata))
  ]
  canonical_projection <- project_keyed_metadata_tibble(
    reordered_metadata,
    key_col = "GEM_well_ID",
    requested_columns = c(
      "GEM_well_run_harmony",
      "GEM_well_multiplex_batch"
    ),
    strict = TRUE,
    canonical = TRUE
  )
  canonical_projection_from_original <- project_keyed_metadata_tibble(
    aggregation_metadata,
    key_col = "GEM_well_ID",
    requested_columns = c(
      "GEM_well_multiplex_batch",
      "GEM_well_run_harmony"
    ),
    strict = TRUE,
    canonical = TRUE
  )

  testthat::expect_identical(
    canonical_projection,
    canonical_projection_from_original
  )
  testthat::expect_identical(
    canonical_projection$GEM_well_ID,
    c("well_A", "well_B")
  )
  testthat::expect_identical(
    names(canonical_projection),
    c(
      "GEM_well_ID",
      "GEM_well_multiplex_batch",
      "GEM_well_run_harmony"
    )
  )
})

testthat::test_that("changes to unused rows do not affect the aggregation subset", {
  canonical <- make_GEM_well_metadata_fixture()
  aggregation_metadata <- make_aggregation_metadata_fixture()
  canonical$GEM_well_notes[canonical$GEM_well_ID == "well_unused"] <- "changed"
  changed_unused_subset <- subset_keyed_metadata_tibble(
    canonical,
    "GEM_well_ID",
    c("well_A", "well_B"),
    "synthetic canonical metadata"
  )

  testthat::expect_identical(aggregation_metadata, changed_unused_subset)
})

testthat::test_that("metadata views respond only to their projected columns", {
  aggregation_metadata <- make_aggregation_metadata_fixture()
  harmony <- project_keyed_metadata_tibble(
    aggregation_metadata,
    "GEM_well_ID",
    "GEM_well_run_harmony"
  )
  analysis <- project_keyed_metadata_tibble(
    aggregation_metadata,
    "GEM_well_ID",
    "GEM_well_multiplex_batch"
  )
  annotation <- get_GEM_well_annotation_metadata_tibble(aggregation_metadata)
  changed_note <- aggregation_metadata
  changed_note$GEM_well_notes[changed_note$GEM_well_ID == "well_A"] <- "changed"

  testthat::expect_identical(
    harmony,
    project_keyed_metadata_tibble(
      changed_note,
      "GEM_well_ID",
      "GEM_well_run_harmony"
    )
  )
  testthat::expect_identical(
    analysis,
    project_keyed_metadata_tibble(
      changed_note,
      "GEM_well_ID",
      "GEM_well_multiplex_batch"
    )
  )
  testthat::expect_false(
    identical(annotation, get_GEM_well_annotation_metadata_tibble(changed_note))
  )
})

testthat::test_that("batch changes invalidate analysis but not Harmony metadata", {
  aggregation_metadata <- make_aggregation_metadata_fixture()
  harmony <- project_keyed_metadata_tibble(
    aggregation_metadata,
    "GEM_well_ID",
    "GEM_well_run_harmony"
  )
  analysis <- project_keyed_metadata_tibble(
    aggregation_metadata,
    "GEM_well_ID",
    "GEM_well_multiplex_batch"
  )
  changed_batch <- aggregation_metadata
  changed_batch$GEM_well_multiplex_batch[
    changed_batch$GEM_well_ID == "well_A"
  ] <- "changed"

  testthat::expect_identical(
    harmony,
    project_keyed_metadata_tibble(
      changed_batch,
      "GEM_well_ID",
      "GEM_well_run_harmony"
    )
  )
  testthat::expect_false(
    identical(
      analysis,
      project_keyed_metadata_tibble(
        changed_batch,
        "GEM_well_ID",
        "GEM_well_multiplex_batch"
      )
    )
  )
})
