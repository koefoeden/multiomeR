source("R/bootstrap_helpers.R")
load_project_runtime(force = TRUE)

canonical <- tibble::tibble(
  GEM_well_ID = c("well_B", "well_A", "well_unused"),
  GEM_well_dataset = c("study", "study", "other"),
  GEM_well_cellranger_arc_count_dir = c("B", "A", "unused"),
  GEM_well_QC_exclude_list = c("nCount_RNA < 10", "", ""),
  GEM_well_multiplex_batch = c("batch_2", "batch_1", "batch_unused"),
  GEM_well_run_harmony = c(TRUE, TRUE, FALSE),
  GEM_well_notes = c("note B", "note A", "unused note")
)

aggregation_metadata <- subset_keyed_metadata_tibble(
  canonical,
  key_col = "GEM_well_ID",
  keys = c("well_A", "well_B"),
  source_label = "synthetic canonical metadata"
)
stopifnot(identical(aggregation_metadata$GEM_well_ID, c("well_A", "well_B")))

harmony <- project_keyed_metadata_tibble(
  aggregation_metadata,
  key_col = "GEM_well_ID",
  requested_columns = "GEM_well_run_harmony"
)
analysis <- project_keyed_metadata_tibble(
  aggregation_metadata,
  key_col = "GEM_well_ID",
  requested_columns = "GEM_well_multiplex_batch"
)
annotation <- get_GEM_well_annotation_metadata_tibble(aggregation_metadata)

changed_unused_row <- canonical
changed_unused_row$GEM_well_notes[changed_unused_row$GEM_well_ID == "well_unused"] <- "changed"
changed_unused_subset <- subset_keyed_metadata_tibble(
  changed_unused_row,
  "GEM_well_ID",
  c("well_A", "well_B"),
  "synthetic canonical metadata"
)
stopifnot(identical(aggregation_metadata, changed_unused_subset))

changed_note <- aggregation_metadata
changed_note$GEM_well_notes[changed_note$GEM_well_ID == "well_A"] <- "changed"
stopifnot(
  identical(
    harmony,
    project_keyed_metadata_tibble(
      changed_note,
      "GEM_well_ID",
      "GEM_well_run_harmony"
    )
  ),
  identical(
    analysis,
    project_keyed_metadata_tibble(
      changed_note,
      "GEM_well_ID",
      "GEM_well_multiplex_batch"
    )
  ),
  !identical(annotation, get_GEM_well_annotation_metadata_tibble(changed_note))
)

changed_batch <- aggregation_metadata
changed_batch$GEM_well_multiplex_batch[changed_batch$GEM_well_ID == "well_A"] <- "changed"
stopifnot(
  identical(
    harmony,
    project_keyed_metadata_tibble(
      changed_batch,
      "GEM_well_ID",
      "GEM_well_run_harmony"
    )
  ),
  !identical(
    analysis,
    project_keyed_metadata_tibble(
      changed_batch,
      "GEM_well_ID",
      "GEM_well_multiplex_batch"
    )
  )
)

message("GEM-well metadata view validation passed.")
