# Both profiles exercise their complete configured graph.
targets::tar_make()
local({
  metadata <- targets::tar_meta(fields = c(name, type, children, error), complete_only = FALSE)
  current_names <- targets::tar_manifest(fields = name)$name
  patterns <- metadata[metadata$name %in% current_names & metadata$type == "pattern", ]
  current_names <- union(current_names, unlist(patterns$children, use.names = FALSE))
  stopifnot(all(current_names %in% metadata$name))
  current_metadata <- metadata[metadata$name %in% current_names, ]
  errors <- current_metadata[!is.na(current_metadata$error) & nzchar(current_metadata$error), ]
  if (nrow(errors)) stop(paste(errors$name, errors$error, collapse = "\n"))
  config <- read_aggregation_config_tibble()
  config <- config[vapply(config$is_active, isTRUE, logical(1)), ]
  stopifnot(nrow(config) > 0L)
  summary <- lapply(seq_len(nrow(config)), function(i) {
    aggregation <- config$aggregation[[i]]
    object <- targets::tar_read_raw(paste0("multimodal_Seurat_object.8_multimodal_QC.", aggregation))
    stopifnot(inherits(object, "Seurat"), ncol(object) > 0L)
    stopifnot(all(c("RNA", "ATAC") %in% names(object@assays)))
    stopifnot(setequal(unique(object$GEM_well_ID), config$aggregation_GEM_well_IDs[[i]]))
    data.frame(aggregation = aggregation, cells = ncol(object), GEM_wells = length(unique(object$GEM_well_ID)))
  })
  run <- Sys.getenv("MULTIOMER_VALIDATION_RUN")
  private_checks <- file.path(run, "profile", "validate.R")
  if (file.exists(private_checks)) source(private_checks, local = TRUE)
  stopifnot(length(targets::tar_outdated()) == 0L)
  write.table(do.call(rbind, summary), file.path(run, "summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines("Complete configured graph and output checks passed", file.path(run, "validated"))
  print(do.call(rbind, summary))
})
