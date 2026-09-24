# Verify the full validation graph after tar_make(), without invalidating targets.
local({
  stopifnot(length(validation_aggregations()) == 3L)
  commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
  metadata <- targets::tar_meta(fields = c(name, type, children, error), complete_only = FALSE)
  current_names <- targets::tar_manifest(fields = name)$name
  patterns <- metadata[metadata$name %in% current_names & metadata$type == 'pattern', ]
  current_names <- union(current_names, unlist(patterns$children, use.names = FALSE))
  stopifnot(all(current_names %in% metadata$name))
  current <- metadata[metadata$name %in% current_names, ]
  errors <- current[!is.na(current$error) & nzchar(current$error), ]
  if (nrow(errors)) stop(paste(errors$name, errors$error, collapse = '\n'))
  config <- read_aggregation_config_tibble()
  config <- config[vapply(config$is_active, isTRUE, logical(1)), ]
  stopifnot(setequal(config$aggregation, validation_aggregations()))
  summary <- lapply(seq_len(nrow(config)), function(i) {
    aggregation <- config$aggregation[[i]]
    object <- targets::tar_read_raw(paste0('multimodal_Seurat_object.8_multimodal_QC.', aggregation))
    stopifnot(inherits(object, 'Seurat'), ncol(object) > 0L)
    stopifnot(all(c('RNA', 'ATAC') %in% names(object@assays)))
    stopifnot(setequal(unique(object$GEM_well_ID), config$aggregation_GEM_well_IDs[[i]]))
    data.frame(commit = commit, aggregation = aggregation,
               cells = ncol(object), GEM_wells = length(unique(object$GEM_well_ID)))
  })
  source('validation/check_differential.R', local = TRUE)
  stopifnot(length(targets::tar_outdated()) == 0L)
  report <- Sys.getenv('MULTIOMER_VALIDATION_REPORT')
  stopifnot(nzchar(report))
  dir.create(report, recursive = TRUE, showWarnings = FALSE)
  write.table(do.call(rbind, summary), file.path(report, 'summary.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
  print(do.call(rbind, summary))
  writeLines('PASS', file.path(report, 'PASS'))
  cat('PASS: current graph and output checks; incremental target reuse allowed.\nReport:', report, '\n')
})
