# Explicit local validation; reuse the configured store and all valid targets.
Sys.setenv(MULTIOMER_VALIDATION = '1')
local({
  pid <- targets::tar_pid()
  if (length(pid) == 1L && !is.na(pid)) {
    live <- suppressWarnings(system2('ps', c('-p', pid, '-o', 'pid='), stdout = TRUE))
    if (length(live)) stop('A process still owns the recorded targets PID: ', pid)
  }
  wells <- build_active_GEM_well_tibble(build_GEM_well_tibble())
  inputs <- c(wells$GEM_well_cellranger_arc_count_dir, wells$GEM_well_cellbender_h5_file)
  missing <- inputs[!is.na(inputs) & !file.exists(inputs)]
  if (length(missing)) stop('Missing validation inputs:\n', paste(missing, collapse = '\n'), '\nSee validation/README.md.')
  report <- file.path(targets::tar_config_get('store'), 'validation', format(Sys.time(), '%Y%m%dT%H%M%S'))
  dir.create(report, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(MULTIOMER_VALIDATION_REPORT = normalizePath(report))
  writeLines(system2('git', c('rev-parse', 'HEAD'), stdout = TRUE), file.path(report, 'commit.txt'))
  writeLines(system2('git', c('status', '--porcelain'), stdout = TRUE), file.path(report, 'working-tree.txt'))
  writeLines(capture.output(utils::sessionInfo()), file.path(report, 'session.txt'))
  source('validation/preflight.R', local = TRUE)
  checks <- run_validation_preflight()
  write.table(checks, file.path(report, 'preflight.tsv'), sep = '\t', row.names = FALSE, quote = TRUE)
  print(checks)
  if (any(checks$result != 'PASS')) stop('Preflight failed; pipeline was not started. Report: ', report)
})
targets::tar_make()
source('validation/check.R')
