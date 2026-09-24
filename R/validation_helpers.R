# Explicit selection prevents unrelated inactive examples from entering validation.
validation_aggregations <- function() {
  mode <- Sys.getenv('MULTIOMER_VALIDATION', unset = '0')
  if (!mode %in% c('0', '1')) stop('MULTIOMER_VALIDATION must be 0 or 1.')
  if (mode == '1') c('immune_human_2x', 'brain_mouse', 'mixed_human_31x') else character()
}
