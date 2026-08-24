skip_w_dummy_file_if <- function(condition_expr, targets_list, skip_expr = character(0)) {
  # capture (quote) both pieces without evaluating now
  condition_lang <- substitute(condition_expr)
  skip_lang <- substitute(skip_expr)

  # build the hook: if (<condition>) <skip_expr> else .x
  hook_lang <- bquote(if (.(condition_lang)) .(skip_lang) else .x)

  # hand the *quoted* hook to the raw variant
  tarchetypes::tar_hook_outer_raw(targets_list, hook = hook_lang)
}

# Called by tar_load()
w_def <- function(packages) {
  targets::tar_option_get("packages") %>%
    c(packages) %>%
    unique()
}
