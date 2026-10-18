.suffix_stripped_target_name <- function(name) {
  stringr::str_remove(stringr::str_remove(name, "\\..*"), "_[a-f0-9]{16}$")
}

.target_commands <- function(names) {
  capture.output(manifest <- targets::tar_manifest(fields = c("name", "command")))
  manifest$command[match(stringr::str_remove(names, "_[a-f0-9]{16}$"), manifest$name)]
}

.target_traceback <- function(name, store = targets::tar_config_get("store")) {
  tryCatch(
    targets:::workspace_read(name = name, path_store = store)$target$metrics$traceback,
    error = function(e) NULL
  )
}

load_CFG <- function(name) {
  build_GEM_well_tibble() |>
    build_dataset_config_tibble() |>
    dplyr::filter(.data$dataset == .env$name) |>
    as.list() |>
    purrr::flatten() |>
    list2env(envir = .GlobalEnv)
}

#' Inspect a saved target workspace
#'
#' Print the target's command, stored error, warnings, and traceback, and the
#' top-level structure of each dependency reconstructed from the workspace.
#'
#' @param name Full target name, including the dynamic branch suffix when relevant.
#' @param envir Environment used by `targets` internals to materialize dependency frames.
#' @param store Targets store path.
#' @param print Logical; if TRUE, print the summary.
#' @return Invisibly returns a list with target metadata, command, traceback,
#'   and an environment holding the dependencies.
#' @keywords internal
inspect_target_workspace <- function(name, envir = parent.frame(), store = targets::tar_config_get("store"), print = TRUE) {
  workspace <- targets:::workspace_read(name = name, path_store = store)
  targets:::workspace_populate(workspace)
  dependencies <- targets:::frames_get_envir(targets:::frames_produce(envir, workspace$target, workspace$subpipeline))
  meta <- targets::tar_meta(
    names = tidyselect::all_of(name),
    fields = c("type", "parent", "error", "warnings"),
    store = store
  )
  out <- list(
    name = name,
    meta = meta,
    command = .target_commands(name),
    traceback = workspace$target$metrics$traceback,
    dependencies = dependencies
  )
  if (!isTRUE(print)) {
    return(invisible(out))
  }

  cat("Target: ", name, "\nType: ", meta$type[1], "\n", sep = "")
  if (!is.na(meta$parent[1])) {
    cat("Parent: ", meta$parent[1], "\n", sep = "")
  }
  cat("Command:\n", out$command, "\n", sep = "")
  for (field in c("error", "warnings")) {
    if (!is.na(meta[[field]][1])) {
      cat("\nStored ", field, ":\n", stringr::str_trunc(meta[[field]][1], width = 1000), "\n", sep = "")
    }
  }
  if (length(out$traceback)) {
    cat("\nTraceback:\n\t", paste(stringr::str_trunc(out$traceback, width = 1000), collapse = "\n\t"), "\n", sep = "")
  }
  cat("\nDependencies:\n")
  for (dependency in sort(ls(dependencies, all.names = TRUE))) {
    cat("- ", dependency, ": ", sep = "")
    utils::str(dependencies[[dependency]], max.level = 1, list.len = 16, give.attr = FALSE)
  }
  invisible(out)
}

#' List errored targets of the current graph
#'
#' Stored errors last until a target succeeds, whichever run raised them, so
#' errors of removed targets are dropped by matching each target, or a branch's
#' pattern, against the manifest.
#'
#' @param target_name_pattern Regex used to filter target names.
#' @return Tibble of errored targets with their stored error, warnings and command.
#' @keywords internal
list_errored_targets <- function(target_name_pattern = ".") {
  capture.output(manifest <- targets::tar_manifest(fields = c("name", "command")))
  targets::tar_meta(fields = c("name", "parent", "error", "warnings")) |>
    dplyr::mutate(command = manifest$command[match(dplyr::coalesce(.data$parent, .data$name), manifest$name)]) |>
    dplyr::filter(
      !is.na(.data$error),
      !is.na(.data$command),
      stringr::str_detect(.data$name, target_name_pattern)
    ) |>
    dplyr::arrange(.data$error) |>
    dplyr::select(dplyr::all_of(c("name", "error", "warnings", "command")))
}

#' List selected targets that the latest run did not finish
#'
#' A target is finished when its metadata holds no error and the latest run's
#' progress reports it completed or skipped. Neither the exit status of
#' `tar_make()`, which returns normally when targets error under
#' `error = "trim"`, nor its log shows this: a worker that crashes beyond
#' `crashes_max` aborts the run without storing an error. Every run start
#' empties the progress file, so check right after the run.
#'
#' @param target_names Names of the selected stems and patterns.
#' @param store Targets store path.
#' @return Tibble of unfinished targets with their stored error and progress;
#'   `NA` progress means the run never reached the target.
#' @keywords internal
list_unfinished_targets <- function(target_names, store = targets::tar_config_get("store")) {
  tibble::tibble(name = target_names) |>
    dplyr::left_join(
      targets::tar_meta(names = tidyselect::any_of(target_names), fields = "error", store = store),
      by = "name"
    ) |>
    dplyr::left_join(
      targets::tar_progress(names = tidyselect::any_of(target_names), store = store),
      by = "name"
    ) |>
    dplyr::filter(!is.na(.data$error) | !.data$progress %in% c("completed", "skipped"))
}

#' Print one errored target per distinct error message
#'
#' @param target_name_pattern Regex used to filter target names.
#' @param tracebacks Logical; if TRUE, also print each target's command and traceback.
#' @return Invisibly returns the distinct errored-target tibble.
#' @keywords internal
list_distinct_errored_targets <- function(target_name_pattern = ".", tracebacks = FALSE) {
  errors <- list_errored_targets(target_name_pattern) |>
    dplyr::distinct(.data$error, .keep_all = TRUE)
  width <- if (!tracebacks) 500 else if (identical(target_name_pattern, ".")) 200 else 1000
  for (i in seq_len(nrow(errors))) {
    name <- errors$name[[i]]
    cat(
      "- Full, tar_map-resolved target name: ", name, "\n",
      "- Suffix-stripped target name: ", .suffix_stripped_target_name(name), "\n",
      sep = ""
    )
    if (tracebacks) {
      cat("- Command: ", errors$command[[i]], "\n", sep = "")
    }
    cat("- Error: ", stringr::str_trunc(errors$error[[i]], width = width), "\n", sep = "")
    if (tracebacks) {
      traceback <- stringr::str_trunc(.target_traceback(name), width = width)
      cat("- Traceback:\n\t", paste(traceback, collapse = "\n\t"), "\n", sep = "")
    }
    cat("- Warnings: ", stringr::str_trunc(errors$warnings[[i]], width = width), "\n\n", sep = "")
  }
  invisible(errors)
}

list_distinct_errored_targets_w_tracebacks <- function(target_name_pattern = ".") {
  list_distinct_errored_targets(target_name_pattern, tracebacks = TRUE)
}
