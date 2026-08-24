.suffix_stripped_target_name <- function(name) {
  name |>
    stringr::str_remove("\\..*") |>
    stringr::str_remove("_[a-f0-9]{16}$")
}

load_CFG <- function(name) {
  build_GEM_well_tibble() |>
    build_dataset_config_tibble() |>
    dplyr::filter(.data$dataset == .env$name) |>
    as.list() |>
    purrr::flatten() |>
    list2env(envir = .GlobalEnv)
}

.target_traceback <- function(name, store = targets::tar_config_get("store")) {
  if (!is.character(name) || length(name) != 1L) {
    stop("`name` must be a length-one character vector.", call. = FALSE)
  }

  workspace <- tryCatch(
    targets:::workspace_read(name = name, path_store = store),
    error = function(e) NULL
  )
  if (is.null(workspace)) {
    return(character(0))
  }

  out <- workspace$target$metrics$traceback
  if (is.null(out)) {
    return(character(0))
  }
  out
}


.target_manifest_command <- function(name) {
  name_wo_dynamic_hash <- stringr::str_remove(name, "_[a-f0-9]{16}$")
  capture.output(command_tibble <- targets::tar_manifest())
  command_tibble |>
    dplyr::filter(.data$name == name_wo_dynamic_hash) |>
    dplyr::pull(.data$command) |>
    dplyr::first(default = NA_character_)
}


.target_workspace_dependency_env <- function(workspace, envir) {
  targets:::workspace_populate(workspace)
  frames <- targets:::frames_produce(envir, workspace$target, workspace$subpipeline)
  targets:::frames_get_envir(frames)
}


.collapse_head <- function(x, n = 8) {
  if (is.null(x)) {
    return(NA_character_)
  }
  x <- utils::head(as.character(x), n)
  if (length(x) == 0) {
    return("<none>")
  }
  paste(x, collapse = ", ")
}


.object_dim_chr <- function(object) {
  out <- dim(object)
  if (is.null(out)) {
    return(NA_character_)
  }
  paste(out, collapse = " x ")
}


.metadata_value <- function(meta, column) {
  if (!column %in% colnames(meta) || nrow(meta) == 0) {
    return(NA_character_)
  }
  out <- meta[[column]][[1]]
  if (is.null(out)) {
    return(NA_character_)
  }
  out
}


.target_workspace_object_details <- function(object) {
  details <- character()
  if (is.data.frame(object)) {
    list_cols <- names(object)[vapply(object, is.list, logical(1))]
    if (length(list_cols) > 0) {
      details <- c(details, paste0("list columns: ", paste(list_cols, collapse = ", ")))
    }
  }
  if (is.list(object) && !is.data.frame(object) && !is.null(object$design)) {
    details <- c(details, paste0("design: ", .object_dim_chr(object$design)))
    details <- c(details, paste0("design columns: ", .collapse_head(colnames(object$design), n = 16)))
  }
  if (is.list(object) && !is.data.frame(object) && is.data.frame(object$samples)) {
    details <- c(details, paste0("samples: ", .object_dim_chr(object$samples)))
    details <- c(details, paste0("sample columns: ", .collapse_head(colnames(object$samples), n = 16)))
  }
  if (length(details) == 0) {
    return(NA_character_)
  }
  paste(details, collapse = "; ")
}


.target_workspace_object_summary <- function(name, object) {
  dim_chr <- .object_dim_chr(object)
  list(
    name = name,
    class = paste(class(object), collapse = ", "),
    type = typeof(object),
    length = length(object),
    dim = dim_chr,
    names = .collapse_head(names(object)),
    rownames = .collapse_head(rownames(object)),
    colnames = .collapse_head(colnames(object)),
    columns = if (is.data.frame(object)) .collapse_head(colnames(object), n = 16) else NA_character_,
    details = .target_workspace_object_details(object),
    note = .target_workspace_object_note(name, object, dim_chr)
  )
}


.target_workspace_object_note <- function(name, object, dim_chr) {
  notes <- character()
  object_dim <- dim(object)

  if (!is.null(object_dim) && any(object_dim == 0)) {
    notes <- c(notes, paste0("zero dimension: ", dim_chr))
  }
  if (is.data.frame(object) && nrow(object) == 0) {
    notes <- c(notes, "zero-row data frame")
  }
  if (is.list(object) && length(object) == 0) {
    notes <- c(notes, "empty list")
  }
  if (is.atomic(object) && length(object) == 0) {
    notes <- c(notes, "empty vector")
  }
  if (is.atomic(object) && anyNA(object)) {
    notes <- c(notes, "contains NA")
  }

  if (length(notes) == 0) {
    return(NA_character_)
  }
  paste(unique(notes), collapse = "; ")
}


.print_target_workspace_inspection <- function(out) {
  cat("Target: ", out$name, "\n", sep = "")
  cat("Type: ", .metadata_value(out$meta, "type"), "\n", sep = "")
  parent <- .metadata_value(out$meta, "parent")
  if (!is.na(parent)) {
    cat("Parent: ", parent, "\n", sep = "")
  }
  cat("Command:\n", out$command, "\n", sep = "")

  error <- .metadata_value(out$meta, "error")
  if (!is.na(error)) {
    cat("\nStored error:\n", stringr::str_trunc(error, width = 1000), "\n", sep = "")
  }
  warnings <- .metadata_value(out$meta, "warnings")
  if (!is.na(warnings)) {
    cat("\nStored warnings:\n", stringr::str_trunc(warnings, width = 1000), "\n", sep = "")
  }
  if (length(out$traceback) > 0) {
    traceback <- stringr::str_trunc(out$traceback, width = 1000)
    cat("\nTraceback:\n\t", paste(traceback, collapse = "\n\t"), "\n", sep = "")
  }

  cat("\nDependencies:\n")
  purrr::pwalk(
    out$dependencies,
    function(name, class, type, length, dim, names, rownames, colnames, columns, details, note, ...) {
      cat("- ", name, "\n", sep = "")
      cat("  class: ", class, "\n", sep = "")
      cat("  type/length: ", type, " / ", length, "\n", sep = "")
      if (!is.na(dim)) {
        cat("  dim: ", dim, "\n", sep = "")
      }
      if (!is.na(columns)) {
        cat("  columns: ", columns, "\n", sep = "")
      } else if (!is.na(colnames)) {
        cat("  colnames: ", colnames, "\n", sep = "")
      }
      if (!is.na(rownames)) {
        cat("  rownames: ", rownames, "\n", sep = "")
      } else if (!is.na(names)) {
        cat("  names: ", names, "\n", sep = "")
      }
      if (!is.na(details)) {
        cat("  details: ", details, "\n", sep = "")
      }
      if (!is.na(note)) {
        cat("  note: ", note, "\n", sep = "")
      }
    }
  )

  invisible(out)
}


#' Inspect a saved target workspace
#'
#' @param name Full target name, including the dynamic branch suffix when relevant.
#' @param envir Environment used by `targets` internals to materialize dependency frames.
#' @param store Targets store path.
#' @param print Logical; if TRUE, print a compact target and dependency summary.
#' @return Invisibly returns a list with target metadata, command, traceback, and
#'   dependency summaries.
#' @keywords internal

inspect_target_workspace <- function(name, envir = parent.frame(), store = targets::tar_config_get("store"), print = TRUE) {
  if (!is.character(name) || length(name) != 1L) {
    stop("`name` must be a length-one character vector.", call. = FALSE)
  }

  workspace <- targets:::workspace_read(name = name, path_store = store)
  dependency_env <- .target_workspace_dependency_env(workspace, envir)
  dependency_names <- sort(ls(dependency_env, all.names = TRUE))
  dependency_summaries <- purrr::map(
    dependency_names,
    \(dependency_name) {
      object <- get(dependency_name, envir = dependency_env, inherits = FALSE)
      .target_workspace_object_summary(dependency_name, object)
    }
  ) |>
    dplyr::bind_rows()

  meta <- targets::tar_meta(fields = tidyselect::any_of(c("name", "type", "parent", "error", "warnings", "seconds"))) |>
    dplyr::filter(.data$name == .env$name) |>
    dplyr::slice_head(n = 1)
  traceback <- workspace$target$metrics$traceback
  if (is.null(traceback)) {
    traceback <- character()
  }

  out <- list(
    name = name,
    meta = meta,
    command = .target_manifest_command(name),
    traceback = traceback,
    dependencies = dependency_summaries
  )

  if (isTRUE(print)) {
    .print_target_workspace_inspection(out)
  }

  invisible(out)
}


#' List currently errored target metadata
#'
#' @param target_name_pattern Regex used to filter target names.
#' @return Tibble of errored targets.
#' @keywords internal

list_errored_targets <- function(target_name_pattern = ".") {
  targets::tar_meta(fields = c("name", "error", "warnings")) |>
    dplyr::filter(
      !is.na(.data$error),
      .data$name %in% targets::tar_errored(),
      stringr::str_detect(.data$name, target_name_pattern)
    ) |>
    dplyr::arrange(.data$error) |>
    dplyr::select(dplyr::all_of(c("name", "error", "warnings")))
}

#' Print distinct errored targets
#'
#' @return Invisibly returns the distinct errored-target tibble.
#' @keywords internal

list_distinct_errored_targets <- function() {
  error_tibble <- list_errored_targets() |>
    dplyr::distinct(.data$error, .keep_all = TRUE)

  purrr::pwalk(
    error_tibble,
    function(name, error, warnings, ...) {
      cat(
        "- Full, tar_map-resolved target name: ",
        name,
        "\n",
        "- Suffix-stripped target name: ",
        .suffix_stripped_target_name(name),
        "\n",
        "- Error: ",
        stringr::str_trunc(error, width = 500),
        "\n",
        "- Warnings: ",
        stringr::str_trunc(warnings, width = 500),
        "\n\n",
        sep = ""
      )
    }
  )

  invisible(error_tibble)
}

#' Print distinct errored targets with tracebacks
#'
#' @param target_name_pattern Regex used to filter target names.
#' @return Invisibly returns the joined error, traceback, and command tibble.
#' @keywords internal

list_distinct_errored_targets_w_tracebacks <- function(target_name_pattern = ".") {
  trunc_width <- if (identical(target_name_pattern, ".")) 200 else 1000
  error_w_tracebacks_tibble <- list_errored_targets(target_name_pattern = target_name_pattern) |>
    dplyr::mutate(
      traceback_str = purrr::map_chr(
        .data$name,
        \(name) {
          .target_traceback(name) |>
            stringr::str_trunc(width = trunc_width) |>
            paste(collapse = "\n\t")
        }
      ),
      name_wo_dynamic_hash = stringr::str_remove(.data$name, "_[a-f0-9]{16}$")
    )

  capture.output(command_tibble <- targets::tar_manifest())

  combined_tibble <- error_w_tracebacks_tibble |>
    dplyr::left_join(command_tibble, by = c("name_wo_dynamic_hash" = "name"))

  purrr::pwalk(
    combined_tibble,
    function(name, error, warnings, command, traceback_str, ...) {
      cat(
        "- Full, tar_map-resolved target name: ",
        name,
        "\n",
        "- Suffix-stripped target name: ",
        .suffix_stripped_target_name(name),
        "\n",
        "- Command: ",
        command,
        "\n",
        "- Error: ",
        stringr::str_trunc(error, width = trunc_width),
        "\n",
        "- Traceback:\n\t",
        traceback_str,
        "\n",
        "- Warnings: ",
        stringr::str_trunc(warnings, width = trunc_width),
        "\n\n",
        sep = ""
      )
    }
  )

  invisible(combined_tibble)
}
