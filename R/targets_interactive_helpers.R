.suffix_stripped_target_name <- function(name) {
  name |>
    stringr::str_remove("\\..*") |>
    stringr::str_remove("_[a-f0-9]{16}$")
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
