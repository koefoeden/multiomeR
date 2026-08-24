assert_command_on_path <- function(command) {
  if (file.exists(command) || nzchar(Sys.which(command))) {
    return(invisible(TRUE))
  }

  stop(
    "Command '", command, "' was not found on PATH. ",
    "Run through `pixi run ...` or add the tool to pixi.toml.",
    call. = FALSE
  )
}

run_w_error_check <- function(
  command_string,
  arguments_chr = character(),
  ...
) {
  assert_command_on_path(command_string)

  result_list <- processx::run(
    command = command_string,
    args = arguments_chr,
    error_on_status = FALSE,
    stdout_line_callback = function(line, proc) {
      cat(line, "\n", file = stdout())
      flush(stdout())
    },
    stderr_line_callback = function(line, proc) {
      cat(line, "\n", file = stderr())
      flush(stderr())
    },
    echo_cmd = TRUE,
    ...
  )

  if (result_list$status != 0) {
    format_process_stream <- function(x) {
      if (is.null(x) || !nzchar(x)) "<empty>" else x
    }

    base::stop(
      "Command failed (exit code ", result_list$status, "):\n",
      "Command:\n",
      stringr::str_flatten(c(command_string, shQuote(arguments_chr)), " "),
      "\n\nstderr:\n",
      format_process_stream(result_list$stderr),
      "\n\nstdout:\n",
      format_process_stream(result_list$stdout),
      call. = FALSE
    )
  }

  utils::tail(arguments_chr, 1)
}
