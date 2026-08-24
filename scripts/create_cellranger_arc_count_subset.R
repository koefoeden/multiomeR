#!/usr/bin/env Rscript

repo_root <- rprojroot::find_root(rprojroot::has_file("_targets.R"))
setwd(repo_root)

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) {
  cat(paste0(
    "Create a deterministic, barcode-consistent Cell Ranger ARC smoke input.\n\n",
    "Usage:\n",
    "  create_cellranger_arc_count_subset.R --input=<count-dir> [options]\n\n",
    "Options:\n",
    "  --output=<dir>   Output directory (default: <count-dir>_subset).\n",
    "  --cells=<n>      Number of high-fragment Cell Ranger cells (default: 300).\n",
    "  --include-bams   Also filter and index the ATAC and GEX BAMs by CB tag.\n",
    "  --help           Show this help.\n"
  ))
  quit(save = "no", status = 0L)
}

known_args <- grepl("^--(input|output|cells)=", args) | args == "--include-bams"
if (any(!known_args)) {
  stop("Unknown argument(s): ", paste(args[!known_args], collapse = ", "), call. = FALSE)
}
argument_value <- function(name, required = FALSE, default = NULL) {
  matches <- args[startsWith(args, paste0("--", name, "="))]
  if (length(matches) > 1L || (required && length(matches) != 1L)) {
    stop("Pass exactly one --", name, "=<value> argument.", call. = FALSE)
  }
  if (!length(matches)) {
    return(default)
  }
  sub(paste0("^--", name, "="), "", matches)
}

input_dir <- argument_value("input", required = TRUE)
output_dir <- argument_value("output")
n_cells <- suppressWarnings(as.integer(argument_value("cells", default = "300")))
if (is.na(n_cells)) {
  stop("--cells must be an integer.", call. = FALSE)
}

source("R/cellranger_arc_subset_helpers.R")
create_cellranger_arc_count_subset(
  input_dir = input_dir,
  output_dir = output_dir,
  n_cells = n_cells,
  include_bams = "--include-bams" %in% args
)
