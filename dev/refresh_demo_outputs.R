# Refresh the printed output under the R blocks of website/demo_outputs.md, and
# its UMAP snapshot, from the built public demo. Run it from the root of a
# checkout whose targets store holds the demo: it evaluates the blocks there and
# writes into the repository that contains this script.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
website <- file.path(dirname(dirname(normalizePath(script))), "website")
page <- paste(readLines(file.path(website, "demo_outputs.md")), collapse = "\n")
block_pattern <- '(?s)```\\{\\.r filename="R"\\}\n(.*?)\n```\n\n\\{\\{< include (data/demo_outputs/\\w+\\.md) >\\}\\}'
blocks <- regmatches(page, gregexpr(block_pattern, page, perl = TRUE))[[1]]
stopifnot(length(blocks) > 0)

options(width = 70, cli.num_colors = 1, cli.unicode = FALSE)
env <- new.env()
for (block in blocks) {
  parts <- regmatches(block, regexec(block_pattern, block, perl = TRUE))[[1]]
  output <- capture.output(for (expr in parse(text = parts[[2]])) {
    result <- withVisible(eval(expr, env))
    if (result$visible) print(result$value)
  })
  output <- sub("\\s+$", "", gsub(paste0(normalizePath(getwd()), "/"), "", output, fixed = TRUE))
  writeLines(c('```{.default filename="Output"}', output, "```"), file.path(website, parts[[3]]))
}

UMAPs <- unlist(targets::tar_read(categorical.UMAPs.8_multimodal_QC.immune_human_2x))
invisible(file.copy(
  grep("/WNN_harmony_SNN_cluster_cell_type[.]png$", UMAPs, value = TRUE),
  file.path(website, "figures", "demo_WNN_cell_type_UMAP.png"),
  overwrite = TRUE
))
