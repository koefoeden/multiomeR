targets::tar_glimpse(
  names = tidyselect::matches("immune_human_2x"),
) |>
  targets_graph_mermaid_lines() |>
  writeLines("full.mmd")
