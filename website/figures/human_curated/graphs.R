# targets::tar_mermaid(
#   names = tidyselect::matches("embryonic_brain_mouse"),
# ) |>
#   targets_graph_mermaid_lines() |>
#   writeLines("parallel_reactions.mmd")

tar_glimpse(
  names = tidyselect::matches("GEX_Seurat_object.brain_mouse"),
  exclude = matches("embryonic_brain_mouse")
) |>
  targets_graph_mermaid_lines() |>
  writeLines("GEX.mmd")

targets::tar_glimpse(
  names = tidyselect::matches("multimodal_Seurat_object.brain_mouse"),
  exclude = matches("embryonic_brain_mouse|GEX"),
) |>
  targets_graph_mermaid_lines() |>
  writeLines("multimodal.mmd")


targets::tar_glimpse(
  names = tidyselect::matches("multimodal_Seurat_object.brain_mouse"),
) |>
  targets_graph_mermaid_lines() |>
  writeLines("full.mmd")

fs::file_copy(
  "full.mmd",
  "full_pruned.mmd",
  overwrite = TRUE
)
