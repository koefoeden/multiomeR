targets::tar_glimpse(
  names = tidyselect::matches("multimodal_Seurat_object.brain_mouse"),
) |>
  targets_graph_mermaid_lines() |>
  writeLines("full.mmd")
