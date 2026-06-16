manifest <- targets::tar_manifest(callr_function = NULL)
graph_ids <- targets_graph_part_of_graph_ids(manifest)
label_suffixes <- targets_graph_default_label_suffixes()
output_dir <- file.path("website", "figures", "human_curated")

for (graph_id in graph_ids) {
  output_file <- file.path(output_dir, paste0(graph_id, "_v2.mmd"))
  targets_graph_write_part_of_graph_mermaid(
    graph_id = graph_id,
    output_file = output_file,
    manifest = manifest,
    label_suffixes = label_suffixes
  )
  message("wrote ", output_file)
}
