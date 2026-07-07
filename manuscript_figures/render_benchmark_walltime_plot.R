repo_root <- rprojroot::find_root(rprojroot::has_file("_targets.R"))
setwd(repo_root)

source("R/bootstrap_helpers.R")
load_project_runtime(force = TRUE)

aggregations <- c("benchmark_5x", "benchmark_10x", "benchmark_20x", "muscle", "multi_tissue")
store <- targets::tar_config_get("store")
output_dir <- file.path("manuscript_figures", "outputs", "benchmark")

# Set to character() to include optional AMULET and cellsnp-lite/vireo runtimes.
exclude_target_regex <- c(
  "^amulet_metrics_tibble(\\.|$)",
  "^(cellsnp_dir|vireo_donor_ids_tibble)(\\.|$)"
)

# Keep this TRUE when the final Seurat export target has not completed yet.
skip_endpoint_runtime <- TRUE

# Set TRUE to recompute from target metadata and raw CellRanger H5 files.
force_cache <- FALSE

target_step_labels <- c(
  Other = "Other critical-path steps",
  peaks_per_cluster_narrowPeaks.peaks.ATAC = "MACS3 peak-calling (ATAC)",
  fragments_per_peak_calling_cluster_discovery.fragments.ATAC = "MACS3 fragments-preparation (ATAC)",
  LSI_BPCells.ATAC = "LSI (ATAC)",
  consensus_peak_BPCells_matrix_dir.ATAC = "Consensus peak matrix (ATAC)",
  WNN_results = "Weighted nearest-neighbor graph (multimodal)",
  WNN_results_raw = "Raw WNN graph (multimodal)",
  metadata_w_QC_tibble.ATAC = "QC metadata (ATAC)",
  fragments_w_prefix_bpcells_dir = "BPCells fragment dir creation (ATAC)",
  UMAP_embeddings_tibble.GEX = "UMAP embeddings (GEX)",
  UMAP_embeddings_tibble.ATAC = "UMAP embeddings (ATAC)"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

benchmark_results <- cache_multimodal_seurat_walltime(
  aggregations = aggregations,
  cache_file = file.path(output_dir, "multimodal_seurat_walltime.rds"),
  force = force_cache,
  store = store,
  exclude_target_regex = exclude_target_regex,
  skip_endpoint_runtime = skip_endpoint_runtime
)

summary_tibble <- benchmark_results[, c(
  "aggregation",
  "reaction_count",
  "cellranger_input_nuclei",
  "critical_path_minutes",
  "serial_sum_minutes",
  "excluded_runtime_minutes",
  "endpoint_runtime_skipped",
  "endpoint_seconds_before_skip"
)]

print(summary_tibble)

plot <- plot_multimodal_seurat_walltime(
  benchmark_results,
  target_step_labels = target_step_labels,
  color_critical_path_steps = FALSE
)

ggplot2::ggsave(
  filename = file.path(output_dir, "multimodal_seurat_walltime.png"),
  plot = plot,
  width = 8.5,
  height = 5.2,
  dpi = 320
)

saveRDS(
  plot,
  file.path(output_dir, "multimodal_seurat_walltime_plot.rds")
)

readr::write_tsv(
  summary_tibble,
  file.path(output_dir, "multimodal_seurat_walltime.tsv")
)

cat("cache ", file.path(output_dir, "multimodal_seurat_walltime.rds"), "\n", sep = "")
cat("saved ", file.path(output_dir, "multimodal_seurat_walltime.png"), "\n", sep = "")
cat("saved ", file.path(output_dir, "multimodal_seurat_walltime_plot.rds"), "\n", sep = "")
cat("saved ", file.path(output_dir, "multimodal_seurat_walltime.tsv"), "\n", sep = "")
