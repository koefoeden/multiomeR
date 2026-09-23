# Recompute the parked UCell annotation perturbation diagnostics from stored
# targets. Run from the repository root:
#   Rscript dev/cluster_annotation_diagnostics/run_diagnostics.R <aggregation> [GEX|ATAC|WNN]
args <- commandArgs(trailingOnly = TRUE)
aggregation <- args[[1]]
modality <- if (length(args) > 1L) args[[2]] else "GEX"
cluster_column <- c(GEX = "PCA_harmony_SNN_cluster", ATAC = "LSI_harmony_SNN_cluster",
  WNN = "WNN_harmony_SNN_cluster")[[modality]]
source("dev/cluster_annotation_diagnostics/helpers.R")

read_target <- function(name) targets::tar_read_raw(paste(name, aggregation, sep = "."))
control <- read_target("cluster_UCell_controls.GEX")
min_advantage <- read_aggregation_config_tibble() |>
  dplyr::filter(.data$aggregation == .env$aggregation) |>
  dplyr::pull("aggregation_cluster_annotation_min_advantage") |>
  unlist()
summaries <- summarize_cluster_UCell_diagnostic_counts(
  read_target("aggregated_counts_BPCells_matrix.GEX"),
  read_target(paste0("metadata_w_clusters_tibble.", modality)),
  control,
  cluster_column
)
diagnostics <- evaluate_cluster_UCell_diagnostics(
  score_cluster_UCell_diagnostic_summaries(summaries, control),
  min_advantage
)

output_dir <- file.path(targets::tar_config_get("store"), "files", aggregation,
  "cluster_UCell_perturbation_diagnostics", modality)
fs::dir_create(output_dir)
readr::write_tsv(diagnostics$decisions, file.path(output_dir, "clusters.tsv"))
if (nrow(diagnostics$sample_decisions)) {
  readr::write_tsv(diagnostics$sample_decisions, file.path(output_dir, "GEM_well_agreement.tsv"))
}
cat("wrote ", output_dir, "\n", sep = "")
