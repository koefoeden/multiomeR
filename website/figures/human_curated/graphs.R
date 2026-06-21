# Main pipeline
# Generates important processing intermediates, without being contamianted with downstream plots
# Can subsequently be used to subset into 4 parts (parallel, GEX, ATAC and WNN)
# using website/figures/human_curated/ATAC_keep_patterns.txt and similar for parallel, GEX and WNN
targets::tar_glimpse(
  names = tidyselect::matches("multimodal_Seurat_object.immune_human_2x"),
) |>
  targets_graph_mermaid_lines() |>
  writeLines("website/figures/human_curated/multimodal_Seurat_object.immune_human_2x.mmd")

# Genetic enrichment - across genes
targets::tar_glimpse(
  names = tidyselect::matches("^(ZScores_tibble[.]summarized[.]single_nucleus|TRS_summary_tibble).*immune_human_2x")
) |>
  targets_graph_mermaid_lines() |>
  writeLines("website/figures/human_curated/genetic_enrichment_across_genes.mmd")

# Genetic enrichment - per genes
targets::tar_glimpse(
  names = tidyselect::matches("^(ZScores_heatmap|coefficient_ranges|QC_plots|results_tibble|volcano_plot_files|locus_tracks_plot_files)[.](trait_level|gene_level)[.]pseudobulk.*immune_human_2x")
) |>
  targets_graph_mermaid_lines() |>
  writeLines("website/figures/human_curated/psbulk_GWAS.immune_human_2x.mmd")


# Differential analyses
# Needs
targets::tar_glimpse(
  names = tidyselect::matches("immune_human_2x") &
    targets::tar_described_as(tidyselect::contains("checkpoint:differential_analyses"))
) |>
  targets_graph_mermaid_lines() |>
  writeLines("website/figures/human_curated/checkpoint.differential_analyses.immune_human_2x.mmd")
