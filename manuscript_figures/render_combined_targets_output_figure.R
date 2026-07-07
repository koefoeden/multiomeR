repo_root <- rprojroot::find_root(rprojroot::has_file("_targets.R"))
setwd(repo_root)

source("R/bootstrap_helpers.R")
load_project_runtime()

library(patchwork)

figure_2_UMAP_modifier <- NA

output_dir <- "manuscript_figures/outputs/combined_targets_output_figures"

panel_specs <- tibble::tribble(
  ~figure_number , ~tag , ~plot_call                                                                                                                  , ~pick_regex                          , ~title        , ~caption                                                                                                                              , ~plot_modifier                              ,
  "S1"           , "A"  , quote(targets::tar_read(categorical.UMAPs.WNN.PBMC_human_6x))                                                               , "WNN_harmony_SNN_cluster_cell_type"  , NA_character_ , "UMAP embedding derived from native RNA+ATAC weighted-nearest-neighbor graph. Nuclei are colored by cluster-level cell type labels."  , NA                                          ,
  "S1"           , "B"  , quote(targets::tar_read(markers_dot_plot.GEX.PBMC_human_6x))                                                                , NA_character_                        , NA_character_ , "Expression of canonical marker genes across cell types."                                                                             , quote(plot + ggplot2::labs(y = "celltype")) ,
  "S1"           , "C"  , quote(targets::tar_read(coverage_tracks_plots.ATAC.PBMC_human_6x))                                                          , "MS4A1"                              , NA_character_ , "Chromatin accessibility track at the MS4A1 locus."                                                                                   , NA                                          ,
  "2"            , "A"  , quote(targets::tar_read(harmony.categorical.UMAPs.GEX.PBMC_human_6x))                                                       , "PCA_harmony_SNN_cluster_cell_type"  , NA_character_ , "RNA UMAP colored by GEX-derived cluster-level cell type labels."                                                                     , figure_2_UMAP_modifier                      ,
  "2"            , "B"  , quote(targets::tar_read(categorical.UMAPs.ATAC.PBMC_human_6x))                                                              , "LSI_harmony_SNN_cluster_cell_type"  , NA_character_ , "ATAC UMAP colored by ATAC-derived cluster-level cell type labels."                                                                   , figure_2_UMAP_modifier                      ,
  "2"            , "C"  , quote(targets::tar_read(categorical.UMAPs.WNN.PBMC_human_6x))                                                               , "WNN_harmony_SNN_cluster_cell_type"  , NA_character_ , "WNN UMAP colored by multimodal cluster-level cell type labels."                                                                      , figure_2_UMAP_modifier                      ,
  "2"            , "D"  , quote(targets::tar_read(TF_activity_dot_plot.ATAC.PBMC_human_6x))                                                           , NA_character_                        , NA_character_ , "ChromVAR transcription-factor activity scores for configured PBMC marker TFs."                                                        , NA                                          ,
  "2"            , "E"  , quote(targets::tar_read(peak_gene_correlation_top_link_ATAC_tracks_plots.peak_gene_correlation.WNN.PBMC_human_6x))           , "B_rank001_ODF2[.]AS1_chr9"          , NA_character_ , "Example top peak-gene link with ATAC accessibility and loop track."                                                                  , NA                                          ,
  "2"            , "F"  , quote(targets::tar_read(TRS_UMAPs.WNN_harmony_SNN.SCAVENGE.single_nucleus.genetic_enrichment.PBMC_human_6x))                , "MultipleSclerosis_IMSGC2019"        , NA_character_ , "SCAVENGE trait relevance scores for multiple sclerosis projected onto the WNN UMAP."                                                 , figure_2_UMAP_modifier                      ,
  "2"            , "G"  , quote(targets::tar_read(ZScores_heatmap.trait_level.pseudobulk.genetic_enrichment.PBMC_human_6x))                           , NA_character_                        , NA_character_ , "Pseudobulk GWAS chromVAR Z-scores across PBMC cell types and GWAS traits."                                                           , NA                                          ,
)

figure_specs <- tibble::tribble(
  ~figure_number , ~height , ~layout                                , ~subtitle                                                                            ,
  "S1"           ,       9 , quote((A / B / C))                     , "Visualization of the observed cell types and expected markers in the PBMC dataset." ,
  "2"            ,       7 , quote((A | B | C) / (D | E) / (F | G)) , NA_character_                                                                        ,
)

expected_pngs <- file.path(output_dir, paste0(figure_specs$figure_number, ".png"))
stale_pngs <- setdiff(list.files(output_dir, pattern = "[.]png$", full.names = TRUE), expected_pngs)
unlink(stale_pngs)

shared_render_figures(
  panel_specs = panel_specs,
  figure_specs = figure_specs,
  output_dir = output_dir
)
