```{.default filename="Output"}
save_plots_structured(plot_UMAP_from_metadata(metadata_w_cell_types_analysis_tibble.WNN.immune_human_2x,
     variable = categorical_UMAP_var.WNN.immune_human_2x, umap_cols = c("WNN_UMAP_1",
         "WNN_UMAP_2")), dyn_suffix_in_subdir = TRUE, override_suffix = stringr::str_replace_all(categorical_UMAP_var.WNN.immune_human_2x,
     "[/\\\\]", "_"))
```
