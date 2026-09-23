# Inspect the demo results

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE
)
```

The demo saved its results in the targets store, the `outputs/` folder of the clone (set by `store` in `_targets.yaml`; later pages write `<store>`). The store holds four kinds of output:

- `objects/`: R objects, such as tables and the Seurat/Signac object.
- `files/`: data files, such as matrix folders and TSV tables.
- `plots/`: PNG plot images.
- `plot_objects/`: an editable R copy of each plot.

This page opens one example of each kind in the R session. The [output reference](review_outputs.md#output-folders) describes the folders in detail.

## Objects

Read a stored object by its target name with `targets::tar_read()`, here the cell metadata and the multimodal Seurat/Signac object:

```{.r filename="R"}
cell_metadata <- targets::tar_read(
  metadata_w_cell_types_tibble.WNN.immune_human_2x
)
dim(cell_metadata)
head(cell_metadata)

demo_object <- targets::tar_read(
  multimodal_Seurat_object.8_multimodal_QC.immune_human_2x
)
demo_object
```

The metadata table has one row per retained nucleus, with its QC metrics, GEX, ATAC and WNN clusters, cell-type labels and UMAP coordinates. The Seurat/Signac object is a convenience export for exploring the results with Seurat and Signac. It reads its count matrices and ATAC fragments from files in the store and `example_data/`, so keep those folders in place.

## Files

For a file target, `tar_read()` returns the path of the saved file or folder. Open it with a suitable reader:

```{.r filename="R"}
targets::tar_read(cellranger_barcodes_tsv.healthy_PBMC_human)
#> [1] "outputs/files/healthy_PBMC_human/cellranger_barcodes_tsv.tsv"

GEX_counts <- BPCells::open_matrix_dir(
  targets::tar_read(aggregated_GEX_BPCells_matrix_dir.GEX.immune_human_2x)
)
GEX_counts
```

The first file belongs to one GEM well, the second to the aggregation. Paths follow the target name from right to left: `aggregated_GEX_BPCells_matrix_dir.GEX.immune_human_2x` is saved in `outputs/files/immune_human_2x/GEX/aggregated_GEX_BPCells_matrix_dir/`. `consensus_peak_BPCells_matrix_dir.ATAC.immune_human_2x` holds the ATAC peak counts.

## Plots

Plot targets return image paths in the same way, under `outputs/plots/`. The demo built one categorical UMAP per variable:

```{.r filename="R"}
targets::tar_read(categorical.UMAPs.8_multimodal_QC.immune_human_2x)
```

Open `WNN_harmony_SNN_cluster_cell_type.png` from that list to see the WNN clusters labeled by cell type. It should resemble this snapshot:

![WNN UMAP of the two demo GEM wells, colored by cluster and cell type](figures/demo_WNN_cell_type_UMAP.png){width="70%"}

Plot subtitles and captions explain how to read each plot.

## Plot objects

Each image has an R copy under `plot_objects/`, at the same relative path with `.rds` instead of `.png`. Reopen a plot without rerunning the analysis:

```{.r filename="R"}
plot_file <- file.path(
  targets::tar_config_get("store"),
  "plot_objects/immune_human_2x/8_multimodal_QC/UMAPs/categorical",
  "WNN_harmony_SNN_cluster_cell_type.rds"
)
p <- readRDS(plot_file)
p
```

Edit a ggplot object with the usual ggplot2 functions and save your copy outside the store, where a rerun cannot overwrite it:

```{.r filename="R"}
p <- p + ggplot2::labs(title = "My integrated cell types")
ggplot2::ggsave("my_cell_types.png", p, width = 10, height = 8)
```

Some plots are composites rather than single ggplot objects and need their own editing methods.

## Possible next steps

- Build the rest of the demo aggregation, including every checkpoint plot, with `targets::tar_make()`. The demo configuration has no other active aggregation and no optional modules.
- Browse the [output gallery](gallery.md) for an example of every plot the pipeline and its optional modules save.
- Start your own analysis with [Plan your analysis](main_overview.md).
