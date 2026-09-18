# Inspect the demo results

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE
)
source("../packages/multiomeRCore/R/null_default.R")
source("../R/output_gallery_helpers.R")
gallery_items <- check_output_gallery_assets(
  manifest_file = "output_gallery.yaml",
  gallery_root = "."
)
```

The pipeline saves four main kinds of output in its targets store, normally `outputs/`. The `store` setting in `_targets.yaml` selects this folder.

- Serialized R objects in `objects/`, managed by targets.
- Data files in `files/`, grouped by target and analysis.
- Plot images in `plots/`, normally in PNG format.
- Editable plot objects in `plot_objects/`, saved as RDS files alongside the corresponding image hierarchy.

## Objects

Most intermediate and final result objects are saved automatically by `targets` and can be loaded into any repository-root R-session using `targets::tar_read()`:

``` {.r filename="R"}
cell_metadata <- targets::tar_read(
  metadata_w_cell_types_tibble.WNN.immune_human_2x
)
dim(cell_metadata)
head(cell_metadata)

demo_object <- targets::tar_read(multimodal_Seurat_object.8_multimodal_QC.immune_human_2x)
demo_object
```

The metadata table describes the retained nuclei and their annotations. WNN means *weighted nearest neighbors*: the integrated representation uses information from both RNA and ATAC. The Seurat/Signac object is a convenient export for further exploration; the pipeline also retains its matrices in BPCells format on disk.

## Files

File targets also load with `targets::tar_read()`, but their value is a path rather than an in-memory result, so you will have to load them yourself using the appropriate tool.

``` {.r filename="R"}
targets::tar_read(cellranger_barcodes_tsv.healthy_PBMC_human)
targets::tar_read(aggregated_GEX_BPCells_matrix_dir.GEX.immune_human_2x)
targets::tar_read(consensus_peak_BPCells_matrix_dir.ATAC.immune_human_2x)
```

As you can see from the output above, pipeline-generated files generally follow a folder hierarchy derived from their target names; `tar_read()` gives their actual paths.

## Plots

The demo command also built `categorical.UMAPs.8_multimodal_QC.immune_human_2x`. Plots use the same scope-based layout as files, but under `outputs/plots/` instead. Note that a target might produce multiple files, as seen in the example below:

``` {.r filename="R"}
targets::tar_read(categorical.UMAPs.8_multimodal_QC.immune_human_2x)
```

Open `WNN_harmony_SNN_cluster_cell_type.png` there to see the integrated clusters and cell-type labels from your own run. It should resemble this documentation snapshot:

```{r, echo = FALSE, eval = TRUE, results = "asis"}
render_gallery_grid(gallery_items[gallery_items$id == "wnn-umap", ])
```

## Plot objects

Each plot saved by the standard plotting helper also has an `.rds` copy under `plot_objects/`, unless plot-object saving was disabled. This lets you reopen a plot in R without repeating the analysis. For example:

``` {.r filename="R"}
plot_file <- file.path(
  targets::tar_config_get("store"),
  "plot_objects/immune_human_2x/8_multimodal_QC/UMAPs/categorical",
  "WNN_harmony_SNN_cluster_cell_type.rds"
)
p <- readRDS(plot_file)
p
```

For a ggplot object, edit it with the usual ggplot2 functions and save a separate copy:

``` {.r filename="R"}
p <- p + ggplot2::labs(title = "My integrated cell types")
ggplot2::ggsave("my_cell_types.png", p, width = 10, height = 8)
```

Some outputs are composite plots rather than ordinary ggplot objects and need their own editing methods. Keep custom exports separate from pipeline outputs, which can be overwritten on a rerun.

## Possible next steps

- To continue with the demo-aggregation, and explore other outputs, run `targets::tar_make()` without `names`.
- To get started with your own data, please continue at [Plan your analysis](main_overview.md).
- To diagnose a failed or unexpectedly stale target, use [Troubleshooting](troubleshooting.md).