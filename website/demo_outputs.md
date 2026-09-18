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

The pipeline produces 4 main kinds of outputs in the targets-store, defaulted to `outputs/` (can be configured via the `store` parameter in `_targets.yaml)`

- Serialized R objects (qs2-format) in a flat file hierarchy managed by targets itself inside `objects/`.

- Various file types inside `files/` in a structured folder hierarchy

- Plots in .png-format inside `plots/` in a structured folder hierarchy

- Serialized ggplot2-objects inside `plots_objects/` in a structured folder hierarchy

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

As you can see from the output above, file-targets are always saved in a location derived directly from their target-name.

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

\<WIP\>\
The [Main pipeline gallery](gallery_main.md) shows the other plot families the workflow produces. Those are saved snapshots and do not reflect the state of your analysis. To build one more of them, follow [Request an additional result](main_running.md#after-the-steps); to build all review plots for a stage, run its step in [Run your own analysis](main_running.md#steps).

## Possible next steps

- To continue with the demo-aggregation, and explore other outputs, try to run all targets in the pipeline by leaving out the `names` parameter in the `tar_make()`-call on the previous page.
- To adopt the workflow, continue with [Plan your analysis](main_overview.md) and the steps in [Run your own analysis](main_running.md#steps).
- To diagnose a failed or unexpectedly stale target, use [Troubleshooting](troubleshooting.md).