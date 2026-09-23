# Run the demo

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE
)
```

The supplied configuration has one active aggregation, `immune_human_2x`, which combines the two demo GEM wells. In the R session opened during installation, request two of its results: the multimodal Seurat/Signac object and UMAPs of the integrated WNN embedding, colored by cluster, cell type and other categorical metadata. WNN (weighted nearest neighbors) combines the GEX and ATAC data of each nucleus.

```{.r filename="R"}
demo_targets <- c(
  "multimodal_Seurat_object.8_multimodal_QC.immune_human_2x",
  "categorical.UMAPs.8_multimodal_QC.immune_human_2x"
)

targets::tar_make(names = tidyselect::all_of(demo_targets))
```

`tar_make()` builds these two targets and every result they depend on, from per-GEM-well quality control to WNN integration. It skips results they do not need, such as most checkpoint plots in the [output gallery](gallery.md).

With 16 threads the run takes about 30 minutes and writes about 6 GB. Keep the R session open until it finishes; progress messages show each target as it is dispatched and completed.

## Confirm success

This command should return `character(0)`, meaning that the requested results and their dependencies are up to date:

```{.r filename="R"}
targets::tar_outdated(
  names = tidyselect::all_of(demo_targets),
  callr_function = NULL
)
```

If target names are returned, run the same `tar_make()` command again; completed targets are reused. If the run stopped with an error, first find and fix the cause with [Troubleshooting](troubleshooting.md).

Continue to [Inspect the demo results](demo_outputs.md).
