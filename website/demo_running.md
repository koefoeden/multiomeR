# Run the demo

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE
)
```

In the R session opened during installation, run the command below to process the `immune_human_2x` aggregation. It combines the two demo GEM wells, produces a Seurat/Signac object containing the multimodal results, and draws the integrated WNN UMAPs colored by cluster, cell type, and the other categorical metadata.

`names` selects these two targets by their exact names using `all_of()`. `tar_make()` also builds the dependencies needed for them, but does not build every plot in the [output gallery](gallery.md).

``` {.r filename="R"}
demo_targets <- c(
  "multimodal_Seurat_object.8_multimodal_QC.immune_human_2x",
  "categorical.UMAPs.8_multimodal_QC.immune_human_2x"
)

targets::tar_make(names = tidyselect::all_of(demo_targets))
```

Keep the R session open until the command finishes. Progress messages report targets being dispatched, completed, or skipped because they are already up to date. Using 16 threads, this should take \~ 30 minutes, writing about 6 GB to disk.

## Confirm success

After the run, the following command should return `character(0)`, meaning the requested results and their dependencies are up to date:

``` {.r filename="R"}
targets::tar_outdated(
  names = tidyselect::all_of(demo_targets),
  callr_function = NULL
)
```

If names are returned, those results still need building. If the run failed, follow [Troubleshooting](troubleshooting.md), fix the reported cause, and run the same `tar_make()` command again. Completed results can be reused.

Continue to [Inspect the demo results](demo_outputs.md) to read the object and find the associated files.