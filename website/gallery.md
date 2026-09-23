# Output gallery

```{r, include = FALSE}
source("../packages/multiomeRCore/R/null_default.R")
source("../R/output_gallery_helpers.R")
```

This page shows one saved plot from each plot target of the `mixed_human_31x` aggregation: 31 public GEM wells from 10x Genomics and ENCODE covering heart, blood, pancreas, liver, colon, lung and cerebellum, with all three optional modules enabled. The public default configuration includes this aggregation as an inactive example; its raw data are public, but the reprocessed `cellranger-arc count` outputs are not supplied.

Each card shows that plot, the target's description and its name without the `.mixed_human_31x` suffix. Select a preview to enlarge it. The [running guide](main_running.md#steps) explains when to review each checkpoint.

```{r, echo = FALSE, results = "asis"}
render_output_gallery()
```
