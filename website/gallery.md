# Output gallery

```{r, include = FALSE}
source("../packages/multiomeRCore/R/null_default.R")
source("../R/output_gallery_helpers.R")
```

Each plot target has one example here, taken from `mixed_human_31x`: 31 public GEM wells from 10x Genomics and ENCODE covering heart, blood, pancreas, liver, colon, lung and cerebellum, with all optional analyses enabled. A card shows one of the files its target saves, the target's description and its name; select a preview to enlarge it. [Run your own analysis](main_running.md#steps) explains when to review each checkpoint.

```{r, echo = FALSE, results = "asis"}
render_output_gallery()
```
