# Peak–gene correlation

```{r setup, include = FALSE}
pipeline_name <- "processing_and_aggregation"
source("helpers/_setup.R")
```

`module_peak_gene_correlation/targets.R` maps only opted-in aggregations, binds the primary-module inputs they consume, and maps the target files in `module_peak_gene_correlation/`. Parameters use the `peak_gene_correlation` manifest scope, targets end in `.peak_gene_correlation.<aggregation>`, and plot checkpoint tags use `peak_gene_correlation`, outside the numbered QC selections.

The methods and settings are in [Peak–gene correlation methods](methods_peak_gene_correlation.md), and the prerequisites and module selector in [Peak–gene correlation](../downstream_peak_gene_correlation.html).

## Candidate pairs, pseudobulks and tests

This view covers the TSS table, candidate peak–gene pairs, the broad WNN cell groups, donor–state pseudobulk construction, measurement-support filtering, the per-chromosome analysis branches and the diagnostic and top-link outputs.

```{r, echo = FALSE, eval = TRUE, results = "asis"}
emit_mermaid("website/figures/human_curated/peak_gene_correlation_v2.mmd")
```