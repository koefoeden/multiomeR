# Peak–gene correlation

```{r setup, include = FALSE}
pipeline_name <- "processing_and_aggregation"
source("helpers/_setup.R")
```

`module_peak_gene_correlation/targets.R` maps only opted-in aggregations and binds their existing WNN metadata, GEX and ATAC matrices, ATAC embeddings, fragments and reference annotations. `correlation_targets.R` owns the candidate-pair construction, donor–state pseudobulking, filtering, the conditional and hierarchical analyses, SuSiE prioritization, exports and plots.

Parameters use the `peak_gene_correlation` manifest scope and matching module YAML rows. Targets end in `.peak_gene_correlation.<aggregation>`, with `.WNN` before that suffix for intermediate results. Plot checkpoint tags use `peak_gene_correlation`, keeping this analysis outside the numbered QC selections.

The fixed thresholds and the configurable settings of this module are listed in [Peak–gene correlation methods](methods_peak_gene_correlation.md). The user-facing prerequisites and module selector are documented in [Peak–gene correlation](../downstream_peak_gene_correlation.html).

## Candidate pairs, pseudobulks and tests

This view covers the TSS table, candidate peak–gene pairs, the broad WNN cell groups, donor–state pseudobulk construction, measurement-support filtering, the per-chromosome analysis branches and the diagnostic and top-link outputs.

```{r, echo = FALSE, eval = TRUE, results = "asis"}
emit_mermaid("website/figures/human_curated/peak_gene_correlation_v2.mmd")
```