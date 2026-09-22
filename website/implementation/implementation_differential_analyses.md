# Differential analyses

```{r setup, include = FALSE}
pipeline_name <- "differential_analyses"
source("helpers/_setup.R")
```

`module_differential_analyses/targets.R` filters aggregations that enabled the module, joins their module config, attaches symbols for accepted WNN metadata and pseudobulk inputs, and maps the composition, pseudobulk, gene-set enrichment, and cross-modality target fragments. The generic pseudobulk model family is instantiated for gene expression, chromatin accessibility, motif-family accessibility, and expression-derived CollecTRI activity (transcription-factor activity). transcription-factor activity first converts filtered, normalized GEX pseudobulks to signed ULM scores and then reuses the same model and contrast machinery.

The cross-modality fragment creates a CollecTRI-to-JASPAR family crosswalk, a detailed regulator-level table containing transcription-factor activity, motif-family accessibility, and TF-expression results, a family-level comparison table, a contrast-level concordance summary, and its plot. CollecTRI complexes remain intact in transcription-factor activity; complex-member mappings are introduced only by the comparison crosswalk.

The graph below is an orientation view. Inspect `setup_and_cell_type_composition_targets.R`, `pseudobulk_differential_targets.R`, `gene_set_enrichment_targets.R`, and `cross_modality_targets.R` for the complete model and plotting commands. The methods, model routes and every fixed or configurable setting are listed in [Differential analyses methods](methods_differential_analyses.md). The user-facing prerequisites and module selector are documented in [Differential analyses](../downstream_differential_analyses.html).

```{r, echo = FALSE, eval = TRUE, results = "asis"}
emit_mermaid("website/figures/human_curated/differential_analyses_v2.mmd")
```
