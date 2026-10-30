# Peak–gene correlation methods

```{r setup, include = FALSE}
pipeline_name <- "processing_and_aggregation"
source("helpers/_setup.R")
```

This page describes the models of the optional peak–gene correlation module. [Peak–gene correlation](downstream_peak_gene_correlation.md) describes its prerequisites, configuration and outputs, and its settings and defaults are listed in the [parameter browser](parameters.html#workflow=peak_gene_correlation).

## Pseudobulks and support filters {#pseudobulks-and-support-filters}

{{< include _shared_methods/peak_gene_pseudobulks.md >}}

## Models and candidate links {#models-and-candidate-links}

{{< include _shared_methods/peak_gene_models.md >}}

The model fits are compared with `lme4` and `pbkrtest` in [Algorithm validation](algorithm_validation.md). The ABC contact exponent and the 5 kb and 5 Mb distances are the defaults of the [ABC pipeline](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction). Its current version also sets self-promoter scores to 1 after the division. This module keeps the original ratio, so a gene's scores sum to at most 1.

## How links are checked {#link-checks}

No check in the module validates an individual link. Each addresses a different way the links could mislead:

- **Model fits.** The hierarchical fits agree with `lme4` and `pbkrtest` to within 1e-6, as [Algorithm validation](algorithm_validation.md) lists.
- **Shared chromatin states.** The accessibility components and each gene's background of distant peaks remove correlation that a gene shares with many peaks. In `distance_correlation_plot.png`, the median calibrated −log10 P of pairs far from the TSS should stay near the background's 0.30. A curve that stays higher far from the gene points to shared variation that the calibration misses.
- **Power.** `significant_pairs_vs_technical_features_plot.png` compares each cell type's count of FDR-significant pairs with its pseudobulks, depth, detection and residual degrees of freedom. Counts that follow these features show power, not more regulation.
- **Intragenic associations.** A peak in the target gene's own body can vary with the gene's transcription, so its association need not mark an enhancer. These links have the `link_type` `intragenic_association` and are counted apart in `support_counts_plot.png`. Only the others are candidate enhancer links.
- **Independent evidence.** ABC scores rank the same pairs from accessibility and distance alone, without donor variation. With the genetic-enrichment module, S-LDSC and fine-mapped variants test whether the top peaks of the ABC and correlation rankings hold more GWAS signal than their territory and than rankings by activity, accessibility share or distance ([S-LDSC methods](methods_genetic_enrichment.md#heritability-enrichment)). A ranking that beats these baselines supports the ranking as a whole, not any one link.

The module does not compare links with CRISPR perturbations or fine-mapped eQTLs, so they remain candidates for such tests.

## Source and target graph

`extra_targets/module_peak_gene_correlation/targets.R` maps only opted-in aggregations and binds the primary-module inputs they consume. Parameters use the `peak_gene_correlation` manifest scope, targets end in `.peak_gene_correlation.<aggregation>`, and plot checkpoint tags use `peak_gene_correlation`, outside the numbered QC selections.

**Source:** [`extra_targets/module_peak_gene_correlation/`](https://github.com/koefoeden/multiomeR/tree/main/extra_targets/module_peak_gene_correlation), [`R/peak_gene_correlation_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/peak_gene_correlation_helpers.R), [`R/peak_gene_filter_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/peak_gene_filter_helpers.R), [`R/peak_gene_hierarchical_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/peak_gene_hierarchical_helpers.R), [`R/peak_gene_ABC_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/peak_gene_ABC_helpers.R), [`src/peak_gene_REML.cpp`](https://github.com/koefoeden/multiomeR/blob/main/src/peak_gene_REML.cpp).

This view covers the TSS table, candidate peak–gene pairs, the WNN cell groups, donor–state pseudobulk construction, measurement-support filtering, the per-chromosome analysis branches and the diagnostic and top-link outputs.

```{r, echo = FALSE, eval = TRUE, results = "asis"}
emit_mermaid("website/figures/human_curated/peak_gene_correlation_v2.mmd")
```
