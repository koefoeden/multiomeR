# Peak–gene correlation

Run this optional module after accepting the final WNN cell set. It relates ATAC accessibility to RNA expression within broad GEX-derived cell types, using the retained WNN nuclei; it is no longer part of checkpoint 8.

## Configure

Add `peak_gene_correlation` to the aggregation's existing [`modules`](parameters.html#modules) list in `cfg_aggregations.yaml`, and add a matching row in `cfg_module_peak_gene_correlation.yaml` in the selected configuration directory:

```yaml
my_aggregation:
  peak_gene_correlation_top_links_per_cell_group: 3
```

The top-link count controls the number of detail figures per cell type. [`peak_gene_correlation_filter`](parameters.html#peak_gene_correlation_filter) selects `lenient` (default), `moderate`, or `strict` measurement-support filtering. Disabled aggregations contribute no module targets.

## Run

```r
targets::tar_make(
  names = tidyselect::ends_with(".peak_gene_correlation.my_aggregation")
)
```

## Review

All module targets have description tag `[checkpoint:peak_gene_correlation]`. Their paths are `<store>/plots/my_aggregation/peak_gene_correlation/`; file exports use the corresponding `files` directory, with any modality suffixes as deeper subdirectories. For example, read selected links with:

```r
targets::tar_read(
  peak_gene_correlation_links_tibble.WNN.peak_gene_correlation.my_aggregation
)
```

## Model and scope

Peak–gene links are candidate regulatory relationships. Cells are aggregated within donor and ATAC-defined state using WNN cell-type annotations, without reusing a cell across aggregates. Measurement-support filters select hypotheses before fitting; the default requires shared support from at least two donors, and the strict preset requires three. Review `filter_retention_plot` first.

The hierarchical analysis fits a mean peak effect with donor-specific slope variation, donor intercepts and RNA/ATAC depth adjustment. It reports Kenward–Roger p-values, BH FDR and numerical reliability diagnostics. Top-link figures rank positive, estimable nonpromoter slopes by p-value without a significance cutoff, so appearing in a figure is not evidence of significance. They combine focal-cell-type coverage, gene context, hierarchical evidence and adjusted aggregate scatterplots.

The existing HC3 correlation summaries remain a separate conditional analysis. Neither analysis establishes causal regulation; numerical reference parity does not establish statistical calibration across datasets. See the [module methods and filtering reference](https://github.com/koefoeden/multiomeR/blob/main/module_peak_gene_correlation/README.md) for inference limits, support thresholds and output details.
