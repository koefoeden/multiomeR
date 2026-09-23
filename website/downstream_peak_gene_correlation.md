# Peak–gene correlation

## When to use

Use this module to nominate candidate regulatory links between accessible regions and nearby genes within each WNN cell type. It pairs consensus peaks with nearby gene transcription start sites and tests whether accessibility and expression vary together across **donor–state pseudobulks**: nuclei of one cell type summed by donor and ATAC-defined state, with each nucleus in at most one pseudobulk.

A hierarchical model adjusts for donor and sequencing depth and lets the peak–gene slope vary between donors. A separate conditional correlation analysis (HC3) provides the correlation summary plots and a SuSiE prioritization of peaks per gene. Links are hypotheses: neither analysis establishes causal regulation, and the hierarchical tests are approximate and have not been broadly calibrated.

## Prerequisites

Before enabling the module, confirm that:

- you have reviewed the aggregation through [checkpoint 8](main_running.md#checkpoint-8) and accept its final WNN cell set and cell-type annotations; and
- the cell types you want to study contain nuclei from several donors.

Cell types and donors with too few nuclei are skipped. A cell type from a single donor yields diagnostics but no tests, and the `strict` support filter requires more shared donors than the other presets.

## Outputs

| Question | Output |
|---|---|
| How many candidate pairs does each support filter retain? | Filter-retention plot |
| Which cell types or chromosome branches were skipped, and why? | Diagnostics plot and table |
| Which peaks are associated with a gene's expression? | Hierarchical results table and top-link figures |
| How do conditional correlations vary by cell type and distance? | HC3 summary plots and link table |
| Which peaks best explain a linked gene? | SuSiE prioritization table |

## Configure

Add `peak_gene_correlation` to the aggregation's [`modules`](parameters.html#modules) list, keeping its other settings:

```{.yaml filename="cfg_aggregations.yaml"}
my_aggregation:
  modules: [peak_gene_correlation]
```

Then add an entry for `my_aggregation` to `cfg_module_peak_gene_correlation.yaml` in the [selected configuration directory](main_overview.md#configuration-directory):

```{.yaml filename="cfg_module_peak_gene_correlation.yaml"}
my_aggregation:
  peak_gene_correlation_top_links_per_cell_group: 3
  peak_gene_correlation_filter: lenient
```

[`peak_gene_correlation_top_links_per_cell_group`](parameters.html#peak_gene_correlation_top_links_per_cell_group) sets the number of top-link figures per cell type; it does not change which pairs are tested. [`peak_gene_correlation_filter`](parameters.html#peak_gene_correlation_filter) selects the `lenient` (default), `moderate` or `strict` measurement-support filter, which removes pairs without enough expression, accessibility and shared donor support before testing. An empty entry, `my_aggregation: {}`, keeps both defaults.

## Run

Preview the targets tagged for this module:

```{.r filename="R"}
targets::tar_manifest(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:peak_gene_correlation")
  ) & tidyselect::ends_with(".my_aggregation"),
  callr_function = NULL
)[, c("name", "description")]
```

Then run the same selection:

```{.r filename="R"}
targets::tar_make(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:peak_gene_correlation")
  ) & tidyselect::ends_with(".my_aggregation")
)
```

## Review

Review the plots ([examples](gallery.md#peak-gene-correlation)); each plot's subtitle and caption say what to look for.

```text
<store>/plots/my_aggregation/peak_gene_correlation/
├── filter_retention_plot.png
├── diagnostics_plot.png
├── top_link_aggregate_scatter_plots/
├── support_counts_plot.png
├── correlation_histogram_plot.png
├── distance_correlation_plot.png
└── significant_pairs_vs_technical_features_plot.png
```

Start with `filter_retention_plot.png`, which compares the three support filters by cell type and marks the active one. Top-link figures rank positive, estimable slopes outside self-promoter peaks by hierarchical p-value, without a significance cutoff, so appearing in a figure is not evidence of significance. `diagnostics_plot.png` shows skipped branches, and the other plots summarize the HC3 analysis.

Read the tables in R, for example the hierarchical results with BH FDR within each cell type:

```{.r filename="R"}
targets::tar_read(
  peak_gene_correlation_hierarchical_results_tibble.WNN.peak_gene_correlation.my_aggregation
)
```

The HC3 links and SuSiE prioritization are in `peak_gene_correlation_links_tibble.WNN` and `peak_gene_correlation_finemapped_links_tibble.WNN`, with the same `.peak_gene_correlation.my_aggregation` suffix.

## Parameter reference

[Open the searchable parameter browser](parameters.html#workflow=peak_gene_correlation). The [peak–gene correlation implementation page](implementation/methods_peak_gene_correlation.html) describes the methods and their key fixed values, including the support-filter presets, links to the source files and shows the target structure.
