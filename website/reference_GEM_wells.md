# GEM well table

```{r setup, include = FALSE}
pipeline_name <- "processing_and_aggregation"
source("helpers/_setup.R")
```

`cfg_GEM_wells.tsv` in the [selected configuration directory](main_overview.md#configuration-directory) has one row per GEM well, that is, one `cellranger-arc count` output. You fill it in at [step 1 of Configure your data](main_running.md#steps), and aggregations select wells by their `GEM_well_ID`.

## Columns

Copy a row of the [public example](#public-example), give it a unique `GEM_well_ID`, and fill in these columns:

```{r, echo = FALSE, eval = TRUE, results = "asis"}
emit_GEM_well_column_table("website/data/GEM_well_columns.tsv")
```

Add library or batch annotations as extra columns prefixed with `GEM_well_`, such as `GEM_well_multiplex_batch`. They become cell metadata: select them for plots with [`aggregation_categorical_vars`](parameters.html#aggregation_categorical_vars) or [`aggregation_continuous_vars`](parameters.html#aggregation_continuous_vars), and for batch correction with [`aggregation_harmony_correction_metadata_col_names`](parameters.html#aggregation_harmony_correction_metadata_col_names). Keep donor phenotypes in the [donor metadata table](reference_donor_metadata.md).

## Cell Ranger inputs {#cellranger-inputs}

Each `GEM_well_cellranger_arc_count_dir` must contain:

``` text
outs/
├── summary.csv
├── filtered_feature_bc_matrix.h5
├── atac_fragments.tsv.gz
├── atac_fragments.tsv.gz.tbi
└── per_barcode_metrics.csv
```

The pipeline is tested with the outputs of Cell Ranger ARC 2.0.0 and 2.1.0.

A [pooled GEM well](#pooled-wells) also needs `atac_possorted_bam.bam` in `outs/`. For CellBender counts, run [CellBender remove-background](https://cellbender.readthedocs.io/en/latest/usage/) first and give its H5 file in `GEM_well_cellbender_h5_file`.

The pipeline identifies each well's Cell Ranger ARC reference by matching the FASTA and GTF hashes in the `atac_fragments.tsv.gz` header to a `reference.json` in the repository's `reference_metadata/` folder, so keep that header intact. JSON files for the GRCh38 2020-A, GRCh38 2024-A and mm10 2020-A references are included; for another reference, copy its `reference.json` into a new subdirectory there.

## Pooled GEM wells {#pooled-wells}

A GEM well that pools several donors needs a donor-genotype VCF in `GEM_well_donors_VCF_file`, the number of donors in `GEM_well_n_donors`, and `NA` in `GEM_well_donor_id`; the [Vireo documentation](https://vireosnp.readthedocs.io/en/stable/manual.html) describes the VCF. The pipeline genotypes each nucleus with cellsnp-lite and assigns it with Vireo. Every nucleus receives its most likely single donor, named by the VCF's sample name, and `vireo_type` records whether Vireo called it a singlet, a doublet or unassigned; remove doublets and unassigned nuclei with a [QC filter](#qc-filters) such as `vireo_type != "singlet"`.

Without a VCF, no demultiplexing runs and every nucleus receives `GEM_well_donor_id`, so a pooled well needs its VCF to obtain donors.

## Rules checked {#rules}

When the pipeline starts, it checks that:

- `GEM_well_ID` values are unique;
- every GEM well that an active aggregation selects exists and is active; and
- each QC filter is a valid R expression.

When the targets of a well run, exactly one `reference.json` must match its fragments, and all wells of an aggregation must share the same reference. Apart from the keys, no column name may appear both here and in the [donor metadata table](reference_donor_metadata.md#rules).

## Set QC filters after the first run {#qc-filters}

Leave `GEM_well_QC_exclude_list` as `NA` for the first run. After reviewing the QC distributions at [checkpoint 1](main_running.md#checkpoint-1), enter filters and rerun checkpoint 1. Each filter is a complete R expression over the per-nucleus metadata; separate filters with `;;`:

``` {.text filename="cfg_GEM_wells.tsv"}
TSS.enrichment < 4 ;; nucleosome_signal > 4 ;; nCount_RNA < 250
```

A nucleus for which any expression is `TRUE` is excluded, and each expression is reported as a separate exclusion reason in the checkpoint 1 UpSet and retention plots. Filters can use the metrics that [`QC_metric_manifest.tsv`](https://github.com/koefoeden/multiomeR/blob/main/QC_metric_manifest.tsv) lists as available from checkpoint 1, and other per-nucleus columns such as `vireo_type`. The checkpoint 1 comparison plots show each metric's distribution and draw simple cutoffs such as these; they are examples, not recommendations for your tissue.

## Public example {#public-example}

The table below shows the two public demo wells with the required columns in bold and one optional annotation, `GEM_well_cell_sorting`. Scroll horizontally, and focus or hover over a column's **i** button for its meaning. The [committed example](https://github.com/koefoeden/multiomeR/blob/main/configuration/cfg_GEM_wells.tsv) contains further inactive rows and annotation columns.

```{r, echo = FALSE, eval = TRUE, results = "asis"}
emit_GEM_well_demo_table(
  GEM_well_config_file = "website/data/demo_GEM_wells.tsv",
  dictionary_file = "website/data/GEM_well_columns.tsv"
)
```
