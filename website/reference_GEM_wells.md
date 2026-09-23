# GEM well table

```{r setup, include = FALSE}
pipeline_name <- "processing_and_aggregation"
source("helpers/_setup.R")
```

`cfg_GEM_wells.tsv` in the [selected configuration directory](main_overview.md#configuration-directory) has one row per GEM well, that is, one `cellranger-arc count` output. Aggregations select wells by their `GEM_well_ID`.

## Fill in a row

Copy an example row, give it a unique `GEM_well_ID`, and fill in these columns:

| Column | What to enter |
|---|---|
| `GEM_well_ID` | A unique identifier used in target names and output folders. |
| `GEM_well_dataset` | A label for the dataset or study. |
| `GEM_well_cellranger_arc_count_dir` | The directory containing `outs/`, not `outs/` itself. |
| `GEM_well_n_donors` | Number of donors in the well. |
| `GEM_well_donor_id` | For a single-donor well, an ID matching the [donor metadata table](reference_donor_metadata.md); otherwise `NA`. |
| `GEM_well_donors_VCF_file` | For a multiplexed well, the donor-genotype VCF used for demultiplexing; otherwise `NA`. |
| `GEM_well_add_cellbender` | `TRUE` to use externally generated CellBender counts; otherwise `FALSE`. |
| `GEM_well_cellbender_h5_file` | The CellBender H5 file when `GEM_well_add_cellbender` is `TRUE`; otherwise `NA`. |
| `GEM_well_QC_exclude_list` | `NA` for the first run; later, [QC filters](#qc-filters). |
| `GEM_well_is_active` | `TRUE` to process the well; `FALSE` for unused rows, including unused demo wells. Every well selected by an active aggregation must be active. |

Add library or batch annotations as extra columns prefixed with `GEM_well_`, such as `GEM_well_multiplex_batch`. They become cell metadata for batch correction and plots. Keep donor phenotypes in the donor metadata table; apart from their keys, the two tables must not share column names.

## Required inputs

Each `GEM_well_cellranger_arc_count_dir` must contain:

``` text
outs/
├── summary.csv
├── filtered_feature_bc_matrix.h5
├── atac_fragments.tsv.gz
├── atac_fragments.tsv.gz.tbi
└── per_barcode_metrics.csv
```

Genotype demultiplexing also requires `atac_possorted_bam.bam` in `outs/`. Prepare optional VCF and CellBender inputs as described in [Plan your analysis](main_overview.md).

The pipeline identifies each well's Cell Ranger ARC reference by matching the FASTA and GTF hashes in the `atac_fragments.tsv.gz` header to a `reference.json` under `reference_metadata/`, so keep that header intact. JSON files for the GRCh38 2020-A, GRCh38 2024-A and mm10 2020-A references are included; for another reference, copy its `reference.json` into a new subdirectory there. Exactly one JSON must match, and all wells in an aggregation must share the same reference.

## Set QC filters after the first run {#qc-filters}

Leave `GEM_well_QC_exclude_list` as `NA` for the first run. After reviewing the QC distributions at [checkpoint 1](main_running.md#checkpoint-1), enter filters and rerun checkpoint 1. Each filter is a complete R expression over the per-nucleus QC metrics; separate filters with `;;`:

``` {.text filename="cfg_GEM_wells.tsv"}
TSS.enrichment < 4 ;; nucleosome_signal > 4 ;; nCount_RNA < 250
```

Nuclei for which an expression is `TRUE` are excluded, and each expression is reported as a separate exclusion reason. The checkpoint 1 GEM well comparison plots draw simple cutoffs such as these on each metric's distribution. These cutoffs are examples, not recommendations for your tissue.

## Example rows

The table below shows the two public demo wells with the required columns in bold and one optional annotation, `GEM_well_cell_sorting`. Scroll horizontally, and focus or hover over a column's **i** button for its meaning. The [committed example](https://github.com/koefoeden/multiomeR/blob/main/configuration/cfg_GEM_wells.tsv) contains further inactive rows and annotation columns.

```{r, echo = FALSE, eval = TRUE, results = "asis"}
emit_GEM_well_demo_table(
  GEM_well_config_file = "website/data/demo_GEM_wells.tsv",
  dictionary_file = "website/data/GEM_well_columns.tsv"
)
```
