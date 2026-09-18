# GEM well table

```{r setup, include = FALSE}
pipeline_name <- "processing_and_aggregation"
source("helpers/_setup.R")
```

`cfg_GEM_wells.tsv` has one row per GEM well (one `cellranger-arc count` output). Copy an example row, give it a unique `GEM_well_ID`, and fill in the settings below. Aggregations select wells by these IDs. See [Run your own analysis](main_running.md#steps) for the first QC run.

## Fill in a row

| Setting | What to enter |
|---|---|
| `GEM_well_ID` | A unique identifier used in target names and output folders. |
| `GEM_well_dataset` | A label for the dataset or study. |
| `GEM_well_cellranger_arc_count_dir` | The directory containing `outs/`, not `outs/` itself. |
| `GEM_well_n_donors` | Number of donors in the well. |
| `GEM_well_donor_id` | For a single-donor well, an ID matching the [donor metadata table](reference_donor_metadata.md). |
| `GEM_well_donors_VCF_file` | For a multiplexed well, the donor-genotype VCF used for demultiplexing; otherwise `NA`. |
| `GEM_well_add_cellbender` | `TRUE` to use externally generated CellBender counts; otherwise `FALSE`. |
| `GEM_well_cellbender_h5_file` | Path to the CellBender H5 file when enabled; otherwise `NA`. |
| `GEM_well_QC_exclude_list` | `NA` for the first run; then exclusion expressions separated by `;;`. |
| `GEM_well_is_active` | `TRUE` for wells you want to process. Every well selected by an active aggregation must be active. |

Add library or batch annotations as extra columns prefixed with `GEM_well_`, such as `GEM_well_multiplex_batch`. These become cell metadata and can be used for batch correction or plots. Keep donor phenotypes in the donor table and avoid duplicate column names between the two tables, apart from their keys.

## Required inputs

``` text
<GEM_well_cellranger_arc_count_dir>/outs/
├── summary.csv
├── filtered_feature_bc_matrix.h5
├── atac_fragments.tsv.gz
├── atac_fragments.tsv.gz.tbi
└── per_barcode_metrics.csv
```

Genotype demultiplexing additionally requires `atac_possorted_bam.bam` in `outs/`. Prepare optional VCF and CellBender inputs as described in [Plan your analysis](main_overview.md).

The public pipeline identifies the Cell Ranger reference from the fragment-file header and matches it to a `reference.json` under `reference_metadata/`. Keep that header intact and add the corresponding JSON for a new reference. All wells in an aggregation must share the same reference. If your checkout requires an explicit `GEM_well_cellranger_arc_reference_json` column, follow its README.

## Set QC filters after the first run

Enter complete R expressions in `GEM_well_QC_exclude_list`, separated by `;;`:

``` {.text filename="cfg_GEM_wells.tsv"}
TSS.enrichment < 4 ;; nucleosome_signal > 4 ;; nCount_RNA < 250
```

Each expression identifies cells to exclude and is recorded as a separate exclusion reason. These cutoffs are examples, not recommendations for your tissue. Start with `NA`, inspect the first-step plots, and then choose filters. Set unused rows, including unused demo wells, to `GEM_well_is_active = FALSE` before running the whole pipeline.

## Column dictionary

Use the [committed example](https://github.com/koefoeden/multiomeR/blob/main/configuration/cfg_GEM_wells.tsv) to check the full set of columns.

The table below is a documentation snapshot of the two public demo wells, showing the core columns and one optional annotation. Bold columns must be present in the TSV; some allow an NA value. Scroll horizontally and focus or hover over a column's **i** button for its meaning. The other inactive rows and metadata columns in the public configuration remain available as examples.

```{r, echo = FALSE, eval = TRUE, results = "asis"}
emit_GEM_well_demo_table(
  GEM_well_config_file = "website/data/demo_GEM_wells.tsv",
  dictionary_file = "website/data/GEM_well_columns.tsv"
)
```
