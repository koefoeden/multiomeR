# Plan your analysis

Before configuring anything, check that the inputs below are available. A
**GEM well** is one 10x Multiome library with its `cellranger-arc count`
output. An **aggregation** is a joint analysis of one or more GEM wells; its
settings are chosen step by step in [Run your own
analysis](main_running.md#steps).

## Per GEM well

- The `cellranger-arc count` output directory containing `outs/` with
  `summary.csv`, `filtered_feature_bc_matrix.h5`, `atac_fragments.tsv.gz`, its
  `.tbi` index, and `per_barcode_metrics.csv`.
- The Cell Ranger ARC reference must be one of the references known under
  `reference_metadata/`; it is identified automatically from the fragment
  file. All GEM wells analyzed together must share the same reference.
- For a well with several donors: a VCF with the donors' genotypes and the
  `atac_possorted_bam.bam` file from the same `outs/` directory.
- For ambient RNA removal: a CellBender H5 file produced outside multiomeR.

## Per aggregation

- A donor metadata TSV with one row per donor and the phenotypes or
  covariates you will use; see the [donor metadata
  table](reference_donor_metadata.md).
- Marker genes for the cell types you expect in the tissue, as gene symbols
  matching the reference annotation. Marker transcription factors for ATAC are
  optional.
- For the optional modules: [differential
  analyses](downstream_differential_analyses.md) require biological
  replication across donors, and [genetic
  enrichment](downstream_genetic_enrichment.md) requires a human aggregation
  and network access to Open Targets.

## Compute

- A Linux machine or scheduler meeting the [system
  requirements](demo_installation.md#system-requirements), with
  `crew_controllers.R` adjusted as described in [Choose where the analysis
  runs](performance_distributed_computing.md).
- Disk space for the store: the two-GEM-well demo writes about 6 GB.
- Time: the [recorded runtimes](performance_overview.md) for two and six GEM
  wells indicate what to expect; each step is rerun several times while
  settings are tuned.

## How the data move through the workflow

| Stage | What it does |
|---|---|
| GEM well processing | Reads Cell Ranger matrices and fragments, calculates QC metrics, and prefixes barcodes with the GEM well identifier. |
| GEX (gene expression) | Combines selected GEM wells, reduces dimensions, clusters nuclei, and annotates cell types using marker genes. |
| ATAC (chromatin accessibility) | Combines fragments, defines peaks, builds the peak-count matrix, and summarizes accessibility and motif signals. |
| WNN (weighted nearest neighbors) | Combines RNA and ATAC representations to produce integrated clusters, metadata, and a Seurat/Signac export. |

The [Main pipeline gallery](gallery_main.md) shows representative outputs;
the [implementation graph](implementation/implementation_main.html) provides
the detailed computational dependencies when you need them.

Continue to [Run your own analysis](main_running.md#steps).
