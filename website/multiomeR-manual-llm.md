# multiomeR Website LLM Export

This file is generated from the Quarto book outlines and resolves Quarto include shortcodes.
Hidden setup chunks, generated helper chunks, Mermaid graph bodies, and verbose image metadata are omitted by default.


# Book: multiomeR Manual


## Part: Start here


<!-- source: website/index.md -->

[Image omitted; source: `figures/multiomeR-logo.svg`; alt: Image]

# Start here

## What is multiomeR?
multiomeR is a [targets](https://books.ropensci.org/targets/)-based workflow for processing and analyzing single-nucleus 10x Genomics Multiome ATAC + Gene Expression datasets. It is meant to be adapted to your own data, compute setup, and biological questions.

The workflow starts from `cellranger-arc count` outputs, processes gene-expression (GEX) and ATAC data, builds multimodal aggregations, and supports optional downstream modules for differential analyses and genetic enrichment for human datasets.

::: {.callout-warning title="Beta software"}
multiomeR is in beta and may introduce breaking changes between releases. The
steps in [Run your own analysis](main_running.md#steps)
show what to inspect before each stage, but acceptance criteria still depend on
the tissue, study design, and intended use. Report problems or questions
through [GitHub issues](https://github.com/koefoeden/multiomeR/issues).
:::

## Your first analysis

Start with the public demo: two human GEM wells with supplied configuration.
You will install the software, run one joint analysis, and read its cell
metadata and multimodal object. This gives you a working example before you
choose settings for your own study.

You need basic R skills, a Linux terminal, and a machine with sufficient
[memory and disk space](demo_installation.md#system-requirements). You do not
need to know how to write a `targets` pipeline. Commands labelled **Bash** run
in the terminal; commands labelled **R** run in the R session opened during
installation. Run both from the repository folder unless stated otherwise.

## Find what you need

| If you want to... | Start here |
|---|---|
| See what the workflow produces | Browse the [Main pipeline gallery](gallery_main.md). |
| Try multiomeR on public data | Follow [Install and prepare the demo](demo_installation.md), [Run the demo](demo_running.md), then [Inspect the demo results](demo_outputs.md). |
| Analyze your own data | Check the inputs in [Plan your analysis](main_overview.md), then follow the steps in [Run your own analysis](main_running.md). |
| Look up a configuration column or parameter | Open the [GEM well table](reference_GEM_wells.md), [Donor metadata table](reference_donor_metadata.md), or [Aggregation configuration](reference_aggregations.md) reference. |
| Understand a review plot or table | Look it up in [Output files and metadata](review_outputs.md). |
| Add a downstream analysis | Check the prerequisites for [Differential analyses](downstream_differential_analyses.md) or [Genetic enrichment](downstream_genetic_enrichment.md). |
| Run on a cluster or a smaller machine | Read [Choose where the analysis runs](performance_distributed_computing.md). |
| Fix a failed or stale run | Start with [Troubleshooting](troubleshooting.md). |
| Understand or modify the internals | Use the separate [implementation book](implementation/). |

## Terms used in this manual

In this manual, a **GEM well** is one configured 10x library and output
directory, an **aggregation** is a joint analysis of one or more GEM wells, and
a **donor** is the individual identified by `donor_id`. One GEM well may contain
multiple donors.

A **target** is a named result, such as a metadata table, matrix directory, or
plot. Its **dependencies** are the inputs and earlier results needed to build
it. You request the result you want; `targets` works out the order and reuses
results that are up to date. The **store** is the folder where it keeps results
and the records needed for reruns. For a small worked introduction, see the
[targets walkthrough](https://books.ropensci.org/targets/walkthrough.html).

## Workflow at a glance

The **main pipeline** processes each GEM well, aggregates selected GEM wells, and builds multimodal RNA/ATAC outputs for clustering, cell typing, and WNN integration. Two optional modules extend completed aggregations with differential analyses or genetic enrichment.

[Image omitted; source: `figures/multiomeR_overview_simplified.drawio.svg`; alt: multiomeR workflow from Cell Ranger ARC GEM well outputs through per GEM well processing, aggregation-level GEX and A...]

Continue to [Install and prepare the demo](demo_installation.md).


## Part: Output gallery


<!-- source: website/gallery_main.md -->

# Main pipeline gallery



These documentation snapshots show representative outputs from the public
`immune_human_2x` configuration with its two active GEM wells. The cards follow
the steps in [Run your own analysis](main_running.md#steps), where each step
also shows its own cards. Each card names the target that
generated the displayed demo result; click an image to open it at full
resolution.

To reproduce these plot families, follow [Install and prepare the
demo](demo_installation.md) and [Run the demo](demo_running.md), then run the
steps in [Run your own analysis](main_running.md#steps). The
demo command builds the final object and the WNN UMAPs only; [Inspect the demo
results](demo_outputs.md) explains the distinction.

[Generated Quarto chunk omitted: `render_gallery_section( gallery_items, "Main pipeline", subsection_descriptions = c( "Pre-aggregation QC" = "Unfilter...`]


<!-- source: website/gallery_differential_analyses.md -->

# Differential analyses gallery



These cards are a curated subset from the public `immune_human_2x` configuration. They illustrate diagnostics, not acceptable effect sizes or significance patterns for another study. See [Differential analyses](downstream_differential_analyses.md) for prerequisites, models, and the module run command.

The full module additionally produces expression-derived CollecTRI activity results and a CollecTRI-JASPAR concordance plot. They are not shown below until stable public example assets are available.

Named abundance models now produce donor-proportion plots and contrast plots with 95% Wald intervals. The previous coefficient and pooled-baseline-change examples have been retired; updated public abundance examples are pending.

[Generated Quarto chunk omitted: `render_gallery_section( gallery_items, "Differential analyses module", subsection_descriptions = c( "Gene expression"...`]


<!-- source: website/gallery_genetic_enrichment.md -->

# Genetic enrichment gallery



These curated outputs use the larger `PBMC_human_6x` aggregation, not the quickstart with two GEM wells, and show selected SCAVENGE/WNN results rather than every attribution output. See [Genetic enrichment](downstream_genetic_enrichment.md) for prerequisites, module target selection, and interpretation guidance.

[Generated Quarto chunk omitted: `render_gallery_section( gallery_items, "Genetic enrichment module", subsection_descriptions = c( "Single-nucleus chro...`]


## Part: Try the public demo


<!-- source: website/demo_installation.md -->

# Install and prepare the demo



## System requirements

These instructions use the public repository and its demo configuration.
An institutional checkout may supply different input paths, a different
output folder, and cluster controllers. Use its local setup instructions
before running the demo commands.

- Linux with `git` and `curl`, plus HTTPS access to GitHub, Pixi, and 10x Genomics downloads.
- At least 60 GB of RAM. This is enough for one heavy target at a time; machines near the minimum should reduce concurrent workers in `crew_controllers.R`.
- At least 30 GB of free disk space for the public inputs, pixi environment, temporary files, and approximately 6 GB of demo outputs.
- Multiple CPU cores are strongly recommended. The timing quoted in the next chapter was measured with 16 logical threads.

::: {.callout-tip title="Machines with less than 256 GB of RAM"}
The committed `crew_controllers.R` is sized for a 16-CPU, 256-GB workstation:
four light workers and two heavy workers that may each use 60 GB. On a machine
near the 60-GB minimum, edit the two `workers` values in `crew_controllers.R`
before running the demo so that only one heavy target runs at a time:

```{.r filename="crew_controllers.R"}
controller_list <- list(
  crew::crew_controller_local(
    name = "local-light",
    workers = 2
  ),
  crew::crew_controller_local(
    name = "local-heavy",
    workers = 1,
    crashes_max = 1
  )
)
```

The `RAM_GB` values in the same file describe routing capacity, not enforced
limits, so the workers that can run at once must fit in physical memory. To
run on a SLURM or other scheduler instead, see [Choose where the analysis
runs](performance_distributed_computing.md).
:::


## Set up the demo

Run this block from the directory where you want to clone multiomeR. The single
`pixi run` setup command installs the locked environment before its
`setup-demo` task downloads the two configured public inputs and installs the
pinned GitHub-only R packages.

```{.bash filename="Bash"}
# Clone the repository and enter its root directory.
git clone https://github.com/koefoeden/multiomeR.git
cd multiomeR

# Skip these two lines when pixi is already available on PATH.
curl -fsSL https://pixi.sh/install.sh | sh
export PATH="$HOME/.pixi/bin:$PATH"

# Install the locked environment, download 3.9 GB of demo inputs, and install
# the pinned GitHub versions of BPCells, Signac, and betterChromVAR.
pixi run --use-environment-activation-cache --locked --run-post-link-scripts setup-demo

# Start R in the configured environment for the commands in the next chapter.
pixi run --use-environment-activation-cache --locked R
```

The download task is restart-safe: non-empty files already present under
`example_data` are skipped. The repository includes the small `reference.json`
from the exact `refdata-cellranger-arc-GRCh38-2020-A-2.0.0` reference used for
both public outputs, so the full Cell Ranger ARC reference is not required.

Continue to [Run the demo](demo_running.md) from the R prompt.


<!-- source: website/demo_running.md -->

# Run the demo



In the R session opened during installation, run the command below to process
`immune_human_2x`. It combines the two demo GEM wells, produces a
Seurat/Signac object containing the multimodal results, and draws the
integrated WNN UMAPs coloured by cluster, cell type, and the other categorical
metadata.

`names` selects these two results by their exact names using `all_of()`.
`tar_make()` also builds the dependencies needed for them, but does
not build every plot in the gallery.

```{.r filename="R"}
demo_results <- c(
  "multimodal_Seurat_object.8_multimodal_QC.immune_human_2x",
  "categorical.UMAPs.8_multimodal_QC.immune_human_2x"
)

targets::tar_make(names = tidyselect::all_of(demo_results))
```

Keep the R session open until the command finishes. Progress messages report
targets being dispatched, completed, or skipped because they are already
up to date. A previous demo run on an AMD EPYC 7543 system with 16 logical
threads took under 25 minutes and wrote about 6 GB; your runtime may differ.

## Confirm success

After the run, the following command should return `character(0)`, meaning the requested results and their dependencies are up to date:

```{.r filename="R"}
targets::tar_outdated(
  names = tidyselect::all_of(demo_results),
  callr_function = NULL
)
```

If names are returned, those results still need building. If the run failed,
follow [Troubleshooting](troubleshooting.md), fix the reported cause, and run
the same `tar_make()` command again. Completed results can be reused.

Continue to [Inspect the demo results](demo_outputs.md) to read the object
and find the associated files.

## Local validation

The full local validation command also enables the inactive mouse and public
ENCODE examples, checks external resources, and verifies analysis outputs:

```sh
pixi run --use-environment-activation-cache validate-local
```

It reuses the existing targets store. The extra inputs and validation scope are
described in the repository's `validation/README.md`; the normal demo remains
limited to `immune_human_2x`.


<!-- source: website/demo_outputs.md -->

# Inspect the demo results



The public demo uses `outputs/` as its results store. Other checkouts may
use a different folder, recorded under `store` in `_targets.yaml`.
The R commands below use that configuration automatically.

The store contains serialized R objects in `objects/`, file artifacts in `files/`, and requested review figures in `plots/`.

## Objects

Most intermediate and final result objects are saved automatically by `targets` under `outputs/objects` during a pipeline run. Read them with `targets::tar_read()` from a repository-root R session after the demo has completed.

```{.r filename="R"}
cell_metadata <- targets::tar_read(
  metadata_w_cell_types_tibble.WNN.immune_human_2x
)
dim(cell_metadata)
head(cell_metadata)

demo_object <- targets::tar_read(multimodal_Seurat_object.8_multimodal_QC.immune_human_2x)
demo_object
```

The metadata table describes the retained nuclei and their annotations.
WNN means *weighted nearest neighbors*: the integrated representation uses
information from both RNA and ATAC. The Seurat/Signac object is a convenient
export for further exploration; the pipeline also retains its matrices in
BPCells format on disk.

## Files

File targets also load with `targets::tar_read()`, but their value is a path rather than an in-memory result.

```{.r filename="R"}
targets::tar_read(cellranger_barcodes_tsv.healthy_PBMC_human)
targets::tar_read(aggregated_GEX_BPCells_matrix_dir.GEX.immune_human_2x)
targets::tar_read(consensus_peak_BPCells_matrix_dir.ATAC.immune_human_2x)
```

For example, `aggregated_GEX_BPCells_matrix_dir.GEX.immune_human_2x` is placed under:

```text
outputs/files/immune_human_2x/GEX/
```

Other files are grouped under `outputs/files/<scope>/`, where the scope is a
GEM well, an internal pre-aggregation QC group, or an aggregation.

## Plots

The demo command also built `categorical.UMAPs.8_multimodal_QC.immune_human_2x`.
Plots use the same scope-based layout under `outputs/plots/`, so its files are
in:

```text
outputs/plots/immune_human_2x/8_multimodal_QC/UMAPs/categorical/
```

Open `WNN_harmony_SNN_cluster_cell_type.png` there to see the integrated
clusters and cell-type labels from your own run. It should resemble this
documentation snapshot:

[Generated Quarto chunk omitted: `render_gallery_grid(gallery_items[gallery_items$id == "wnn-umap", ])`]

The [Main pipeline gallery](gallery_main.md) shows the other plot families
the workflow produces. Those are saved snapshots and do not reflect the state of
your analysis. To build one more of them, follow [Request an additional
result](main_running.md#after-the-steps); to build all review
plots for a stage, run its step in [Run your own
analysis](main_running.md#steps).

## Next steps

- To adopt the workflow, continue with [Plan your analysis](main_overview.md) and the steps in [Run your own analysis](main_running.md#steps).
- To request more results or rerun after a change, use [Run your own analysis](main_running.md).
- To diagnose a failed or unexpectedly stale target, use [Troubleshooting](troubleshooting.md).


## Part: Analyze your own data


<!-- source: website/main_overview.md -->

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


<!-- source: website/main_running.md -->

# Run your own analysis

Start from the working demo and replace its inputs and settings one step at a
time. Each step has three parts: **Configure** the relevant settings, **Run**
the step, and **Review** its plots. Every plot carries a subtitle that states
how to read it and what to look for, so this chapter only lists which plots to
open. The trees show the main review outputs; configured variables and optional
analyses determine the additional files. Revise the settings and rerun a step until its plots are acceptable, then
continue to the next.

Conventions used below:

- Commands run in a repository-root R session, opened as in [Install and
  prepare the demo](demo_installation.md).
- Replace `my_GEM_well` and `my_aggregation` with your identifiers.
- `<store>` is the results folder configured in `_targets.yaml`, `outputs/` in
  the demo.
- File formats are documented in the [GEM well table](reference_GEM_wells.md),
  [donor metadata table](reference_donor_metadata.md), and [aggregation
  configuration](reference_aggregations.md) references.
- `crew_controllers.R` must describe your machine or scheduler; see [Choose
  where the analysis runs](performance_distributed_computing.md).

Settings default to `configuration/`. To use your own directory, copy it and
write the new directory path into the ignored root `configuration.local` file.
Keep all settings together; disabled modules may omit their files. See the
[configuration guide](https://github.com/koefoeden/multiomeR/blob/main/configuration/README.md).

These instructions describe the public defaults. For installation-specific
setup and reference selection, follow the README in your checkout.

## How a step is run {#how-to-run}

Target names and plot folders contain the name of the step's **checkpoint**.
For example, `harmony.categorical.UMAPs.3_GEX_QC.my_aggregation` is saved
under `<store>/plots/my_aggregation/3_GEX_QC/UMAPs/categorical/harmony/`.
Selecting by checkpoint name and scope suffix builds the step's plots; the
tables they depend on are built automatically:

```{.r filename="R"}
targets::tar_make(
  names = tidyselect::contains("3_GEX_QC") & tidyselect::ends_with(".my_aggregation")
)
```

To preview a selection, pass the same `names` to `targets::tar_manifest()`
with `callr_function = NULL`. An empty result means a wrong suffix. After
editing configuration, `targets::tar_manifest(callr_function = NULL)` must
build without errors; see [Troubleshooting](troubleshooting.md) if it does not.

In the trees below, a trailing `/` marks a folder with one plot per variable,
metric, or component. How the plots and the tables behind them are constructed
is documented in [Output files and metadata](review_outputs.md).

## 1. Process the GEM wells {#steps}

Reads each GEM well's Cell Ranger output and reviews per-well exclusions.
Compare distributions across wells after defining an aggregation in step 2.

**Configure** one row per GEM well in `cfg_GEM_wells.tsv`
([GEM well table](reference_GEM_wells.md)):

- a unique `GEM_well_ID` and an appropriate `GEM_well_dataset` label;
- `GEM_well_cellranger_arc_count_dir`;
- donors: `GEM_well_n_donors` and `GEM_well_donor_id` for a single-donor
  well, or `GEM_well_donors_VCF_file` to demultiplex several donors by
  genotype;
- `GEM_well_add_cellbender` and `GEM_well_cellbender_h5_file` to replace the
  Cell Ranger counts with CellBender output;
- further `GEM_well_` columns for library-level variables you will need as
  batch covariates or plot variables later;
- `GEM_well_QC_exclude_list`: `NA` for the first run, then filter expressions
  such as `TSS.enrichment < 4 ;; nCount_RNA < 250`;
- `GEM_well_is_active`: `TRUE`.

**Run**

```{.r filename="R"}
targets::tar_make(
  names = tidyselect::contains("1_pre_aggregation_QC") & tidyselect::ends_with(".my_GEM_well")
)
```

**Review**

```{.text}
<store>/plots/my_GEM_well/1_pre_aggregation_QC/
├── excluded_barcodes_by_type_upset.png
└── excluded_cellranger_only_barcodes_by_type_upset.png
```

## 2. Aggregate the GEM wells

Applies the per-GEM-well filters, combines the selected wells, and checks that
they share one reference and identical gene definitions.

**Configure**

- a [donor metadata table](reference_donor_metadata.md) with one row per
  `donor_id` holding the donor-level phenotypes and covariates;
- an entry in `cfg_aggregations.yaml` ([aggregation
  configuration](reference_aggregations.md)) with:
  - `aggregation_GEM_well_IDs`: the wells to analyze together. The group
    should reflect the biological comparison you intend to make, keep donor,
    preparation, and batch distinguishable, and share one reference;
  - `aggregation_donor_id_metadata_tsv`;
  - `aggregation_GEX_marker_genes`: valid positive marker genes for the expected cell types, refined in step 4;
  - `is_active: true`.

**Run**

```{.r filename="R"}
targets::tar_make(
  names = tidyselect::contains("1_pre_aggregation_QC") & tidyselect::ends_with(".my_aggregation")
)
```

**Review**

```{.text}
<store>/plots/my_aggregation/1_pre_aggregation_QC/
├── per_aggregation_GEM_well_QC_comparisons/
├── aggregation_excluded_cellranger_only_barcodes_by_type_upset.png
├── aggregation_excluded_barcodes_by_type_upset.png
├── nuclei_per_donor_id_bars.png
└── cell_retention_flow_plot.png
```

## 3. Build the GEX embedding

Normalizes the GEX counts, selects variable genes, and computes PCA with
optional Harmony correction.

**Configure**

- `aggregation_GEX_data_PCs`: candidate PCA dimensions;
- `aggregation_harmony_correction_metadata_col_names`: metadata columns to
  correct with Harmony, such as `GEM_well_ID` or a batch column. Omit it for
  the first run and add it only if this step's association plots show batch
  structure;
- normalization and variable-gene settings under the GEX topic of the
  [parameter reference](reference_aggregations.md#parameter-reference).

**Run**

```{.r filename="R"}
targets::tar_make(
  names = tidyselect::contains("2_GEX_PCA_QC") & tidyselect::ends_with(".my_aggregation")
)
```

**Review**

```{.text}
<store>/plots/my_aggregation/2_GEX_PCA_QC/
├── variable_feature_plot.png
├── VizDimLoadings_plots/
├── PCA_singular_values_elbow_plot.png
├── PCA_embedding_sdev_plot.png
└── PCA_metadata_association_barplots/
```

## 4. Cluster and label cell types with GEX

Builds the neighbour graph, clusters, annotates cell types from the marker
genes, and scores doublets. The accepted cells and labels guide peak calling.

**Configure**

- `aggregation_GEX_marker_genes`: marker genes per expected cell type, as
  symbols in the reference annotation; a `-` suffix marks a gene that should
  be absent;
- neighbour and clustering-resolution settings;
- `aggregation_scDblFinder_GEX_remove_called_doublets: false` and
  `aggregation_scDblFinder_GEX_max_doublet_fraction_per_cluster: null` for the
  first run, then the chosen doublet policy;
- categorical and continuous metadata variables for the review plots.

**Run**

```{.r filename="R"}
targets::tar_make(
  names = tidyselect::contains("3_GEX_QC") & tidyselect::ends_with(".my_aggregation")
)
```

**Review**

```{.text}
<store>/plots/my_aggregation/3_GEX_QC/
├── UMAPs/
│   ├── categorical/{harmony,non_harmony}/
│   ├── continuous/{harmony,non_harmony}/
│   └── cross/
├── categorical_by_cell_type_bars_plots/
├── categorical_by_cluster_bars_plots/
├── continuous_by_cluster_violin_plot/
├── continuous_by_cell_type_violin_plot/
├── markers_by_cluster_dot_plot.png
├── markers_by_cell_type_dot_plot.png
├── module_scores_by_cluster_dot_plot.png
├── module_scores_by_cell_type_dot_plot.png
├── cluster_UCell_advantage_plots/
├── cluster_marker_volcano_plots/
├── cell_type_marker_volcano_plots.png
└── cell_retention_flow_plot.png
```

## 5. Call peaks and inspect ATAC quality

Calls peaks per accepted GEX group, builds the consensus peak matrix, and
computes peak-based QC metrics before any ATAC filter is applied.

**Configure**

- `aggregation_QC_exclude_list_combined_object`: omit for the first run, then
  filter expressions such as:

```{.yaml filename="cfg_aggregations.yaml"}
aggregation_QC_exclude_list_combined_object:
  - nCount_ATAC < 1000
  - atac_peak_counts_frac < 0.1
  - atac_peak_counts_blacklist_frac > 0.01
```

**Run**

```{.r filename="R"}
targets::tar_make(
  names = tidyselect::contains("4_peak_QC") & tidyselect::ends_with(".my_aggregation")
)
```

**Review**

```{.text}
<store>/plots/my_aggregation/4_peak_QC/
├── peaks_QC_violins_plot/
└── peaks_similarity_tiles_plot.png
```

## 6. Filter nuclei on ATAC quality

Applies the expressions from step 5 and shows what they removed.

**Configure** nothing new.

**Run**

```{.r filename="R"}
targets::tar_make(
  names = tidyselect::contains("5_pre_LSI_QC") & tidyselect::ends_with(".my_aggregation")
)
```

**Review**

```{.text}
<store>/plots/my_aggregation/5_pre_LSI_QC/
├── QC_excluded_upset_plot.png
└── cell_retention_flow_plot.png
```

## 7. Build the ATAC embedding

Computes LSI on the retained nuclei with optional Harmony correction.

**Configure**

- `aggregation_ATAC_data_PCs`: candidate LSI components;
- `aggregation_extra_harmony_covars_ATAC`: extra columns to correct in ATAC
  only, added to the GEX Harmony columns; add them only if this step's
  association plots show batch structure.

**Run**

```{.r filename="R"}
targets::tar_make(
  names = tidyselect::contains("6_ATAC_LSI_QC") & tidyselect::ends_with(".my_aggregation")
)
```

**Review**

```{.text}
<store>/plots/my_aggregation/6_ATAC_LSI_QC/
├── LSI_singular_values_elbow_plot.png
├── LSI_embedding_sdev_plot.png
├── VizDimLoadings_plots/
└── LSI_metadata_association_barplots/
```

## 8. Cluster and label cell types with ATAC

Clusters the ATAC embedding, transfers cell-type labels from the GEX marker
scores, scores doublets, and summarizes motif-family accessibility and gene
activity.

**Configure**

- neighbour and clustering-resolution settings;
- `aggregation_ATAC_marker_TFs`: marker transcription factors per cell type;
- `aggregation_scDblFinder_ATAC_remove_called_doublets: false` and
  `aggregation_scDblFinder_ATAC_max_doublet_fraction_per_cluster: null` for the
  first run, then the chosen doublet policy.

**Run**

```{.r filename="R"}
targets::tar_make(
  names = tidyselect::contains("7_ATAC_QC") & tidyselect::ends_with(".my_aggregation")
)
```

**Review**

```{.text}
<store>/plots/my_aggregation/7_ATAC_QC/
├── UMAPs/{categorical,continuous,cross}/
├── categorical_by_cell_type_bars_plots/
├── categorical_by_cluster_bars_plots/
├── marker_gene_activity_dot_plot.png
├── motif_family_accessibility_by_ATAC_cluster_heatmap.png
├── motif_family_accessibility_by_GEX_cluster_heatmap.png
├── motif_family_accessibility_by_GEX_cell_type_heatmap.png
├── motif_family_accessibility_marker_volcano_plots.png
├── coverage_tracks_plots/
├── cluster_UCell_advantage_plots/
├── confusion_matrices_plots.png
└── cell_retention_flow_plot.png
```

## 9. Integrate GEX and ATAC

Combines the accepted GEX and ATAC representations with WNN, clusters and
labels the integrated embedding, and exports the Seurat/Signac object.

**Configure** the WNN neighbour, resolution, and UMAP settings.

**Run**

```{.r filename="R"}
targets::tar_make(
  names = tidyselect::contains("8_multimodal_QC") & tidyselect::ends_with(".my_aggregation")
)
```

**Review**

```{.text}
<store>/plots/my_aggregation/8_multimodal_QC/
├── UMAPs/{categorical,continuous,cross}/
├── categorical_by_cell_type_bars_plots/
├── categorical_by_cluster_bars_plots/
├── UMAPs/cluster_named_dim_tri_plot.png
├── UMAPs/cluster_cell_type_dim_tri_plot.png
├── markers_by_cluster_dot_plot.png
├── module_scores_by_cluster_dot_plot.png
├── continuous_by_cell_type_violin_plot/
├── continuous_by_cluster_violin_plot/
├── WNN_weight_metadata_associations_plot.png
├── confusion_matrices_plots/
├── cluster_UCell_advantage_plots/
└── cell_retention_flow_plot.png
```

The final object is `multimodal_Seurat_object.8_multimodal_QC.my_aggregation`;
[Inspect the demo results](demo_outputs.md) shows how to read it. Continue to
[Differential analyses](downstream_differential_analyses.md) or [Genetic
enrichment](downstream_genetic_enrichment.md) once the aggregation is
accepted.

## After the steps

**Request one result.** Use its exact name; `all_of()` errors if the name is
absent.

```{.r filename="R"}
targets::tar_make(
  names = tidyselect::all_of("categorical.UMAPs.8_multimodal_QC.my_aggregation")
)
```

**Rerun after a change.** Reuse the step's selection. `targets::tar_outdated()`
with the same `names` and `callr_function = NULL` lists what will rebuild.

**Build everything active.** `targets::tar_make()` without `names` builds
every active GEM well, aggregation, and enabled module. Set `GEM_well_is_active`
and `is_active` to `FALSE` for entries you do not want built, including unused
demo entries, before running it.


## Part: Add an optional analysis


<!-- source: website/downstream_differential_analyses.md -->

# Differential analyses



## When to use this module

Use this module to ask how cell-type proportions, gene expression, or chromatin accessibility differ with a condition or donor phenotype. It requires biological replication across donors or samples. Molecular measurements are combined into **pseudobulks**: counts summarized for each cell type within a donor or sample.

The module does not create biological replication. The donor structure, covariates, design formula, and contrasts must be defensible for the intended analysis before the workflow is run.

See the [Differential analyses gallery](gallery_differential_analyses.md) for representative diagnostics and the [implementation graph](implementation/implementation_differential_analyses.html) for target structure.

## Prerequisites

Before enabling the module, confirm that:

- you have reviewed the main aggregation and its cell-type annotations;
- WNN cell-type metadata and GEX and ATAC pseudobulk matrices are available;
- the donor metadata contains one unique row per `donor_id` and every variable used in a model;
- model variables are donor- or pseudobulk-sample-level variables, not duplicated cell-level measurements; and
- the number and distribution of donors support the specified design and contrasts.

Use `differential_analyses_extended_donor_id_metadata_tsv` when the modelling table needs variables beyond the aggregation's normal donor metadata. It must retain the same unique `donor_id` key.

## Outputs

Choose the output that matches your question:

| Question | Output family |
|---|---|
| Do cell-type proportions differ? | `cell_type_composition` |
| Which genes change expression? | `gene_expression` |
| Which peaks change accessibility? | `chromatin_accessibility` |
| Which motif families change accessibility? | `motif_family_accessibility` (JASPAR) |
| Which regulators show altered expression-based activity? | `transcription_factor_activity` (CollecTRI) |

The module also produces model diagnostics, comparisons across modalities,
and gene-set tests for Hallmark and Reactome pathways. Motif-family accessibility
summarizes ATAC evidence; transcription-factor activity is inferred from gene
expression using CollecTRI. Interpret each in the context of its measurement.

Plot directories use these descriptive family names below
`plots/<aggregation>/differential_analyses/`. Gene-set plots appear under
`gene_expression/gene_set_enrichment/Hallmark/enrichment_plots/<model>/`
or the corresponding `Reactome` directory. Volcano outputs use
`<family>/volcano_plots/<model>/`; saved plot targets omit redundant `_file`
and `_files` suffixes. Renaming targets creates new cache entries and output
paths on the next run; existing output directories are not migrated.

See the [method details](implementation/implementation_differential_analyses.html#method-details)
for activity inference, motif-family definitions, and gene-set testing.

## Configure

Add `modules` to the existing aggregation entry, keeping its input and marker settings:

```{.yaml filename="cfg_aggregations.yaml"}
your_aggregation:
  modules: [differential_analyses]
```

Then create a matching row directly in
`configuration/cfg_module_differential_analyses.yaml`.

```{.yaml filename="configuration/cfg_module_differential_analyses.yaml"}
your_aggregation:
  differential_analyses_cell_type_composition_models:
    condition_abundance:
      formula: ~ condition
      contrast_specs_vec:
        treated_vs_control: conditiontreated
      plot_phenotype_vars: condition
      color_by: condition
  differential_analyses_pseudobulk_models:
    condition_model:
      cell_type_subset: NULL
      design_matrix_func_name: NULL
      formula: ~ 0 + cluster + condition
      random_effect: NULL
      contrast_specs_vec:
        treated_vs_control: conditiontreated
```

Both branches use named models, donor eligibility checks and named contrasts. Abundance models use `differential_analyses_cell_type_composition_models`; feature models use `differential_analyses_pseudobulk_models`. Omitting abundance models disables that branch. The former aggregation-wide cell-type composition formula, phenotype and colour settings have been replaced by fields inside each named model.

For mixed tissues, set `GEM_well_IDs` inside an abundance model to define its population, for example the six left-ventricle wells. Optional `donor_ids` can further restrict donors in either branch. Donors with missing model metadata or no selected samples are excluded and recorded in model-specific cohort TSVs. Feature cohorts also report retained pseudobulk sample counts and depth-filter exclusions.

By default, abundance models test all cell-type labels observed in their eligible population, including unassigned labels. Optional `cell_types_to_test` restricts the response cell types **without changing the denominator**: every retained nucleus in the selected wells contributes to its donor's total. Zero donor–cell-type counts remain in the analysis. In contrast, feature-model `cell_type_subset` selects the cells represented by the pseudobulks. Feature matrices already pool wells within donors and cell types, so they cannot support a late `GEM_well_IDs` filter; the module rejects that field for feature models.

Abundance models fit a separate fixed-effects beta-binomial logit model per cell type. Use a one-sided predictor formula, `formula: ~ ...`, and named `contrast_specs_vec`. The implementation supplies the fixed response `cbind(n_nuclei, n_other_nuclei)`; two-sided formulas are rejected. Custom design/contrast functions and random-effects formulas are not supported in this branch. Contrast tables report log-odds effects, Wald uncertainty, donor counts and BH FDR across tested cell types within each model/contrast. Failed fits are explicitly marked non-estimable. Counts, cohort tables and plots use the same eligible donors and denominators. Plot-only variables do not exclude donors from the fit.

The module selection below requests both configured abundance and pseudobulk outputs. Formula terms and contrast coefficients must match columns produced by the model matrix. The two branches retain their distinct response construction and fitting methods; sharing configuration does not make their effect estimates interchangeable.

The model example assumes `condition` distinguishes treated and control donors.
Check which group is the reference and what each model coefficient represents
before using `conditiontreated` as a contrast. Replace the example formula
and contrast to match your study.

## Run

Preview the selected module outputs before running them:

```{.r filename="R"}
targets::tar_manifest(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:differential_analyses")
  ) & tidyselect::ends_with(".your_aggregation"),
  callr_function = NULL
)[, c("name", "description")]
```

Then run the same selection:

```{.r filename="R"}
targets::tar_make(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:differential_analyses")
  ) & tidyselect::ends_with(".your_aggregation")
)
```

## Review

Open the configured model outputs listed above and the [differential gallery](gallery_differential_analyses.md). Interpretation and method details are included in the plot subtitles and captions.

Runtime depends on donors, cell types, models, contrasts, and gene-set analyses. Use [Troubleshooting](troubleshooting.md) if a formula, contrast, or metadata join fails.

## Parameter reference

The OLINK and bulk-RNA path fields are reserved optional integration inputs and are not consumed by the current public differential-analysis selection. Leave them `NULL` unless the corresponding integration is implemented in your downstream workflow.

[Generated Quarto chunk omitted: `emit_parameter_overview("differential_analyses")`]

The public demos leave this module disabled. Comparing one healthy PBMC donor
with one lymphoma lymph-node donor cannot separate condition, donor, and tissue
effects. Configure differential analyses for a design with biological replication.


<!-- source: website/downstream_genetic_enrichment.md -->

# Genetic enrichment



## When to use this module

Use this module to ask which cell types or nuclei have accessible regions
overlapping genetic evidence for a human trait. It connects fine-mapped GWAS
variants to ATAC peaks, calculates accessibility-based enrichment, and uses
[`SCAVENGE`](https://github.com/sankaranlab/SCAVENGE) to summarize trait
relevance across related nuclei.

A **credible set** contains candidate causal variants at a GWAS locus, with
probabilities from fine-mapping. Enrichment helps prioritize cellular contexts;
it does not by itself identify a causal cell type, gene, or mechanism.

See the [Genetic enrichment gallery](gallery_genetic_enrichment.md) for representative results and the [implementation graph](implementation/implementation_genetic_enrichment.html) for upstream ATAC and WNN dependencies.

## Prerequisites

Before enabling the module, confirm that:

- the aggregation is human and you have reviewed its main results and cell-type annotations;
- WNN metadata and graph results, consensus peaks, ATAC counts, chromVAR objects, and GEX/ATAC embeddings are available;
- each configured `sourceId` represents the intended trait and population;
- the machine can download and retain the Open Targets study and credible-set Parquet datasets; and
- you have chosen whether to interpret individual nuclei, graph-smoothed scores, or cell-type summaries; these answer related but different questions.

## Outputs

| Result | What to inspect |
|---|---|
| Study and variant-to-peak tables | Which studies, variants, and accessible regions contributed |
| Single-nucleus deviations | Accessibility enrichment for each nucleus |
| SCAVENGE plots | Trait-relevance scores propagated through the cell-neighbor graph |
| Cell-type heatmaps and attribution tables | Enrichment by cell type and the loci or variants contributing to it |

## Configure

Add `modules` to the existing human aggregation entry, keeping its other settings:

```{.yaml filename="cfg_aggregations.yaml"}
your_aggregation:
  modules: [genetic_enrichment]
```

Then create a matching row directly in `configuration/cfg_module_genetic_enrichment.yaml`.

```{.yaml filename="configuration/cfg_module_genetic_enrichment.yaml"}
your_aggregation:
  genetic_enrichment_GWAS_studies:
    lymphocyte_count:
      Category: positive_control
      sourceId: GCST90002388
      finemappingMethod: auto
```

`sourceId` values beginning with `GCST` use the pinned Open Targets datasets.
Every other value is a local Parquet filename, resolved from the project root
and tracked as a file target. Local files must satisfy the schema enforced by
`validate_local_finemapped_GWAS_tibble()`; their study ID, fine-mapping method,
build, credible-set probability, and provenance are read from the file rather
than repeated in YAML.

The root workflow currently pins Open Targets release `26.03`. That release identifier is recorded in downstream metadata and determines the available studies, credible sets, and fine-mapping methods.

For `finemappingMethod: auto`, multiomeR selects the first available supported method in this order: `SuSie`, `SuSiE-inf`, then `PICS`. Specify a method explicitly when the method itself is part of the analysis contract; the workflow fails if that method is unavailable for the study.

## Run

Preview the selected module outputs before running them:

```{.r filename="R"}
targets::tar_manifest(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:genetic_enrichment")
  ) & tidyselect::ends_with(".your_aggregation"),
  callr_function = NULL
)[, c("name", "description")]
```

Then run the same selection:

```{.r filename="R"}
targets::tar_make(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:genetic_enrichment")
  ) & tidyselect::ends_with(".your_aggregation")
)
```

## Review

Open the study-selection summaries, chromVAR summaries, SCAVENGE heatmaps and locus-contribution plots produced for your configured studies. See the [genetic enrichment gallery](gallery_genetic_enrichment.md) for available previews; interpretation belongs to each plot.

Runtime and disk use grow with studies, cells, graph representations, permutations, and attributed loci. The [Genetic enrichment gallery](gallery_genetic_enrichment.md) uses a larger aggregation with six GEM wells and is not produced by the minimal quickstart.

## Parameter reference

[Generated Quarto chunk omitted: `emit_parameter_overview("genetic_enrichment")`]

<details>
<summary>Show the public <code>immune_human_2x</code> example</summary>

[Generated Quarto chunk omitted: `emit_yaml_entry(module_config_file, "immune_human_2x")`]

</details>


<!-- source: website/downstream_peak_gene_correlation.md -->

# Peak–gene correlation

Run this optional module after accepting the final WNN cell set. It relates
ATAC accessibility to RNA expression within broad GEX-derived cell types,
using the retained WNN nuclei; it is no longer part of checkpoint 8.

## Configure

Add `peak_gene_correlation` to the aggregation's existing `modules` list in
`cfg_aggregations.yaml`, and add a matching row in
`cfg_module_peak_gene_correlation.yaml` in the selected configuration directory:

```yaml
my_aggregation:
  peak_gene_correlation_top_links_per_cell_group: 3
```

The top-link count controls the number of detail figures per cell type.
`peak_gene_correlation_filter` selects `lenient` (default), `moderate`, or
`strict` measurement-support filtering. Disabled aggregations contribute no
module targets.

## Run

```r
targets::tar_make(
  names = tidyselect::ends_with(".peak_gene_correlation.my_aggregation")
)
```

## Review

All module targets have description tag `[checkpoint:peak_gene_correlation]`.
Their paths are `<store>/plots/my_aggregation/peak_gene_correlation/`;
file exports use the corresponding `files` directory, with any modality
suffixes as deeper subdirectories. For example, read selected links with:

```r
targets::tar_read(
  peak_gene_correlation_links_tibble.WNN.peak_gene_correlation.my_aggregation
)
```

## Model and scope

Peak–gene links are candidate regulatory relationships. Cells are aggregated
within donor and ATAC-defined state using WNN cell-type annotations, without
reusing a cell across aggregates. Measurement-support filters select hypotheses
before fitting; the default requires shared support from at least two donors,
and the strict preset requires three. Review `filter_retention_plot` first.

The hierarchical analysis fits a mean peak effect with donor-specific slope
variation, donor intercepts and RNA/ATAC depth adjustment. It reports
Kenward–Roger p-values, BH FDR and numerical reliability diagnostics. Top-link
figures rank positive, estimable nonpromoter slopes by p-value without a
significance cutoff, so appearing in a figure is not evidence of significance.
They combine focal-cell-type coverage, gene context, hierarchical evidence and
adjusted aggregate scatterplots.

The existing HC3 correlation summaries remain a separate conditional analysis.
Neither analysis establishes causal regulation; numerical reference parity does
not establish statistical calibration across datasets. See the
[module methods and filtering reference](https://github.com/koefoeden/multiomeR/blob/main/module_peak_gene_correlation/README.md)
for inference limits, support thresholds and output details.


## Part: Operation and scaling


<!-- source: website/performance_overview.md -->

# Performance and scaling

multiomeR gains performance from two complementary mechanisms: `{targets}` exposes independent GEM wells, modalities, and analysis branches for concurrent execution, while BPCells keeps large matrices on disk and streams many operations instead of materializing dense objects.

Actual wall time depends on input nuclei, retained cells, peak counts, enabled plots and modules, hardware, scheduler latency, controller concurrency, and which targets are already up to date. Treat benchmark results as workload descriptions, not promises for another machine or configuration.

## Recorded examples

The recorded benchmark below contains two public aggregations and estimates the critical path to each final `multimodal_Seurat_object.8_multimodal_QC` from recorded `{targets}` runtime metadata.

| Aggregation | GEM wells | Cell Ranger input nuclei | Estimated critical path | Sum if ancestor targets ran serially |
|---|---:|---:|---:|---:|
| `immune_human_2x` | 2 | 17,277 | 17.3 minutes | 31.4 minutes |
| `PBMC_human_6x` | 6 | 51,291 | 23.8 minutes | 52.4 minutes |

The critical path is the longest chain of dependent tasks: it estimates the fastest completion if enough workers are available. This estimate treats dynamic branches as concurrently runnable. The difference between that estimate and the serial sum illustrates available parallelism for these two runs; it does not establish a general scaling law from two observations.

## How to use these numbers

- Compare changes only when the pipeline version, configuration, target endpoint, and estimation method are recorded together.
- Expect scheduler queue time and limited worker capacity to increase observed wall time beyond the graph-only critical-path estimate.
- Inspect target-level runtime and memory on your own representative aggregation before sizing a production run.
- Start with one aggregation and a selected result, then increase worker concurrency only when memory headroom is known.

Configure execution capacity in [Choose where the analysis runs](performance_distributed_computing.md). If a run is unexpectedly slow or repeatedly rebuilds targets, use [Troubleshooting](troubleshooting.md).


<!-- source: website/performance_distributed_computing.md -->

# Choose where the analysis runs



A **worker** is an R process that runs an analysis task. A **controller** starts
and manages those workers, either on your machine or through a cluster
scheduler. multiomeR uses `crew` for this. The committed configuration runs
the demo unchanged on a machine that meets the [system
requirements](demo_installation.md#system-requirements); adjust it as
described below before running on a smaller machine or a scheduler.

Use local workers on a suitable workstation. On a shared cluster, ask your
support team which scheduler, account, and resource limits to use. The
[targets distributed-computing guide](https://books.ropensci.org/targets/crew.html)
explains the general setup; this page covers multiomeR's configuration file.

## Local execution

A fresh clone includes a local `crew_controllers.R` sized for a 16-CPU,
256-GB workstation, with four light workers and two heavy workers. A machine
near the 60-GB minimum should reduce concurrency to one heavy worker and should
not run several memory-intensive targets simultaneously. Edit the worker
counts and resource tiers directly in `crew_controllers.R`, keeping controller
names identical between `controller_list` and
`controller_resources_tibble`.

After changing the file, restart R or reload the project runtime explicitly:

```{.r filename="R"}
load_project_runtime(force = TRUE)
```

Rebuild a narrow manifest selection before starting the data run to validate the controller contract.

## Scheduler execution

For SLURM, PBS, SGE, or LSF, replace the local controllers with the corresponding
`crew.cluster` controllers. The commented SLURM section in
`crew_controllers.R` shows the expected shape.

For every scheduler tier:

1. Match the controller name in both the controller object and resource table.
2. Align scheduler CPU and memory requests with the capacity declared in the table.
3. Set queue, account, wall-time, module, and worker-startup options required by the cluster.
4. Keep GPU tiers separate; GPU controllers are considered only for targets requesting GPUs.
5. Test a small target selection before increasing worker counts.

Scheduler startup failures, resource-routing errors, and target failures are handled separately in [Troubleshooting](troubleshooting.md). Developer-facing details about runtime bootstrap and `get_tar_resources()` are in [Implementation conventions](implementation/implementation_conventions.html#runtime-bootstrap).

## Controller contract

`crew_controllers.R` is sourced during `load_project_runtime()` and must return a named list containing these components:

- `controller_list`: a non-empty list of `crew` controllers with unique controller names.
- `controller_resources_tibble`: a data frame with exactly `controller_name`, `cores`, `RAM_GB`, and `gpus`, in that order.

The resource columns must be numeric, non-missing, and contain one unique row per controller name represented in `controller_list`. The first resource-table row is the default controller. For explicit requests, `get_tar_resources()` selects the first compatible row after applying the requested CPU, RAM, and GPU constraints.

```r
controller_resources_tibble <- tibble::tribble(
  ~controller_name, ~cores, ~RAM_GB, ~gpus,
  "local-light",        1,      16,     0,
  "local-heavy",        6,      60,     0
)
```

The table describes controller capacity for routing. A local controller does not create physical memory: its `workers` value must be low enough that concurrent jobs cannot exhaust the machine.


<!-- source: website/troubleshooting.md -->

# Troubleshooting

Find the target name and first error message in the run output. Fix that
cause, then rerun the same selection: `targets` can reuse completed work.
Keep the store intact, since it also contains the records needed to diagnose
the failure.

| What happened? | Start here |
|---|---|
| Setup fails before any data processing | [The manifest does not build](#the-manifest-does-not-build) |
| A target reports an error | [A run reports errored targets](#a-run-reports-errored-targets) |
| A completed result needs rebuilding | [A target is unexpectedly outdated](#a-target-is-unexpectedly-outdated) |
| Workers do not start or are killed | [Controller and scheduler failures](#controller-and-scheduler-failures) |

Run these commands in the repository's Pixi R session. For general debugging
techniques beyond the project helpers below, see the
[targets debugging guide](https://books.ropensci.org/targets/debugging.html).

## The manifest does not build

Run the project bootstrap and manifest in a repository-root R session:

```{.r filename="R"}
load_project_runtime(force = TRUE)
targets::tar_manifest(callr_function = NULL)
```

Failures at this stage usually indicate:

- an unknown, misspelled, wrongly typed, or required YAML parameter;
- an aggregation referencing a GEM well absent from `cfg_GEM_wells.tsv`;
- an unknown module or a missing module-config row;
- a malformed `crew_controllers.R` return value; or
- an R package or startup problem.

Correct the field or file named in the error, then repeat the manifest check
before running the analysis.

## A run reports errored targets

List the errors with the project's helper. It groups repeated failures so
you can start with their common cause:

```{.r filename="R"}
list_distinct_errored_targets()
```

For commands and stored tracebacks matching a target, aggregation, or module:

```{.r filename="R"}
list_distinct_errored_targets_w_tracebacks(
  target_name_pattern = "your_target_or_aggregation"
)
```

Some targets run separately for multiple groups; these runs are called
**branches**. Copy the full target or branch name from the error listing to
inspect its saved workspace:

```{.r filename="R"}
inspect_target_workspace("full_target_or_branch_name")
```

The inspection reports the target command, stored error and traceback, plus compact summaries of its dependencies. Empty tables, zero-dimensional matrices, missing model columns, or unexpected labels usually point to an upstream data or configuration problem.

## Inspect one upstream value

Read a target by its full stored name when a focused probe is necessary:

```{.r filename="R"}
value <- targets::tar_read_raw("full_target_or_branch_name")
value
```

Prefer `head()`, `dim()`, `names()`, or a small subset over printing a large object.

## A target is unexpectedly outdated

Ask `{targets}` which part of a narrow endpoint needs rebuilding:

```{.r filename="R"}
targets::tar_outdated(
  names = tidyselect::matches("your_target.*your_aggregation"),
  callr_function = NULL
)
```

Code, configuration, input files, controller-independent global objects, or an upstream invalidation can all make downstream targets outdated. Review the returned upstream names before assuming the final target itself is the cause.

## Input and metadata failures

Check the contracts in the [GEM well table](reference_GEM_wells.md), [Donor metadata table](reference_donor_metadata.md), and [Aggregation configuration](reference_aggregations.md) references:

- `GEM_well_cellranger_arc_count_dir` contains the required `outs/` files;
- VCF-backed demultiplexing also has `atac_possorted_bam.bam`;
- donor metadata keys and canonical `cfg_GEM_wells.tsv` keys are present and unique;
- non-key metadata columns belong to only one metadata table; and
- configured donor and GEM well IDs match the metadata values exactly.

## Controller and scheduler failures

If a worker does not start or no controller can satisfy a target request:

1. Validate the names, column order, numeric resource values, and controller membership described in [Choose where the analysis runs](performance_distributed_computing.md).
2. Reload with `load_project_runtime(force = TRUE)` after edits.
3. For scheduler controllers, inspect the scheduler output/error log and confirm queue, account, wall time, memory, CPU, module, and filesystem settings.
4. Reduce concurrency when local workers are being killed for memory pressure.

## Rerun safely

After fixing the cause, rerun the same target selection. Successful upstream
results remain cached.

```{.r filename="R"}
targets::tar_make(
  names = tidyselect::matches("your_target.*your_aggregation")
)
```

Use an unqualified `targets::tar_make()` only when every active aggregation and enabled module is intentionally in scope. Avoid deleting the target store as a debugging step: it removes evidence and forces unrelated recomputation.


## Part: Reference


<!-- source: website/review_outputs.md -->

# Output files and metadata

The [running guide](main_running.md#steps) lists the files to open after each
checkpoint. Read each plot's subtitle and caption for interpretation and
method details.

## Output folders and metric selection {#output-folders}

Plots live under `<store>/plots/<scope>/<checkpoint>/`. The store comes from
`_targets.yaml`; scope is a GEM well or aggregation. Existing files are not a
record of which targets are current: use `targets::tar_outdated()` before
reviewing results after a configuration change.

`QC_metric_manifest.tsv` selects metrics, display labels and plotting quantiles.
`do_plot = FALSE` hides a metric; plotting quantiles change the displayed range,
not the cells retained by the pipeline. Filtering is configured separately in
`cfg_GEM_wells.tsv` and `cfg_aggregations.yaml`.

GEX, ATAC and WNN parameter sweeps belong to checkpoints 3, 7 and 8.
The GEX compatibility object belongs to checkpoint 3 and the multimodal object
to checkpoint 8. Peak–gene correlation is a separate optional module.

## Cell retention tables {#cell-retention}

Append `.my_aggregation` to these target names when reading them with
`targets::tar_read()`:

| Target | Checkpoint |
|---|---|
| `cell_retention_tibble.GEX_input` | 1 |
| `cell_retention_tibble.GEX` | 3 |
| `cell_retention_tibble.ATAC_input` | 5 |
| `cell_retention_tibble.ATAC` | 7 |
| `cell_retention_tibble.WNN` | 8 |

Tables accumulate stages and retain wells with zero surviving cells. Rows identify
an exclusion action and its order, input count, excluded count and retained count.

## Cluster numbering {#cluster-numbering}

Leiden clusters receive size-ranked IDs before minimum-size filtering. Filtering
does not renumber survivors, so IDs can have gaps.

## Cluster annotation {#cluster-annotation}

Cluster labelling uses adjusted-score advantage over matched random controls.
Marker lists accept unsigned genes, positive `+` suffixes and negative `-` suffixes.
Signed UCell scores are clipped at zero per cell before averaging, for observed
and matched-control signatures alike. The method uses the same
GEX control reference for GEX, ATAC and WNN clusters, and returns `Assigned`,
or `Unassigned` annotation status. The adjusted score is the mean UCell score
minus that label's matched-control 95th percentile. The highest adjusted score
nominates the candidate. Set `aggregation_cluster_annotation_min_advantage`
(default **0.05**) to the required lead over both zero background and the
next-best adjusted score. Exact ties and scores at or below background remain
unassigned even at a zero threshold. Raising the threshold only withdraws
assignments; it never switches the nominated label. An unassigned cluster
retains its candidate and a reason (weak background advantage or competition)
in the diagnostics. This status does not remove cells. Mixture detection is
not included. The former module-score labeller and its
`aggregation_allow_multiple_cell_types` setting have been removed. Remove that
setting and `aggregation_cluster_annotation_method` from older configuration
files; there is now one annotation method. Signed marker panels remain supported.

Annotation exports include `clusters.tsv`, `marker_evidence.tsv`,
`control_gene_matching.tsv` and, where eligible, `GEM_well_agreement.tsv`.
The `cluster_UCell_evidence` targets cache scoring separately from the
threshold-dependent annotation. `marker_set_UCell_summary.3_GEX_QC` exports
`marker_sets.tsv` and `method.txt` before GEX doublet filtering.

## Further methods

See [algorithm validation](implementation/algorithm_validation.html) for reference
comparisons and deviations, and the [peak–gene module](downstream_peak_gene_correlation.md)
for its model, support filters and output paths.


<!-- source: website/reference_GEM_wells.md -->

# GEM well table



`cfg_GEM_wells.tsv` has one row per GEM well: one `cellranger-arc count`
output together with its donor assignment, optional inputs, and
pre-aggregation QC filters. Aggregations refer to rows by `GEM_well_ID`. Edit
the committed file directly; the demo rows can stay as worked examples. This
page describes the columns. When to edit them, and how to review the effect,
is step 1 of [Run your own analysis](main_running.md#steps).

## Minimal row

The example below describes one non-multiplexed GEM well from one donor. This
vertical view is a reading aid; the saved TSV has one GEM well per row.

::: {.scrollable-table}

| Column | Example value |
|---|---|
| `GEM_well_ID` | `your_GEM_well` |
| `GEM_well_dataset` | `your_dataset` |
| `GEM_well_donor_id` | `donor_1` |
| `GEM_well_n_donors` | `1` |
| `GEM_well_cellranger_arc_count_dir` | `/path/to/your_GEM_well` |
| `GEM_well_add_cellbender` | `FALSE` |
| `GEM_well_cellbender_h5_file` | `NA` |
| `GEM_well_donors_VCF_file` | `NA` |
| `GEM_well_QC_exclude_list` | `NA` |
| `GEM_well_is_active` | `TRUE` |
| `GEM_well_multiplex_batch` | `batch_1` |

:::

## Identity and grouping

`GEM_well_ID` must be unique and is used in target names, output folders, and
`aggregation_GEM_well_IDs`. `GEM_well_dataset` is a reporting label, not the donor or the biological
design. Cross-well QC comparisons follow the wells selected for an aggregation.

## Input files

`GEM_well_cellranger_arc_count_dir` points to the directory containing `outs/`, not to `outs/` itself. The baseline pipeline requires:

``` {.text filename="Text"}
<GEM_well_cellranger_arc_count_dir>/outs/
|-- summary.csv
|-- filtered_feature_bc_matrix.h5
|-- atac_fragments.tsv.gz
|-- atac_fragments.tsv.gz.tbi
`-- per_barcode_metrics.csv
```

multiomeR identifies the reference from the FASTA/GTF hashes in `outs/atac_fragments.tsv.gz` and selects the matching `reference.json` under `reference_metadata/`. For another reference, add its JSON there. The fragment header must be intact and exactly one JSON must match. All GEM wells in an aggregation must use the same reference. Deployments that
require an explicit `GEM_well_cellranger_arc_reference_json` column must supply
that path for each well; follow the checkout README for its input contract.

`GEM_well_add_cellbender` set to `TRUE` replaces the Cell Ranger GEX counts
with a CellBender H5 file produced outside multiomeR, given in
`GEM_well_cellbender_h5_file` or found in the standard layout.

## Donors

`GEM_well_n_donors` states how many donors the GEM well contains. For a
non-multiplexed well, set it to `1` and put the donor's identifier in
`GEM_well_donor_id`; every called nucleus receives that ID, which must match a
`donor_id` row in the [donor metadata table](reference_donor_metadata.md).

For a multiplexed well, supply the donors' genotypes in
`GEM_well_donors_VCF_file`. `atac_possorted_bam.bam` is then also required in
`outs/` for `cellsnp-lite`, and nuclei are assigned to donors by genotype.

## Pre-aggregation QC filters

`GEM_well_QC_exclude_list` contains zero or more complete R filter expressions
separated by `;;`, for example:

```{.text filename="cfg_GEM_wells.tsv"}
TSS.enrichment < 4 ;; nucleosome_signal > 4 ;; nCount_RNA < 250
```

Expressions are evaluated individually against per-barcode metadata,
preserving their order and their separate exclusion reasons. An empty field
(`NA`) applies no pre-aggregation QC filters. Leave it empty for a new GEM well
until step 1 of [Run your own analysis](main_running.md#steps) has shown its
distributions; do not treat the demo's numerical cutoffs as recommendations for
your tissue.

AMULET doublet detection runs as part of the standard QC calculation; see the
[implementation and validation](implementation/algorithm_validation.html#bpcells-native-amulet)
for details.

## Active flag

`GEM_well_is_active` controls whether processing targets are constructed for
the GEM well. Every GEM well selected by an active aggregation must also be
active. Before an unqualified `targets::tar_make()`, deactivate every GEM well
you are not ready to run, including demo rows you do not need.

## Extra metadata columns

Add library-, run-, or batch-specific variables as further columns with a
`GEM_well_` prefix, such as `GEM_well_multiplex_batch` or
`GEM_well_cell_sorting`. They are joined to the cell metadata and can be used
as Harmony covariates or plotting variables. Donor-level phenotypes belong in
the [donor metadata table](reference_donor_metadata.md) instead. Apart from
their key columns, the two tables must not reuse column names.

## Column dictionary

Use the [committed example](https://github.com/koefoeden/multiomeR/blob/main/configuration/cfg_GEM_wells.tsv)
to check the full set of columns.

The table below is a documentation snapshot of the two public demo wells,
showing the core columns and one optional annotation. Bold columns must be
present in the TSV; some allow an NA value. Scroll horizontally and focus or
hover over a column's **i** button for its meaning. The other inactive rows and
metadata columns in the public configuration remain available as examples.

[Generated Quarto chunk omitted: `emit_GEM_well_demo_table( GEM_well_config_file = "website/data/demo_GEM_wells.tsv", dictionary_file = "website/data/G...`]


<!-- source: website/reference_donor_metadata.md -->

# Donor metadata table

The donor metadata table is a TSV with one unique row per `donor_id`. Each
aggregation points to one such file through
`aggregation_donor_id_metadata_tsv` in the [aggregation
configuration](reference_aggregations.md). It is created in step 2 of [Run
your own analysis](main_running.md#steps).

## Minimal table

``` {.text filename="donor_metadata.tsv"}
donor_id	condition
donor_1	control
```

## Matching donors to nuclei

Every nucleus receives a `donor_id` from its GEM well: the configured
`GEM_well_donor_id` for a non-multiplexed well, or a genotype-based assignment
for a well with a configured VCF. Each of those IDs must appear exactly once in
this table; see the [GEM well table](reference_GEM_wells.md#donors).

## Which variables belong here

Put donor-specific phenotypes and covariates in this table, for example
condition, age, or sex. Put library-, run-, or batch-specific variables in the
GEM well table with a `GEM_well_` prefix. Apart from their key columns, the two
tables must not reuse column names. Metadata file contents are validated when
their targets run, not when the manifest is built.

## Extended table for differential analyses

The [differential analyses](downstream_differential_analyses.md) module can
read additional donor-level model variables from a second table given in
`differential_analyses_extended_donor_id_metadata_tsv`. It must keep the same
unique `donor_id` key. If it is not set, the module inherits the aggregation's
donor table.


<!-- source: website/reference_aggregations.md -->

# Aggregation configuration



`cfg_aggregations.yaml` has one top-level entry per aggregation: a joint GEX,
ATAC, and WNN analysis of one or more GEM wells. The committed file enables
the two human GEM wells in `immune_human_2x`, with optional modules disabled.
The mouse and ENCODE validation examples are inactive by default. Edit the file
directly; the demo entries can stay as worked examples. This page describes the
entry structure and lists every parameter. When to set each parameter, and how
to review the effect, is given step by step in [Run your own
analysis](main_running.md#steps).

## Minimal entry

``` {.yaml filename="cfg_aggregations.yaml"}
your_aggregation:
  aggregation_GEM_well_IDs: [your_GEM_well]
  aggregation_donor_id_metadata_tsv: /path/to/donor_metadata.tsv
  aggregation_GEX_marker_genes:
    Cell_type_A: [GENE1, GENE2]
    Cell_type_B: [GENE3, GENE4]
  is_active: true
```

Every other parameter has a default from `cfg_pipeline_parameters.tsv`, listed
in the [parameter reference](#parameter-reference) below. Add a parameter to
the entry only when you want to change its default.

## Required keys

- `aggregation_GEM_well_IDs`: the `GEM_well_ID` values to combine. Each must be
  an active row of the [GEM well table](reference_GEM_wells.md), and all must
  use the same Cell Ranger reference.
- `aggregation_donor_id_metadata_tsv`: the [donor metadata
  table](reference_donor_metadata.md) for these GEM wells.
- `aggregation_GEX_marker_genes`: a named list of marker genes per expected
  cell type, used for cluster annotation and marker plots.
- `is_active`: whether targets are constructed for the aggregation. Deactivate
  aggregations you are not ready to run before an unqualified
  `targets::tar_make()`.

## Marker genes and transcription factors

Replace the placeholder genes with symbols appropriate for the tissue and
reference. A gene listed without a suffix or with a `+` suffix is a positive
marker; a `-` suffix marks a gene that should be absent. The optional
`aggregation_ATAC_marker_TFs` list names transcription factors per cell type for
the motif-activity plots. How the annotation uses these lists is described in
[Output files and metadata](review_outputs.md#cluster-annotation).

## QC filters after peak calling

`aggregation_QC_exclude_list_combined_object` lists dplyr filter expressions
applied to the peak-based ATAC metrics of the combined object, for example:

```{.yaml filename="cfg_aggregations.yaml"}
aggregation_QC_exclude_list_combined_object:
  - nCount_ATAC < 1000
  - atac_peak_counts_frac < 0.1
  - atac_peak_counts_blacklist_frac > 0.01
```

Omit it or set it to `null` until step 5 of [Run your own
analysis](main_running.md#steps) has shown the distributions.

## Optional modules

Omit `modules` for the first run. After reviewing the main results, enable an
optional analysis by listing its name and adding a matching entry for the
aggregation in the module's own configuration file:

```{.yaml filename="cfg_aggregations.yaml"}
your_aggregation:
  modules: [differential_analyses]
```

See [Differential analyses](downstream_differential_analyses.md) and
[Genetic enrichment](downstream_genetic_enrichment.md) for the module
entries and their parameters.

## Parameter reference {#parameter-reference}

The searchable overview below is generated from `cfg_pipeline_parameters.tsv`,
using a shared snapshot of the public runtime defaults and validation schema. Search by name or
purpose, or choose a topic. Defaults are visible beside each parameter; open a
row for its type and example. See the [committed
example](https://github.com/koefoeden/multiomeR/blob/main/configuration/cfg_aggregations.yaml)
for a complete configuration.

[Generated Quarto chunk omitted: `emit_parameter_overview("aggregation")`]

<details>

<summary>Show the public <code>immune_human_2x</code> example</summary>

[Generated Quarto chunk omitted: `emit_yaml_entry(aggregations_config_file, "immune_human_2x")`]

</details>


# Book: multiomeR Implementation


<!-- source: website/implementation/index.md -->

# Introduction

Use this book to trace a result back to its code or change how multiomeR
works. For installation, configuration, execution, and output inspection,
start with the [user manual](../). You do not need to read this book to run
the demo.

Use this book when you need to trace a configuration value into mapped targets, understand how the simplified graph views relate to the real `{targets}` graph, or decide where an implementation change belongs.

## Where to start

For a first implementation pass:

1. Read [Reading the graph views](graph_methodology.md) and follow its configuration-to-target trace.
2. Open the [main pipeline](implementation_main.md) graph for the modality or checkpoint you plan to change.
3. Use [Implementation conventions](implementation_conventions.md) to understand the relevant manifest, mapping, symbol, tag, and runtime contracts.
4. Read [Background and design philosophy](background_philosophy.md) when you need the rationale for the editable-workflow design.

The [differential analyses](implementation_differential_analyses.md) and [genetic enrichment](implementation_genetic_enrichment.md) chapters cover the optional module graphs.

## Common entry points

| Change | Start with |
|---|---|
| Add or revise a YAML parameter | `cfg_pipeline_parameters.tsv`, then the owning config reader or target. |
| Change GEM well preprocessing | `_targets.R` mapping plus `extra_targets/per_GEM_well_targets.R`. |
| Change aggregation GEX, ATAC, or WNN processing | The corresponding graph section and `extra_targets/*_targets.R` file. |
| Inspect existing review selections | `[checkpoint:<name>]` description tags and the steps in [Run your own analysis](../main_running.html#steps). |
| Add a graph-visible target | Existing `[part_of_graph:<graph_id>]` tags and graph-pruning rules. |
| Change resource routing | `crew_controllers.R`, `packages/multiomeRCore/R/resource_helpers.R`, and the runtime bootstrap convention. |

If you are trying to run multiomeR rather than modify it, start with the [main manual](../).


## Part: Orientation


<!-- source: website/implementation/graph_methodology.md -->

# Reading the graph views



The graph chapters collect simplified views of the real `{targets}` dependency graph. They are meant to make the workflow easier to reason about before reading the target code directly.

The diagrams are generated from tagged target metadata and the real dependency graph, then simplified by pruning or bypassing lower-level nodes that would make each view harder to read. They keep real target names and preserve the dependency structure where practical, while staying compact enough to build intuition about the main control points.

The following chapters cover the main pipeline, the differential analyses module, and the genetic enrichment module.

## Trace one configured aggregation

The quickest way to understand the implementation is to follow one value across the graph:

1. `immune_human_2x` is a key in `cfg_aggregations.yaml`.
2. `read_aggregation_config_tibble()` resolves manifest defaults and inheritance into one aggregation row.
3. `build_aggregation_tibble()` filters active rows and adds symbols for GEM-well upstream targets, including those used for aggregation-level QC summaries.
4. The root `_targets.R` passes that row through `tar_map(names = aggregation, delimiter = ".")`.
5. A base target such as `multimodal_Seurat_object.8_multimodal_QC` becomes `multimodal_Seurat_object.8_multimodal_QC.immune_human_2x`.
6. Description tags make selected targets discoverable as checkpoints or graph nodes, while structured file helpers derive output paths from the active target name.

Inspect the exact target command and description without running it:

```r
targets::tar_manifest(
  names = tidyselect::matches(
    "^multimodal_Seurat_object[.]8_multimodal_QC[.]immune_human_2x$"
  ),
  fields = c(name, command, description),
  callr_function = NULL
)
```

This trace connects the [parameter manifest](implementation_conventions.md#parameter-manifest), [mapping tibbles](implementation_conventions.md#mapping-tibbles), and [target-symbol columns](implementation_conventions.md#target-symbol-columns) before the larger diagrams introduce many nodes at once.

## What the diagrams omit

The curated graph views are orientation aids, not alternate target definitions. A node can be absent because it was pruned as a lower-level implementation detail, bypassed to preserve a useful dependency path, or omitted because it lacks the graph-membership tag for that view. Use `tar_manifest()` or the source target files when exact completeness matters.

[Mermaid graph omitted; source: `website/figures/standard_node_color_legend.mmd`]


<!-- source: website/implementation/implementation_conventions.md -->

# Implementation conventions

This chapter explains how configuration becomes target definitions: the
parameter manifest supplies defaults and validation rules, mapping tables
define repeated analyses, and target symbols connect their dependencies.
Description tags support result selection and graph views. The final section
covers project startup and resource configuration.

## Target metadata tags

multiomeR stores lightweight target metadata in the `description` argument of `targets::tar_target()` and `tarchetypes::tar_file()` calls. The descriptions should remain readable prose, with bracketed tags appended when a target needs to be discoverable from the manifest.

```r
targets::tar_target(
  name = harmony_embeddings_matrix.GEX,
  description = "Harmony-corrected SCTransform GEX PCA embeddings [part_of_graph:GEX] [part_of_graph:WNN]",
  command = ...
)
```

The currently meaningful tag families are:

```text
[checkpoint:<name>]             review or execution checkpoint
[part_of_graph:<graph_id>]      curated membership in an implementation graph
[resource_observation:<note>]   compact empirical resource note
```

`[checkpoint:<name>]` marks targets selectable with `targets::tar_described_as()`.
The eight numbered main-pipeline groups are listed in
`QC_checkpoint_manifest.tsv`; optional module groups remain unnumbered.
UMAP parameter sweeps belong to their modality's numbered checkpoint, and
compatibility objects belong to GEX checkpoint 3 or multimodal checkpoint 8.
Selection matches description substrings; include the
closing `]` to match a complete checkpoint tag. Dependencies still come from
the target commands. [Run your own analysis](../main_running.html#steps)
explains each boundary; acceptance criteria depend on the study.

Numbered checkpoint plot targets end in the checkpoint name with hyphens
replaced by underscores, before the mapped dataset or aggregation suffix.
For example, `VizDimLoadings_plots.2_GEX_PCA_QC.my_aggregation` writes beneath
`<store>/plots/my_aggregation/2_GEX_PCA_QC/`. Only plot targets use this naming
convention; computational and metadata targets retain their modality suffixes.

`[part_of_graph:<graph_id>]` marks targets that should stay visible in a named implementation graph after graph-pruning helpers remove less informative intermediate nodes. This is the strictest tag family: `graph_id` must contain only letters, numbers, and underscores, and helper code parses these tags directly from target descriptions. A target may belong to several graph views.

```r
description = paste(
  "Build the lightweight BPCells-backed chromVAR RSE",
  "[part_of_graph:ATAC]",
  "[part_of_graph:seurat_export]",
  "[part_of_graph:genetic_enrichment_single_nucleus]"
)
```

`[resource_observation:<note>]` is currently best treated as provisional documentation. It is useful when a target has a compact empirical runtime or memory observation worth keeping near the target definition, but it is not yet a structured resource-estimation system. Keep these notes short, dated when relevant, and self-explanatory.

Use tags only when they create a durable handle for readers, graph helpers, or checkpoint commands. Ordinary internal dependencies can stay untagged.

## Parameter manifest

`cfg_pipeline_parameters.tsv` is the schema for YAML-backed pipeline configuration. Each row defines one parameter for one scope:

```text
aggregation
differential_analyses
genetic_enrichment
```

For each parameter, the manifest records its name, type, cardinality, default value, missing-value rule, allowed values, example values, topic, graph/module ownership, and human description. The YAML files then only need to specify values that differ from the manifest defaults, plus values that are required because their resolved value may not be missing.

At read time, the pipeline loads the manifest for a scope and parses each `default_value` as YAML. This allows defaults to be literal scalars, `NULL`, YAML lists, or evaluated YAML expressions such as `!expr 1:30`. Manifest defaults seed every config row before inheritance and row-specific overrides are applied.

The resolution order is:

1. Start with manifest defaults for the requested scope.
2. Resolve each parent listed in `inherits`.
3. Overlay parent values onto the defaults.
4. Overlay the child row onto the inherited values.
5. Validate the fully resolved row.

```yaml
immune_human_2x:
  aggregation_GEM_well_IDs: [healthy_PBMC_human, lymphoma_lymph_human]
  aggregation_GEX_marker_genes:
    B: [MS4A1, CD79A]
    T: [TRAC, CD3D]

PBMC_human_6x:
  inherits: immune_human_2x
  aggregation_GEM_well_IDs:
    - healthy_PBMC_human
    - pbmc_10k_chromium_controller
  modules: [genetic_enrichment]
```

Validation is manifest-driven and happens before target construction. Unknown YAML parameters fail early. Resolved values are then checked for missingness, cardinality, type, and allowed values.

```text
scalar      one non-list value
vector      atomic vector
list        list
named_list  list with non-empty names
```

The `data_type` column checks the R type after YAML parsing. `path` and `regex` are currently character-like schema labels; the validator does not check file existence or compile regular expressions. `allowed_values` is a comma-separated allow-list checked after coercing resolved values to character.

The website renders parameter tables from a generated snapshot of the public runtime manifest. Refresh that snapshot with the shared documentation inputs when public defaults change; private configuration is never used for these tables.

Module-specific YAML uses the same mechanism. Aggregations opt into modules through the aggregation config, and each enabled aggregation must have a matching module config row. Some cross-scope fallbacks are still implemented by target code rather than by manifest inheritance; for example, a module parameter may intentionally allow `NULL` and then fall back to an aggregation-level path during module setup.

## Mapping tibbles

The root `_targets.R` builds the target graph from mapping tibbles. Each mapping tibble is a row-wise contract: one row becomes one set of mapped target instances, and columns in that row become local symbols inside the corresponding `tarchetypes::tar_map()` block.

The core mapping flow is:

Configuration readers use `configuration_path()` to resolve a basename within
`configuration/` or the directory named by the ignored `configuration.local`.
Relative selections are anchored at the repository root. This selection does not
change data-path interpretation or `_targets.yaml`. The selected GEM-well path
is a graph global consumed by the file target, so changing directories also
changes its dependency. Aggregation and enabled-module settings are resolved
during graph construction. Disabled modules do not read their configuration.

1. `GEM_well_tibble_all` reads only the pre-aggregation processing columns from every row in the canonical `cfg_GEM_wells.tsv`.
2. `aggregation_tibble_all_from_yaml` is read from `cfg_aggregations.yaml`.
3. `aggregation_tibble` keeps active aggregations, validates their GEM well references against the complete view, and adds upstream target-symbol columns.
4. `GEM_well_tibble` keeps GEM wells whose `GEM_well_is_active` value is true.
5. `_targets.R` expands active GEM wells and aggregations with `tar_map()`, then appends module target files. Cross-GEM-well QC summaries use the aggregation's selected wells.

Within each aggregation, `GEM_well_metadata_tibble` reads the same canonical
file, subsets it to `aggregation_GEM_well_IDs`, and preserves that order. Cheap
keyed projection targets then expose only the columns requested for SCT,
Harmony or configured analyses. Complete non-processing
annotations are joined only for explicit export objects. These projection
targets are cache boundaries: a newly added or edited online column can update
the canonical table without changing expensive consumers whose selected view
is identical.

```r
tarchetypes::tar_map(
  values = GEM_well_tibble,
  names = GEM_well_ID,
  delimiter = ".",
  source("extra_targets/per_GEM_well_targets.R")$value
)
```

With `GEM_well_ID = "healthy_PBMC_human"`, a target named `cellranger_summary_file` becomes `cellranger_summary_file.healthy_PBMC_human`. The same dot-delimited suffix convention is used for aggregations, module targets, and nested module maps.

```text
active cfg_GEM_wells.tsv row -> GEM_well_tibble row    -> per GEM well targets
cfg_aggregations.yaml key   -> aggregation_tibble row -> per-aggregation targets
```

Aggregation rows may opt into optional modules through `modules`. `_targets.R` validates module names against the known module list, and module target files then filter `aggregation_tibble` to the active aggregations that requested that module. Each opted-in aggregation must have a matching module config row.

The naming convention is therefore compositional:

```text
<target>.<GEM_well_ID>
<target>.<dataset_name>
<target>.<aggregation_name>
<module_target>.<module_name>.<aggregation_name>
<nested_module_target>.<nested_suffix>.<module_name>.<aggregation_name>
```

Because these suffixes become target names and cache identity, config keys should be stable, human-readable, and free of unnecessary punctuation. In particular, avoid dots in GEM well, dataset, aggregation, and module IDs unless there is a compelling reason.

## Target-symbol columns

Mapped target tables sometimes need to carry references to other mapped targets. multiomeR represents those references as columns of `rlang` symbols. Each row stores the upstream target symbols that should be spliced into downstream target commands generated for that row.

The compact constructor is `target_sym_col()`. It records a base target name, the source column containing suffixes, the separator, and an optional transform. `add_target_sym_cols()` then turns those specifications into list-columns of `rlang::syms()`.

```r
aggregation_tibble |>
  add_target_sym_cols(
    aggregation_GEX_counts_BPCells_matrix_syms =
      target_sym_col("GEX_counts_BPCells_matrix", "aggregation_GEM_well_IDs")
  )
```

For an aggregation whose `aggregation_GEM_well_IDs` are `c("rx1", "rx2")`, this creates a row value equivalent to:

```r
rlang::syms(c(
  "GEX_counts_BPCells_matrix.rx1",
  "GEX_counts_BPCells_matrix.rx2"
))
```

The aggregation target can then consume the row-local symbol list directly:

```r
combined_counts_matrix <- purrr::reduce(
  aggregation_GEX_counts_BPCells_matrix_syms,
  cbind
)
```

Column names should describe the downstream scope, the upstream target, and the fact that the value is a symbol list. The `aggregation_*_syms` columns, such as `aggregation_GEX_counts_BPCells_matrix_syms`, splice per GEM well targets into aggregation-level targets.

Module target files also need aggregation-specific references to main-pipeline targets. For this, `add_aggregation_target_syms()` creates one symbol per row, suffixed by the aggregation name. These columns are named like the target they replace rather than with `*_syms`, because each cell is a single symbol rather than a list.

```r
differential_analyses_tibble |>
  add_aggregation_target_syms(c(
    "metadata_w_cell_types_tibble.WNN",
    "pseudobulk_counts_matrix.GEX",
    "organism_chr"
  ))
```

For aggregation `PBMC`, the column `metadata_w_cell_types_tibble.WNN` contains the symbol `metadata_w_cell_types_tibble.WNN.PBMC`. Inside a module target, the command can be written against the unsuffixed local name; `tar_map()` resolves it to the aggregation-specific upstream target for that row.

## Runtime bootstrap

multiomeR assumes that the repository runtime is bootstrapped before the target graph is inspected or run. The root `.Rprofile` is intentionally minimal: it sources `R/bootstrap_helpers.R` and calls `load_project_runtime()`.

`load_project_runtime()` is the single entry point for:

1. loading core workflow packages and conflict preferences,
2. sourcing generally reusable helpers from `packages/multiomeRCore/R`,
3. sourcing pipeline-specific helpers from the root `R/` directory,
4. applying global plotting and `{targets}` options,
5. sourcing `crew_controllers.R` and installing controller resources.

The nested `multiomeRCore` directory is both ordinary editable pipeline source
and an installable package boundary for standalone repositories. multiomeR does
not install or attach that package itself: `targets::tar_source()` loads the
same implementation files before the root helpers. Keep domain-specific code
under `R/`, but do not duplicate the generally reusable implementations there.

For commands that intentionally bypass startup side effects, source the bootstrap helper directly and then load the runtime:

```r
source("R/bootstrap_helpers.R")
load_project_runtime(force = TRUE)
targets::tar_manifest(callr_function = NULL)
```

Bootstrap state is cached in `bootstrap_state_env`. This avoids reloading packages, re-sourcing helpers, reapplying target options, reassigning patches, and reloading controllers on every call. Use `force = TRUE` when the current R session may be stale, such as after changing helper files, switching checkout roots, editing `crew_controllers.R`, or reusing a long-lived interactive session.

Project-root detection walks upward from the current working directory until it finds `pixi.toml`. Bootstrap commands should therefore be run from inside the multiomeR checkout.

Controller loading is part of the runtime contract, not a later execution detail. `crew_controllers.R` must return a named list with `controller_resources_tibble` and `controller_list`. The bootstrap validates that shape, installs a grouped `crew` controller into `{targets}`, and stores the resource table for `get_tar_resources()`.

```r
targets::tar_target(
  example_target,
  example_function(),
  resources = get_tar_resources(cores_req = 6, RAM_GB_req = 60)
)
```

If `get_tar_resources()` is called before controller resources are loaded, it fails deliberately with an instruction to call `load_project_runtime()` first. Scheduler-specific examples belong in the main manual's [Choose where the analysis runs](../performance_distributed_computing.html) page; the implementation contract is that target code can request resources declaratively once the runtime has been loaded.

## How to read the rest of the implementation book

These conventions are the connective tissue behind the graph chapters. The parameter manifest explains why config rows can be compact. Mapping tibbles explain why target names have stable suffixes. Target-symbol columns explain how mapped targets pass sets of upstream targets across graph levels. Target metadata tags explain why some nodes remain visible in curated graph views. The bootstrap contract explains why helper functions, controller resources, and target options are available before `_targets.R` is evaluated.

When modifying the implementation, preserve these contracts unless the change is explicitly meant to replace one of them.

### Peak–gene correlation module

`module_peak_gene_correlation/targets.R` maps only opted-in aggregations and
binds their existing WNN metadata, GEX/ATAC matrices, ATAC embeddings,
fragments and reference annotations. `correlation_targets.R` owns the analysis,
SuSiE prioritization, exports and plots. Parameters use the
`peak_gene_correlation` manifest scope and matching module YAML rows.
Targets end in `.peak_gene_correlation.<aggregation>` (with `.WNN` before
that suffix for intermediate results). Plot checkpoint tags use
`peak_gene_correlation`, keeping this analysis outside numbered QC selections.
Renamed targets rebuild on the first module run; existing core target names
and numerical analysis defaults are unchanged.


<!-- source: website/implementation/algorithm_validation.md -->

# Algorithmic implementations, deviations and validation

multiomeR reimplements a small number of reference algorithms so they can operate on the workflow's native matrices and graph state. This page is the maintained record of what those implementations preserve, where they deliberately differ, and what the executable validation establishes.

The evidence labels used below are intentionally narrow:

- **Reference-parity tested** means the repository and named reference implementation run on the same deterministic fixture and their returned values are compared directly.
- **Reference-similarity tested** means exact equality is not an appropriate contract, so predefined similarity thresholds are checked against the named reference implementation.
- **Algorithmically derived** means the implementation is checked against an independent mathematical result, not against another software implementation.
- **Internally checked** means a deterministic repository fixture exercises an internal contract without establishing reference parity.

Passing these fixtures does not validate every dataset, parameter regime, approximate-neighbor realization, biological interpretation, or downstream target. The test suite distinguishes fast unit tests for isolated data contracts from slower parity and integration tests that load the project runtime or compare external reference implementations. The [CI workflow](https://github.com/koefoeden/multiomeR/blob/main/.github/workflows/algorithm-validation.yaml) runs the complete suite when tests, relevant helpers, or the Pixi environment change.

Run the complete suite with `pixi run --use-environment-activation-cache test`. The narrower `pixi run --use-environment-activation-cache test-algorithm-validation` task runs only the slow UCell, AMULET, WNN, and SCAVENGE parity and acceptance tests.

| Implementation | Evidence status | Maintained reference | Current fixed-fixture result |
|---|---|---|---|
| BPCells-native UCell | Reference-parity tested | UCell 2.14.0 | Exact values, dimensions, and dimnames |
| BPCells-native AMULET | Reference-parity tested | scDblFinder 1.24.0 | Exact metrics and multi-chromosome loci, including order |
| Native WNN | Reference-similarity tested | Seurat 5.5.0 | Small-SNN pilot: weight Spearman 0.989; mean neighbor overlap 0.996 |
| Sparse SCAVENGE propagation | Algorithmically derived and reference-parity tested | SCAVENGE 1.0.2 at `8ee8b173d965` | Closed-form delta 8.61e-13; pinned-reference propagation delta 1.11e-16; exact streamed exceedance counts and significant-cell calls |

## BPCells-native UCell scoring

**Reference algorithm.** [`UCell::ScoreSignatures_UCell()`](https://bioconductor.org/packages/release/bioc/html/UCell.html) calculates per-cell signature scores from descending feature ranks, caps ranks at `maxRank`, combines positive and negative signatures, and clips negative combined scores to zero. The maintained comparison also covers `UCell::AddModuleScore_UCell()`.

**Reason for reimplementation.** The workflow keeps gene-by-cell counts in BPCells-backed matrices. Materializing the complete matrix in memory or building a Seurat object solely for marker scoring would discard that storage contract, so multiomeR ranks bounded cell chunks and returns metadata-ready scores directly.

**Behavior preserved.** The implementation preserves UCell signature syntax (`+` and `-` suffixes), descending per-cell ranks, configurable tie handling, `maxRank` capping, impute/skip behavior for missing genes, negative-signature weighting, lower-bound clipping, signature names, and cell order.

**Deliberate deviations and consequences.** Matrix materialization is limited to one cell chunk at a time, and optional fork workers operate across chunks. This changes memory and execution behavior but not the tested score values. The helper returns a data frame instead of mutating a Seurat object. The target-level marker validator rejects configured genes missing from the Cell Ranger reference before normal pipeline scoring, whereas the lower-level helper still exposes UCell's impute/skip modes for explicit use.

**Implementation and wiring.** The scorer is [`calculate_BPCells_UCell_scores_from_matrix()` in `R/processing_GEX_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/processing_GEX_helpers.R). [`extra_targets/general_aggregation_targets.R`](https://github.com/koefoeden/multiomeR/blob/main/extra_targets/general_aggregation_targets.R) validates `UCell_GEX_marker_genes_list`. The matched-control annotation method in `R/cluster_annotation_helpers.R` ranks bounded GEX count chunks once per modality's cluster partition, retaining exact sufficient statistics for positive signatures and per-cell scores for GEX metadata. For signed signatures, it scores and clips each cell before aggregation, preserving positive/negative membership in observed signatures, matched controls and marker-deletion variants. This costs more computation than the positive-only rank-summary shortcut; cell chunks bound temporary memory.

**Matched-control cluster annotation.** The method freezes 999 random gene mappings
using up to 50 reference cells per GEM well. Controls exclude all marker genes,
match normalized abundance and detection in standardized transformed coordinates,
and preserve signature overlap through one mapping per marker and replicate.
Each mapping samples without replacement from the nearest 50 eligible candidates.
Reference and matching diagnostics remain available for review.

For each label, the adjusted score is its observed cluster mean UCell score
minus its matched-control 95th percentile. The highest adjusted score nominates
the candidate, with no preliminary tail-score or marker-detection filter.
Its advantage is `best - max(0, second_best)`. The sole assignment threshold is
`aggregation_cluster_annotation_min_advantage`, default 0.05. Assignment requires
an advantage at least this large, a positive best score and no exact tie. Otherwise
the result is `Unassigned`, with the candidate and abstention reason retained.
Candidate ordering is independent of this threshold, so increasing stringency
can only withdraw assignments. Scores are cached separately from decisions.

Marker-deletion diagnostics remove one candidate marker and its matched control
gene from each signature, then compare the resulting adjusted score against the
unchanged competing labels and zero background. Their agreement fraction never
vetoes assignment. Single-marker signatures have no deletion diagnostic.
Ten deterministic deletion blocks stratified by cluster and GEM well assess
the fraction retaining the same assigned candidate at the selected threshold.
Fewer than two assessable deletions yield an undefined diagnostic without
overriding the label. Detection counts (positive markers detected in at least 10% of cells)
and GEM-well subgroup agreement (groups of at least 25 cells) are also diagnostic
only. No mixture subdivision is performed.

Faceted advantage plots use shared label order and y limits across pages,
display non-leading negative scores at zero, and mark zero background and the facet-specific cutoff
`best - min_advantage`. Thus the plotted decision geometry matches the numerical
rule, including explicit handling of exact ties at a zero threshold. Negative
leading scores and cutoff lines remain visible; display clipping never alters
the scores used for assignment.

The GEX module dot plots reuse this pre-doublet-filtering evidence directly.
Both views order marker sets by Euclidean distance between unstandardized
cluster-adjusted profiles and `hclust(method = "ward.D2")`, weighting clusters
equally. Rows follow the assigned cell types in that order, leaving unassigned
groups last. Cluster colours equal the cached adjusted scores; cell-type colours
are cell-count-weighted averages of those scores, rather than a recalculated
background for pooled cells. Both views share symmetric colour limits centred
at zero. Dot area reports the fraction of the same pre-filtering cells with raw
UCell above zero. Marker-expression plots also use pre-filtering cells but keep
their existing configured marker-set order.

These are shared technical defaults, not learned identity probabilities or
validated biological error rates. The control comparison can remain imperfect
for unusually abundant or rare markers. Scoring-parity tests cover exact
positive-signature aggregation, chunk/worker agreement, singleton reference
sampling, threshold boundaries, ties, monotonic abstention and separate
doublet-detection groups. The
signed cluster tests compare cell scores, observed cluster means and matched controls with UCell 2.14.0, including marker deletions and chunk/worker agreement.

**Validation.** [`tests/testthat/test-scoring-parity.R`](https://github.com/koefoeden/multiomeR/blob/main/tests/testthat/test-scoring-parity.R) creates a deterministic 500-gene by 37-cell matrix, writes the project input as BPCells, and compares signed signatures with imputed and skipped missing genes. It requires `identical()` values, dimensions, and dimnames against UCell 2.14.0. It also retains exact Seurat `AddModuleScore` and cell-cycle helper checks plus a metadata-join contract.

The UCell 2.14.0 reference call with imputed missing genes can emit non-fatal R stack-imbalance warnings under the repository's R 4.5 environment. A narrow diagnostic reproduced them in the UCell reference call but not in the repository scorer. The validation therefore runs the reference call in a disposable `callr` process and compares its returned matrix in the clean parent session; this isolates the package warning without weakening the equality assertion or changing the runtime.

**Status and rerun.** Reference-parity tested against UCell 2.14.0 with exact equality; passing in the current locked environment.

```bash
pixi run --use-environment-activation-cache test-scoring-parity
```

## BPCells-native AMULET

**Reference algorithm.** [`scDblFinder::amulet()`](https://bioconductor.org/packages/release/bioc/html/scDblFinder.html) detects likely scATAC-seq doublets from the number of genomic loci covered by more than two fragments. Its underlying `getFragmentOverlaps()` implementation filters fragment sizes and excluded regions, calculates per-barcode fragment and overlap counts, removes loci recurrently covered across many cells, and derives Poisson p-values with Benjamini-Hochberg correction.

**Reason for reimplementation.** The per-GEM-well workflow already stores Cell Ranger ATAC fragments as compressed BPCells directories. Passing the original fragment TSV to scDblFinder materializes chromosome-scale `GRanges` objects and previously requested six cores and 60 GB. The local implementation streams the existing BPCells fragment target and retains only one chromosome's selected fragments while calculating coverage runs.

**Behavior preserved.** [`calculate_amulet_metrics_BPCells()`](https://github.com/koefoeden/multiomeR/blob/main/R/amulet_BPCells_helpers.R) preserves barcode selection, minimum-fragment thresholds, maximum fragment size, excluded regions, `nFrags`, `uniqFrags`, `nAbove2`, `total.nAbove2`, p-values, q-values, and high-overlap-site removal. The lower-level loci return also preserves scDblFinder's cell-major, chromosome, and coordinate ordering. Cell Ranger's inclusive end-insertion convention is shifted back by one base before calculation so the BPCells representation matches scDblFinder's BED import.

**Deliberate deviations and consequences.** Only unique-fragment operation is supported. BPCells fragment objects do not retain Cell Ranger's PCR-duplicate count column, so requesting non-unique expansion fails explicitly instead of silently changing `nFrags`. BPCells does not export its fragment iterator header; the native helper therefore mirrors that private C++ interface, verifies the exact project-pinned BPCells commit `28759cdd5125` before use, and compiles a small shared library in each worker's temporary directory. A BPCells upgrade must revalidate this interface and the exact parity fixture before updating the pin. The target is single-threaded and requests the standard 16-GB worker tier. A native-only probe on the stored `healthy_PBMC_human` fragments processed 2,711 selected cells in 15.4 seconds with 0.64 GB peak RSS, including R startup and native compilation; this supports the reduced allocation for the public fixture but is not a memory guarantee for larger datasets.

**Implementation and wiring.** Native iteration and coverage-run calculation are implemented in [`src/amulet_bpcells.cpp`](https://github.com/koefoeden/multiomeR/blob/main/src/amulet_bpcells.cpp). The R wrapper, ABI check, high-overlap filtering, and AMULET statistics are in [`R/amulet_BPCells_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/amulet_BPCells_helpers.R). [`amulet_metrics_tibble` in `extra_targets/per_GEM_well_targets.R`](https://github.com/koefoeden/multiomeR/blob/main/extra_targets/per_GEM_well_targets.R) consumes the existing prefixed BPCells fragments, restores unprefixed barcode keys, and preserves the existing downstream metrics shape.

**Validation.** [`tests/testthat/test-amulet-parity.R`](https://github.com/koefoeden/multiomeR/blob/main/tests/testthat/test-amulet-parity.R) compares against scDblFinder 1.24.0 in disposable `callr` processes. It requires `identical()` results for the bundled fragment-file metrics, prefixed pipeline barcodes, a deterministic 12,000-fragment multi-chromosome loci fixture, and the corresponding full AMULET metrics. It also requires an explicit error for unsupported PCR-duplicate expansion. The disposable reference processes isolate stack-imbalance warnings emitted by the current scDblFinder reference under R 4.5 without weakening the returned-object comparison.

**Status and rerun.** Reference-parity tested against scDblFinder 1.24.0 with exact equality; passing in the current locked environment.

```bash
pixi run --use-environment-activation-cache test-amulet-parity
```

## BPCells-backed ATAC scDblFinder feature aggregation


**Reference algorithm.** [`scDblFinder::scDblFinder()`](https://bioconductor.org/packages/release/bioc/html/scDblFinder.html) can aggregate a high-dimensional ATAC count matrix before artificial-doublet classification. With `aggregateFeatures = TRUE`, scDblFinder performs its own TF-IDF-based feature clustering and sums peaks into the requested number of feature groups.

**Reason for adaptation.** Materializing and transforming every peak within each GEM well can exhaust worker memory before scDblFinder reaches classification. The pipeline instead derives 50 feature groups once from the aggregation's global LSI loadings, sums the disk-backed peak matrix with BPCells, and passes the compact matrix to scDblFinder with `aggregateFeatures = FALSE`.

**Behavior preserved.** This is not a reimplementation of scDblFinder's artificial-doublet classifier. Both paths use scDblFinder 1.24.0 with the same GEM-well barcodes, supplied biological cluster labels, `dbr.sd = 1`, 50 aggregated features, `processing = "normFeatures"`, serial BiocParallel execution, and returned score/class columns. Only the upstream construction of the 50-feature matrix changes.

**Deliberate deviations and consequences.** Global LSI-derived groups replace scDblFinder's per-GEM-well TF-IDF feature groups. The global groups are reusable across GEM wells and let BPCells aggregate before sparse-matrix materialization, but the resulting feature matrix is not expected to equal scDblFinder's internal aggregation. Doublet scores and calls can therefore differ; agreement must be evaluated at the final score and call level; exact parity is not expected.

**Implementation and wiring.** [`get_feature_groups_from_LSI_loadings()` and `aggregate_BPCells_rows_by_group()`](https://github.com/koefoeden/multiomeR/blob/main/R/processing_GEX_helpers.R) construct the compact feature matrix. [`extra_targets/ATAC_targets.R`](https://github.com/koefoeden/multiomeR/blob/main/extra_targets/ATAC_targets.R) computes the groups and aggregation once, then maps the unchanged scDblFinder classifier over GEM wells. The GEX path also uses the per-GEM-well wrapper but calls `scDblFinder::scDblFinder()` directly on GEX counts; it is a memory-bounding wrapper, not another algorithm reimplementation.

**Validation scope.** The public reference-parity fixtures do not establish equivalence for this alternative feature-aggregation path.

## Native weighted nearest neighbors

**Reference algorithm.** [`Seurat::FindMultiModalNeighbors()`](https://satijalab.org/seurat/reference/findmultimodalneighbors) constructs cell-specific modality weights from within- and cross-modality neighborhood prediction, collects candidate neighbors across modalities, and selects a weighted multimodal neighbor set.

**Reason for reimplementation.** The pipeline already has aligned RNA PCA/Harmony and ATAC LSI/Harmony matrices and needs reusable neighbor indices, distances, and modality weights without creating a Seurat object. Native graph state also feeds UMAP, Leiden clustering, SCAVENGE, and the optional Seurat/Signac export.

**Behavior preserved.** [`weighted_nearest_neighbors_BPCells()`](https://github.com/koefoeden/multiomeR/blob/main/R/processing_multimodal_helpers.R) preserves row-wise L2 normalization, per-modality nearest neighbors, nearest-distance correction, Seurat's small-SNN far-neighbour kernel bandwidth, within/cross prediction kernels, capped modality affinity ratios, normalized cell-specific modality weights, candidate-set union, weighted neighbor ranking, and Seurat's transformation from weighted affinity to neighbor distance.

**Deliberate deviations and consequences.** BPCells HNSW replaces Seurat's Annoy search, so approximate candidate sets need not be identical. The helper does not expose Seurat's optional smoothing or cross-constant list, and BPCells builds downstream SNN state rather than storing Seurat `Neighbor` and `Graph` objects. Its current `seed` argument is not consulted by the HNSW calls, so it must not be interpreted as controlling neighbor-search randomness. These choices can change weights, selected neighbors, SNN edges, clusters, and UMAP coordinates; correlation and overlap are therefore the validation contract, not exact equality.

**Implementation and wiring.** The implementation and graph consumers are in [`R/processing_multimodal_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/processing_multimodal_helpers.R); the tracked project-owned small-SNN kernel is in [`src/wnn_snn_bandwidth.cpp`](https://github.com/koefoeden/multiomeR/blob/main/src/wnn_snn_bandwidth.cpp). [`extra_targets/WNN_targets.R`](https://github.com/koefoeden/multiomeR/blob/main/extra_targets/WNN_targets.R) aligns modality embeddings, creates `WNN_results_raw`, filters small clusters, optionally recomputes `WNN_results`, and wires that state into WNN UMAP, clustering, metadata, and cell-type targets.

**Validation.** [`tests/testthat/test-wnn-parity.R`](https://github.com/koefoeden/multiomeR/blob/main/tests/testthat/test-wnn-parity.R) compares two deterministic RNA/ATAC fixtures with Seurat 5.5.0. The headline production-like fixture uses 400 cells, 12 dimensions, `k = 30`, and candidate range 200, matching the pipeline's configured neighbor count and native WNN candidate range. The pre-migration small-SNN pilot gave modality-weight Spearman 0.989 and mean neighbour-set overlap 0.996 on this fixture. On the 26,667-cell production object, it gave RNA/ATAC weight Spearman 0.988/0.988 and mean neighbour overlap 0.979. These figures supported the migration, but the maintained validation was deliberately not rerun as part of the production commit. The existing thresholds remain the acceptance contract. The current locked BPCells 0.3.1 build is pinned at `28759cdd5125`.

**Status and rerun.** The small-SNN default is based on the reference-similarity pilot above; post-migration validation remains pending. This does not assert exact equality of selected neighbours, SNN weights, clustering, or UMAP.

```bash
pixi run --use-environment-activation-cache test-algorithm-validation
```

## Sparse SCAVENGE propagation and significance

**Reference algorithm.** [SCAVENGE 1.0.2 at commit `8ee8b173d965`](https://github.com/sankaranlab/SCAVENGE/tree/8ee8b173d965009a696b2a590d5b17b28b7cf851) selects high chromVAR Z-score seed cells, constructs a binary mutual-nearest-neighbor adjacency graph, performs a column-normalized random walk with restart, caps and rescales the propagation score into a trait relevance score (TRS), and uses degree-matched seed permutations to identify significant cells.

**Reason for reimplementation.** The reference package's last commit and dependency stack predate the pipeline's current R/Bioconductor environment. multiomeR needs sparse propagation over native RNA PCA, ATAC LSI, and multimodal WNN SNN matrices and must avoid materializing a cell-by-permutation score matrix for large cell sets.

**Behavior preserved.** The implementation converts nonzero graph support to binary adjacency before analysis and preserves the one-sided Z-score seed threshold and top-percent cap, column-normalized transition matrix, equal seed restart mass, iterative random walk, 0.95 propagation-score cap, min-max scaling, Z-score scale factor, sequential base-R degree-matched seed sampling, and strict per-cell comparison with permuted propagation scores. The sampled seed-index lists are retained, but the cell-by-permutation score matrix is not: a native worker streams random walks and accumulates only per-cell exceedance counts and the cluster medians needed downstream. Random walks, rather than random-number generation, are parallelized, so the sampled null is invariant to the requested core count.

**Deliberate deviations and consequences.** The reference workflow constructs a mutual-kNN graph, whereas multiomeR uses the binary support of its BPCells-derived PCA, LSI, or WNN SNN graph; edge weights are discarded, but graph topology can still differ. Seed and scale-factor helpers guarantee at least one selected cell for small inputs. The degree sampler also handles a one-cell candidate stratum explicitly, avoiding base R's special interpretation of `sample(x, 1)` when `x` is one positive integer. The random walk validates graph inputs and has a maximum-iteration guard. Cell-level empirical P-values and significance calls follow the reference exceedance fraction and threshold. Cluster-level permutation medians, add-one P-values, and Benjamini--Hochberg adjustment within each grouping column are pipeline extensions.

**Implementation and wiring.** Seed selection, sparse random walk, streaming degree-matched permutations, TRS construction, and cluster-level null statistics are in [`R/SCAVENGE_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/SCAVENGE_helpers.R). [`module_genetic_enrichment/SCAVENGE_graph_targets.R`](https://github.com/koefoeden/multiomeR/blob/main/module_genetic_enrichment/SCAVENGE_graph_targets.R) constructs each graph and maps chromVAR Z-score records into `SCAVENGE_result_records`, from which cell-level TRS and cluster-level summaries are extracted; [`SCAVENGE_group_targets.R`](https://github.com/koefoeden/multiomeR/blob/main/module_genetic_enrichment/SCAVENGE_group_targets.R) combines summaries and plots.

**Validation.** [`tests/testthat/test-scavenge-parity.R`](https://github.com/koefoeden/multiomeR/blob/main/tests/testthat/test-scavenge-parity.R) uses a deterministic 60-cell fixture with repeated heterogeneous-degree graph blocks, nonuniform input edge weights, and three enriched seeds. The production helper receives the weighted graph, so the fixture also tests conversion to binary adjacency. First, the iterative sparse random walk is compared with the closed-form solution

\[
s = r\left(I - (1-r)P\right)^{-1}p_0,
\]

with a maximum absolute tolerance of 1e-10; the current delta is 8.61e-13. Second, compact local reference functions reproduce the relevant SCAVENGE 1.0.2 code at the pinned commit without installing its historical dependency stack. The random-walk delta against that reference is 1.11e-16, the transformed-score delta is 3.33e-16, and all 199 fixed-RNG degree-matched seed samples, streamed per-cell exceedance counts, empirical P-values, and significant-cell calls are identical. One- and two-core native results are also identical.

**Status and rerun.** Random-walk propagation is algorithmically derived against the closed form. Seed selection, binary propagation, transformed scores, sequential permutation sampling, streamed exceedance counts, empirical P-values, and significant-cell calls are reference-parity tested against the pinned source calculation; cluster summaries are documented pipeline extensions. This does not establish parity of mutual-kNN versus pipeline graph construction, chromVAR inputs, or biological interpretation.

```bash
pixi run --use-environment-activation-cache test-algorithm-validation
```

Run the complete maintained test suite, including the fast helper contracts, with:

```bash
pixi run --use-environment-activation-cache test
```


## Part: Background


<!-- source: website/implementation/background_philosophy.md -->

<!-- begin include: website/background_philosophy.md -->

# Why an editable workflow?

multiomeR keeps the analysis steps in an editable repository. Configuration
covers common choices such as inputs, markers, dimensions, and models; R
helpers and target definitions are available when a study needs a change
beyond those settings. This flexibility also means that users must review
which methods and assumptions fit their data.

## Reuse completed work

`targets` records dependencies between results so that a change can rebuild
the affected parts of an analysis. Independent tasks can run concurrently
when worker capacity permits. This is useful when processing several GEM
wells or repeating analyses with revised settings. The
[targets manual](https://books.ropensci.org/targets/) explains the execution
model and its limits.

## Keep large matrices on disk

BPCells provides disk-backed matrices and streaming operations that can
reduce the need to hold full matrices in memory. Some analysis steps still
need substantial RAM, and performance depends on the data, storage, and
available workers. See the [BPCells documentation](https://bnprks.github.io/BPCells/)
for its matrix operations and [Performance and scaling](../performance_overview.html)
for multiomeR examples.

## Keep the analysis inspectable

Separate targets make intermediate tables, matrices, and files available for
inspection. Seurat/Signac exports provide another way to explore completed
results. The [implementation conventions](implementation_conventions.md)
explain where to change parameters, helpers, and target definitions; the
[user manual](../) covers running an existing configuration.

<!-- end include: website/background_philosophy.md -->


## Part: Target graph views


<!-- source: website/implementation/implementation_main.md -->

<!-- begin include: website/implementation_main.md -->

# Main pipeline



The root `_targets.R` creates GEM-well and aggregation
mapping rows, then maps target fragments from `extra_targets/`. Use the diagrams
to find the relevant stage, then inspect the corresponding source file for the
complete command and resource declaration.

| Stage | Primary source |
|---|---|
| GEM well preprocessing | `extra_targets/per_GEM_well_targets.R` |
| Aggregation setup and shared QC | `extra_targets/general_aggregation_targets.R` |
| GEX | `extra_targets/GEX_merge_and_dim_reduc_targets.R`, `extra_targets/GEX_graph_and_cluster_targets.R` |
| ATAC | `extra_targets/ATAC_targets.R` |
| WNN | `extra_targets/WNN_targets.R` |
| Compatibility export | `extra_targets/Seurat_Signac_export_targets.R` |

## Parallel pre-processing

This view covers per GEM well processing and QC, including optional ambient RNA correction, donor demultiplexing, doublet detection, barcode filtering, and handoffs into aggregation-level GEX and ATAC objects.

[Mermaid graph omitted; source: `website/figures/human_curated/parallel_v2.mmd`]

## GEX processing

This view covers merged RNA processing, clustering, marker detection, cell type annotation, and GEX review outputs.

[Mermaid graph omitted; source: `website/figures/human_curated/GEX_v2.mmd`]

## ATAC processing

This view covers ATAC QC, peak calling, consensus peak construction, chromatin accessibility processing, chromVAR scoring, and coverage tracks. Peak–gene correlation is a separate optional module.

[Mermaid graph omitted; source: `website/figures/human_curated/ATAC_v2.mmd`]

## WNN integration

This view covers GEX and ATAC embedding handoffs, WNN integration, modality weights, cluster comparison, and integrated metadata outputs.

The native implementation, its differences from Seurat, and the maintained similarity thresholds are recorded in [Algorithmic implementations, deviations and validation](algorithm_validation.md#native-weighted-nearest-neighbors).

[Mermaid graph omitted; source: `website/figures/human_curated/WNN_v2.mmd`]

<!-- end include: website/implementation_main.md -->


<!-- source: website/implementation/implementation_differential_analyses.md -->

<!-- begin include: website/implementation_differential_analyses.md -->

# Differential analyses



`module_differential_analyses/targets.R` filters aggregations that enabled the module, joins their module config, attaches symbols for accepted WNN metadata and pseudobulk inputs, and maps the composition, pseudobulk, gene-set enrichment, and cross-modality target fragments. The generic pseudobulk model family is instantiated for gene expression, chromatin accessibility, motif-family accessibility, and expression-derived CollecTRI activity (transcription-factor activity). transcription-factor activity first converts filtered, normalized GEX pseudobulks to signed ULM scores and then reuses the same model and contrast machinery.

motif-family accessibility uses the 233 official JASPAR2026 CORE vertebrate familial root motifs as its complete feature universe. The pipeline scans those family-level profiles directly, rather than scanning individual motifs and taking the union of their peak matches.

The cross-modality fragment creates a CollecTRI-to-JASPAR family crosswalk, a detailed regulator-level table containing transcription-factor activity, motif-family accessibility, and TF-expression results, a family-level comparison table, a contrast-level concordance summary, and its plot. CollecTRI complexes remain intact in transcription-factor activity; complex-member mappings are introduced only by the comparison crosswalk.

The graph below is an orientation view. Inspect `setup_and_cell_type_composition_targets.R`, `pseudobulk_differential_targets.R`, `gene_set_enrichment_targets.R`, and `cross_modality_targets.R` for the complete model and plotting commands. The user-facing prerequisites and module selector are documented in [Differential analyses](../downstream_differential_analyses.html).

[Mermaid graph omitted; source: `website/figures/human_curated/differential_analyses_v2.mmd`]

## Method details

The transcription-factor activity branch infers signed TF or TF-complex activity from normalized GEX pseudobulks with CollecTRI regulons and the `decoupleR` univariate linear model (ULM). Genes are filtered for expression across cell-type pseudobulks, and each retained regulator must have at least five measured targets. Its inferred activities then use the same configured donor-level models and contrasts as gene expression, chromatin accessibility, and motif-family accessibility. The published human CollecTRI network is downloaded from the OmniPath rescue archive and accepted only when it matches the pipeline's pinned SHA-256 checksum.

motif-family accessibility tests the 233 sequence-similarity families in the official JASPAR2026 CORE vertebrate clustering rather than individual TF motifs. Each family is represented by its published root motif, which is scanned directly against the consensus peaks; individual member motifs are used only as family metadata. The same family-level accessibility matrix supports marker plots and the Seurat compatibility export. The CollecTRI-JASPAR comparison maps individual CollecTRI regulators to these JASPAR families and compares model t-statistics, not raw activity scales. AP1 and NFKB remain intact as complex regulons during activity inference; their canonical members are used only to associate the complexes with motif families for comparison. Detailed source-level results retain TF expression as a third reference, while family-level summaries use the median CollecTRI regulator t-statistic and report whether any mapped source is FDR-significant.

Each gene-set collection is tested independently with `cameraPR`, `inter.gene.cor = 0.01`, and a minimum of 10 genes represented in the contrast-specific universe. A significant set is more strongly associated with the contrast than the remaining tested genes, rather than merely showing any collective change. Open Targets evidence annotation is optional.

## Metadata used by the models

For each aggregation, the module retains its own file and full-tibble targets, then projects a canonical analysis view containing only donors in that aggregation and columns required by its configured models and composition plots. Rows are ordered by `donor_id` and non-key columns by name. Changes to unused columns, out-of-aggregation donors, or source row and column order therefore stop at this inexpensive projection boundary.

<!-- end include: website/implementation_differential_analyses.md -->


<!-- source: website/implementation/implementation_genetic_enrichment.md -->

<!-- begin include: website/implementation_genetic_enrichment.md -->

# Genetic enrichment



`module_genetic_enrichment/targets.R` filters enabled human aggregations, resolves one configured Open Targets study set per aggregation, and attaches symbols for WNN metadata, graphs, embeddings, consensus peaks, chromVAR state, and ATAC fragments.

The main target fragments live in `setup_targets.R`, `gchromVAR_targets.R`, `SCAVENGE_graph_targets.R`, `SCAVENGE_group_targets.R`, `GWAS_chromVAR_cell_type_targets.R`, and `GWAS_chromVAR_contribution_targets.R`. The user-facing release, method-selection, and interpretation contracts are documented in [Genetic enrichment](../downstream_genetic_enrichment.html).

## Single-nucleus and graph-based enrichment

This view covers the configured GWAS inputs, single-nucleus enrichment state, graph propagation, and downstream trait summaries. Additional cell-type contribution and locus-attribution branches may be pruned from this compact orientation view; use the manifest for the complete graph.

The sparse SCAVENGE reimplementation, deliberate graph and permutation differences, and validation evidence are recorded in [Algorithmic implementations, deviations and validation](algorithm_validation.md#sparse-scavenge-propagation-and-significance).

[Mermaid graph omitted; source: `website/figures/human_curated/genetic_enrichment_single_nucleus_v2.mmd`]

<!-- end include: website/implementation_genetic_enrichment.md -->


# Orphaned Markdown Pages

Tracked Markdown files not reached from the Quarto book graph or include graph.

- `website/data/README.md`
- `website/figures/human_curated/README.md`
- `website/helpers/_targets_graph_snippet.md`
