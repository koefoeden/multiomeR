# multiomeR Website LLM Export

This file is generated from the Quarto book outlines and resolves Quarto include shortcodes.
Hidden setup chunks, generated helper chunks, Mermaid graph bodies, and verbose image metadata are omitted by default.


# Book: multiomeR Manual


## Part: Overview


<!-- source: website/index.qmd -->

[Image omitted; source: `figures/multiomeR-logo.svg`; alt: Image]

# Introduction

## What is multiomeR?
multiomeR is a [targets](https://books.ropensci.org/targets/)-based workflow for processing and analyzing single-nucleus 10x Genomics Epi Multiome ATAC + Gene Expression datasets. It is meant to be adapted to your own data, compute setup, and biological questions.

The workflow starts from `cellranger-arc count` outputs, processes gene-expression (GEX) and ATAC data, builds multimodal aggregations, and supports optional downstream modules for differential analyses and genetic enrichment for human datasets.

## Development status
multiomeR is currently in beta. Expect breaking changes as the workflow is maturing. Report issues, questions, and feature requests on GitHub.

## Why use it?
- **Adaptable framework**: The repository is meant to be modified through config files, target selections, helper functions, and optional modules.
- **Transparent execution**: Each processing and analysis step is represented as an explicit `targets` target.
- **Efficient reruns**: `targets` reruns only the parts of the workflow affected by changed inputs, code, or config.
- **Multimodal analysis**: The main pipeline covers RNA/ATAC processing, clustering, cell typing, WNN integration, and optional downstream modules.

The framework assumes some familiarity with R and `targets`, but the public example workflow is intended to provide a concrete starting point.


## Workflow at a glance

multiomeR is organized around one main pipeline and two optional downstream modules.

The **main pipeline** starts from `cellranger-arc count` outputs, processes each 10x reaction, aggregates selected reactions, and builds multimodal RNA/ATAC outputs for clustering, cell typing, and WNN integration.

The **modules** are optional branches that run from configured aggregations. Enable them only when the relevant input data and biological question are present.

[Image omitted; source: `figures/multiomeR_overview_simplified.drawio.svg`; alt: Flowchart of the sample processing pipeline]


<!-- source: website/intro_manual_usage.qmd -->

# Using this manual

## Where to start

Use the [public example quickstart](demo_installation.qmd) if you want to install multiomeR, run the focused `immune_human_2x` example, and inspect representative outputs before reading the full pipeline walkthroughs.

Then use the workflow pages as references for the parts of the analysis you plan to run:

- **[Running your own data](main_inputs.qmd)**: Main-pipeline inputs, checkpoint-based execution, implementation diagrams, per-reaction QC, dataset-level merging, BPCells-native RNA/ATAC processing, clustering, cell typing, and WNN integration.
- **[Differential analyses module](downstream_differential_analyses.qmd)**: Optional aggregation module for datasets with experimental or phenotypic variables that should be modelled across cell type composition, pseudobulk expression, chromatin accessibility, and TF activity.
- **[Genetic enrichment module](downstream_genetic_enrichment.qmd)**: Optional human-data module for GWAS-informed chromatin accessibility enrichment and trait-relevance analyses.
- **[Implementation book](implementation/)**: Design philosophy, target graph methodology, and developer-facing details for readers modifying the workflow internals.

## Terminology
The following nomenclature is used in the documentation and code:

- **'10x reaction' (or possibly just 'reaction')**: The experimental/computational outcome of a single lane on the 10x Genomics Controller chip, which might contain nuclei from multiple samples and/or multiple donors/individuals. One 10x reaction thus corresponds to a single output directory from [cellranger-arc count](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/analysis/running-pipelines/single-library-analysis).
- **donor**: A single individual, subject, patient, participant, etc. that has a unique donor_id and specific genotype.
- **dataset**: A group of one or more reactions that share dataset-level processing settings, such as reference genome and reaction-level QC filters.
- **aggregation**: A configured joint analysis of one or more reactions, including merged RNA/ATAC processing, clustering, cell typing, and optional downstream modules.


## Part: Quickstart: Public PBMC data


<!-- source: website/demo_installation.qmd -->

# Installation



## System requirements

- Linux system with at least 60 GB of RAM, preferably equipped with a job scheduler supported by the `crew.cluster` package: SLURM, PBS, SGE, or LSF.


## Clone and enter the repository

```{.bash filename="Bash"}
git clone https://github.com/koefoeden/multiomeR.git
cd multiomeR
```

## Download public inputs

Download the two public human `cellranger-arc count` outputs used by the `immune_human_2x` example aggregation from the repository root (3.9 GB).

```{.bash filename="Bash"}
pbmc_dir="example_data/healthy_PBMC_human/outs"
lymphoma_dir="example_data/lymphoma_lymph_human/outs"

mkdir -p $pbmc_dir $lymphoma_dir

pbmc_prefix="https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k"
curl -fL -o "${pbmc_dir}/web_summary.html" "${pbmc_prefix}_web_summary.html"
curl -fL -o "${pbmc_dir}/summary.csv" "${pbmc_prefix}_summary.csv"
curl -fL -o "${pbmc_dir}/per_barcode_metrics.csv" "${pbmc_prefix}_per_barcode_metrics.csv"
curl -fL -o "${pbmc_dir}/filtered_feature_bc_matrix.h5" "${pbmc_prefix}_filtered_feature_bc_matrix.h5"
curl -fL -o "${pbmc_dir}/raw_feature_bc_matrix.h5" "${pbmc_prefix}_raw_feature_bc_matrix.h5"
curl -fL -o "${pbmc_dir}/atac_fragments.tsv.gz" "${pbmc_prefix}_atac_fragments.tsv.gz"
curl -fL -o "${pbmc_dir}/atac_fragments.tsv.gz.tbi" "${pbmc_prefix}_atac_fragments.tsv.gz.tbi"

lymphoma_prefix="https://cf.10xgenomics.com/samples/cell-arc/2.0.0/lymph_node_lymphoma_14k/lymph_node_lymphoma_14k"
curl -fL -o "${lymphoma_dir}/web_summary.html" "${lymphoma_prefix}_web_summary.html"
curl -fL -o "${lymphoma_dir}/summary.csv" "${lymphoma_prefix}_summary.csv"
curl -fL -o "${lymphoma_dir}/per_barcode_metrics.csv" "${lymphoma_prefix}_per_barcode_metrics.csv"
curl -fL -o "${lymphoma_dir}/filtered_feature_bc_matrix.h5" "${lymphoma_prefix}_filtered_feature_bc_matrix.h5"
curl -fL -o "${lymphoma_dir}/raw_feature_bc_matrix.h5" "${lymphoma_prefix}_raw_feature_bc_matrix.h5"
curl -fL -o "${lymphoma_dir}/atac_fragments.tsv.gz" "${lymphoma_prefix}_atac_fragments.tsv.gz"
curl -fL -o "${lymphoma_dir}/atac_fragments.tsv.gz.tbi" "${lymphoma_prefix}_atac_fragments.tsv.gz.tbi"
```

Then download and unpack the matching human Cell Ranger ARC reference (14.9 GB). If you already have this reference, symlink it to `example_data/refdata-cellranger-arc-GRCh38-2024-A` instead.

```{.bash filename="Bash"}
refdata="refdata-cellranger-arc-GRCh38-2024-A"
refdata_url="https://cf.10xgenomics.com/supp/cell-arc"

curl -fL -o "example_data/${refdata}.tar.gz" "${refdata_url}/${refdata}.tar.gz"
tar -xzf "example_data/${refdata}.tar.gz" -C example_data
```

## Create and launch the pixi environment

If `pixi` is not installed yet, install it first:

```{.bash filename="Bash"}
curl -fsSL https://pixi.sh/install.sh | sh
export PATH="$HOME/.pixi/bin:$PATH"
```

Install the locked pixi environment from the repository root (<5 mins):

```{.bash filename="Bash"}
pixi install --locked --run-post-link-scripts
pixi run install-r-github-packages
```

```{.bash filename="Bash"}
pixi shell
R
```

Use this repository-root R session as the evaluation environment for subsequent R commands in this manual.

::: {.callout-note collapse="true"}
## Why pixi?

multiomeR depends on both R packages and command-line genomics tools. Pixi keeps those dependencies in one locked environment, so R, Bioconductor packages, and external tools such as `macs3`, `tabix`, and `cellsnp-lite` are resolved together instead of being installed manually across separate package managers.
:::

::: {.callout-note collapse="true"}
## Why is `--run-post-link-scripts` needed?

The `--locked` flag keeps installation aligned with `pixi.lock`. The `--run-post-link-scripts` flag is required for Bioconductor data packages such as the BSgenome genome packages. Pixi disables conda post-link scripts by default, but these packages use Bioconda's post-link scripts to materialize the actual R package directories in the pixi environment.
:::

::: {.callout-note collapse="true"}
## Why is `pixi run install-r-github-packages` needed?

A small number of R packages still need GitHub versions, since they are not currently available on conda-forge/bioconda. The `pixi run install-r-github-packages` task installs the pinned GitHub versions of BPCells, Signac, and betterChromVAR into the pixi R library. The pixi environment also preinstalls system libraries needed by these source installs, including `libhwy` for BPCells.
:::


<!-- source: website/demo_running.qmd -->

# Run the main pipeline



Run the command below to process the example `immune_human_2x` aggregation and produce a multimodal Seurat object. This object is primarily a convenience export: the pipeline does not use Seurat for the core processing, and instead keeps BPCells-backed matrices as the primary state. Using 16 logical threads on an AMD EPYC 7543, this target selection takes <25 minutes to run and writes about 6 GB of data to disk.

```{.r filename="R"}
targets::tar_make(
  tidyselect::matches("multimodal_Seurat_object.immune_human_2x")
  )
```


<!-- source: website/demo_outputs.qmd -->

# Understanding the outputs



The public `immune_human_2x` demo writes its `targets` store to `outputs`, as configured in `_targets.yaml`. That store contains three kinds of output that are useful to inspect after the run:

- `objects/`: serialized R objects managed by `targets`, returned by `targets::tar_target()` targets
- `files/`: file and directory artifacts saved by `get_structured_file_path()` inside `tarchetypes::tar_file()` targets
- `plots/`: PNG and SVG review plots saved by `save_plots_structured()` inside `tarchetypes::tar_file()` targets

## Objects

Most intermediate and final result objects are saved automatically by `targets` under `outputs/objects` during a pipeline run. Read them with `targets::tar_read()` from a repository-root R session after the demo has completed.

```{.r filename="R"}
targets::tar_read(UMAP_embeddings_tibble.WNN.immune_human_2x)
targets::tar_read(metadata_w_cell_types_tibble.WNN.immune_human_2x)
targets::tar_read(cell_type_marker_tibbles.GEX.immune_human_2x)
```

See the official [`tar_read()` documentation](https://docs.ropensci.org/targets/reference/tar_read.html) for details on reading target values from storage.

## Files

Other targets in the pipeline explicitly write files to disk inside `tarchetypes::tar_file()` targets. These targets can still be loaded by `targets::tar_read()`, but their value is simply the file path rather than an in-memory analysis object.

```{.r filename="R"}
targets::tar_read(cellranger_barcodes_tsv.healthy_PBMC_human)
targets::tar_read(aggregated_GEX_BPCells_matrix_dir.GEX.immune_human_2x)
targets::tar_read(consensus_peak_BED.ATAC.immune_human_2x)
```

Internally, these paths are created with `get_structured_file_path()`. The helper derives a stable location from the active target name and writes it under `outputs/files`. Target name parts after dots become parent directories in reverse order, so a target like `aggregated_GEX_BPCells_matrix_dir.GEX.immune_human_2x` is placed under:

```text
outputs/files/immune_human_2x/GEX/
```

After the public demo run, you should therefore find produced files ordered separated by individual reactions, datasets, and the aggregation name:

```text
outputs/files/
|-- healthy_PBMC_human/ (reaction 1 specific files)
|-- lymphoma_lymph_human/ (reaction 2 specific files)
`-- immune_human_2x/ (aggregation-level files)
```

See the official [`targets` local data documentation](https://books.ropensci.org/targets/data.html#external-files) for more information on file targets.

## Plots

Plots are saved in a very similar manner to files. The `save_plots_structured()` helper function derives paths from the target name, like `get_structured_file_path()`, but writes under `outputs/plots` instead:
```text
outputs/plots/
|-- healthy_PBMC_human/ (reaction 1 specific plots)
|-- lymphoma_lymph_human/ (reaction 2 specific plots)
|-- immune_human_dataset/ (dataset-level plots)
`-- immune_human_2x/ (aggregation-level plots)
```

See the [main pipeline output gallery](gallery_main.qmd) for selected plots.


## Part: Configure and run your own data


<!-- source: website/main_overview.qmd -->

# Introduction to the main pipeline



## Purpose

This workflow is defined in `_targets.R` and represents the main analysis pipeline after cellranger-arc count output generation. It processes individual cellranger-arc count outputs (reactions) in parallel, combines selected reactions into specified aggregations, and runs several multimodal workflows.

Representative outputs from the public example workflow are collected in the [main pipeline output gallery](gallery_main.qmd).

## Features

- **Genotype-aware demultiplexing and metadata integration**: per-reaction barcode extraction, optional [`cellsnp-lite`](https://cellsnp-lite.readthedocs.io/) and [`vireo`](https://vireosnp.readthedocs.io/) donor assignment, donor/reaction metadata validation, and aggregation-level metadata joins.
- **Extensive configurable QC**: RNA count/feature/mitochondrial/novelty metrics, ATAC TSS enrichment and nucleosome signal, Cell Ranger-only and full-metadata exclusion sets, [`AMULET`](https://github.com/UcarLab/AMULET) ATAC doublet metrics, GEX/ATAC [`scDblFinder`](https://www.bioconductor.org/packages/release/bioc/html/scDblFinder.html), and interpretable UpSet/violin/checkpoint plots across reaction, dataset, and aggregation levels.
- **BPCells-native GEX processing**: optional CellBender input, BPCells-backed count matrices, SCTransform/PCA, [`Harmony`](https://portals.broadinstitute.org/harmony/) correction, UMAP, Leiden clustering, marker-gene module scoring, cell-type labelling, cluster marker testing, and marker/metadata review plots.
- **BPCells-native ATAC processing**: fragment merging, genome blacklist handling, cluster-aware peak calling with MACS3 or BPCells tiling, consensus peak construction, peak annotation, BPCells peak matrices, TF-IDF/LSI, Harmony correction, UMAP, Leiden clustering, differential accessibility, motif/TF activity summaries, and coverage-track plots.
- **Multimodal integration**: native WNN construction from aligned GEX PCA and ATAC LSI embeddings, WNN clustering and UMAPs, RNA/ATAC modality-weight summaries, cross-modality cluster comparison plots, and final multimodal metadata.
- **Regulatory feature outputs**: consensus peak tables, marker TF activity, differential accessibility BEDs, ArchR-style peak-gene correlation candidates/results, and pseudobulk GEX/ATAC matrices for downstream modules.
- **Compatibility exports and focused reprocessing**: GEX-only review objects and final [`Seurat`](https://satijalab.org/seurat/)/[`Signac`](https://stuartlab.org/signac/) multimodal objects, while keeping BPCells-native targets as the primary pipeline state; optional subgroup reprocessing reruns GEX, ATAC, and WNN analyses within large parent cell groups.


<!-- source: website/main_inputs.qmd -->

# Configuration and inputs



## Three configuration layers

The main pipeline is configured in three layers: 
- **datasets** define shared genome and QC settings for a number of reactions
- **reactions** point to individual cellranger-arc count outputs
- **aggregations** choose which reactions to process together.

Defaults, missing-value rules, data types, and allowed values for YAML-backed configuration live in `cfg_pipeline_parameters.tsv`. Dataset, aggregation, and module YAML files contain concrete config rows and optional row-to-row inheritance with `inherits:`.

## `cfg_datasets.yaml`

Dataset-level settings live in `cfg_datasets.yaml`. A dataset is the shared context for one or more 10x reactions, typically a genome/reference combination plus dataset-wide QC rules. Reactions in `cfg_reactions.tsv` refer to these dataset keys, so define datasets before adding reaction rows.

In the public repository, `cfg_datasets.yaml`, `cfg_reactions.tsv`, and `cfg_aggregations.yaml` point to template files. Keep those public example entries as a reference, and use project-specific copies or downstream private config files for study-specific paths, metadata, and analysis settings.

Defaults and validation rules for dataset, aggregation, and module YAML parameters live in `cfg_pipeline_parameters.tsv`. A YAML row only needs to list values that differ from the manifest defaults, plus any parameters where `allow_missing_after_inheritance` is `FALSE` and the manifest default is `NULL`.

The manifest-backed overview below is the reference for dataset parameters. Open a topic to inspect defaults, data types, allowed values, and examples.

[Generated Quarto chunk omitted: `emit_parameter_overview("dataset")`]

The public human example uses the following dataset entry, rendered directly from `cfg_datasets.yaml`:

<details>
<summary>Show <code>immune_human_dataset</code></summary>

[Generated Quarto chunk omitted: `emit_yaml_template_entry(datasets_config_file, "immune_human_dataset")`]

</details>

## `cfg_reactions.tsv`

Reaction-level configuration lives in `cfg_reactions.tsv`. Add one row per 10x reaction with a stable `reaction_ID`, the matching `dataset`, the `reaction_cellranger_count_dir`, and any donor demultiplexing fields relevant to that reaction.

Each reaction should point to a [`cellranger-arc count`](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/analysis/running-pipelines/single-library-analysis) output directory. The pipeline expects the configured `reaction_cellranger_count_dir` to contain the standard `outs/` files, including:

```{.text filename="Text"}
<reaction_cellranger_count_dir>/outs/
|-- raw_feature_bc_matrix.h5
|-- filtered_feature_bc_matrix.h5
|-- atac_fragments.tsv.gz
|-- atac_fragments.tsv.gz.tbi
|-- atac_possorted_bam.bam
|-- per_barcode_metrics.csv
```

[Generated Quarto chunk omitted: `reaction_template <- readr::read_tsv("cfg_reactions.tsv", show_col_types = FALSE) cat('<div class="scrollable-table">...`]

The reaction table columns are:

| Column | Meaning |
|---|---|
| `reaction_ID` | Stable reaction identifier used throughout target names and downstream config. |
| `dataset` | Dataset key used to join dataset-level settings from `cfg_datasets.yaml`. |
| `reaction_donor_id` | Donor ID assigned to all nuclei in a non-multiplexed reaction. Use `NA` for multiplexed reactions. |
| `reaction_n_donors` | Number of donors expected in the reaction. |
| `reaction_cellranger_count_dir` | Path to the `cellranger-arc count` output directory. |
| `reaction_cellbender_h5_file` | Optional CellBender-corrected GEX H5 file. Use `NA` when not used. |
| `reaction_donors_VCF_file` | Optional VCF for donor demultiplexing. Use `NA` when not used. |

## `cfg_aggregations.yaml`

Aggregation-level settings live in `cfg_aggregations.yaml`. An aggregation combines one or more configured reactions and controls merged RNA/ATAC processing, clustering, annotation, plotting, and optional downstream modules.

The manifest-backed overview below is the reference for aggregation parameters. Required rows must resolve to a non-missing value after inheritance and defaults are applied; optional rows can remain `NULL` when the corresponding behavior should be disabled.

[Generated Quarto chunk omitted: `emit_parameter_overview("aggregation")`]

Each aggregation also points to metadata files:

- **`aggregation_donor_id_metadata_tsv`**: donor-level metadata with a `donor_id` column matching the donor IDs assigned during per-reaction processing.
- **`aggregation_reaction_ID_metadata_tsv`**: reaction-level metadata with a `TENX_reaction_ID` column matching the configured reaction IDs.

These metadata tables are joined to the per-cell metadata before downstream correction, plotting, and modelling.

The public `immune_human_2x` aggregation below is rendered directly from `cfg_aggregations.yaml`. It shows the required aggregation inputs plus the example-specific optional overrides; standard optional settings come from `cfg_pipeline_parameters.tsv`.

<details>
<summary>Show <code>immune_human_2x</code></summary>

[Generated Quarto chunk omitted: `emit_yaml_template_entry(aggregations_config_file, "immune_human_2x")`]

</details>


<!-- source: website/main_running.qmd -->

# Running the workflow

From the repository-root R session described in [Installation](demo_installation.qmd), run the main workflow in checkpoint-sized pieces. Replace `your_aggregation` with the aggregation name from `cfg_aggregations.yaml`, for example `immune_human_2x`.

## GEX checkpoint

Run the RNA-processing checkpoint first:

```{.r filename="R"}
targets::tar_make(
  names = targets::tar_described_as(tidyselect::contains("checkpoint:GEX")) &
    tidyselect::matches("your_aggregation")
)
```

At this checkpoint, inspect the aggregation-level barcode exclusion summaries, RNA-only QC plots, marker summaries, and `GEX_Seurat_object.your_aggregation`. If the GEX filtering, clustering, and annotation look acceptable, continue to the ATAC checkpoint.

## ATAC checkpoint

Run the ATAC-processing checkpoint after the GEX review:

```{.r filename="R"}
targets::tar_make(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:ATAC")
    ) &
    tidyselect::matches("your_aggregation")
)
```

At this checkpoint, inspect ATAC QC plots, LSI diagnostics, peak and motif summaries, coverage tracks, and TF or gene activity plots. If the ATAC processing looks acceptable, continue to the multimodal checkpoint.

## Multimodal checkpoint

```{.r filename="R"}
targets::tar_make(
  names = targets::tar_described_as(tidyselect::contains("checkpoint:multimodal")) &
    tidyselect::matches("your_aggregation")
)
```

At this checkpoint, inspect WNN UMAPs, modality-weight plots, multimodal cluster comparisons, and `multimodal_Seurat_object.your_aggregation`. Optional downstream modules are described in their own walkthroughs.

To run all configured aggregations for a checkpoint, omit the `& tidyselect::matches("your_aggregation")` part. To run the full workflow after checkpoint review, run:

```{.r filename="R"}
targets::tar_make()
```

The full main pipeline can be computationally expensive. Start with the public example or a narrow target selection, then scale to larger datasets after local paths, metadata files, and controller settings are working.

The workflow output directory can be changed with the `store` field in [`_targets.yaml`](https://github.com/koefoeden/multiomeR/tree/main/_targets.yaml).


## Part: Optional modules


<!-- source: website/downstream_differential_analyses.qmd -->

# Differential Analyses



## Purpose

This optional module is relevant when an aggregation has experimental, clinical, phenotypic, or technical variables that should be modelled across processed Multiome data. It is defined inside the root `_targets.R` workflow, with module-specific target fragments and configuration in `module_differential_analyses/`.

Representative outputs from this module are shown in the [differential analyses output gallery](gallery_differential_analyses.qmd).

## Features

The differential analyses module extends a completed aggregation with configurable statistical models for composition and pseudobulk feature testing.

- **Differential modelling**
  - Cell-type composition models for phenotype, clinical, technical, or experimental variables.
  - Pseudobulk differential gene expression, chromatin accessibility, and transcription-factor activity.
  - Model fitting with [`edgeR`](https://bioconductor.org/packages/release/bioc/html/edgeR.html) and [`limma`](https://bioconductor.org/packages/release/bioc/html/limma.html).

- **Configuration-driven contrasts**
  - Dynamic model branches defined from module configuration.
  - Shared donor metadata and modality-specific filtering.
  - Reusable model specifications across GEX, ATAC, and TF-activity analyses.

- **Result interpretation**
  - Pseudobulk depth summaries and P-value density diagnostics.
  - Volcano plots and significant-element summaries.
  - Cross-modality summaries comparing significant signals across data layers.

- **Biological annotation**
  - GSEA on DGE results using selected [`MSigDB`](https://www.gsea-msigdb.org/gsea/msigdb) collections and optional custom gene sets.
  - Optional links between top features and [`Open Targets`](https://platform.opentargets.org/) GWAS evidence.

## Configuration

Enable the module for an aggregation in `cfg_aggregations.yaml`, then configure the module-specific model settings in `module_differential_analyses/cfg.yaml`. Only edit module rows for aggregations that enable this module.

For example, this aggregation setting enables the module:

```{.yaml filename="YAML"}
modules: [differential_analyses]
```

The module config controls which variables are modelled, which metadata columns are used, and which contrasts are run. The manifest-backed overview below is the reference for differential-analysis module parameters.

[Generated Quarto chunk omitted: `emit_parameter_overview("differential_analyses")`]

The public example below is rendered directly from `module_differential_analyses/cfg_template.yaml`:

<details>
<summary>Show <code>immune_human_2x</code></summary>

[Generated Quarto chunk omitted: `emit_yaml_template_entry(module_cfg_template_file, "immune_human_2x")`]

</details>

Edit `module_differential_analyses/cfg.yaml`, keeping disabled or example rows available for reference until the corresponding aggregation is ready.

## Running the pipeline

From the repository-root R session described in [Installation](demo_installation.qmd), first complete the main pipeline checkpoints for the aggregation. Then run the module checkpoint. Replace `your_aggregation` with an aggregation that enables the module in `cfg_aggregations.yaml`.

### Differential analyses checkpoint

```{.r filename="R"}
targets::tar_make(
  names = targets::tar_described_as(tidyselect::contains("checkpoint:differential_analyses")) &
    tidyselect::matches("your_aggregation")
)
```

This checkpoint writes the differential cell-type composition plots, pseudobulk model diagnostics, volcano plots, GSEA summaries, and cross-modality summaries configured for the aggregation. Runtime depends on the number of modelled variables, cell types, donors, and pseudobulk contrasts.


<!-- source: website/downstream_genetic_enrichment.qmd -->

# Genetic Enrichment



## Purpose

This optional module is relevant for human-data aggregations where GWAS or fine-mapped trait evidence should be integrated with chromatin accessibility data. It is defined inside the root `_targets.R` workflow, with module-specific target fragments and configuration in `module_genetic_enrichment/`.

Representative outputs from this module are shown in the [genetic enrichment output gallery](gallery_genetic_enrichment.qmd).

## Features

The genetic enrichment module links processed human multiome aggregations to fine-mapped GWAS and Open Targets evidence.

- **GWAS input handling**
  - Configured study metadata, credible sets, ancestry/sample metadata, and annotation tracks.
  - Posterior-probability weighting of credible-set variants.
  - Mapping of genetic evidence onto consensus ATAC peaks.

- **Single-nucleus enrichment**
  - [`chromVAR`](https://greenleaflab.github.io/chromVAR/articles/Introduction.html) / [`gchromVAR`](https://caleblareau.github.io/gchromVAR/reference/index.html) Z-scores for configured traits.
  - Summaries by clusters, donors, and cell types.
  - Heatmaps, boxplots, and metadata-linked review plots.

- **Graph-based enrichment**
  - [`SCAVENGE`](https://github.com/sankaranlab/SCAVENGE) trait-score propagation across PCA, LSI, and WNN graph representations.
  - TRS UMAPs, enrichment heatmaps, cluster-level boxplots, and significant-cell proportion summaries.

- **Cell-type enrichment and attribution**
  - GWAS chromVAR deviations from ATAC counts summed over each cell type, with nuclei support shown alongside the heatmap.
  - Exact decomposition of each deviation into peak, credible-variant, and locus contributions.
  - Relative-deviation heatmaps, locus waterfalls, and variant-level accessibility tracks.

## Configuration

Enable the module for a human-data aggregation in `cfg_aggregations.yaml`, then configure Open Targets studies and module-specific settings in `module_genetic_enrichment/cfg.yaml`. Only edit module rows for aggregations that enable this module.

For example, this aggregation setting enables the module:

```{.yaml filename="YAML"}
modules: [genetic_enrichment]
```

The module config controls which Open Targets studies are used, how enrichment inputs are selected, and which trait-relevance summaries are produced. The manifest-backed overview below is the reference for genetic-enrichment module parameters.

[Generated Quarto chunk omitted: `emit_parameter_overview("genetic_enrichment")`]

The public example below is rendered directly from `module_genetic_enrichment/cfg_template.yaml`:

<details>
<summary>Show <code>immune_human_2x</code></summary>

[Generated Quarto chunk omitted: `emit_yaml_template_entry(module_cfg_template_file, "immune_human_2x")`]

</details>

Edit `module_genetic_enrichment/cfg.yaml`, keeping disabled or example rows available for reference until the corresponding human-data aggregation is ready.

## Running the pipeline

From the repository-root R session described in [Installation](demo_installation.qmd), first complete the main pipeline checkpoints for the human-data aggregation. Then run the module checkpoint. Replace `your_aggregation` with an aggregation that enables the module in `cfg_aggregations.yaml`.

### Genetic enrichment checkpoint

```{.r filename="R"}
targets::tar_make(
  names = targets::tar_described_as(tidyselect::contains("checkpoint:genetic_enrichment")) &
    tidyselect::matches("your_aggregation")
)
```

This checkpoint writes the configured cell-type, attribution, and single-nucleus trait-relevance summaries, including GWAS chromVAR and SCAVENGE outputs. Runtime depends on the number of configured GWAS studies, enrichment methods, cells, and attributed loci.


## Part: Performance


<!-- source: website/performance_overview.qmd -->

# Performance overview

## Performance results

The initial example benchmark shows sub-linear scaling across the current 2- and 6-reaction example aggregations: the larger CellRanger input increases estimated wall time by only ~1.4-fold. In practice, this means that larger projects should benefit substantially from the pipeline's reaction-level and branch-level parallelism when sufficient compute resources are available, rather than scaling like a fully serial workflow.

[Image omitted; source: `../outputs/benchmark/multimodal_seurat_walltime.png`; alt: Estimated wall time to each aggregation's final `multimodal_Seurat_object` target by total CellRanger input nuclei. S...]


<!-- source: website/performance_distributed_computing.qmd -->

# Distributed computing



multiomeR runs targets through the [targets](https://books.ropensci.org/targets/crew.html) and `crew` backend. The public controller template, `crew_controllers_template.R`, defines a local-controller setup that is sufficient for the public demo.

For larger datasets, adapt your active `crew_controllers.R` to match your scheduler, queue names, resource limits, and worker startup requirements.

## Controller contract

`crew_controllers.R` must return a named list with:

- `controller_list`: the list of `crew` controllers.
- `controller_resources_tibble`: a table with `controller_name`, `cores`, `RAM_GB`, and `gpus`.

Targets request resources by cores, RAM, and optionally GPUs through `get_tar_resources()`. The first row of `controller_resources_tibble` is the default controller, and row order controls which matching controller is preferred. GPU controllers are excluded unless `gpus_req > 0`.

## Local default

The default template uses a local `crew` controller:

```{.bash filename="Bash"}
crew_controllers.R -> crew_controllers_template.R
```

Increase local workers, cores, and memory only after the example workflow succeeds on your machine.

## Cluster execution

For cluster execution, replace the local controller with a scheduler-backed controller such as `crew.cluster::crew_controller_slurm()`. Keep the `controller_resources_tibble` rows aligned with the controller names, because that table controls how target-level resource requests are routed.


## Part: Output gallery


<!-- source: website/gallery_main.qmd -->

# Main Pipeline



This gallery shows representative outputs produced by multiomeR on the public `immune_human_2x` demo. It is meant to give a quick sense of the QC summaries and multimodal embeddings that the main pipeline can produce before you run it on your own data.

Each card includes the corresponding `targets` target name so you can connect the displayed output back to the workflow.

[Generated Quarto chunk omitted: `render_gallery_section(gallery_items, "Main pipeline")`]


<!-- source: website/gallery_differential_analyses.qmd -->

# Differential Analyses



[Generated Quarto chunk omitted: `render_gallery_section(gallery_items, "Differential analyses module")`]


<!-- source: website/gallery_genetic_enrichment.qmd -->

# Genetic Enrichment



[Generated Quarto chunk omitted: `render_gallery_section( gallery_items, "Genetic enrichment module", subsection_descriptions = c( "Single-nucleus chro...`]


## Part: Reference


# Book: multiomeR Implementation


<!-- source: website/implementation/index.qmd -->

# Introduction

This book collects developer-facing details for readers who want to understand or modify the internals of multiomeR. It complements the [main user manual](../), which covers installation, configuration, checkpoint-based execution, and output inspection.

Use this book when you need to reason about why the workflow is structured the way it is, how the simplified target graph views relate to the real `{targets}` graph, or where a change should fit into the main pipeline and optional downstream modules.

## Where to start

- **[Background and design philosophy](background_philosophy.qmd)**: Why multiomeR is organized as a configurable, editable workflow rather than a conventional high-level package interface.
- **[Implementation conventions](implementation_conventions.qmd)**: The metadata, configuration, mapping, symbol, and runtime contracts used throughout the target graph.
- **[Reading the graph views](graph_methodology.qmd)**: How the implementation diagrams are generated, simplified, and intended to be interpreted.
- **[Main pipeline](implementation_main.qmd)**: Graph views for reaction-level preprocessing, GEX, ATAC, WNN integration, and subgroup reprocessing.
- **[Differential analyses](implementation_differential_analyses.qmd)**: Graph view for the optional differential-analysis module.
- **[Genetic enrichment](implementation_genetic_enrichment.qmd)**: Graph views for the optional genetic-enrichment module.

If you are trying to run multiomeR rather than modify it, start with the [main manual](../).


## Part: Background


<!-- source: website/implementation/background_philosophy.qmd -->

<!-- begin include: website/background_philosophy.qmd -->

# Background and design philosophy

## Philosophy

multiomeR challenges the commonly used single-cell & single-nucleus processing workflows, where most processing steps are run sequentially and are implemented in hard-to-customize packaged functions. Widely used examples include brilliant packages such as Seurat and Signac. However, while offering superb ease-of-use, this traditional workflow has two important downsides:

1) The built-in processing steps (usually functions) are only customizable to levels chosen by the package developers, which often falls short in the face of the enormous variance and complexity of real-life single-cell and single-nucleus sequencing datasets
2) The traditional workflows leave a lot of performance on the table, since they rarely take advantage of the inherent parallelism when processing these large sequencing datasets. This includes commonly found multiplicities of independent experimental samples (donors, reactions), data modalities, data representations, downstream analyses, etc.

Until recently, this simplified and easy-to-use workflow has been necessary to make the tools practically usable in a busy research environment. However, with the recent advent of powerful, AI coding agents, which can quickly ingest large complex codebases, this inherent tug-of-war between performance/customizability and ease-of-use/code interpretability has started to move the needle in favor of complex, performant, and customizable approaches. This is especially relevant because single-cell and single-nucleus datasets continue to grow faster than routine compute capacity. These observations, together with a lack of an established framework for the processing of Multiome data, constitute the main motivation for the development of this pipeline.

## Achieving customizability

By not packaging the pipeline into an ordinary R package, users retain straightforward access to the codebase and can modify any part of the workflow to fit their needs. The pipeline thus functions as a template of suggested best practices, which can be used as-is, modified, or extended as needed.

## Achieving performance

multiomeR achieves its high performance through two key design choices:

1) The workflow is implemented as a `targets` pipeline, which explicitly represents each processing and analysis step as a node, and the dependencies between them as edges, in a directed acyclic graph (DAG). This allows `targets` to automatically determine which parts of the workflow can be run in parallel, and which parts need to run sequentially - and whether they need to run at all through the use of cached data and intelligent invalidation.
2) multiomeR builds heavily on the BPCells R package, which is designed from the ground up for computational efficiency through lazy, streaming data processing and on-disk storage. Please see the BPCells documentation for more details.

<!-- end include: website/background_philosophy.qmd -->


<!-- source: website/implementation/implementation_conventions.qmd -->

# Implementation conventions

The implementation book assumes a few repository conventions that are easy to miss when reading individual target files in isolation. These conventions are not separate framework features; they are the small contracts that make the `{targets}` graph configurable, inspectable, and readable across reaction-level preprocessing, aggregation-level analysis, and optional downstream modules.

This chapter records the conventions that currently carry the most structural weight: target description tags, the parameter manifest, mapping tibbles, target-symbol columns, and the runtime bootstrap contract.

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

`[checkpoint:<name>]` marks targets that users can request as review boundaries with `targets::tar_described_as()`. These are practical handles for running or inspecting major checkpoints such as GEX, ATAC, multimodal, differential analyses, or genetic enrichment. They are not dependency groups and they are not parsed with a strict tag API; the user-facing commands currently match them by description substring.

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
dataset
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
  aggregation_reaction_IDs: [healthy_PBMC_human, lymphoma_lymph_human]
  aggregation_GEX_marker_genes:
    B: [MS4A1, CD79A]
    T: [TRAC, CD3D]

PBMC_human_6x:
  inherits: immune_human_2x
  aggregation_reaction_IDs:
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

The website renders parameter tables from the same manifest rather than maintaining a second documentation schema. This keeps the user-facing configuration reference tied to the runtime validation contract.

Module-specific YAML uses the same mechanism. Aggregations opt into modules through the aggregation config, and each enabled aggregation must have a matching module config row. Some cross-scope fallbacks are still implemented by target code rather than by manifest inheritance; for example, a module parameter may intentionally allow `NULL` and then fall back to an aggregation-level path during module setup.

## Mapping tibbles

The root `_targets.R` builds the target graph from mapping tibbles. Each mapping tibble is a row-wise contract: one row becomes one set of mapped target instances, and columns in that row become local symbols inside the corresponding `tarchetypes::tar_map()` block.

The core mapping flow is:

1. `dataset_tibble_from_yaml` is read from `cfg_datasets.yaml`.
2. `reaction_tibble` is read from `cfg_reactions.tsv` and joined to dataset-level config.
3. `aggregation_tibble_all_from_yaml` is read from `cfg_aggregations.yaml`.
4. `aggregation_tibble` keeps active aggregations and adds upstream target-symbol columns.
5. `dataset_tibble` collapses reaction rows back to one row per dataset.
6. `_targets.R` expands reactions, datasets, and aggregations with `tar_map()`, then appends module target files.

```r
tarchetypes::tar_map(
  values = reaction_tibble,
  names = reaction_ID,
  delimiter = ".",
  source("extra_targets/per_reaction_targets.R")$value
)
```

With `reaction_ID = "healthy_PBMC_human"`, a target named `cellranger_summary_file` becomes `cellranger_summary_file.healthy_PBMC_human`. The same dot-delimited suffix convention is used for datasets, aggregations, module targets, and nested module maps.

```text
cfg_reactions.tsv row       -> reaction_tibble row    -> per-reaction targets
cfg_datasets.yaml key       -> dataset_tibble row     -> per-dataset targets
cfg_aggregations.yaml key   -> aggregation_tibble row -> per-aggregation targets
```

Aggregation rows may opt into optional modules through `modules`. `_targets.R` validates module names against the known module list, and module target files then filter `aggregation_tibble` to the active aggregations that requested that module. Each opted-in aggregation must have a matching module config row.

The naming convention is therefore compositional:

```text
<target>.<reaction_ID>
<target>.<dataset_name>
<target>.<aggregation_name>
<module_target>.<module_name>.<aggregation_name>
<nested_module_target>.<nested_suffix>.<module_name>.<aggregation_name>
```

Because these suffixes become target names and cache identity, config keys should be stable, human-readable, and free of unnecessary punctuation. In particular, avoid dots in reaction, dataset, aggregation, and module IDs unless there is a compelling reason.

## Target-symbol columns

Mapped target tables sometimes need to carry references to other mapped targets. multiomeR represents those references as columns of `rlang` symbols. Each row stores the upstream target symbols that should be spliced into downstream target commands generated for that row.

The compact constructor is `target_sym_col()`. It records a base target name, the source column containing suffixes, the separator, and an optional transform. `add_target_sym_cols()` then turns those specifications into list-columns of `rlang::syms()`.

```r
aggregation_tibble |>
  add_target_sym_cols(
    aggregation_GEX_counts_BPCells_matrix_syms =
      target_sym_col("GEX_counts_BPCells_matrix", "aggregation_reaction_IDs")
  )
```

For an aggregation whose `aggregation_reaction_IDs` are `c("rx1", "rx2")`, this creates a row value equivalent to:

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

Column names should describe the downstream scope, the upstream target, and the fact that the value is a symbol list. Existing symbol-list columns use the `*_syms` suffix, such as `aggregation_GEX_counts_BPCells_matrix_syms`, `dataset_unfiltered_cells_n_vecs_syms`, and `per_dataset_QC_violins_syms`.

The shared pattern appears at three boundaries:

- `aggregation_*_syms` columns splice per-reaction targets into aggregation-level targets.
- `dataset_*_syms` columns splice per-reaction targets into dataset-level targets.
- `per_dataset_*_syms` columns splice per-dataset targets into aggregation-level summary targets.

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
2. sourcing project helpers with `targets::tar_source("R")`,
3. applying global plotting and `{targets}` options,
4. installing the current `{targets}` patches,
5. sourcing `crew_controllers.R` and installing controller resources.

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

If `get_tar_resources()` is called before controller resources are loaded, it fails deliberately with an instruction to call `load_project_runtime()` first. Scheduler-specific examples belong in the main manual's [Distributed computing](../performance_distributed_computing.html) page; the implementation contract is that target code can request resources declaratively once the runtime has been loaded.

## How to read the rest of the implementation book

These conventions are the connective tissue behind the graph chapters. The parameter manifest explains why config rows can be compact. Mapping tibbles explain why target names have stable suffixes. Target-symbol columns explain how mapped targets pass sets of upstream targets across graph levels. Target metadata tags explain why some nodes remain visible in curated graph views. The bootstrap contract explains why helper functions, controller resources, and target options are available before `_targets.R` is evaluated.

When modifying the implementation, preserve these contracts unless the change is explicitly meant to replace one of them.


## Part: Target graph methodology


<!-- source: website/implementation/graph_methodology.qmd -->

# Reading the graph views



The graph chapters collect simplified views of the real `{targets}` dependency graph. They are meant to make the workflow easier to reason about before reading the target code directly.

The diagrams are generated from tagged target metadata and the real dependency graph, then simplified by pruning or bypassing lower-level nodes that would make each view harder to read. They keep real target names and preserve the dependency structure where practical, while staying compact enough to build intuition about the main control points.

The following chapters cover the main pipeline, the differential analyses module, and the genetic enrichment module.

[Mermaid graph omitted; source: `website/figures/standard_node_color_legend.mmd`]


<!-- source: website/implementation/implementation_main.qmd -->

<!-- begin include: website/implementation_main.qmd -->

# Main pipeline



## Parallel pre-processing

This view covers per-reaction processing and QC, including optional ambient RNA correction, donor demultiplexing, doublet detection, barcode filtering, and handoffs into aggregation-level GEX and ATAC objects.

[Mermaid graph omitted; source: `website/figures/human_curated/parallel_v2.mmd`]

## GEX processing

This view covers merged RNA processing, clustering, marker detection, cell type annotation, and GEX review outputs.

[Mermaid graph omitted; source: `website/figures/human_curated/GEX_v2.mmd`]

## ATAC processing

This view covers ATAC QC, peak calling, consensus peak construction, chromatin accessibility processing, chromVAR scoring, coverage tracks, and peak-gene links.

[Mermaid graph omitted; source: `website/figures/human_curated/ATAC_v2.mmd`]

## WNN integration

This view covers GEX and ATAC embedding handoffs, WNN integration, modality weights, cluster comparison, and integrated metadata outputs.

[Mermaid graph omitted; source: `website/figures/human_curated/WNN_v2.mmd`]

## Subgroup reprocessing

This view covers optional subgroup-native GEX, ATAC, and WNN reprocessing for sufficiently large parent groups.

[Mermaid graph omitted; source: `website/figures/human_curated/full_subgroups_v2.mmd`]

<!-- end include: website/implementation_main.qmd -->


<!-- source: website/implementation/implementation_differential_analyses.qmd -->

<!-- begin include: website/implementation_differential_analyses.qmd -->

# Differential analyses



[Mermaid graph omitted; source: `website/figures/human_curated/differential_analyses_v2.mmd`]

<!-- end include: website/implementation_differential_analyses.qmd -->


<!-- source: website/implementation/implementation_genetic_enrichment.qmd -->

<!-- begin include: website/implementation_genetic_enrichment.qmd -->

# Genetic enrichment



## Single-nucleus and graph-based enrichment

This view covers across-gene enrichment outputs and downstream trait summaries.

[Mermaid graph omitted; source: `website/figures/human_curated/genetic_enrichment_single_nucleus_v2.mmd`]

<!-- end include: website/implementation_genetic_enrichment.qmd -->


# Orphaned QMD Pages

Tracked QMD files not reached from the Quarto book graph or include graph.

- `website/helpers/_targets_graph_snippet.qmd`
- `website/implementation_overview.qmd`
