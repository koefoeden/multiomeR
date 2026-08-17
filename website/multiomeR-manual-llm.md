# multiomeR Website LLM Export

This file is generated from the Quarto book outlines and resolves Quarto include shortcodes.
Hidden setup chunks, generated helper chunks, Mermaid graph bodies, and verbose image metadata are omitted by default.


# Book: multiomeR Manual


## Part: Overview


<!-- source: website/index.qmd -->

[Image omitted; source: `figures/multiomeR-logo.svg`; alt: Image]

# Introduction

## What is multiomeR?
multiomeR is a [targets](https://books.ropensci.org/targets/)-based workflow for processing and analyzing single-nucleus 10x Genomics Multiome ATAC + Gene Expression datasets. It is meant to be adapted to your own data, compute setup, and biological questions.

The workflow starts from `cellranger-arc count` outputs, processes gene-expression (GEX) and ATAC data, builds multimodal aggregations, and supports optional downstream modules for differential analyses and genetic enrichment for human datasets.

## Development status
multiomeR is currently in beta. Expect breaking changes as the workflow is maturing. Report issues, questions, and feature requests on GitHub.

## Why use it?
- **Adaptable framework**: The repository is meant to be modified through config files, target selections, helper functions, and optional modules.
- **Transparent execution**: Each processing and analysis step is represented as an explicit `targets` target.
- **Efficient reruns**: `targets` reruns only the parts of the workflow affected by changed inputs, code, or config.
- **Multimodal analysis**: The main pipeline covers RNA/ATAC processing, clustering, cell typing, WNN integration, and optional downstream modules.

## Choose your route

| If you want to... | Start here |
|---|---|
| See what the workflow produces | Browse the [main pipeline output gallery](gallery_main.qmd). |
| Try multiomeR on public data | Follow the three-part quickstart: [install](demo_installation.qmd), [run](demo_running.qmd), then [inspect the outputs](demo_outputs.qmd). |
| Configure your own data | Read the [main-pipeline overview](main_overview.qmd), prepare the [configuration and inputs](main_inputs.qmd), then [run the checkpoints](main_running.qmd). |
| Add a downstream analysis | Check the prerequisites for [differential analyses](downstream_differential_analyses.qmd) or [genetic enrichment](downstream_genetic_enrichment.qmd). |
| Understand or modify the internals | Use the separate [implementation book](implementation/). |

Before starting, you need Linux and basic familiarity with R, tabular metadata, and `targets`. The demo supplies a working configuration; your own analysis also needs complete `cellranger-arc count` outputs, the matching reference metadata JSON, and keyed GEM well and donor metadata. Tutorial pages provide a safe route through the workflow; searchable parameter tables and the implementation book are references to consult as needed.

In this manual, a **GEM well** is one configured 10x library and output directory, a **dataset** groups GEM wells that share reference and QC settings, and an **aggregation** is a joint analysis of one or more GEM wells. A **donor** is the individual identified by `donor_id`; one GEM well may contain multiple donors.

## Workflow at a glance

The **main pipeline** processes each GEM well, aggregates selected GEM wells, and builds multimodal RNA/ATAC outputs for clustering, cell typing, and WNN integration. Two optional modules extend completed aggregations with differential analyses or genetic enrichment.

[Image omitted; source: `figures/multiomeR_overview_simplified.drawio.svg`; alt: multiomeR workflow from Cell Ranger ARC GEM well outputs through per GEM well processing, aggregation-level GEX and A...]


## Part: Quickstart: Public demo with two GEM wells


<!-- source: website/demo_installation.qmd -->

# Install and prepare the demo



## System requirements

- Linux with `git` and `curl`, plus HTTPS access to GitHub, Pixi, and 10x Genomics downloads.
- At least 60 GB of RAM. This is enough for one heavy target at a time; machines near the minimum should reduce concurrent workers in `crew_controllers.R`.
- At least 30 GB of free disk space for the public inputs, pixi environment, temporary files, and approximately 6 GB of demo outputs.
- Multiple CPU cores are strongly recommended. The timing quoted in the next chapter was measured with 16 logical threads.

The committed `crew_controllers.R` symlink points to the local controller template. That template assumes a 16-CPU, 256-GB workstation and can run several workers concurrently. Review [Distributed computing](performance_distributed_computing.qmd) before running on a smaller machine or a scheduler.


## Clone and enter the repository

```{.bash filename="Bash"}
git clone https://github.com/koefoeden/multiomeR.git
cd multiomeR
```

## Download public inputs

Download the subset of two public human `cellranger-arc count` outputs required by the `immune_human_2x` example aggregation from the repository root (3.9 GB). For your own data, retain the complete `outs/` directories.

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

The repository already includes the small `reference.json` from the exact
`refdata-cellranger-arc-GRCh38-2020-A-2.0.0` reference used to create both
public outputs. The demo configuration points to that file, so the full Cell
Ranger ARC reference is not required.

## Create and launch the pixi environment

If `pixi` is not installed yet, install it first:

```{.bash filename="Bash"}
curl -fsSL https://pixi.sh/install.sh | sh
export PATH="$HOME/.pixi/bin:$PATH"
```

Install the locked pixi environment from the repository root. The duration depends on the package cache, network connection, and the source builds in the following step.

```{.bash filename="Bash"}
pixi install --locked --run-post-link-scripts
pixi run install-r-github-packages
```

```{.bash filename="Bash"}
pixi shell
R
```

Use this repository-root R session as the evaluation environment for subsequent R commands in this manual.

## Validate the setup

Before starting a compute-intensive run, build the target manifest for the exact demo endpoint:

```{.r filename="R"}
targets::tar_manifest(
  names = tidyselect::matches(
    "^multimodal_Seurat_object[.]immune_human_2x$"
  ),
  callr_function = NULL
)[, c("name", "description")]
```

The command should return one row named `multimodal_Seurat_object.immune_human_2x`. Constructing the manifest also evaluates the configuration and target graph, so configuration, controller, or startup errors appear before the data run begins. Continue to [Run the demo](demo_running.qmd) when this check succeeds.

::: {.callout-note collapse="true"}
## Why are both installation commands needed?

Pixi resolves R, Bioconductor, and command-line tools in one locked environment. `--run-post-link-scripts` materializes Bioconductor data packages such as BSgenome, while `pixi run install-r-github-packages` installs the pinned GitHub versions of BPCells, Signac, and betterChromVAR that are not supplied by the locked conda environment.
:::


<!-- source: website/demo_running.qmd -->

# Run the demo



Run the command below to process the example `immune_human_2x` aggregation and produce its final multimodal Seurat/Signac compatibility object. The pipeline keeps BPCells-backed matrices as its primary state; the Seurat object is a convenience export and a clear end point for the quickstart.

This selection builds the upstream objects and file targets needed by the final export. It does not request every sibling review plot in the GEX, ATAC, and WNN checkpoints. On one AMD EPYC 7543 system using 16 logical threads, a previous run of this example completed in under 25 minutes and wrote about 6 GB. Treat that as orientation, not a runtime guarantee.

```{.r filename="R"}
targets::tar_make(
  names = tidyselect::matches(
    "^multimodal_Seurat_object[.]immune_human_2x$"
  )
)
```

## Confirm success

After the run, the following command should return `character(0)`, meaning the requested endpoint and its dependencies are up to date:

```{.r filename="R"}
targets::tar_outdated(
  names = tidyselect::matches(
    "^multimodal_Seurat_object[.]immune_human_2x$"
  ),
  callr_function = NULL
)
```

Confirm that the final object can be read:

```{.r filename="R"}
demo_object <- targets::tar_read(
  multimodal_Seurat_object.immune_human_2x
)
demo_object
```

Continue to [Verify and inspect the outputs](demo_outputs.qmd). If the run stops with an error, use [Troubleshooting](troubleshooting.qmd) before broadening the target selection.


<!-- source: website/demo_outputs.qmd -->

# Verify and inspect the outputs



The public `immune_human_2x` demo writes its `targets` store to `outputs`, as configured in `_targets.yaml`. The quickstart selection creates the objects and file targets required by the final compatibility export. Plot targets are separate review outputs and are not all dependencies of that endpoint.

The store contains serialized R objects in `objects/`, file artifacts in `files/`, and requested review figures in `plots/`.

## Objects

Most intermediate and final result objects are saved automatically by `targets` under `outputs/objects` during a pipeline run. Read them with `targets::tar_read()` from a repository-root R session after the demo has completed.

```{.r filename="R"}
targets::tar_read(UMAP_embeddings_tibble.WNN.immune_human_2x)
targets::tar_read(metadata_w_cell_types_tibble.WNN.immune_human_2x)
targets::tar_read(multimodal_Seurat_object.immune_human_2x)
```

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

Other files are grouped under `outputs/files/<scope>/`, where the scope is a GEM well, dataset, or aggregation.

## Plots

Plots are saved in a similar structure under `outputs/plots`, but the final-object quickstart selection does not request every checkpoint plot. To build and review those outputs, follow the GEX, ATAC, and multimodal selections in [Running the workflow](main_running.qmd).

Requested plots use the same scope-based layout under `outputs/plots/`.

The [main pipeline output gallery](gallery_main.qmd) contains curated snapshots showing what representative checkpoint outputs look like. Gallery assets are documentation snapshots; they are not evidence that the corresponding target was built in your local store.

## Next steps

- To adopt the workflow, continue with the [main-pipeline overview](main_overview.qmd) and [configuration walkthrough](main_inputs.qmd).
- To generate all checkpoint review outputs for an aggregation, use [Running the workflow](main_running.qmd).
- To diagnose a failed or unexpectedly stale target, use [Troubleshooting](troubleshooting.qmd).


## Part: Configure and run your own data


<!-- source: website/main_overview.qmd -->

# Main-pipeline overview

## Purpose

This page maps the stable stages of the main processing branch after `cellranger-arc count`. For required files and metadata, continue to [Configuration and inputs](main_inputs.qmd); for representative results, see the [main pipeline gallery](gallery_main.qmd).

## Processing stages

1. **GEM well preprocessing** reads Cell Ranger matrices and fragments, calculates RNA and ATAC QC, and assigns stable GEM well-prefixed barcodes.
2. **Dataset-level review** compares GEM wells that share a reference genome and QC policy before they enter an aggregation.
3. **GEX processing** merges selected GEM wells into BPCells-backed matrices, performs dimension reduction and clustering, assigns cell-type labels, and creates the first review checkpoint.
4. **ATAC processing** merges fragments, calls or tiles peaks, creates the consensus peak matrix, performs LSI and clustering, and summarizes motif and regulatory activity.
5. **Multimodal integration** combines aligned GEX and ATAC embeddings with WNN, produces integrated metadata and review plots, and optionally exports a [`Seurat`](https://satijalab.org/seurat/)/[`Signac`](https://stuartlab.org/signac/) compatibility object.

Configuration can add genotype-aware donor assignment, CellBender inputs, Harmony correction, chromHMM annotation, peak-gene correlation, coverage tracks, and subgroup reprocessing.

The [implementation graph](implementation/implementation_main.html) expands these stages to target level. After configuring the inputs, [run and review the GEX, ATAC, and multimodal checkpoints](main_running.qmd) in order.


<!-- source: website/main_inputs.qmd -->

# Configuration and inputs



multiomeR uses three linked configuration layers. A **dataset** defines a shared reference and GEM well-level QC policy, a **GEM well** points to one `cellranger-arc count` output, and an **aggregation** selects GEM wells for joint GEX, ATAC, and WNN analysis.

| Layer | Configuration | Key relationship |
|---|---|---|
| Dataset | `cfg_datasets.yaml` | GEM wells refer to the dataset key. |
| GEM well | `cfg_GEM_wells.tsv` | Aggregations refer to one or more `GEM_well_ID` values. |
| Aggregation | `cfg_aggregations.yaml` | Metadata tables must cover the configured donor and GEM well IDs. |

The committed `cfg_datasets.yaml`, `cfg_GEM_wells.tsv`, and `cfg_aggregations.yaml` files are symlinks to public templates. Keep those templates as examples. For a real project, replace the symlinks with project-specific files or provide them from a downstream private repository.

## Migration to GEM well terminology

This release makes a breaking terminology transition to the 10x Genomics term **GEM well**. There are no compatibility aliases: update configuration, metadata, and any code that calls the renamed helpers together. Existing identifier values do not need to change.

| Previous name | Current name |
|---|---|
| `cfg_reactions.tsv` / `cfg_reactions_template.tsv` | `cfg_GEM_wells.tsv` / `cfg_GEM_wells_template.tsv` |
| `reaction_ID` or `TENX_reaction_ID` | `GEM_well_ID` |
| `aggregation_reaction_IDs` | `aggregation_GEM_well_IDs` |
| `aggregation_reaction_ID_metadata_tsv` | `aggregation_GEM_well_metadata_tsv` |
| `reaction_donor_id`, `reaction_n_donors`, `reaction_cellbender_h5_file`, `reaction_donors_VCF_file` | The corresponding `GEM_well_*` columns |
| `reaction_cellranger_count_dir` | `GEM_well_cellranger_arc_count_dir` |
| `per_reaction_targets.R` | `per_GEM_well_targets.R` |

The reference input is also narrowed from a directory to the metadata file the
pipeline actually reads. Replace `dataset_cellranger_arc_refdata_dir` with
`dataset_cellranger_arc_reference_json` and append `/reference.json` to existing
reference-directory values.

## Start with one explicitly scoped analysis

`is_active` defaults to `true` for every aggregation. Before an unqualified `targets::tar_make()`, either remove unavailable template rows or add `is_active: false` to every aggregation you are not ready to run. During setup, always select one aggregation explicitly.

The following minimum example shows the relationship between the three layers. Replace the paths, identifiers, and marker genes with values appropriate for your study.

### 1. Define one dataset

```{.yaml filename="cfg_datasets.yaml"}
your_dataset:
  dataset_cellranger_arc_reference_json: /path/to/reference.json
  dataset_QC_exclude_list_per_GEM_well:
    - TSS.enrichment < 4
    - nucleosome_signal > 4
    - nCount_RNA < 250
  dataset_run_amulet: false
```

`dataset_cellranger_arc_reference_json` is required and must point to the
`reference.json` from the exact Cell Ranger ARC reference used to create the
dataset outputs. The repository includes metadata for the public human and
mouse 2020-A references and the human 2024-A reference under
`reference_metadata/`; custom references can supply their own JSON without
retaining the full reference directory. multiomeR checks that the JSON genome
matches the feature HDF5 and rejects aggregations whose GEM wells use different
references.

The QC expressions are evaluated against per-barcode metadata. Starting with `dataset_run_amulet: false` avoids adding AMULET to the first data run; enable it deliberately after the basic input contract works. When enabled, multiomeR runs its [validated BPCells-native AMULET implementation](implementation/algorithm_validation.html#bpcells-native-amulet) on the existing per-GEM-well fragment directory. The implementation preserves scDblFinder's tested unique-fragment metrics while avoiding full fragment materialization; PCR-duplicate expansion (`uniqueFrags = FALSE`) is not supported.

### 2. Define one GEM well

```{.text filename="cfg_GEM_wells.tsv"}
GEM_well_ID	dataset	GEM_well_donor_id	GEM_well_n_donors	GEM_well_cellranger_arc_count_dir	GEM_well_cellbender_h5_file	GEM_well_donors_VCF_file
your_GEM_well	your_dataset	donor_1	1	/path/to/your_GEM_well	NA	NA
```

`GEM_well_cellranger_arc_count_dir` points to the directory containing `outs/`, not to `outs/` itself. The baseline pipeline requires:

```{.text filename="Text"}
<GEM_well_cellranger_arc_count_dir>/outs/
|-- summary.csv
|-- filtered_feature_bc_matrix.h5
|-- atac_fragments.tsv.gz
|-- atac_fragments.tsv.gz.tbi
`-- per_barcode_metrics.csv
```

If `GEM_well_donors_VCF_file` is configured, `atac_possorted_bam.bam` is also required for `cellsnp-lite`. Without a VCF, the pipeline skips genotype demultiplexing and assigns `GEM_well_donor_id` to every called nucleus. That donor ID must match the donor metadata table.

### 3. Create keyed metadata

The donor metadata table must contain one unique row per `donor_id`:

```{.text filename="donor_metadata.tsv"}
donor_id	condition
donor_1	control
```

The GEM well metadata table must contain one unique row per `GEM_well_ID`:

```{.text filename="GEM_well_metadata.tsv"}
GEM_well_ID	batch
your_GEM_well	batch_1
```

Apart from their key columns, the donor and GEM well tables must not reuse column names. Put donor-specific phenotypes and covariates in the donor table; put library-, run-, or batch-specific variables in the GEM well table.

The GEM well table may contain rows for multiple aggregations. Each aggregation
subsets it to `aggregation_GEM_well_IDs` before expensive targets consume it,
so metadata changes for unrelated GEM wells do not propagate downstream.

### 4. Define one aggregation

```{.yaml filename="cfg_aggregations.yaml"}
your_aggregation:
  aggregation_GEM_well_IDs: [your_GEM_well]
  aggregation_donor_id_metadata_tsv: /path/to/donor_metadata.tsv
  aggregation_GEM_well_metadata_tsv: /path/to/GEM_well_metadata.tsv
  aggregation_GEX_marker_genes:
    Cell_type_A: [GENE1, GENE2]
    Cell_type_B: [GENE3, GENE4]
  is_active: true
```

Omit `modules` for the first main-pipeline run. Add optional module names and matching module-config rows only after the main aggregation has passed its review checkpoints.

### 5. Validate before running

From the repository-root R session, construct the graph and inspect the targets created for the aggregation:

```{.r filename="R"}
manifest <- targets::tar_manifest(callr_function = NULL)

manifest |>
  dplyr::filter(stringr::str_ends(name, ".your_aggregation")) |>
  dplyr::select(name, description)
```

Manifest construction validates the YAML parameter schema, aggregation references to GEM wells, module names and rows, and controller setup. Metadata file contents are validated when their targets run. Fix manifest-time errors before calling `tar_make()`, then continue to [Running the workflow](main_running.qmd).

## Configuration reference

The searchable overviews below are generated from `cfg_pipeline_parameters.tsv`, the same manifest used for runtime defaults and validation. Use them to change a default after the minimum configuration works.

### Dataset parameters

Dataset-level settings live in `cfg_datasets.yaml`.

[Generated Quarto chunk omitted: `emit_parameter_overview("dataset")`]

<details>
<summary>Show the public <code>immune_human_dataset</code> example</summary>

[Generated Quarto chunk omitted: `emit_yaml_template_entry(datasets_config_file, "immune_human_dataset")`]

</details>

### GEM well columns

GEM well-level configuration lives in `cfg_GEM_wells.tsv`.

| Column | Meaning |
|---|---|
| `GEM_well_ID` | Stable GEM well identifier used in target names and aggregation config. |
| `dataset` | Dataset key from `cfg_datasets.yaml`. |
| `GEM_well_donor_id` | Donor assigned to a non-multiplexed GEM well; use `NA` when donor identities are resolved from a VCF. |
| `GEM_well_n_donors` | Expected donor count. |
| `GEM_well_cellranger_arc_count_dir` | Path to the directory containing `outs/`. |
| `GEM_well_cellbender_h5_file` | Optional CellBender-corrected GEX H5 file; otherwise `NA`. |
| `GEM_well_donors_VCF_file` | Optional donor-genotype VCF for `cellsnp-lite` and `vireo`; otherwise `NA`. |

<details>
<summary>Show the public GEM well table</summary>

[Generated Quarto chunk omitted: `GEM_well_template <- readr::read_tsv("cfg_GEM_wells.tsv", show_col_types = FALSE) cat('<div class="scrollable-table">...`]

</details>

### Aggregation parameters

Aggregation-level settings live in `cfg_aggregations.yaml`. Required parameters must resolve to a non-missing value after defaults and optional `inherits:` parents are applied. Optional parameters can remain `NULL` to disable the corresponding behavior.

[Generated Quarto chunk omitted: `emit_parameter_overview("aggregation")`]

<details>
<summary>Show the public <code>immune_human_2x</code> example</summary>

[Generated Quarto chunk omitted: `emit_yaml_template_entry(aggregations_config_file, "immune_human_2x")`]

</details>


<!-- source: website/main_running.qmd -->

# Run and review the workflow

From the repository-root R session described in [Install and prepare the demo](demo_installation.qmd), run one aggregation in checkpoint-sized pieces. Replace `your_aggregation` with a configured aggregation name such as `immune_human_2x`.

Before running, confirm that `crew_controllers.R` matches the machine or scheduler you intend to use. See [Distributed computing](performance_distributed_computing.qmd) for the controller contract.

## Preview a checkpoint selection

Use `tar_manifest()` with the same selector you plan to pass to `tar_make()`. This shows the requested checkpoint targets without executing them:

```{.r filename="R"}
targets::tar_manifest(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:GEX")
  ) & tidyselect::matches("your_aggregation"),
  callr_function = NULL
)[, c("name", "description")]
```

If the result is empty, check the aggregation spelling. If manifest construction fails, resolve the configuration or controller error before starting compute work.

## GEX checkpoint

```{.r filename="R"}
targets::tar_make(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:GEX")
  ) & tidyselect::matches("your_aggregation")
)
```

Review before continuing:

- GEM well- and dataset-level barcode exclusion summaries for unexpected sample loss;
- PCA, Harmony, UMAP, and metadata-association diagnostics for technical structure;
- cluster markers and marker-module scores for coherent cell-type labels; and
- `GEX_Seurat_object.your_aggregation` plus `metadata_w_cell_types_tibble.GEX.your_aggregation` for the cells entering ATAC processing.

The [main gallery](gallery_main.qmd) shows representative GEX review plots. Change the dataset QC, aggregation dimensions, clustering, or marker configuration and rerun this checkpoint if the result is not biologically and technically credible.

## ATAC checkpoint

Run this only after accepting the GEX checkpoint:

```{.r filename="R"}
targets::tar_make(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:ATAC")
  ) & tidyselect::matches("your_aggregation")
)
```

Review before continuing:

- peak-based QC exclusions and retained cell counts;
- LSI dimensions, metadata associations, UMAPs, and cluster stability;
- consensus peaks, motif matches, TF activity, and marker-gene activity; and
- coverage tracks and differential-accessibility summaries where configured.

Key objects include `consensus_peak_tibble.ATAC.your_aggregation`, `metadata_w_cell_types_tibble.ATAC.your_aggregation`, and the BPCells peak and TF-activity matrices. Revisit the peak-calling, ATAC dimensions, QC, or marker-TF settings if the review fails.

## Multimodal checkpoint

```{.r filename="R"}
targets::tar_make(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:multimodal")
  ) & tidyselect::matches("your_aggregation")
)
```

Review before continuing:

- WNN UMAPs and cluster identities against the accepted GEX and ATAC results;
- RNA and ATAC modality weights for unexpected modality dominance;
- cross-modality cluster comparisons and final integrated metadata; and
- `multimodal_Seurat_object.your_aggregation` if a compatibility export is needed.

Use [Verify and inspect the outputs](demo_outputs.qmd) for `tar_read()` and output-path examples.

## Finish the intended scope

To build the final main-pipeline endpoint for one aggregation after checkpoint review:

```{.r filename="R"}
targets::tar_make(
  names = tidyselect::matches(
    "^multimodal_Seurat_object[.]your_aggregation$"
  )
)
```

An unqualified `targets::tar_make()` runs every active aggregation and every enabled optional module in the root workflow. Use it only when that is the intended scope. Because `is_active` defaults to `true`, disable unavailable template aggregations before broad execution.

Continue to the [differential analyses](downstream_differential_analyses.qmd) or [genetic enrichment](downstream_genetic_enrichment.qmd) module only after the main aggregation is accepted. If a target fails, follow [Troubleshooting](troubleshooting.qmd) and rerun the narrowest affected selection.


## Part: Optional modules


<!-- source: website/downstream_differential_analyses.qmd -->

# Differential analyses



## When to use this module

Use this optional module when a completed aggregation contains replicated donor- or sample-level experimental, clinical, phenotypic, or technical variables that should be modelled across cell-type composition or pseudobulk molecular features.

The module does not create biological replication. The donor structure, covariates, design formula, and contrasts must be defensible for the intended analysis before the workflow is run.

See the [output gallery](gallery_differential_analyses.qmd) for representative diagnostics and the [implementation graph](implementation/implementation_differential_analyses.html) for target structure.

## Prerequisites

Before enabling the module, confirm that:

- the aggregation has passed the GEX, ATAC, and multimodal checkpoints;
- WNN cell-type metadata and GEX, ATAC, and TF-activity pseudobulk matrices are available;
- the donor metadata contains one unique row per `donor_id` and every variable used in a model;
- model variables are donor- or pseudobulk-sample-level variables, not duplicated cell-level measurements; and
- the number and distribution of donors support the specified design and contrasts.

Use `differential_analyses_extended_donor_id_metadata_tsv` when the modelling table needs variables beyond the aggregation's normal donor metadata. It must retain the same unique `donor_id` key.

## Outputs

It produces cell-type-composition models; pseudobulk differential gene expression, accessibility, and TF activity; diagnostics and cross-modality summaries; and competitive gene-set tests against the MSigDB Hallmark and Reactome collections. Each collection is tested independently with `cameraPR`, `inter.gene.cor = 0.01`, and a minimum of 10 genes represented in the contrast-specific universe. A significant set is more strongly associated with the contrast than the remaining tested genes, rather than merely showing any collective change. Open Targets evidence annotation is optional.

## Configure the module

First opt the aggregation into the module:

```{.yaml filename="cfg_aggregations.yaml"}
your_aggregation:
  modules: [differential_analyses]
```

Then create a matching row in `module_differential_analyses/cfg.yaml`. The committed file is a symlink to the public template; preserve the template as a reference and use a project-specific file for a real analysis.

```{.yaml filename="module_differential_analyses/cfg.yaml"}
your_aggregation:
  differential_analyses_DCTC_plot_phenotype_vars: condition
  differential_analyses_DCTC_formula_chr: >-
    cbind(n_nuclei, n_other_nuclei) ~ 0 + condition
  differential_analyses_psbulk_DX_models:
    condition_model:
      cell_type_subset: NULL
      design_matrix_func_name: NULL
      formula: ~ 0 + cluster + condition
      random_effect: NULL
      contrast_specs_vec:
        treated_vs_control: conditiontreated
```

The full module checkpoint requests both composition and pseudobulk review outputs. Configure the DCTC phenotype/formula and at least one pseudobulk model before using that broad selector. Formula terms and contrast coefficients must match columns produced by the model matrix.

### Parameter reference

The OLINK and bulk-RNA path fields are reserved optional integration inputs and are not consumed by the current public differential-analysis checkpoint. Leave them `NULL` unless the corresponding integration is implemented in your downstream workflow.

[Generated Quarto chunk omitted: `emit_parameter_overview("differential_analyses")`]

<details>
<summary>Show the public <code>immune_human_2x</code> example</summary>

[Generated Quarto chunk omitted: `emit_yaml_template_entry(module_cfg_template_file, "immune_human_2x")`]

</details>

## Run and review

Preview the selected targets first:

```{.r filename="R"}
targets::tar_manifest(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:differential_analyses")
  ) & tidyselect::matches("your_aggregation"),
  callr_function = NULL
)[, c("name", "description")]
```

Then run the checkpoint:

```{.r filename="R"}
targets::tar_make(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:differential_analyses")
  ) & tidyselect::matches("your_aggregation")
)
```

Review pseudobulk depths and retained donor counts before interpreting coefficients. Check model-matrix terms, P-value distributions, effect directions, and agreement or disagreement across DGE, DCA, and DTFA. Treat the [gallery](gallery_differential_analyses.qmd) as a visual reference, not as a statistical acceptance threshold.

Runtime depends on donors, cell types, models, contrasts, and GSEA branches. Use [Troubleshooting](troubleshooting.qmd) if a formula, contrast, or metadata join fails.


<!-- source: website/downstream_genetic_enrichment.qmd -->

# Genetic enrichment



## When to use this module

Use this optional module for a completed human aggregation when fine-mapped GWAS evidence should be integrated with chromatin accessibility. It maps Open Targets credible-set variants to consensus peaks, calculates GWAS chromVAR-style deviations with the pipeline's chromVAR/betterChromVAR machinery, propagates trait relevance with [`SCAVENGE`](https://github.com/sankaranlab/SCAVENGE), and attributes cell-type deviations back to loci and variants.

See the [output gallery](gallery_genetic_enrichment.qmd) for representative results and the [implementation graph](implementation/implementation_genetic_enrichment.html) for upstream ATAC and WNN dependencies.

## Prerequisites

Before enabling the module, confirm that:

- the aggregation is human and has passed the GEX, ATAC, and multimodal checkpoints;
- WNN metadata and graph results, consensus peaks, ATAC counts, chromVAR objects, and GEX/ATAC embeddings are available;
- each configured Open Targets `studyId` represents the intended trait and population;
- the machine can download and retain the Open Targets study and credible-set Parquet datasets; and
- the biological question justifies cell-level, graph-propagated, or cell-type-pseudobulk enrichment rather than treating those summaries as interchangeable.

## Outputs

It produces Open Targets and fine-mapping metadata, posterior-weighted variant-to-peak mappings, single-nucleus GWAS deviations, SCAVENGE trait-relevance summaries, cell-type heatmaps, and peak-, locus-, and variant-level attribution outputs.

## Configure the module

First opt a human aggregation into the module:

```{.yaml filename="cfg_aggregations.yaml"}
your_aggregation:
  modules: [genetic_enrichment]
```

Then create a matching row in `module_genetic_enrichment/cfg.yaml`. The committed file is a symlink to the public template; keep that template as a reference and use a project-specific file for real studies.

```{.yaml filename="module_genetic_enrichment/cfg.yaml"}
your_aggregation:
  genetic_enrichment_open_targets_studies:
    lymphocyte_count:
      Category: positive_control
      studyId: GCST90002388
      finemappingMethod: auto
```

The root workflow currently pins Open Targets release `26.03`. That release identifier is recorded in downstream metadata and determines the available studies, credible sets, and fine-mapping methods.

For `finemappingMethod: auto`, multiomeR selects the first available supported method in this order: `SuSie`, `SuSiE-inf`, then `PICS`. Specify a method explicitly when the method itself is part of the analysis contract; the workflow fails if that method is unavailable for the study.

### Parameter reference

[Generated Quarto chunk omitted: `emit_parameter_overview("genetic_enrichment")`]

<details>
<summary>Show the public <code>immune_human_2x</code> example</summary>

[Generated Quarto chunk omitted: `emit_yaml_template_entry(module_cfg_template_file, "immune_human_2x")`]

</details>

## Run and review

Preview the selected checkpoint targets:

```{.r filename="R"}
targets::tar_manifest(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:genetic_enrichment")
  ) & tidyselect::matches("your_aggregation"),
  callr_function = NULL
)[, c("name", "description")]
```

Then run the checkpoint:

```{.r filename="R"}
targets::tar_make(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:genetic_enrichment")
  ) & tidyselect::matches("your_aggregation")
)
```

Before interpreting trait scores, verify the resolved Open Targets release and fine-mapping method, the number of credible-set loci and variants retained, and the overlap with consensus peaks. Then compare direct deviation summaries with SCAVENGE-propagated scores and inspect whether apparent cell-type enrichment is supported by enough nuclei and loci.

Runtime and disk use grow with studies, cells, graph representations, permutations, and attributed loci. The [gallery](gallery_genetic_enrichment.qmd) uses a larger aggregation with six GEM wells and is not produced by the minimal quickstart.


## Part: Output gallery


<!-- source: website/gallery_main.qmd -->

# Main pipeline gallery



These documentation snapshots show representative outputs from the public `immune_human_2x` configuration. The first cards summarize its shared dataset; GEX, ATAC, and WNN cards are aggregation-level. Each card names its target.

To reproduce these plot families, [install](demo_installation.qmd) and [run](demo_running.qmd) the demo, then request the broader [review checkpoints](main_running.qmd). The endpoint-only quickstart does not build every gallery plot; [inspect its local outputs](demo_outputs.qmd) to see the distinction.

[Generated Quarto chunk omitted: `render_gallery_section( gallery_items, "Main pipeline", subsection_descriptions = c( "Parallel GEM well pre-processin...`]


<!-- source: website/gallery_differential_analyses.qmd -->

# Differential analyses gallery



These cards are a curated subset from the public `immune_human_2x` configuration. They illustrate diagnostics, not acceptable effect sizes or significance patterns for another study. See [Differential analyses](downstream_differential_analyses.qmd) for prerequisites, models, and the checkpoint command.

[Generated Quarto chunk omitted: `render_gallery_section( gallery_items, "Differential analyses module", subsection_descriptions = c( "Gene expression"...`]


<!-- source: website/gallery_genetic_enrichment.qmd -->

# Genetic enrichment gallery



These curated outputs use the larger `PBMC_human_6x` aggregation, not the quickstart with two GEM wells, and show selected SCAVENGE/WNN results rather than every attribution output. See [Genetic enrichment](downstream_genetic_enrichment.qmd) for prerequisites, checkpoint targets, and interpretation guidance.

[Generated Quarto chunk omitted: `render_gallery_section( gallery_items, "Genetic enrichment module", subsection_descriptions = c( "Single-nucleus chro...`]


## Part: Operation and scaling


<!-- source: website/performance_overview.qmd -->

# Performance and scaling

multiomeR gains performance from two complementary mechanisms: `{targets}` exposes independent GEM wells, modalities, and analysis branches for concurrent execution, while BPCells keeps large matrices on disk and streams many operations instead of materializing dense objects.

Actual wall time depends on input nuclei, retained cells, peak counts, enabled plots and modules, hardware, scheduler latency, controller concurrency, and which targets are already up to date. Treat benchmark results as workload descriptions, not promises for another machine or configuration.

## Current example snapshot

The current stored benchmark contains two public aggregations and estimates the critical path to each final `multimodal_Seurat_object` from recorded `{targets}` runtime metadata.

| Aggregation | GEM wells | Cell Ranger input nuclei | Estimated critical path | Sum if ancestor targets ran serially |
|---|---:|---:|---:|---:|
| `immune_human_2x` | 2 | 17,277 | 17.3 minutes | 31.4 minutes |
| `PBMC_human_6x` | 6 | 51,291 | 23.8 minutes | 52.4 minutes |

The critical-path estimate treats dynamic branches as concurrently runnable. The difference between that estimate and the serial sum illustrates available parallelism for these two runs; it does not establish a general scaling law from two observations.

## How to use these numbers

- Compare changes only when the pipeline version, configuration, target endpoint, and estimation method are recorded together.
- Expect scheduler queue time and limited worker capacity to increase observed wall time beyond the graph-only critical-path estimate.
- Inspect target-level runtime and memory on your own representative aggregation before sizing a production run.
- Start with one narrow aggregation and checkpoint, then increase worker concurrency only when memory headroom is known.

Configure execution capacity in [Distributed computing](performance_distributed_computing.qmd). If a run is unexpectedly slow or repeatedly rebuilds targets, use [Troubleshooting](troubleshooting.qmd).


<!-- source: website/performance_distributed_computing.qmd -->

# Distributed computing



multiomeR runs `{targets}` through the [`crew` backend](https://books.ropensci.org/targets/crew.html). Check the active controller before the [first demo run](demo_running.qmd) and again before moving from a workstation to a scheduler.

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

## Local execution

A fresh clone points `crew_controllers.R` to `crew_controllers_template.R`. The template is sized for a 16-CPU, 256-GB workstation, with four light workers and two heavy workers. A machine near the 60-GB minimum should reduce concurrency to one heavy worker and should not run several memory-intensive targets simultaneously.

To create a machine-specific active file without editing the template in place:

```{.bash filename="Bash"}
cp --remove-destination crew_controllers_template.R crew_controllers.R
```

Then edit the worker counts and resource tiers in `crew_controllers.R`. Keep controller names identical between `controller_list` and `controller_resources_tibble`.

After changing the file, restart R or reload the project runtime explicitly:

```{.r filename="R"}
load_project_runtime(force = TRUE)
```

Rebuild a narrow manifest selection before starting the data run to validate the controller contract.

## Scheduler execution

For SLURM, PBS, SGE, or LSF, replace the local controllers with the corresponding `crew.cluster` controllers. The commented SLURM section in `crew_controllers_template.R` shows the expected shape.

For every scheduler tier:

1. Match the controller name in both the controller object and resource table.
2. Align scheduler CPU and memory requests with the capacity declared in the table.
3. Set queue, account, wall-time, module, and worker-startup options required by the cluster.
4. Keep GPU tiers separate; GPU controllers are considered only for targets requesting GPUs.
5. Test one small checkpoint before increasing worker counts.

Scheduler startup failures, resource-routing errors, and target failures are handled separately in [Troubleshooting](troubleshooting.qmd). Developer-facing details about runtime bootstrap and `get_tar_resources()` are in [Implementation conventions](implementation/implementation_conventions.html#runtime-bootstrap).


<!-- source: website/troubleshooting.qmd -->

# Troubleshooting

Start with the narrowest failing target or checkpoint. `{targets}` preserves successful work, so fixing the cause and rerunning the same selection is normally safer and faster than deleting `outputs/` or starting a broad workflow again.

## The manifest does not build

Run the project bootstrap and manifest in a repository-root R session:

```{.r filename="R"}
load_project_runtime(force = TRUE)
targets::tar_manifest(callr_function = NULL)
```

Manifest-time failures usually indicate one of these contracts:

- an unknown, misspelled, wrongly typed, or required YAML parameter;
- an aggregation referencing a GEM well absent from `cfg_GEM_wells.tsv`;
- an unknown module or a missing module-config row;
- a malformed `crew_controllers.R` return value; or
- an R package or startup problem.

Fix the named configuration or controller boundary first. Do not launch `tar_make()` until the manifest builds.

## A run reports errored targets

The project runtime provides helpers that collapse repeated branch failures into a readable first triage:

```{.r filename="R"}
list_distinct_errored_targets()
```

For commands and stored tracebacks matching a target, aggregation, or module:

```{.r filename="R"}
list_distinct_errored_targets_w_tracebacks(
  target_name_pattern = "your_target_or_aggregation"
)
```

Use the full target or dynamic-branch name from that output. Inspect the saved target workspace before loading large dependencies into the interactive session:

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

Check the contracts in [Configuration and inputs](main_inputs.qmd):

- `GEM_well_cellranger_arc_count_dir` contains the required `outs/` files;
- VCF-backed demultiplexing also has `atac_possorted_bam.bam`;
- donor and GEM well metadata keys are present and unique;
- non-key metadata columns belong to only one metadata table; and
- configured donor and GEM well IDs match the metadata values exactly.

## Controller and scheduler failures

If a worker does not start or no controller can satisfy a target request:

1. Validate the names, column order, numeric resource values, and controller membership described in [Distributed computing](performance_distributed_computing.qmd).
2. Reload with `load_project_runtime(force = TRUE)` after edits.
3. For scheduler controllers, inspect the scheduler output/error log and confirm queue, account, wall time, memory, CPU, module, and filesystem settings.
4. Reduce concurrency when local workers are being killed for memory pressure.

## Rerun safely

After fixing the cause, rerun the same narrow target or checkpoint selector. Successful upstream targets remain cached.

```{.r filename="R"}
targets::tar_make(
  names = tidyselect::matches("your_target.*your_aggregation")
)
```

Use an unqualified `targets::tar_make()` only when every active aggregation and enabled module is intentionally in scope. Avoid deleting the target store as a debugging step: it removes evidence and forces unrelated recomputation.


## Part: Reference


# Book: multiomeR Implementation


<!-- source: website/implementation/index.qmd -->

# Introduction

This book collects developer-facing details for readers who want to understand or modify the internals of multiomeR. It complements the [main user manual](../), which covers installation, configuration, checkpoint-based execution, and output inspection.

Use this book when you need to trace a configuration value into mapped targets, understand how the simplified graph views relate to the real `{targets}` graph, or decide where an implementation change belongs.

## Where to start

For a first implementation pass:

1. Read [Reading the graph views](graph_methodology.qmd) and follow its configuration-to-target trace.
2. Open the [main pipeline](implementation_main.qmd) graph for the modality or checkpoint you plan to change.
3. Use [Implementation conventions](implementation_conventions.qmd) to understand the relevant manifest, mapping, symbol, tag, and runtime contracts.
4. Read [Background and design philosophy](background_philosophy.qmd) when you need the rationale for the editable-workflow design.

The [differential analyses](implementation_differential_analyses.qmd) and [genetic enrichment](implementation_genetic_enrichment.qmd) chapters cover the optional module graphs.

## Common entry points

| Change | Start with |
|---|---|
| Add or revise a YAML parameter | `cfg_pipeline_parameters.tsv`, then the owning config reader or target. |
| Change GEM well preprocessing | `_targets.R` mapping plus `extra_targets/per_GEM_well_targets.R`. |
| Change aggregation GEX, ATAC, or WNN processing | The corresponding graph section and `extra_targets/*_targets.R` file. |
| Add a review boundary | Existing `[checkpoint:<name>]` description tags and the user-manual selector contract. |
| Add a graph-visible target | Existing `[part_of_graph:<graph_id>]` tags and graph-pruning rules. |
| Change resource routing | `crew_controllers.R`, `R/resource_helpers.R`, and the runtime bootstrap convention. |

If you are trying to run multiomeR rather than modify it, start with the [main manual](../).


## Part: Orientation


<!-- source: website/implementation/graph_methodology.qmd -->

# Reading the graph views



The graph chapters collect simplified views of the real `{targets}` dependency graph. They are meant to make the workflow easier to reason about before reading the target code directly.

The diagrams are generated from tagged target metadata and the real dependency graph, then simplified by pruning or bypassing lower-level nodes that would make each view harder to read. They keep real target names and preserve the dependency structure where practical, while staying compact enough to build intuition about the main control points.

The following chapters cover the main pipeline, the differential analyses module, and the genetic enrichment module.

## Trace one configured aggregation

The quickest way to understand the implementation is to follow one value across the graph:

1. `immune_human_2x` is a key in `cfg_aggregations.yaml`.
2. `read_aggregation_config_tibble()` resolves manifest defaults and inheritance into one aggregation row.
3. `build_aggregation_tibble()` filters active rows and adds symbols for GEM well-, dataset-, and aggregation-level upstream targets.
4. The root `_targets.R` passes that row through `tar_map(names = aggregation, delimiter = ".")`.
5. A base target such as `multimodal_Seurat_object` becomes `multimodal_Seurat_object.immune_human_2x`.
6. Description tags make selected targets discoverable as checkpoints or graph nodes, while structured file helpers derive output paths from the active target name.

Inspect the exact target command and description without running it:

```r
targets::tar_manifest(
  names = tidyselect::matches(
    "^multimodal_Seurat_object[.]immune_human_2x$"
  ),
  fields = c(name, command, description),
  callr_function = NULL
)
```

This trace connects the [parameter manifest](implementation_conventions.qmd#parameter-manifest), [mapping tibbles](implementation_conventions.qmd#mapping-tibbles), and [target-symbol columns](implementation_conventions.qmd#target-symbol-columns) before the larger diagrams introduce many nodes at once.

## What the diagrams omit

The curated graph views are orientation aids, not alternate target definitions. A node can be absent because it was pruned as a lower-level implementation detail, bypassed to preserve a useful dependency path, or omitted because it lacks the graph-membership tag for that view. Use `tar_manifest()` or the source target files when exact completeness matters.

[Mermaid graph omitted; source: `website/figures/standard_node_color_legend.mmd`]


<!-- source: website/implementation/implementation_conventions.qmd -->

# Implementation conventions

The implementation book assumes a few repository conventions that are easy to miss when reading individual target files in isolation. These conventions are not separate framework features; they are the small contracts that make the `{targets}` graph configurable, inspectable, and readable across GEM well-level preprocessing, aggregation-level analysis, and optional downstream modules.

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

The website renders parameter tables from the same manifest rather than maintaining a second documentation schema. This keeps the user-facing configuration reference tied to the runtime validation contract.

Module-specific YAML uses the same mechanism. Aggregations opt into modules through the aggregation config, and each enabled aggregation must have a matching module config row. Some cross-scope fallbacks are still implemented by target code rather than by manifest inheritance; for example, a module parameter may intentionally allow `NULL` and then fall back to an aggregation-level path during module setup.

## Mapping tibbles

The root `_targets.R` builds the target graph from mapping tibbles. Each mapping tibble is a row-wise contract: one row becomes one set of mapped target instances, and columns in that row become local symbols inside the corresponding `tarchetypes::tar_map()` block.

The core mapping flow is:

1. `dataset_tibble_from_yaml` is read from `cfg_datasets.yaml`.
2. `GEM_well_tibble` is read from `cfg_GEM_wells.tsv` and joined to dataset-level config.
3. `aggregation_tibble_all_from_yaml` is read from `cfg_aggregations.yaml`.
4. `aggregation_tibble` keeps active aggregations and adds upstream target-symbol columns.
5. `dataset_tibble` collapses GEM well rows back to one row per dataset.
6. `_targets.R` expands GEM wells, datasets, and aggregations with `tar_map()`, then appends module target files.

```r
tarchetypes::tar_map(
  values = GEM_well_tibble,
  names = GEM_well_ID,
  delimiter = ".",
  source("extra_targets/per_GEM_well_targets.R")$value
)
```

With `GEM_well_ID = "healthy_PBMC_human"`, a target named `cellranger_summary_file` becomes `cellranger_summary_file.healthy_PBMC_human`. The same dot-delimited suffix convention is used for datasets, aggregations, module targets, and nested module maps.

```text
cfg_GEM_wells.tsv row       -> GEM_well_tibble row    -> per GEM well targets
cfg_datasets.yaml key       -> dataset_tibble row     -> per-dataset targets
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

Column names should describe the downstream scope, the upstream target, and the fact that the value is a symbol list. Existing symbol-list columns use the `*_syms` suffix, such as `aggregation_GEX_counts_BPCells_matrix_syms`, `dataset_unfiltered_cells_n_vecs_syms`, and `per_dataset_QC_violins_syms`.

The shared pattern appears at three boundaries:

- `aggregation_*_syms` columns splice per GEM well targets into aggregation-level targets.
- `dataset_*_syms` columns splice per GEM well targets into dataset-level targets.
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


<!-- source: website/implementation/algorithm_validation.qmd -->

# Algorithmic implementations, deviations and validation

multiomeR reimplements a small number of reference algorithms so they can operate on the workflow's native matrices and graph state. This page is the maintained record of what those implementations preserve, where they deliberately differ, and what the executable validation establishes.

The evidence labels used below are intentionally narrow:

- **Reference-parity tested** means the repository and named reference implementation run on the same deterministic fixture and their returned values are compared directly.
- **Reference-similarity tested** means exact equality is not an appropriate contract, so predefined similarity thresholds are checked against the named reference implementation.
- **Algorithmically derived** means the implementation is checked against an independent mathematical result, not against another software implementation.
- **Internally checked** means a deterministic repository fixture exercises an internal contract without establishing reference parity.

Passing these fixtures does not validate every dataset, parameter regime, approximate-neighbor realization, biological interpretation, or downstream target. The [CI workflow](https://github.com/koefoeden/multiomeR/blob/main/.github/workflows/algorithm-validation.yaml) reruns the checks when the validation scripts, relevant helpers, or Pixi environment change.

| Implementation | Evidence status | Maintained reference | Current fixed-fixture result |
|---|---|---|---|
| BPCells-native UCell | Reference-parity tested | UCell 2.14.0 | Exact values, dimensions, and dimnames |
| BPCells-native AMULET | Reference-parity tested | scDblFinder 1.24.0 | Exact metrics and multi-chromosome loci, including order |
| Native WNN | Reference-similarity tested | Seurat 5.5.0 | Production-like fixture: weight Spearman 0.925; mean neighbor overlap 0.952 |
| Sparse SCAVENGE propagation | Algorithmically derived and reference-similarity tested | SCAVENGE 1.0.2 at `8ee8b173d965` | Closed-form delta 3.55e-13; score-rank Spearman 0.9966; significant-cell Jaccard 1.0 |

## BPCells-native UCell scoring

**Reference algorithm.** [`UCell::ScoreSignatures_UCell()`](https://bioconductor.org/packages/release/bioc/html/UCell.html) calculates per-cell signature scores from descending feature ranks, caps ranks at `maxRank`, combines positive and negative signatures, and clips negative combined scores to zero. The maintained comparison also covers `UCell::AddModuleScore_UCell()`.

**Reason for reimplementation.** The workflow keeps gene-by-cell counts in BPCells-backed matrices. Materializing the complete matrix in memory or building a Seurat object solely for marker scoring would discard that storage contract, so multiomeR ranks bounded cell chunks and returns metadata-ready scores directly.

**Behavior preserved.** The implementation preserves UCell signature syntax (`+` and `-` suffixes), descending per-cell ranks, configurable tie handling, `maxRank` capping, impute/skip behavior for missing genes, negative-signature weighting, lower-bound clipping, signature names, and cell order.

**Deliberate deviations and consequences.** Matrix materialization is limited to one cell chunk at a time, and optional fork workers operate across chunks. This changes memory and execution behavior but not the tested score values. The helper returns a data frame instead of mutating a Seurat object. The target-level marker validator rejects configured genes missing from the Cell Ranger reference before normal pipeline scoring, whereas the lower-level helper still exposes UCell's impute/skip modes for explicit use.

**Implementation and wiring.** The scorer is [`calculate_BPCells_UCell_scores_from_matrix()` in `R/processing_GEX_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/processing_GEX_helpers.R). [`extra_targets/general_aggregation_targets.R`](https://github.com/koefoeden/multiomeR/blob/main/extra_targets/general_aggregation_targets.R) validates `UCell_GEX_marker_genes_list`, and [`extra_targets/GEX_graph_and_cluster_targets.R`](https://github.com/koefoeden/multiomeR/blob/main/extra_targets/GEX_graph_and_cluster_targets.R) computes cell-level scores before GEX cell-type assignment. ATAC and WNN cell-type targets reuse those metadata scores.

**Validation.** [`scripts/validate_scoring_parity.R`](https://github.com/koefoeden/multiomeR/blob/main/scripts/validate_scoring_parity.R) creates a deterministic 500-gene by 37-cell matrix, writes the project input as BPCells, and compares signed signatures with imputed and skipped missing genes. It requires `identical()` values, dimensions, and dimnames against UCell 2.14.0. It also retains exact Seurat `AddModuleScore` and cell-cycle helper checks plus a metadata-join contract.

The UCell 2.14.0 reference call with imputed missing genes can emit non-fatal R stack-imbalance warnings under the repository's R 4.5 environment. A narrow diagnostic reproduced them in the UCell reference call but not in the repository scorer. The validation therefore runs the reference call in a disposable `callr` process and compares its returned matrix in the clean parent session; this isolates the package warning without weakening the equality assertion or changing the runtime.

**Status and rerun.** Reference-parity tested against UCell 2.14.0 with exact equality; passing in the current locked environment.

```bash
pixi run test-scoring-parity
```

## BPCells-native AMULET

**Reference algorithm.** [`scDblFinder::amulet()`](https://bioconductor.org/packages/release/bioc/html/scDblFinder.html) detects likely scATAC-seq doublets from the number of genomic loci covered by more than two fragments. Its underlying `getFragmentOverlaps()` implementation filters fragment sizes and excluded regions, calculates per-barcode fragment and overlap counts, removes loci recurrently covered across many cells, and derives Poisson p-values with Benjamini-Hochberg correction.

**Reason for reimplementation.** The per-GEM-well workflow already stores Cell Ranger ATAC fragments as compressed BPCells directories. Passing the original fragment TSV to scDblFinder materializes chromosome-scale `GRanges` objects and previously requested six cores and 60 GB. The local implementation streams the existing BPCells fragment target and retains only one chromosome's selected fragments while calculating coverage runs.

**Behavior preserved.** [`calculate_amulet_metrics_BPCells()`](https://github.com/koefoeden/multiomeR/blob/main/R/amulet_BPCells_helpers.R) preserves barcode selection, minimum-fragment thresholds, maximum fragment size, excluded regions, `nFrags`, `uniqFrags`, `nAbove2`, `total.nAbove2`, p-values, q-values, and high-overlap-site removal. The lower-level loci return also preserves scDblFinder's cell-major, chromosome, and coordinate ordering. Cell Ranger's inclusive end-insertion convention is shifted back by one base before calculation so the BPCells representation matches scDblFinder's BED import.

**Deliberate deviations and consequences.** Only unique-fragment operation is supported. BPCells fragment objects do not retain Cell Ranger's PCR-duplicate count column, so requesting non-unique expansion fails explicitly instead of silently changing `nFrags`. BPCells does not export its fragment iterator header; the native helper therefore mirrors that private C++ interface, verifies the exact project-pinned BPCells commit `28759cdd5125` before use, and compiles a small shared library in each worker's temporary directory. A BPCells upgrade must revalidate this interface and the exact parity fixture before updating the pin. The target is single-threaded and requests the standard 16-GB worker tier. A native-only probe on the stored `healthy_PBMC_human` fragments processed 2,711 selected cells in 15.4 seconds with 0.64 GB peak RSS, including R startup and native compilation; this supports the reduced allocation for the public fixture but is not a memory guarantee for larger datasets.

**Implementation and wiring.** Native iteration and coverage-run calculation are implemented in [`src/amulet_bpcells.cpp`](https://github.com/koefoeden/multiomeR/blob/main/src/amulet_bpcells.cpp). The R wrapper, ABI check, high-overlap filtering, and AMULET statistics are in [`R/amulet_BPCells_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/amulet_BPCells_helpers.R). [`amulet_metrics_tibble` in `extra_targets/per_GEM_well_targets.R`](https://github.com/koefoeden/multiomeR/blob/main/extra_targets/per_GEM_well_targets.R) consumes the existing prefixed BPCells fragments, restores unprefixed barcode keys, and preserves the existing downstream metrics shape.

**Validation.** [`scripts/validate_amulet_bpcells_parity.R`](https://github.com/koefoeden/multiomeR/blob/main/scripts/validate_amulet_bpcells_parity.R) compares against scDblFinder 1.24.0 in disposable `callr` processes. It requires `identical()` results for the bundled fragment-file metrics, prefixed pipeline barcodes, a deterministic 12,000-fragment multi-chromosome loci fixture, and the corresponding full AMULET metrics. It also requires an explicit error for unsupported PCR-duplicate expansion. The disposable reference processes isolate stack-imbalance warnings emitted by the current scDblFinder reference under R 4.5 without weakening the returned-object comparison.

**Status and rerun.** Reference-parity tested against scDblFinder 1.24.0 with exact equality; passing in the current locked environment.

```bash
pixi run test-amulet-parity
```

## Native weighted nearest neighbors

**Reference algorithm.** [`Seurat::FindMultiModalNeighbors()`](https://satijalab.org/seurat/reference/findmultimodalneighbors) constructs cell-specific modality weights from within- and cross-modality neighborhood prediction, collects candidate neighbors across modalities, and selects a weighted multimodal neighbor set.

**Reason for reimplementation.** The pipeline already has aligned RNA PCA/Harmony and ATAC LSI/Harmony matrices and needs reusable neighbor indices, distances, and modality weights without creating a Seurat object. Native graph state also feeds UMAP, Leiden clustering, SCAVENGE, and the optional Seurat/Signac export.

**Behavior preserved.** [`weighted_nearest_neighbors_BPCells()`](https://github.com/koefoeden/multiomeR/blob/main/R/processing_multimodal_helpers.R) preserves row-wise L2 normalization, per-modality nearest neighbors, nearest-distance correction, within/cross prediction kernels, capped modality affinity ratios, normalized cell-specific modality weights, candidate-set union, weighted neighbor ranking, and Seurat's transformation from weighted affinity to neighbor distance.

**Deliberate deviations and consequences.** BPCells HNSW replaces Seurat's Annoy search, so approximate candidate sets need not be identical. Adaptive bandwidth uses the corrected distance to the configured `k`th neighbor rather than Seurat 5.5.0's default SNN-derived far-neighbor bandwidth. The helper does not expose Seurat's optional smoothing or cross-constant list, and BPCells builds downstream SNN state rather than storing Seurat `Neighbor` and `Graph` objects. Its current `seed` argument is not consulted by the HNSW calls, so it must not be interpreted as controlling neighbor-search randomness. These choices can change weights, selected neighbors, SNN edges, clusters, and UMAP coordinates; correlation and overlap are therefore the validation contract, not exact equality.

**Implementation and wiring.** The implementation and graph consumers are in [`R/processing_multimodal_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/processing_multimodal_helpers.R). [`extra_targets/WNN_targets.R`](https://github.com/koefoeden/multiomeR/blob/main/extra_targets/WNN_targets.R) aligns modality embeddings, creates `WNN_results_raw`, filters small clusters, optionally recomputes `WNN_results`, and wires that state into WNN UMAP, clustering, metadata, and cell-type targets.

**Validation.** [`scripts/validate_algorithm_deviations.R`](https://github.com/koefoeden/multiomeR/blob/main/scripts/validate_algorithm_deviations.R) compares two deterministic RNA/ATAC fixtures with Seurat 5.5.0. The headline production-like fixture uses 400 cells, 12 dimensions, `k = 30`, and candidate range 200, matching the pipeline's configured neighbor count and native WNN candidate range. It requires Spearman correlation at least 0.90 for each modality weight, mean neighbor-set overlap at least 0.94, and first-quartile overlap at least 0.90; observed values are 0.925, 0.952, and 0.933. A smaller 160-cell stress fixture retains opposing per-cell noise gradients at `k = 15` and candidate range 50. Its deliberately looser thresholds are 0.80, 0.85, and 0.75; observed values are 0.851, 0.924, and 0.867. The current locked BPCells 0.3.1 build is pinned at `28759cdd5125`.

**Status and rerun.** Reference-similarity tested against Seurat 5.5.0; passing. This does not assert parity of SNN weights, clustering, or UMAP.

```bash
pixi run Rscript scripts/validate_algorithm_deviations.R
```

## Sparse SCAVENGE propagation and significance

**Reference algorithm.** [SCAVENGE 1.0.2 at commit `8ee8b173d965`](https://github.com/sankaranlab/SCAVENGE/tree/8ee8b173d965009a696b2a590d5b17b28b7cf851) selects high chromVAR Z-score seed cells, performs a column-normalized random walk with restart, caps and rescales the propagation score into a trait relevance score (TRS), and uses degree-matched seed permutations to identify significant cells.

**Reason for reimplementation.** The reference package's last commit and dependency stack predate the pipeline's current R/Bioconductor environment. multiomeR needs sparse propagation over native RNA, ATAC, and WNN SNN matrices and must avoid materializing a cell-by-permutation score matrix for large cell sets.

**Behavior preserved.** The implementation preserves the one-sided Z-score seed threshold and top-percent cap, column-normalized transition matrix, equal seed restart mass, iterative random walk, 0.95 propagation-score cap, min-max scaling, Z-score scale factor, degree-matched seed sampling, and per-cell comparison with permuted propagation scores.

**Deliberate deviations and consequences.** The reference workflow normally uses a mutual-kNN graph; multiomeR propagates over BPCells-derived SNN graphs from RNA, ATAC, or WNN state. Seed and scale-factor helpers guarantee at least one selected cell for small inputs and break boundary ties deterministically. The random walk validates graph inputs and has a maximum-iteration guard. Permutations accumulate exceedance counts rather than retaining every score vector, and report add-one empirical p-values, `(exceedances + 1) / (permutations + 1)`, with `p < 0.05`; the reference thresholds raw exceedance counts. These graph, boundary, and p-value choices can alter TRS values and significant-cell calls, especially near thresholds.

**Implementation and wiring.** Seed selection, sparse random walk, streaming degree-matched permutations, and TRS construction are in [`R/gchromvar_SCAVENGE_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/gchromvar_SCAVENGE_helpers.R). [`module_genetic_enrichment/SCAVENGE_graph_targets.R`](https://github.com/koefoeden/multiomeR/blob/main/module_genetic_enrichment/SCAVENGE_graph_targets.R) constructs each graph and maps chromVAR Z-score records into `TRS_tibbles`; [`SCAVENGE_group_targets.R`](https://github.com/koefoeden/multiomeR/blob/main/module_genetic_enrichment/SCAVENGE_group_targets.R) wires group summaries and plots.

**Validation.** The deterministic 60-cell fixture contains three disconnected, regular graph blocks and three enriched seeds. First, the iterative sparse random walk is compared with the closed-form solution

\[
s = r\left(I - (1-r)P\right)^{-1}p_0,
\]

with a maximum absolute tolerance of 1e-10; the current delta is 3.55e-13. Second, compact local reference functions reproduce the relevant SCAVENGE 1.0.2 code at the pinned commit without installing its historical dependency stack. With 199 fixed-RNG permutations, the test requires TRS rank Spearman correlation at least 0.99 and significant-cell Jaccard overlap at least 0.80. Current values are 0.9966 and 1.0.

**Status and rerun.** Random-walk propagation is algorithmically derived against the closed form, and score ranks/significant-cell calls are reference-similarity tested on one small deterministic fixture. This is not full reference-package parity and does not validate graph construction, chromVAR inputs, or genome-scale permutation stability.

```bash
pixi run Rscript scripts/validate_algorithm_deviations.R
```

Run the complete maintained set, including UCell and AMULET, with:

```bash
pixi run test-algorithm-validation
```


## Part: Background


<!-- source: website/implementation/background_philosophy.qmd -->

<!-- begin include: website/background_philosophy.qmd -->

# Background and design philosophy

## Philosophy

multiomeR challenges the commonly used single-cell & single-nucleus processing workflows, where most processing steps are run sequentially and are implemented in hard-to-customize packaged functions. Widely used examples include brilliant packages such as Seurat and Signac. However, while offering superb ease-of-use, this traditional workflow has two important downsides:

1) The built-in processing steps (usually functions) are only customizable to levels chosen by the package developers, which often falls short in the face of the enormous variance and complexity of real-life single-cell and single-nucleus sequencing datasets
2) The traditional workflows leave a lot of performance on the table, since they rarely take advantage of the inherent parallelism when processing these large sequencing datasets. This includes commonly found multiplicities of independent experimental samples (donors, GEM wells), data modalities, data representations, downstream analyses, etc.

Until recently, this simplified and easy-to-use workflow has been necessary to make the tools practically usable in a busy research environment. However, with the recent advent of powerful, AI coding agents, which can quickly ingest large complex codebases, this inherent tug-of-war between performance/customizability and ease-of-use/code interpretability has started to move the needle in favor of complex, performant, and customizable approaches. This is especially relevant because single-cell and single-nucleus datasets continue to grow faster than routine compute capacity. These observations, together with a lack of an established framework for the processing of Multiome data, constitute the main motivation for the development of this pipeline.

## Achieving customizability

By not packaging the pipeline into an ordinary R package, users retain straightforward access to the codebase and can modify any part of the workflow to fit their needs. The pipeline thus functions as a template of suggested best practices, which can be used as-is, modified, or extended as needed.

## Achieving performance

multiomeR achieves its high performance through two key design choices:

1) The workflow is implemented as a `targets` pipeline, which explicitly represents each processing and analysis step as a node, and the dependencies between them as edges, in a directed acyclic graph (DAG). This allows `targets` to automatically determine which parts of the workflow can be run in parallel, and which parts need to run sequentially - and whether they need to run at all through the use of cached data and intelligent invalidation.
2) multiomeR builds heavily on the BPCells R package, which is designed from the ground up for computational efficiency through lazy, streaming data processing and on-disk storage. Please see the BPCells documentation for more details.

<!-- end include: website/background_philosophy.qmd -->


## Part: Target graph views


<!-- source: website/implementation/implementation_main.qmd -->

<!-- begin include: website/implementation_main.qmd -->

# Main pipeline



The root `_targets.R` creates dataset, GEM well, and aggregation mapping rows, then maps target fragments from `extra_targets/`. Use the diagrams to find the relevant stage, then inspect the corresponding source file for the complete command and resource declaration.

| Stage | Primary source |
|---|---|
| GEM well preprocessing | `extra_targets/per_GEM_well_targets.R` |
| Dataset summaries | `extra_targets/per_dataset_targets.R` |
| Aggregation setup and shared QC | `extra_targets/general_aggregation_targets.R` |
| GEX | `extra_targets/GEX_merge_and_dim_reduc_targets.R`, `extra_targets/GEX_graph_and_cluster_targets.R` |
| ATAC | `extra_targets/ATAC_targets.R` |
| WNN | `extra_targets/WNN_targets.R` |
| Compatibility export | `extra_targets/Seurat_Signac_export_targets.R` |
| Subgroups | `extra_targets/subgroups.R` |

## Parallel pre-processing

This view covers per GEM well processing and QC, including optional ambient RNA correction, donor demultiplexing, doublet detection, barcode filtering, and handoffs into aggregation-level GEX and ATAC objects.

[Mermaid graph omitted; source: `website/figures/human_curated/parallel_v2.mmd`]

## GEX processing

This view covers merged RNA processing, clustering, marker detection, cell type annotation, and GEX review outputs.

[Mermaid graph omitted; source: `website/figures/human_curated/GEX_v2.mmd`]

## ATAC processing

This view covers ATAC QC, peak calling, consensus peak construction, chromatin accessibility processing, chromVAR scoring, coverage tracks, and peak-gene links.

[Mermaid graph omitted; source: `website/figures/human_curated/ATAC_v2.mmd`]

## WNN integration

This view covers GEX and ATAC embedding handoffs, WNN integration, modality weights, cluster comparison, and integrated metadata outputs.

The native implementation, its differences from Seurat, and the maintained similarity thresholds are recorded in [Algorithmic implementations, deviations and validation](algorithm_validation.qmd#native-weighted-nearest-neighbors).

[Mermaid graph omitted; source: `website/figures/human_curated/WNN_v2.mmd`]

## Subgroup reprocessing

This view covers optional subgroup-native GEX, ATAC, and WNN reprocessing for sufficiently large parent groups.

[Mermaid graph omitted; source: `website/figures/human_curated/full_subgroups_v2.mmd`]

Subgroup reprocessing is configuration-dependent and should not be treated as part of the minimum main-pipeline path.

<!-- end include: website/implementation_main.qmd -->


<!-- source: website/implementation/implementation_differential_analyses.qmd -->

<!-- begin include: website/implementation_differential_analyses.qmd -->

# Differential analyses



`module_differential_analyses/targets.R` filters aggregations that enabled the module, joins their module config, attaches symbols for accepted WNN metadata and pseudobulk inputs, and maps the composition, pseudobulk, GSEA, and cross-modality target fragments.

The graph below is an orientation view. Inspect `setup_and_DCTC_targets.R`, `psbulk_DX_targets.R`, `GSEA_targets.R`, and `cross_modality_targets.R` for the complete model and plotting commands. The user-facing prerequisites and checkpoint selector are documented in [Differential analyses](../downstream_differential_analyses.html).

[Mermaid graph omitted; source: `website/figures/human_curated/differential_analyses_v2.mmd`]

<!-- end include: website/implementation_differential_analyses.qmd -->


<!-- source: website/implementation/implementation_genetic_enrichment.qmd -->

<!-- begin include: website/implementation_genetic_enrichment.qmd -->

# Genetic enrichment



`module_genetic_enrichment/targets.R` filters enabled human aggregations, resolves one configured Open Targets study set per aggregation, and attaches symbols for WNN metadata, graphs, embeddings, consensus peaks, chromVAR state, and ATAC fragments.

The main target fragments live in `setup_targets.R`, `gchromVAR_targets.R`, `SCAVENGE_graph_targets.R`, `SCAVENGE_group_targets.R`, `GWAS_chromVAR_cell_type_targets.R`, and `GWAS_chromVAR_contribution_targets.R`. The user-facing release, method-selection, and interpretation contracts are documented in [Genetic enrichment](../downstream_genetic_enrichment.html).

## Single-nucleus and graph-based enrichment

This view covers the configured GWAS inputs, single-nucleus enrichment state, graph propagation, and downstream trait summaries. Additional cell-type contribution and locus-attribution branches may be pruned from this compact orientation view; use the manifest for the complete graph.

The sparse SCAVENGE reimplementation, deliberate graph and permutation differences, and validation evidence are recorded in [Algorithmic implementations, deviations and validation](algorithm_validation.qmd#sparse-scavenge-propagation-and-significance).

[Mermaid graph omitted; source: `website/figures/human_curated/genetic_enrichment_single_nucleus_v2.mmd`]

<!-- end include: website/implementation_genetic_enrichment.qmd -->


# Orphaned QMD Pages

Tracked QMD files not reached from the Quarto book graph or include graph.

- `website/helpers/_targets_graph_snippet.qmd`
- `website/implementation_overview.qmd`
