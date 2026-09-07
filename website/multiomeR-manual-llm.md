# multiomeR Website LLM Export

This file is generated from the Quarto book outlines and resolves Quarto include shortcodes.
Hidden setup chunks, generated helper chunks, Mermaid graph bodies, and verbose image metadata are omitted by default.


# Book: multiomeR Manual


## Part: Start here


<!-- source: website/index.qmd -->

[Image omitted; source: `figures/multiomeR-logo.svg`; alt: Image]

# Start here

## What is multiomeR?
multiomeR is a [targets](https://books.ropensci.org/targets/)-based workflow for processing and analyzing single-nucleus 10x Genomics Multiome ATAC + Gene Expression datasets. It is meant to be adapted to your own data, compute setup, and biological questions.

The workflow starts from `cellranger-arc count` outputs, processes gene-expression (GEX) and ATAC data, builds multimodal aggregations, and supports optional downstream modules for differential analyses and genetic enrichment for human datasets.

## Your first analysis

Start with the public demo: two human GEM wells with supplied configuration.
You will install the software, run one joint analysis, and read its cell
metadata and multimodal object. This gives you a working example before you
choose settings for your own study.

You need basic R skills, a Linux terminal, and a machine with sufficient
[memory and disk space](demo_installation.qmd#system-requirements). You do not
need to know how to write a `targets` pipeline. Commands labelled **Bash** run
in the terminal; commands labelled **R** run in the R session opened during
installation. Run both from the repository folder unless stated otherwise.

## Find what you need

| If you want to... | Start here |
|---|---|
| See what the workflow produces | Browse the [main pipeline output gallery](gallery_main.qmd). |
| Try multiomeR on public data | Follow the three-part quickstart: [install](demo_installation.qmd), [run](demo_running.qmd), then [inspect the outputs](demo_outputs.qmd). |
| Configure your own data | Read the [main-pipeline overview](main_overview.qmd), prepare the [configuration and inputs](main_inputs.qmd), then [run one aggregation](main_running.qmd). |
| Add a downstream analysis | Check the prerequisites for [differential analyses](downstream_differential_analyses.qmd) or [genetic enrichment](downstream_genetic_enrichment.qmd). |
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

## Development status

multiomeR is in beta and may introduce breaking changes. The [numbered QC reviews](main_running.qmd#qc-checkpoints)
show what to inspect before each analysis stage. Acceptance criteria still
depend on the tissue, study design, and intended use. Report problems
or questions through [GitHub issues](https://github.com/koefoeden/multiomeR/issues).

Continue to [Install and prepare the demo](demo_installation.qmd).


## Part: 1. Try the public demo


<!-- source: website/demo_installation.qmd -->

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

The committed `crew_controllers.R` provides a local setup for a 16-CPU,
256-GB workstation and can run several workers concurrently. Review
[Distributed computing](performance_distributed_computing.qmd) before running
on a smaller machine or a scheduler.


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

Continue to [Run the demo](demo_running.qmd) from the R prompt.


<!-- source: website/demo_running.qmd -->

# Run the demo



In the R session opened during installation, run the command below to process
`immune_human_2x`. It combines the two demo GEM wells and produces a
Seurat/Signac object containing the multimodal results.

`names` selects the final result by its exact name using `all_of()`.
`tar_make()` also builds the dependencies needed for that result, but does
not build every plot in the gallery.

```{.r filename="R"}
targets::tar_make(
  names = tidyselect::all_of("multimodal_Seurat_object.immune_human_2x")
)
```

Keep the R session open until the command finishes. Progress messages report
targets being dispatched, completed, or skipped because they are already
up to date. A previous demo run on an AMD EPYC 7543 system with 16 logical
threads took under 25 minutes and wrote about 6 GB; your runtime may differ.

## Confirm success

After the run, the following command should return `character(0)`, meaning the requested endpoint and its dependencies are up to date:

```{.r filename="R"}
targets::tar_outdated(
  names = tidyselect::all_of("multimodal_Seurat_object.immune_human_2x"),
  callr_function = NULL
)
```

If names are returned, those results still need building. If the run failed,
follow [Troubleshooting](troubleshooting.qmd), fix the reported cause, and run
the same `tar_make()` command again. Completed results can be reused.

Continue to [Inspect the demo results](demo_outputs.qmd) to read the object
and find the associated files.


<!-- source: website/demo_outputs.qmd -->

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

demo_object <- targets::tar_read(multimodal_Seurat_object.immune_human_2x)
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

Open the [main pipeline gallery](gallery_main.qmd) to see example QC,
RNA, ATAC, and WNN plots. These are saved documentation snapshots, so you can
browse them without running the demo. They do not show the state of your own
analysis.

Requested plots use the same scope-based layout under `outputs/plots/`.

The final-object demo command does not build all these plots. See
[Request an additional result](main_running.qmd#request-an-additional-result)
for how to build a named gallery target, or follow the
[numbered QC reviews](main_running.qmd#qc-checkpoints) for the broader diagnostics.

## Next steps

- To adopt the workflow, continue with the [main-pipeline overview](main_overview.qmd) and [configuration walkthrough](main_inputs.qmd).
- To request more results or rerun after a change, use [Run your analysis](main_running.qmd).
- To diagnose a failed or unexpectedly stale target, use [Troubleshooting](troubleshooting.qmd).


## Part: 2. Analyze your own data


<!-- source: website/main_overview.qmd -->

# Plan your analysis

After trying the demo, choose one group of GEM wells to analyze together.
This group is an **aggregation**. It should reflect the biological comparison
you intend to make; donor identity, sample preparation, and batch information
need to remain distinguishable within it.

## What you need to prepare

- Cell Ranger ARC outputs for each GEM well and the `reference.json` used to
  produce them. All GEM wells in an aggregation must use the same reference.
- A GEM-well table linking stable identifiers to input paths and library
  annotations.
- A donor table with one row per donor, including the phenotypes or covariates
  you will use.
- Marker genes appropriate for the expected cell types, plus initial processing
  settings that you will review against your data.

Start with the main processing workflow. Donor demultiplexing, CellBender,
Harmony batch correction, and downstream modules require additional inputs
or choices; enable them when your study needs them.

## How the data move through the workflow

| Stage | What it does |
|---|---|
| GEM well processing | Reads Cell Ranger matrices and fragments, calculates QC metrics, and prefixes barcodes with the GEM well identifier. |
| GEX (gene expression) | Combines selected GEM wells, reduces dimensions, clusters nuclei, and annotates cell types using marker genes. |
| ATAC (chromatin accessibility) | Combines fragments, defines peaks, builds the peak-count matrix, and summarizes accessibility and motif signals. |
| WNN (weighted nearest neighbors) | Combines RNA and ATAC representations to produce integrated clusters, metadata, and a Seurat/Signac export. |

The pipeline also produces comparisons of related GEM wells before
aggregation. The [main gallery](gallery_main.qmd) shows representative
outputs; the [implementation graph](implementation/implementation_main.html)
provides the detailed computational dependencies when you need them.

Continue to [Configuration and inputs](main_inputs.qmd) to connect your
files, then follow the [numbered QC reviews](main_running.qmd#qc-checkpoints).


<!-- source: website/main_inputs.qmd -->

# Configuration and inputs



multiomeR uses two linked configuration layers. A **GEM well** points to one `cellranger-arc count` output and defines its pre-aggregation processing and QC. An **aggregation** selects GEM wells for joint GEX, ATAC, and WNN analysis.

::: {.scrollable-table}

| Layer | Configuration | Key relationship |
|------------------------|------------------------|------------------------|
| GEM well | `cfg_GEM_wells.tsv` | Aggregations refer to one or more `GEM_well_ID` values. |
| Aggregation | `cfg_aggregations.yaml` | Selects GEM wells and points to donor-level metadata. |

:::

The committed `cfg_GEM_wells.tsv` and `cfg_aggregations.yaml` files are the
active configuration files. They enable only the two GEM wells and the
`immune_human_2x` aggregation used by the public quickstart. Edit these files
directly when configuring another project.

## Start with one explicitly scoped analysis

`GEM_well_is_active` controls per-GEM-well graph construction, while aggregation
`is_active` controls aggregation graph construction. Every active aggregation
must reference active GEM wells. Before an unqualified `targets::tar_make()`,
deactivate every GEM well and aggregation you are not ready to run. The
committed quickstart configuration already follows this rule.

The example below describes one non-multiplexed GEM well from one donor.
Replace the paths, identifiers, and marker genes for your study. It illustrates
how the files connect; one donor is not enough for a replicated comparison.

### 1. Define one GEM well

Add a row to `cfg_GEM_wells.tsv` with these values. This vertical view is a
reading aid; the saved TSV has one GEM well per row.

::: {.scrollable-table}

| Column | Example value |
|---|---|
| `GEM_well_ID` | `your_GEM_well` |
| `GEM_well_dataset` | `your_dataset` |
| `GEM_well_donor_id` | `donor_1` |
| `GEM_well_n_donors` | `1` |
| `GEM_well_cellranger_arc_count_dir` | `/path/to/your_GEM_well` |
| `GEM_well_cellranger_arc_reference_json` | `/path/to/reference.json` |
| `GEM_well_add_cellbender` | `FALSE` |
| `GEM_well_cellbender_h5_file` | `NA` |
| `GEM_well_donors_VCF_file` | `NA` |
| `GEM_well_is_active` | `TRUE` |
| `GEM_well_multiplex_batch` | `batch_1` |

:::

Set `GEM_well_QC_exclude_list` to the exclusions chosen for your data. Follow the
[numbered QC reviews](main_running.qmd#qc-checkpoints); do not treat the demo's
numerical cutoffs as recommendations for your tissue.

`GEM_well_cellranger_arc_reference_json` must point to the `reference.json` from the exact Cell Ranger ARC reference used to create that GEM well's output. multiomeR checks that the JSON genome matches the feature HDF5 and rejects aggregations whose GEM wells use different references.

`GEM_well_QC_exclude_list` contains zero or more complete R filter expressions separated by `;;`. Expressions are evaluated individually against per-barcode metadata, preserving their order and their separate exclusion reasons. An empty field applies no pre-aggregation QC filters.

AMULET doublet detection runs as part of the standard QC calculation; see the
[implementation and validation](implementation/algorithm_validation.html#bpcells-native-amulet)
for details.

`GEM_well_cellranger_arc_count_dir` points to the directory containing `outs/`, not to `outs/` itself. The baseline pipeline requires:

``` {.text filename="Text"}
<GEM_well_cellranger_arc_count_dir>/outs/
|-- summary.csv
|-- filtered_feature_bc_matrix.h5
|-- atac_fragments.tsv.gz
|-- atac_fragments.tsv.gz.tbi
`-- per_barcode_metrics.csv
```

If `GEM_well_donors_VCF_file` is configured, `atac_possorted_bam.bam` is also required for `cellsnp-lite`. Without a VCF, the pipeline skips genotype demultiplexing and assigns `GEM_well_donor_id` to every called nucleus. That donor ID must match the donor metadata table.

### 2. Create keyed donor metadata

The donor metadata table must contain one unique row per `donor_id`:

``` {.text filename="donor_metadata.tsv"}
donor_id	condition
donor_1	control
```

Put donor-specific phenotypes and covariates in the donor table. Put library-, run-, or batch-specific variables directly in `cfg_GEM_wells.tsv`, using a `GEM_well_` prefix. Apart from their key columns, donor and GEM-well metadata must not reuse column names.

### 3. Define one aggregation

``` {.yaml filename="cfg_aggregations.yaml"}
your_aggregation:
  aggregation_GEM_well_IDs: [your_GEM_well]
  aggregation_donor_id_metadata_tsv: /path/to/donor_metadata.tsv
  aggregation_GEX_marker_genes:
    Cell_type_A: [GENE1, GENE2]
    Cell_type_B: [GENE3, GENE4]
  is_active: true
```

Replace `GENE1`–`GENE4` with gene symbols appropriate for the tissue and
reference. Omit `modules` for the first run. After reviewing the main results,
you can enable optional analyses by adding their module names here and a
matching aggregation entry in each module's configuration.

### 4. Validate before running

From the repository-root R session, construct the graph and inspect the targets created for the aggregation:

``` {.r filename="R"}
manifest <- targets::tar_manifest(callr_function = NULL)

manifest |>
  dplyr::filter(stringr::str_ends(name, ".your_aggregation")) |>
  dplyr::select(name, description)
```

Manifest construction validates the YAML parameter schema, aggregation references to GEM wells, module names and rows, and controller setup. Metadata file contents are validated when their targets run. Fix manifest-time errors before calling `tar_make()`, then continue to [Running the workflow](main_running.qmd).

## Configuration reference

The searchable overviews below are generated from `cfg_pipeline_parameters.tsv`, the same manifest used for runtime defaults and validation. Use them to change a default after the minimum configuration works.

### GEM well columns

Use the [GEM-well example](https://github.com/koefoeden/multiomeR/blob/main/cfg_GEM_wells.tsv)
to check the full set of columns.

The table below is a documentation snapshot of the two public demo wells,
showing the core columns and one optional annotation. Bold columns must be
present in the TSV; some allow an NA value. Scroll horizontally and focus or
hover over a column's **i** button for its meaning. The other inactive rows and
metadata columns in the public configuration remain available as examples.

[Generated Quarto chunk omitted: `emit_GEM_well_demo_table( GEM_well_config_file = "website/data/demo_GEM_wells.tsv", dictionary_file = "website/data/G...`]

### Aggregation parameters

See the [aggregation example](https://github.com/koefoeden/multiomeR/blob/main/cfg_aggregations.yaml)
for a complete configuration. Search by name or purpose, or choose a topic.
Defaults are visible beside each parameter; open a row for its type and example.

[Generated Quarto chunk omitted: `emit_parameter_overview("aggregation")`]

<details>

<summary>Show the public <code>immune_human_2x</code> example</summary>

[Generated Quarto chunk omitted: `emit_yaml_entry(aggregations_config_file, "immune_human_2x")`]

</details>


<!-- source: website/main_running.qmd -->

# Adapt and run the demo workflow

Use the working two-GEM-well demo as the starting point, then replace its
inputs and settings one decision at a time. Each step follows the same loop:

1. edit only the configuration needed for the next decision;
2. preview and run the corresponding `targets` checkpoint;
3. inspect the named outputs; and
4. either revise the settings and rerun, or accept the result and continue.

The examples below use `my_dataset` and `my_aggregation`. Replace these with
the identifiers in your configuration. Commands are run from a repository-root
R session after completing [Install and prepare the demo](demo_installation.qmd).
Keep [Configuration and inputs](main_inputs.qmd) open for the complete TSV and
YAML structures.

The numbered checkpoints are review boundaries, not universal acceptance
criteria. Use thresholds justified for your tissue and sampling design. Each
preview must return targets before you run it; an empty selection usually means
the dataset or aggregation suffix is wrong. Descriptive tags select outputs;
`targets` still builds their upstream dependencies.

Before starting, make sure `crew_controllers.R` describes the computer or
scheduler you intend to use. See
[Distributed computing](performance_distributed_computing.qmd) for that
configuration.

## 1. Review pre-aggregation QC {#qc-checkpoints}

### Add GEM wells and inspect distributions

Start by adding one row per `cellranger-arc count` output to
`cfg_GEM_wells.tsv`. For each new row:

1. assign a unique `GEM_well_ID` and a shared `GEM_well_dataset` for GEM wells
   that should be compared during pre-aggregation QC;
2. set the count-output directory, Cell Ranger reference, and donor fields;
3. set `GEM_well_QC_exclude_list` to `NA` so that custom thresholds are not
   applied during the first inspection;
4. add any required `GEM_well_metadata_` columns; and
5. set `GEM_well_is_active` to `TRUE`.

The committed demo rows already contain reviewed example thresholds and remain
ready to run. The `NA` starting point applies when adding or adapting rows for
new data.

::: {.callout-note title="Required background knowledge"}
The pipeline first joins GEX, ATAC, Cell Ranger, AMULET, and optional Vireo
metrics for each barcode. The pre-aggregation QC checkpoint exposes the
Cell Ranger-called metadata and cross-GEM-well QC distributions before
`GEM_well_QC_exclude_list` is applied. The new GEM wells do not need to be
added to an aggregation yet.
:::

Preview the existing targets selected for this dataset, then run them:

```{.r filename="R"}
targets::tar_manifest(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:1_pre-aggregation-QC]")
  ) & tidyselect::ends_with(".my_dataset"),
  callr_function = NULL
)[, c("name", "description")]

targets::tar_make(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:1_pre-aggregation-QC]")
  ) & tidyselect::ends_with(".my_dataset")
)
```

Use the two-GEM-well demo's [RNA count
distributions](gallery_main.qmd#pre-aggregation-qc-rna-counts) and [TSS
enrichment distributions](gallery_main.qmd#pre-aggregation-qc-tss-enrichment)
as examples of the generated `per_dataset_QC_violins.my_dataset` plot family.
When your distributions need closer investigation, read their source metadata:

```{.r filename="R"}
targets::tar_read(per_dataset_cellranger_kept_metadata_tibble.my_dataset)
```

Compare QC distributions between GEM wells and look for sample-specific tails,
missing metrics, or plausible biological populations that a threshold would
remove. Values such as RNA counts, mitochondrial fraction, TSS enrichment, and
nucleosome signal do not have universally appropriate cutoffs.

When the distributions are understood, add complete R exclusion expressions
to each row's `GEM_well_QC_exclude_list`, separated by ` ;; `:

```{.text filename="cfg_GEM_wells.tsv"}
TSS.enrichment < 4 ;; nucleosome_signal > 4 ;; nCount_RNA < 250
```

Continue only when every active GEM well has either a justified filter or an
intentional `NA` value.

### Define and approve the aggregation input

Now add the metadata and aggregation that connect the accepted GEM wells:

1. add one row per donor to a donor metadata TSV, keyed by `donor_id`;
2. add an aggregation entry to `cfg_aggregations.yaml`;
3. set `aggregation_donor_id_metadata_tsv` and list the intended
   `GEM_well_ID` values under `aggregation_GEM_well_IDs`; and
4. set the aggregation's `is_active` field to `true`.

Use the demo aggregation as the template for schema-required settings, but
defer tuning marker genes and analysis parameters until their review steps.

::: {.callout-note title="Required background knowledge"}
This checkpoint applies the per-GEM-well exclusion expressions, combines the
selected GEM wells, and stops before GEX dimensionality reduction. Its UpSet
plots show overlapping exclusion reasons; the retained metadata shows the
cells that would enter GEX.
:::

```{.r filename="R"}
targets::tar_manifest(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:1_pre-aggregation-QC]")
  ) & tidyselect::ends_with(".my_aggregation"),
  callr_function = NULL
)[, c("name", "description")]

targets::tar_make(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:1_pre-aggregation-QC]")
  ) & tidyselect::ends_with(".my_aggregation")
)
```

Review the two-GEM-well demo's [Cell Ranger-called exclusion-overlap
plot](gallery_main.qmd#aggregation-input-qc-exclusions). The corresponding
all-barcode plot is an advanced diagnostic for investigating disagreement
between Cell Ranger and other barcode calls.

```{.r filename="R"}
cells_before_filtering <-
  targets::tar_read(per_dataset_cellranger_kept_metadata_tibble.my_dataset) |>
  dplyr::count(GEM_well_ID, name = "cells_before_filtering")

cells_after_filtering <-
  targets::tar_read(GEX_cellranger_kept_metadata_tibble.my_aggregation) |>
  dplyr::count(GEM_well_ID, name = "retained_cells")

cells_after_filtering |>
  dplyr::left_join(cells_before_filtering, by = "GEM_well_ID") |>
  dplyr::mutate(retained_fraction = retained_cells / cells_before_filtering)
```

Check whether one GEM well or donor loses an unexpected fraction of its cells
and whether exclusion reasons overlap as intended. Revise
`GEM_well_QC_exclude_list` and rerun this checkpoint until the retained input is
credible.

## 2. Review GEX PCA before clustering

Configure the normalization method, variable-gene selection, PCA dimensions,
and any GEX Harmony covariates before constructing the neighbour graph.

```{.r filename="R"}
targets::tar_manifest(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:2_GEX-PCA-QC]")
  ) & tidyselect::ends_with(".my_aggregation"),
  callr_function = NULL
)[, c("name", "description")]

targets::tar_make(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:2_GEX-PCA-QC]")
  ) & tidyselect::ends_with(".my_aggregation")
)
```

Review variable-gene variance, gene loadings, the PCA singular-value elbow,
embedding spread, and associations with biological and technical metadata.
Check for residual batch or QC effects before choosing PCs. Harmony coordinate
spread is not explained variance. Revise the settings and rerun this checkpoint
before committing to the GEX neighbour graph and clustering.

## 3. Review GEX clusters and cell types

Review the GEX settings for the biological system before this run. In
particular, configure:

- `aggregation_GEX_marker_genes`;
- the requested PCA dimensions and neighbour settings;
- GEX Harmony variables and clustering resolution; and
- GEX categorical and continuous variables used in review plots.

For an initial inspection of scDblFinder evidence, configure:

```{.yaml filename="cfg_aggregations.yaml"}
aggregation_scDblFinder_GEX_remove_called_doublets: false
aggregation_scDblFinder_GEX_max_doublet_fraction_per_cluster: null
```

::: {.callout-note title="Required background knowledge"}
This checkpoint uses the reviewed PCA/Harmony representation to construct
the neighbour graph, clusters, marker scores, cell-type annotations, and
scDblFinder analysis. The accepted GEX cell set and cell types
are subsequently used for ATAC peak calling.
:::

```{.r filename="R"}
targets::tar_manifest(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:3_GEX-QC]")
  ) & tidyselect::ends_with(".my_aggregation"),
  callr_function = NULL
)[, c("name", "description")]

targets::tar_make(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:3_GEX-QC]")
  ) & tidyselect::ends_with(".my_aggregation")
)
```

Review:

- PCA, Harmony, metadata-association diagnostics, and the resulting [GEX
  UMAP](gallery_main.qmd#gex-harmony-umap);
- cluster markers, marker-module scores, and the [marker dot
  plot](gallery_main.qmd#gex-marker-dot-plot);
- categorical composition across clusters and GEM wells;
- the [GEX scDblFinder score
  distributions](gallery_main.qmd#gex-scdblfinder-scores); and
- the cell annotations in `metadata_w_cell_types_tibble.GEX.my_aggregation`.

Set the desired cell- and cluster-level GEX scDblFinder policy only after
reviewing its scores, then rerun the same checkpoint. Continue when this is the
GEX cell set and annotation that should guide peak calling.

## 4. Review peak-based ATAC QC

For the first peak-QC inspection, omit
`aggregation_QC_exclude_list_combined_object` or set it to `null`.

::: {.callout-note title="Required background knowledge"}
Peak-based QC is evaluated after GEX because peak calling uses the accepted GEX
cell set and configured grouping. This checkpoint calculates ATAC QC metrics
from the resulting consensus peak matrix and fragments, but stops before the
aggregation-level ATAC exclusion expressions are applied.
:::

```{.r filename="R"}
targets::tar_manifest(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:4_peak-QC]")
  ) & tidyselect::ends_with(".my_aggregation"),
  callr_function = NULL
)[, c("name", "description")]

targets::tar_make(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:4_peak-QC]")
  ) & tidyselect::ends_with(".my_aggregation")
)
```

Use the two-GEM-well demo's [peak-based QC
distributions](gallery_main.qmd#peak-qc-distributions) as the visual reference.
Read your source metadata when individual distributions need closer
investigation:

```{.r filename="R"}
targets::tar_read(metadata_w_QC_tibble.ATAC.my_aggregation)
```

Compare peak counts, fraction of fragments in peaks, blacklist fraction, and
the other configured peak-based metrics across GEM wells. Then add justified
expressions to `aggregation_QC_exclude_list_combined_object`, for example:

```{.yaml filename="cfg_aggregations.yaml"}
aggregation_QC_exclude_list_combined_object:
  - nCount_ATAC < 1000
  - atac_peak_counts_frac < 0.1
  - atac_peak_counts_blacklist_frac > 0.01
```

## 5. Review the filtered ATAC input before LSI

::: {.callout-note title="Required background knowledge"}
This checkpoint applies the configured peak-based exclusion expressions and
exposes both their overlap and the retained metadata before LSI and ATAC
clustering.
:::

```{.r filename="R"}
targets::tar_manifest(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:5_pre-LSI-QC]")
  ) & tidyselect::ends_with(".my_aggregation"),
  callr_function = NULL
)[, c("name", "description")]

targets::tar_make(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:5_pre-LSI-QC]")
  ) & tidyselect::ends_with(".my_aggregation")
)
```

```{.r filename="R"}
cells_before_peak_filtering <-
  targets::tar_read(metadata_w_QC_tibble.ATAC.my_aggregation) |>
  dplyr::count(GEM_well_ID, name = "cells_before_filtering")

cells_after_peak_filtering <-
  targets::tar_read(metadata_filtered_tibble.ATAC.my_aggregation) |>
  dplyr::count(GEM_well_ID, name = "retained_cells")

cells_after_peak_filtering |>
  dplyr::left_join(cells_before_peak_filtering, by = "GEM_well_ID") |>
  dplyr::mutate(retained_fraction = retained_cells / cells_before_filtering)
```

Review the two-GEM-well demo's [peak-QC exclusion-overlap
plot](gallery_main.qmd#pre-lsi-qc-exclusions) alongside the retained fractions
from your run.

Confirm that the overall loss, loss per GEM well, and overlapping exclusion
reasons are reasonable. Revise the aggregation-level filters and rerun this
checkpoint if they are not.

## 6. Review ATAC

Now configure the ATAC analysis settings, including:

- LSI dimensions and neighbours;
- ATAC Harmony variables and clustering resolution;
- marker transcription factors; and
- ATAC scDblFinder removal settings.

As for GEX, an initial run with cell- and cluster-level ATAC doublet removal
disabled lets the score distributions inform the final policy.

```{.yaml filename="cfg_aggregations.yaml"}
aggregation_scDblFinder_ATAC_remove_called_doublets: false
aggregation_scDblFinder_ATAC_max_doublet_fraction_per_cluster: null
```

::: {.callout-note title="Required background knowledge"}
The ATAC checkpoint starts from the accepted peak-QC cell set, performs LSI,
optional Harmony correction, clustering, cell typing, motif-family analysis,
and configured regulatory summaries.
:::

```{.r filename="R"}
targets::tar_manifest(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:6_ATAC-QC]")
  ) & tidyselect::ends_with(".my_aggregation"),
  callr_function = NULL
)[, c("name", "description")]

targets::tar_make(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:6_ATAC-QC]")
  ) & tidyselect::ends_with(".my_aggregation")
)
```

Review LSI diagnostics, metadata associations, the two-GEM-well demo's [ATAC
UMAP](gallery_main.qmd#atac-harmony-umap), cluster stability, and [ATAC
scDblFinder score distributions](gallery_main.qmd#atac-scdblfinder-scores).
Also review motifs, gene activity, and coverage or differential-accessibility
outputs where configured. Continue when the ATAC result is credible
independently of the RNA result.

## 7. Review the multimodal result

Finally, review the WNN neighbour, resolution, and UMAP settings.

::: {.callout-note title="Required background knowledge"}
The multimodal checkpoint combines the accepted GEX and ATAC representations
using WNN, then produces integrated clusters, metadata, review plots, and the
optional Seurat/Signac compatibility object.
:::

```{.r filename="R"}
targets::tar_manifest(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:7_multimodal-QC]")
  ) & tidyselect::ends_with(".my_aggregation"),
  callr_function = NULL
)[, c("name", "description")]

targets::tar_make(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:7_multimodal-QC]")
  ) & tidyselect::ends_with(".my_aggregation")
)
```

Review the two-GEM-well demo's [integrated WNN
UMAP](gallery_main.qmd#wnn-umap) and [RNA/ATAC modality
weights](gallery_main.qmd#wnn-modality-weights), then compare your WNN cell
types and clusters against the accepted single-modality results. The final
compatibility object is `multimodal_Seurat_object.my_aggregation`.

Use [Verify and inspect the outputs](demo_outputs.qmd) for `tar_read()` and
output-path examples. Continue to [differential
analyses](downstream_differential_analyses.qmd) or [genetic
enrichment](downstream_genetic_enrichment.qmd) only after accepting the main
aggregation.

## Request an additional result

Each gallery card names its target. To request only that output and its
dependencies, use its exact name with your aggregation suffix:

```{.r filename="R"}
targets::tar_make(
  names = tidyselect::all_of("categorical.UMAPs.WNN.my_aggregation")
)
```

Use `multimodal_Seurat_object.my_aggregation` for the final compatibility
object. That endpoint does not include every review plot. Preview an exact
selection with `tar_manifest()` and `callr_function = NULL` first;
`all_of()` reports an error when the name is absent.

## Rerun after a change

Reuse the same selection after changing inputs or settings. Replace
`tar_manifest()` with `tar_outdated()` to inspect which selected targets and
their dependencies need rebuilding, keeping `callr_function = NULL`.
Then rerun the checkpoint and review its outputs again.

## Build the complete active scope

An unqualified `targets::tar_make()` constructs every active GEM well, active
aggregation, derived review output, and enabled optional module. Use it only
after the checkpoint-sized runs are accepted and only when that complete scope
is intended:

```{.r filename="R"}
targets::tar_make()
```

Keep unavailable GEM wells, aggregations, and modules inactive before this
broad execution. If a target fails, use [Troubleshooting](troubleshooting.qmd)
and rerun the narrowest affected checkpoint.


## Part: 3. Add an optional analysis


<!-- source: website/downstream_differential_analyses.qmd -->

# Differential analyses



## When to use this module

Use this module to ask how cell-type proportions, gene expression, or chromatin accessibility differ with a condition or donor phenotype. It requires biological replication across donors or samples. Molecular measurements are combined into **pseudobulks**: counts summarized for each cell type within a donor or sample.

The module does not create biological replication. The donor structure, covariates, design formula, and contrasts must be defensible for the intended analysis before the workflow is run.

See the [output gallery](gallery_differential_analyses.qmd) for representative diagnostics and the [implementation graph](implementation/implementation_differential_analyses.html) for target structure.

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
| Do cell-type proportions differ? | Differential cell-type composition (DCTC) |
| Which genes change expression? | Differential gene expression (DGE) |
| Which peaks change accessibility? | Differential chromatin accessibility (DCA) |
| Which motif families change accessibility? | ATAC-derived differential TF activity (DTFA) |
| Which regulators show altered expression-based activity? | CollecTRI-derived differential TF activity (DCTA) |

The module also produces model diagnostics, comparisons across modalities,
and gene-set tests for Hallmark and Reactome pathways. Activity scores are
inferred from accessibility or expression; interpret them in the context of
the measurement used.

See the [method details](implementation/implementation_differential_analyses.html#method-details)
for activity inference, motif-family definitions, and gene-set testing.

## Configure the module

Add `modules` to the existing aggregation entry, keeping its input and marker settings:

```{.yaml filename="cfg_aggregations.yaml"}
your_aggregation:
  modules: [differential_analyses]
```

Then create a matching row directly in
`module_differential_analyses/cfg.yaml`.

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

The module selection below requests both composition and pseudobulk outputs. Configure the DCTC phenotype/formula and at least one pseudobulk model before using that broad selector. Formula terms and contrast coefficients must match columns produced by the model matrix.

The model example assumes `condition` distinguishes treated and control donors.
Check which group is the reference and what each model coefficient represents
before using `conditiontreated` as a contrast. Replace the example formula
and contrast to match your study.

## Run and review

The existing `checkpoint:differential_analyses` tag selects module outputs.
It does not perform QC approval; the [QC checkpoint procedure](main_running.qmd#qc-checkpoints)
is under development. Preview the selected targets first:

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

Review pseudobulk depths and retained donor counts before interpreting coefficients. Check model-matrix terms, P-value distributions, effect directions, and agreement or disagreement across DGE, DCA, DTFA, and DCTA. The CollecTRI-DTFA concordance target summarizes family coverage, rank correlation, directional agreement, and joint FDR support for every configured contrast. These are complementary regulatory readouts: agreement strengthens a shared interpretation, while disagreement can reflect post-transcriptional regulation, motif-family ambiguity, or different evidence carried by expression and accessibility. Treat the [gallery](gallery_differential_analyses.qmd) as a visual reference, not as a statistical acceptance threshold.

Runtime depends on donors, cell types, models, contrasts, and gene-set analyses. Use [Troubleshooting](troubleshooting.qmd) if a formula, contrast, or metadata join fails.

## Parameter reference

The OLINK and bulk-RNA path fields are reserved optional integration inputs and are not consumed by the current public differential-analysis selection. Leave them `NULL` unless the corresponding integration is implemented in your downstream workflow.

[Generated Quarto chunk omitted: `emit_parameter_overview("differential_analyses")`]

<details>
<summary>Show the public <code>immune_human_2x</code> example</summary>

[Generated Quarto chunk omitted: `emit_yaml_entry(module_config_file, "immune_human_2x")`]

</details>


<!-- source: website/downstream_genetic_enrichment.qmd -->

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

See the [output gallery](gallery_genetic_enrichment.qmd) for representative results and the [implementation graph](implementation/implementation_genetic_enrichment.html) for upstream ATAC and WNN dependencies.

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

## Configure the module

Add `modules` to the existing human aggregation entry, keeping its other settings:

```{.yaml filename="cfg_aggregations.yaml"}
your_aggregation:
  modules: [genetic_enrichment]
```

Then create a matching row directly in `module_genetic_enrichment/cfg.yaml`.

```{.yaml filename="module_genetic_enrichment/cfg.yaml"}
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

## Run and review

The existing `checkpoint:genetic_enrichment` tag selects module outputs. It
does not perform QC approval; the [QC checkpoint procedure](main_running.qmd#qc-checkpoints)
is under development. Preview the selected targets:

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

Before interpreting trait scores, verify the resolved source release and fine-mapping method, the number of credible-set loci and variants retained, and the overlap with consensus peaks. Then compare direct deviation summaries with SCAVENGE-propagated scores and use the heatmap glyphs to identify cluster-median enrichment supported by within-grouping BH-adjusted degree-matched permutation P-values.

Runtime and disk use grow with studies, cells, graph representations, permutations, and attributed loci. The [gallery](gallery_genetic_enrichment.qmd) uses a larger aggregation with six GEM wells and is not produced by the minimal quickstart.

## Parameter reference

[Generated Quarto chunk omitted: `emit_parameter_overview("genetic_enrichment")`]

<details>
<summary>Show the public <code>immune_human_2x</code> example</summary>

[Generated Quarto chunk omitted: `emit_yaml_entry(module_config_file, "immune_human_2x")`]

</details>


## Part: Output gallery


<!-- source: website/gallery_main.qmd -->

# Main pipeline gallery



These documentation snapshots show representative outputs from the public
`immune_human_2x` configuration with its two active GEM wells. The cards follow
the numbered reviews in the adaptation guide. Each
card names the target that generated the displayed demo result; click an image
to open it at full resolution.

To reproduce these plot families, [install](demo_installation.qmd) and [run](demo_running.qmd) the demo, then request the broader [review checkpoints](main_running.qmd). The endpoint-only quickstart does not build every gallery plot; [inspect its local outputs](demo_outputs.qmd) to see the distinction.

[Generated Quarto chunk omitted: `render_gallery_section( gallery_items, "Main pipeline", subsection_descriptions = c( "Pre-aggregation QC" = "Unfilter...`]


<!-- source: website/gallery_differential_analyses.qmd -->

# Differential analyses gallery



These cards are a curated subset from the public `immune_human_2x` configuration. They illustrate diagnostics, not acceptable effect sizes or significance patterns for another study. See [Differential analyses](downstream_differential_analyses.qmd) for prerequisites, models, and the module run command.

The full module additionally produces expression-derived CollecTRI activity results and a CollecTRI-DTFA concordance plot. They are not shown below until stable public example assets are available.

[Generated Quarto chunk omitted: `render_gallery_section( gallery_items, "Differential analyses module", subsection_descriptions = c( "Gene expression"...`]


<!-- source: website/gallery_genetic_enrichment.qmd -->

# Genetic enrichment gallery



These curated outputs use the larger `PBMC_human_6x` aggregation, not the quickstart with two GEM wells, and show selected SCAVENGE/WNN results rather than every attribution output. See [Genetic enrichment](downstream_genetic_enrichment.qmd) for prerequisites, module target selection, and interpretation guidance.

[Generated Quarto chunk omitted: `render_gallery_section( gallery_items, "Genetic enrichment module", subsection_descriptions = c( "Single-nucleus chro...`]


## Part: Operation and scaling


<!-- source: website/performance_overview.qmd -->

# Performance and scaling

multiomeR gains performance from two complementary mechanisms: `{targets}` exposes independent GEM wells, modalities, and analysis branches for concurrent execution, while BPCells keeps large matrices on disk and streams many operations instead of materializing dense objects.

Actual wall time depends on input nuclei, retained cells, peak counts, enabled plots and modules, hardware, scheduler latency, controller concurrency, and which targets are already up to date. Treat benchmark results as workload descriptions, not promises for another machine or configuration.

## Recorded examples

The recorded benchmark below contains two public aggregations and estimates the critical path to each final `multimodal_Seurat_object` from recorded `{targets}` runtime metadata.

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

Configure execution capacity in [Distributed computing](performance_distributed_computing.qmd). If a run is unexpectedly slow or repeatedly rebuilds targets, use [Troubleshooting](troubleshooting.qmd).


<!-- source: website/performance_distributed_computing.qmd -->

# Choose where the analysis runs



A **worker** is an R process that runs an analysis task. A **controller** starts
and manages those workers, either on your machine or through a cluster
scheduler. multiomeR uses `crew` for this. Choose the setup below before the
[first demo run](demo_running.qmd).

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

Scheduler startup failures, resource-routing errors, and target failures are handled separately in [Troubleshooting](troubleshooting.qmd). Developer-facing details about runtime bootstrap and `get_tar_resources()` are in [Implementation conventions](implementation/implementation_conventions.html#runtime-bootstrap).

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


<!-- source: website/troubleshooting.qmd -->

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

Check the contracts in [Configuration and inputs](main_inputs.qmd):

- `GEM_well_cellranger_arc_count_dir` contains the required `outs/` files;
- VCF-backed demultiplexing also has `atac_possorted_bam.bam`;
- donor metadata keys and canonical `cfg_GEM_wells.tsv` keys are present and unique;
- non-key metadata columns belong to only one metadata table; and
- configured donor and GEM well IDs match the metadata values exactly.

## Controller and scheduler failures

If a worker does not start or no controller can satisfy a target request:

1. Validate the names, column order, numeric resource values, and controller membership described in [Distributed computing](performance_distributed_computing.qmd).
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


# Book: multiomeR Implementation


<!-- source: website/implementation/index.qmd -->

# Introduction

Use this book to trace a result back to its code or change how multiomeR
works. For installation, configuration, execution, and output inspection,
start with the [user manual](../). You do not need to read this book to run
the demo.

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
| Inspect existing review selections | `[checkpoint:<name>]` description tags and the numbered QC review guide. |
| Add a graph-visible target | Existing `[part_of_graph:<graph_id>]` tags and graph-pruning rules. |
| Change resource routing | `crew_controllers.R`, `packages/multiomeRCore/R/resource_helpers.R`, and the runtime bootstrap convention. |

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
3. `build_aggregation_tibble()` filters active rows and adds symbols for GEM-well, derived QC-group, and aggregation-level upstream targets.
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
The seven numbered main-pipeline groups are listed in
`QC_checkpoint_manifest.tsv`; optional module and supplementary GEX groups
remain unnumbered. Selection matches description substrings; include the
closing `]` to match a complete checkpoint tag. Dependencies still come from
the target commands. The [review guide](../main_running.html#qc-checkpoints)
explains each boundary; acceptance criteria depend on the study.

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

The website renders parameter tables from the same manifest rather than maintaining a second documentation schema. This keeps the user-facing configuration reference tied to the runtime validation contract.

Module-specific YAML uses the same mechanism. Aggregations opt into modules through the aggregation config, and each enabled aggregation must have a matching module config row. Some cross-scope fallbacks are still implemented by target code rather than by manifest inheritance; for example, a module parameter may intentionally allow `NULL` and then fall back to an aggregation-level path during module setup.

## Mapping tibbles

The root `_targets.R` builds the target graph from mapping tibbles. Each mapping tibble is a row-wise contract: one row becomes one set of mapped target instances, and columns in that row become local symbols inside the corresponding `tarchetypes::tar_map()` block.

The core mapping flow is:

1. `GEM_well_tibble_all` reads only the pre-aggregation processing columns from every row in the canonical `cfg_GEM_wells.tsv`.
2. `aggregation_tibble_all_from_yaml` is read from `cfg_aggregations.yaml`.
3. `aggregation_tibble` keeps active aggregations, validates their GEM well references against the complete view, and adds upstream target-symbol columns.
4. `GEM_well_tibble` keeps GEM wells whose `GEM_well_is_active` value is true.
5. `dataset_tibble` derives internal cross-GEM-well QC groups from those active rows.
6. `_targets.R` expands active GEM wells, derived QC summaries, and aggregations with `tar_map()`, then appends module target files.

Within each aggregation, `GEM_well_metadata_tibble` reads the same canonical
file, subsets it to `aggregation_GEM_well_IDs`, and preserves that order. Cheap
keyed projection targets then expose only the columns requested for SCT,
Harmony, subgroup modelling, or configured analyses. Complete non-processing
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

With `GEM_well_ID = "healthy_PBMC_human"`, a target named `cellranger_summary_file` becomes `cellranger_summary_file.healthy_PBMC_human`. The same dot-delimited suffix convention is used for datasets, aggregations, module targets, and nested module maps.

```text
active cfg_GEM_wells.tsv row -> GEM_well_tibble row    -> per GEM well targets
derived dataset group       -> dataset_tibble row     -> per-dataset targets
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

Passing these fixtures does not validate every dataset, parameter regime, approximate-neighbor realization, biological interpretation, or downstream target. The test suite distinguishes fast unit tests for isolated data contracts from slower parity and integration tests that load the project runtime or compare external reference implementations. The [CI workflow](https://github.com/koefoeden/multiomeR/blob/main/.github/workflows/algorithm-validation.yaml) runs the complete suite when tests, relevant helpers, or the Pixi environment change.

Run the complete suite with `pixi run test`. The narrower `pixi run test-algorithm-validation` task runs only the slow UCell, AMULET, WNN, and SCAVENGE parity and acceptance tests.

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

**Implementation and wiring.** The scorer is [`calculate_BPCells_UCell_scores_from_matrix()` in `R/processing_GEX_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/processing_GEX_helpers.R). [`extra_targets/general_aggregation_targets.R`](https://github.com/koefoeden/multiomeR/blob/main/extra_targets/general_aggregation_targets.R) validates `UCell_GEX_marker_genes_list`, and [`extra_targets/GEX_graph_and_cluster_targets.R`](https://github.com/koefoeden/multiomeR/blob/main/extra_targets/GEX_graph_and_cluster_targets.R) computes cell-level scores before GEX cell-type assignment. ATAC and WNN cell-type targets reuse those metadata scores.

**Validation.** [`tests/testthat/test-scoring-parity.R`](https://github.com/koefoeden/multiomeR/blob/main/tests/testthat/test-scoring-parity.R) creates a deterministic 500-gene by 37-cell matrix, writes the project input as BPCells, and compares signed signatures with imputed and skipped missing genes. It requires `identical()` values, dimensions, and dimnames against UCell 2.14.0. It also retains exact Seurat `AddModuleScore` and cell-cycle helper checks plus a metadata-join contract.

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

**Validation.** [`tests/testthat/test-amulet-parity.R`](https://github.com/koefoeden/multiomeR/blob/main/tests/testthat/test-amulet-parity.R) compares against scDblFinder 1.24.0 in disposable `callr` processes. It requires `identical()` results for the bundled fragment-file metrics, prefixed pipeline barcodes, a deterministic 12,000-fragment multi-chromosome loci fixture, and the corresponding full AMULET metrics. It also requires an explicit error for unsupported PCR-duplicate expansion. The disposable reference processes isolate stack-imbalance warnings emitted by the current scDblFinder reference under R 4.5 without weakening the returned-object comparison.

**Status and rerun.** Reference-parity tested against scDblFinder 1.24.0 with exact equality; passing in the current locked environment.

```bash
pixi run test-amulet-parity
```

## Native weighted nearest neighbors

**Reference algorithm.** [`Seurat::FindMultiModalNeighbors()`](https://satijalab.org/seurat/reference/findmultimodalneighbors) constructs cell-specific modality weights from within- and cross-modality neighborhood prediction, collects candidate neighbors across modalities, and selects a weighted multimodal neighbor set.

**Reason for reimplementation.** The pipeline already has aligned RNA PCA/Harmony and ATAC LSI/Harmony matrices and needs reusable neighbor indices, distances, and modality weights without creating a Seurat object. Native graph state also feeds UMAP, Leiden clustering, SCAVENGE, and the optional Seurat/Signac export.

**Behavior preserved.** [`weighted_nearest_neighbors_BPCells()`](https://github.com/koefoeden/multiomeR/blob/main/R/processing_multimodal_helpers.R) preserves row-wise L2 normalization, per-modality nearest neighbors, nearest-distance correction, Seurat's small-SNN far-neighbour kernel bandwidth, within/cross prediction kernels, capped modality affinity ratios, normalized cell-specific modality weights, candidate-set union, weighted neighbor ranking, and Seurat's transformation from weighted affinity to neighbor distance.

**Deliberate deviations and consequences.** BPCells HNSW replaces Seurat's Annoy search, so approximate candidate sets need not be identical. The helper does not expose Seurat's optional smoothing or cross-constant list, and BPCells builds downstream SNN state rather than storing Seurat `Neighbor` and `Graph` objects. Its current `seed` argument is not consulted by the HNSW calls, so it must not be interpreted as controlling neighbor-search randomness. These choices can change weights, selected neighbors, SNN edges, clusters, and UMAP coordinates; correlation and overlap are therefore the validation contract, not exact equality.

**Implementation and wiring.** The implementation and graph consumers are in [`R/processing_multimodal_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/processing_multimodal_helpers.R); the tracked project-owned small-SNN kernel is in [`src/wnn_snn_bandwidth.cpp`](https://github.com/koefoeden/multiomeR/blob/main/src/wnn_snn_bandwidth.cpp). [`extra_targets/WNN_targets.R`](https://github.com/koefoeden/multiomeR/blob/main/extra_targets/WNN_targets.R) aligns modality embeddings, creates `WNN_results_raw`, filters small clusters, optionally recomputes `WNN_results`, and wires that state into WNN UMAP, clustering, metadata, and cell-type targets.

**Validation.** [`tests/testthat/test-wnn-parity.R`](https://github.com/koefoeden/multiomeR/blob/main/tests/testthat/test-wnn-parity.R) compares two deterministic RNA/ATAC fixtures with Seurat 5.5.0. The headline production-like fixture uses 400 cells, 12 dimensions, `k = 30`, and candidate range 200, matching the pipeline's configured neighbor count and native WNN candidate range. The pre-migration small-SNN pilot gave modality-weight Spearman 0.989 and mean neighbour-set overlap 0.996 on this fixture. On the 26,667-cell production object, it gave RNA/ATAC weight Spearman 0.988/0.988 and mean neighbour overlap 0.979. These figures supported the migration, but the maintained validation was deliberately not rerun as part of the production commit. The existing thresholds remain the acceptance contract. The current locked BPCells 0.3.1 build is pinned at `28759cdd5125`.

**Status and rerun.** The small-SNN default is based on the reference-similarity pilot above; post-migration validation remains pending. This does not assert exact equality of selected neighbours, SNN weights, clustering, or UMAP.

```bash
pixi run test-algorithm-validation
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
pixi run test-algorithm-validation
```

Run the complete maintained test suite, including the fast helper contracts, with:

```bash
pixi run test
```


## Part: Background


<!-- source: website/implementation/background_philosophy.qmd -->

<!-- begin include: website/background_philosophy.qmd -->

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
results. The [implementation conventions](implementation_conventions.qmd)
explain where to change parameters, helpers, and target definitions; the
[user manual](../) covers running an existing configuration.

<!-- end include: website/background_philosophy.qmd -->


## Part: Target graph views


<!-- source: website/implementation/implementation_main.qmd -->

<!-- begin include: website/implementation_main.qmd -->

# Main pipeline



The root `_targets.R` creates GEM-well, derived QC-group, and aggregation
mapping rows, then maps target fragments from `extra_targets/`. Use the diagrams
to find the relevant stage, then inspect the corresponding source file for the
complete command and resource declaration.

| Stage | Primary source |
|---|---|
| GEM well preprocessing | `extra_targets/per_GEM_well_targets.R` |
| Cross-GEM-well QC summaries | `extra_targets/per_dataset_targets.R` |
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



`module_differential_analyses/targets.R` filters aggregations that enabled the module, joins their module config, attaches symbols for accepted WNN metadata and pseudobulk inputs, and maps the composition, pseudobulk, GSEA, and cross-modality target fragments. The generic pseudobulk model family is instantiated for DGE, DCA, DTFA, and expression-derived CollecTRI activity (DCTA). DCTA first converts filtered, normalized GEX pseudobulks to signed ULM scores and then reuses the same model and contrast machinery.

DTFA uses the 233 official JASPAR2026 CORE vertebrate familial root motifs as its complete feature universe. The pipeline scans those family-level profiles directly, rather than scanning individual motifs and taking the union of their peak matches.

The cross-modality fragment creates a CollecTRI-to-JASPAR family crosswalk, a detailed regulator-level table containing DCTA, DTFA, and TF-expression results, a family-level comparison table, a contrast-level concordance summary, and its plot. CollecTRI complexes remain intact in DCTA; complex-member mappings are introduced only by the comparison crosswalk.

The graph below is an orientation view. Inspect `setup_and_DCTC_targets.R`, `psbulk_DX_targets.R`, `GSEA_targets.R`, and `cross_modality_targets.R` for the complete model and plotting commands. The user-facing prerequisites and module selector are documented in [Differential analyses](../downstream_differential_analyses.html).

[Mermaid graph omitted; source: `website/figures/human_curated/differential_analyses_v2.mmd`]

## Method details

The DCTA branch infers signed TF or TF-complex activity from normalized GEX pseudobulks with CollecTRI regulons and the `decoupleR` univariate linear model (ULM). Genes are filtered for expression across cell-type pseudobulks, and each retained regulator must have at least five measured targets. Its inferred activities then use the same configured donor-level models and contrasts as DGE, DCA, and DTFA. The published human CollecTRI network is downloaded from the OmniPath rescue archive and accepted only when it matches the pipeline's pinned SHA-256 checksum.

DTFA tests the 233 sequence-similarity families in the official JASPAR2026 CORE vertebrate clustering rather than individual TF motifs. Each family is represented by its published root motif, which is scanned directly against the consensus peaks; individual member motifs are used only as family metadata. The same family-level accessibility matrix supports marker plots and the Seurat compatibility export. The CollecTRI-DTFA comparison maps individual CollecTRI regulators to these JASPAR families and compares model t-statistics, not raw activity scales. AP1 and NFKB remain intact as complex regulons during activity inference; their canonical members are used only to associate the complexes with motif families for comparison. Detailed source-level results retain TF expression as a third reference, while family-level summaries use the median CollecTRI regulator t-statistic and report whether any mapped source is FDR-significant.

Each gene-set collection is tested independently with `cameraPR`, `inter.gene.cor = 0.01`, and a minimum of 10 genes represented in the contrast-specific universe. A significant set is more strongly associated with the contrast than the remaining tested genes, rather than merely showing any collective change. Open Targets evidence annotation is optional.

## Metadata used by the models

For each aggregation, the module retains its own file and full-tibble targets, then projects a canonical analysis view containing only donors in that aggregation and columns required by its configured models and composition plots. Rows are ordered by `donor_id` and non-key columns by name. Changes to unused columns, out-of-aggregation donors, or source row and column order therefore stop at this inexpensive projection boundary.

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
