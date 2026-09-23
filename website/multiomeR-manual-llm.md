# multiomeR Website LLM Export

This file is generated from the Quarto book outlines and resolves Quarto include shortcodes.
Hidden setup chunks, generated helper chunks, Mermaid graph bodies, and verbose image metadata are omitted by default.


# Book: multiomeR Manual


<!-- source: website/index.md -->

[Image omitted; source: `figures/multiomeR-logo.svg`; alt: Image]

# Start here

## What is multiomeR?

multiomeR is a [targets](https://books.ropensci.org/targets/)-based workflow for single-nucleus 10x Genomics Multiome ATAC + Gene Expression data. It starts from `cellranger-arc count` outputs and is meant to be adapted to your own data, compute setup, and biological questions.

::: {.callout-warning title="Beta software"}
multiomeR is in beta and may introduce breaking changes between releases. The [running guide](main_running.md#steps) shows what to inspect at each checkpoint, but acceptance criteria still depend on the tissue, study design, and intended use. Report problems or questions through [GitHub issues](https://github.com/koefoeden/multiomeR/issues).
:::

## Workflow at a glance

The **primary module** processes each GEM well, combines selected GEM wells into an aggregation, clusters and labels cell types in its gene-expression (GEX) and ATAC data, and integrates both modalities with weighted nearest neighbors (WNN). Three optional modules extend a completed aggregation with differential analyses, genetic enrichment for human traits, or peak–gene correlation.

[Image omitted; source: `figures/multiomeR_overview_simplified.drawio.svg`; alt: Image]

## How this manual is organized

- **Try the public demo:** [install multiomeR](demo_installation.md), then run and inspect a small example analysis.
- **Analyze your own data:** [check your inputs](main_overview.md), then configure, run, and review the primary module one checkpoint at a time.
- **Add an optional analysis:** run [differential analyses](downstream_differential_analyses.md), [genetic enrichment](downstream_genetic_enrichment.md), or [peak–gene correlation](downstream_peak_gene_correlation.md) on a completed aggregation.
- **Operation and scaling:** [run locally or on a scheduler](performance_distributed_computing.md), and [troubleshoot](troubleshooting.md) failed or outdated targets.
- **Reference:** browse the [output gallery](gallery.md) for an example of each plot, and look up output files, configuration tables, and methods.

## Your first analysis

Start with the public demo, even if you plan to analyze your own data. It runs one aggregation of two human GEM wells with supplied configuration and shows how to read its cell metadata and multimodal Seurat object.

You need basic R skills, a Linux terminal, and a machine with sufficient [memory and disk space](demo_installation.md#system-requirements). You do not need to know how to write a `targets` pipeline. Run blocks labeled **Bash** in the terminal and blocks labeled **R** in the R session opened during installation, both from the repository folder unless stated otherwise. A block labeled with a file name, such as `cfg_aggregations.yaml`, shows content for that file.

## Terms used in this manual

- A **GEM well** is one 10x Chromium chip channel and its `cellranger-arc count` output directory. It may contain nuclei from several **donors**, the individuals identified by `donor_id`.
- An **aggregation** is a joint analysis of one or more GEM wells.
- A **checkpoint** is one of the eight stages of the primary module, each ending with plots to review before you continue. Its name, such as `8_multimodal_QC`, appears in target names and plot folders.
- A **target** is a named result, such as a metadata table, matrix directory, or plot. You request the targets you want; `targets` builds the earlier results they depend on, in order, and reuses those that are up to date. Most target names end with their GEM well or aggregation, as in `multimodal_Seurat_object.8_multimodal_QC.immune_human_2x`.
- The **store** is the folder where `targets` keeps results and the records needed for reruns. Paths in this manual write it as `<store>`.

For a small worked introduction to `targets`, see the [targets walkthrough](https://books.ropensci.org/targets/walkthrough.html).

Continue to [Install and prepare the demo](demo_installation.md).


## Part: Try the public demo


<!-- source: website/demo_installation.md -->

# Install and prepare the demo



The demo analyzes two public 10x Genomics GEM wells: peripheral blood mononuclear cells from a healthy donor and a lymph node with lymphoma, 17,277 nuclei in total.

## System requirements

- 64-bit x86 Linux with `git` and `curl`.
- At least 60 GB of RAM, enough for one memory-intensive target at a time.
- About 30 GB of free disk space. This covers 3.9 GB of demo inputs, about 6 GB of results and up to 15 GB for the Pixi environment and its package cache, which Pixi keeps in your home directory by default.
- Multiple CPU cores. The run time in the next chapter was measured with 16 logical threads.

The committed `crew_controllers.R` suits a 16-CPU, 256-GB workstation. On a machine near the 60-GB minimum, lower its worker counts before running the demo, starting with the heavy workers; see [Local execution](performance_distributed_computing.md#local-execution).

## Set up the demo

Run these commands in the folder where the clone should be created.

```{.bash filename="Bash"}
# Clone the repository and enter its root folder.
git clone https://github.com/koefoeden/multiomeR.git
cd multiomeR

# Install Pixi. Skip these two lines if `pixi` is already on your PATH.
curl -fsSL https://pixi.sh/install.sh | sh
export PATH="$HOME/.pixi/bin:$PATH"

# Install the locked environment, download the demo inputs into example_data/,
# and install the pinned GitHub versions of BPCells, Signac and betterChromVAR.
pixi run --use-environment-activation-cache --locked --run-post-link-scripts setup-demo

# Start R in that environment.
pixi run --use-environment-activation-cache --locked R
```

If the download is interrupted, run the `setup-demo` command again; it skips files that are already downloaded.

Keep the R session open and continue to [Run the demo](demo_running.md).


<!-- source: website/demo_running.md -->

# Run the demo



The supplied configuration has one active aggregation, `immune_human_2x`, which combines the two demo GEM wells. In the R session opened during installation, request two of its results: the multimodal Seurat/Signac object and UMAPs of the integrated WNN embedding, colored by cluster, cell type and other categorical metadata. WNN (weighted nearest neighbors) combines the GEX and ATAC data of each nucleus.

```{.r filename="R"}
demo_targets <- c(
  "multimodal_Seurat_object.8_multimodal_QC.immune_human_2x",
  "categorical.UMAPs.8_multimodal_QC.immune_human_2x"
)

targets::tar_make(names = tidyselect::all_of(demo_targets))
```

`tar_make()` builds these two targets and every result they depend on, from per-GEM-well quality control to WNN integration. It skips results they do not need, such as most checkpoint plots in the [output gallery](gallery.md).

With 16 threads the run takes about 30 minutes and writes about 6 GB. Keep the R session open until it finishes; progress messages show each target as it is dispatched and completed.

## Confirm success

This command should return `character(0)`, meaning that the requested results and their dependencies are up to date:

```{.r filename="R"}
targets::tar_outdated(
  names = tidyselect::all_of(demo_targets),
  callr_function = NULL
)
```

If target names are returned, run the same `tar_make()` command again; completed targets are reused. If the run stopped with an error, first find and fix the cause with [Troubleshooting](troubleshooting.md).

Continue to [Inspect the demo results](demo_outputs.md).


<!-- source: website/demo_outputs.md -->

# Inspect the demo results



The demo saved its results in the targets store, the `outputs/` folder of the clone (set by `store` in `_targets.yaml`; later pages write `<store>`). The store holds four kinds of output:

- `objects/`: R objects, such as tables and the Seurat/Signac object.
- `files/`: data files, such as matrix folders and TSV tables.
- `plots/`: PNG plot images.
- `plot_objects/`: an editable R copy of each plot.

This page opens one example of each kind in the R session. The [output reference](review_outputs.md#output-folders) describes the folders in detail.

## Objects

Read a stored object by its target name with `targets::tar_read()`, here the cell metadata and the multimodal Seurat/Signac object:

```{.r filename="R"}
cell_metadata <- targets::tar_read(
  metadata_w_cell_types_tibble.WNN.immune_human_2x
)
dim(cell_metadata)
head(cell_metadata)

demo_object <- targets::tar_read(
  multimodal_Seurat_object.8_multimodal_QC.immune_human_2x
)
demo_object
```

The metadata table has one row per retained nucleus, with its QC metrics, GEX, ATAC and WNN clusters, cell-type labels and UMAP coordinates. The Seurat/Signac object is a convenience export for exploring the results with Seurat and Signac. It reads its count matrices and ATAC fragments from files in the store and `example_data/`, so keep those folders in place.

## Files

For a file target, `tar_read()` returns the path of the saved file or folder. Open it with a suitable reader:

```{.r filename="R"}
targets::tar_read(cellranger_barcodes_tsv.healthy_PBMC_human)
#> [1] "outputs/files/healthy_PBMC_human/cellranger_barcodes_tsv.tsv"

GEX_counts <- BPCells::open_matrix_dir(
  targets::tar_read(aggregated_GEX_BPCells_matrix_dir.GEX.immune_human_2x)
)
GEX_counts
```

The first file belongs to one GEM well, the second to the aggregation. Paths follow the target name from right to left: `aggregated_GEX_BPCells_matrix_dir.GEX.immune_human_2x` is saved in `outputs/files/immune_human_2x/GEX/aggregated_GEX_BPCells_matrix_dir/`. `consensus_peak_BPCells_matrix_dir.ATAC.immune_human_2x` holds the ATAC peak counts.

## Plots

Plot targets return image paths in the same way, under `outputs/plots/`. The demo built one categorical UMAP per variable:

```{.r filename="R"}
targets::tar_read(categorical.UMAPs.8_multimodal_QC.immune_human_2x)
```

Open `WNN_harmony_SNN_cluster_cell_type.png` from that list to see the WNN clusters labeled by cell type. It should resemble this snapshot:

[Image omitted; source: `figures/demo_WNN_cell_type_UMAP.png`; alt: WNN UMAP of the two demo GEM wells, colored by cluster and cell type]

Plot subtitles and captions explain how to read each plot.

## Plot objects

Each image has an R copy under `plot_objects/`, at the same relative path with `.rds` instead of `.png`. Reopen a plot without rerunning the analysis:

```{.r filename="R"}
plot_file <- file.path(
  targets::tar_config_get("store"),
  "plot_objects/immune_human_2x/8_multimodal_QC/UMAPs/categorical",
  "WNN_harmony_SNN_cluster_cell_type.rds"
)
p <- readRDS(plot_file)
p
```

Edit a ggplot object with the usual ggplot2 functions and save your copy outside the store, where a rerun cannot overwrite it:

```{.r filename="R"}
p <- p + ggplot2::labs(title = "My integrated cell types")
ggplot2::ggsave("my_cell_types.png", p, width = 10, height = 8)
```

Some plots are composites rather than single ggplot objects and need their own editing methods.

## Possible next steps

- Build the rest of the demo aggregation, including every checkpoint plot, with `targets::tar_make()`. The demo configuration has no other active aggregation and no optional modules.
- Browse the [output gallery](gallery.md) for an example of every plot the pipeline and its optional modules save.
- Start your own analysis with [Plan your analysis](main_overview.md).


## Part: Analyze your own data


<!-- source: website/main_overview.md -->

# Plan your analysis

Before the first run, prepare the inputs below and choose where your settings will live.

## Inputs to prepare

- **Cell Ranger ARC outputs:** one `cellranger-arc count` output directory per GEM well. Its `outs/` folder must contain `summary.csv`, `filtered_feature_bc_matrix.h5`, `atac_fragments.tsv.gz`, `atac_fragments.tsv.gz.tbi` and `per_barcode_metrics.csv`. GEM wells analyzed together must use the same Cell Ranger ARC reference.
- **Donor metadata table:** a TSV with one row per donor and the phenotypes or covariates you want to analyze; see [Donor metadata table](reference_donor_metadata.md).
- **Donor genotypes (optional):** to demultiplex a GEM well that pools several donors, a VCF with their genotypes, as described in the [Vireo genotype-input documentation](https://vireosnp.readthedocs.io/en/stable/manual.html). The GEM well's `outs/` folder must then also contain `atac_possorted_bam.bam`.
- **CellBender output (optional):** to use gene-expression counts corrected for ambient RNA, run [CellBender remove-background](https://cellbender.readthedocs.io/en/latest/usage/) first and give its H5 file in the [GEM well table](reference_GEM_wells.md).
- **Compute resources:** for large datasets, use a compute cluster with a job scheduler; see [Choose where the analysis runs](performance_distributed_computing.md).

## Where the configuration lives {#configuration-directory}

The pipeline reads its settings from one configuration directory:

- `cfg_GEM_wells.tsv`: one row per GEM well; see [GEM well table](reference_GEM_wells.md).
- `cfg_aggregations.yaml`: one entry per aggregation, including the path to its donor metadata table; see [Aggregation configuration](reference_aggregations.md).
- `cfg_module_<module>.yaml`: the settings of one optional module, needed only when that module is enabled.

By default, this is `configuration/`, which holds the public demo and example settings. You can edit it directly, or keep your project's settings in a copy and select that copy in `configuration.local`:

```{.bash filename="Bash"}
cp -r configuration configuration_my_project
echo configuration_my_project > configuration.local
```

`configuration.local` lives in the repository root, is ignored by Git and contains one directory path; a relative path is read from the repository root. While it exists, every configuration file comes from the selected directory, with no fallback to `configuration/`. Delete it to return to the defaults. Do not change the selection while a pipeline is running. Relative paths inside the files are still read from the repository root, and the targets store set in `_targets.yaml` stays the same.

Continue to [Run your own analysis](main_running.md#steps).


<!-- source: website/main_running.md -->

# Run your own analysis

Configure your data once, then work through the eight checkpoints of the primary module in order. At each checkpoint, run its targets, review its plots and revise its settings until you accept the result. Start with the defaults and revise as the plots suggest: each plot's subtitle and caption say what to look for, and the [output gallery](gallery.md) shows an example of every plot.

Run the R commands from the repository root, as in the demo. Replace `my_aggregation` and `my_GEM_well` with your own names. `<store>` is the targets store set in `_targets.yaml`, normally `outputs/` in the repository root. Edit the files in the [selected configuration directory](main_overview.md#configuration-directory); linked parameters belong in your aggregation's entry in `cfg_aggregations.yaml`.

Each command also builds the earlier results it depends on. After changing a setting, rerun the checkpoint: `targets` rebuilds only the results the change affects, and later checkpoints pick up the change when they run. The methods behind the checkpoints are described in [Preprocessing and nucleus QC](implementation/methods_preprocessing_and_QC.html), [GEX, ATAC, batch correction and WNN](implementation/methods_GEX_ATAC_and_WNN.html) and [Cell-type annotation and motif accessibility](implementation/methods_annotation_and_motifs.html).

## Configure your data {#steps}

1. **GEM wells:** add one row per GEM well to `cfg_GEM_wells.tsv`, with `GEM_well_QC_exclude_list` set to `NA` for the first run; see [GEM well table](reference_GEM_wells.md).
2. **Donors:** prepare the donor metadata TSV, with one row for each donor in these GEM wells; see [Donor metadata table](reference_donor_metadata.md).
3. **Aggregation:** in `cfg_aggregations.yaml`, copy the `template_aggregation` entry and rename the copy `my_aggregation`. Set `is_active: true`, list your GEM wells in `aggregation_GEM_well_IDs`, give the path to the donor table in `aggregation_donor_id_metadata_tsv`, and replace the placeholder `aggregation_GEX_marker_genes` with markers for the cell types you expect. Optionally, list the donor or GEM well columns to show in the plots in [`aggregation_categorical_vars`](parameters.html#aggregation_categorical_vars) and [`aggregation_continuous_vars`](parameters.html#aggregation_continuous_vars). See [Aggregation configuration](reference_aggregations.md).
4. **Demo:** set `is_active: false` for the `immune_human_2x` aggregation and `GEM_well_is_active` to `FALSE` for its GEM wells, `healthy_PBMC_human` and `lymphoma_lymph_human`. Change both together, because an active aggregation cannot use an inactive GEM well. A later `targets::tar_make()` without `names` then builds only your data.

## Checkpoint 1: Pre-aggregation QC {#checkpoint-1}

Computes QC metrics and donor assignments for each GEM well, applies the GEM well QC filters and compares the GEM wells of the aggregation.

**Run pipeline:**

```{.r filename="R"}
targets::tar_make(
  names = tidyselect::ends_with("1_pre_aggregation_QC.my_aggregation")
)
```

**Review plots** ([examples](gallery.md#1-pre-aggregation-qc)):

```text
<store>/plots/my_aggregation/1_pre_aggregation_QC/
├── per_aggregation_GEM_well_QC_comparisons/
├── demultiplexing_assignment_bars.png
├── nuclei_per_donor_id_bars.png
├── aggregation_excluded_barcodes_by_type_upset.png
├── aggregation_excluded_cellranger_only_barcodes_by_type_upset.png
└── cell_retention_flow_plot.png
```

For the exclusion overlaps within one GEM well, run its checkpoint 1 targets:

```{.r filename="R"}
targets::tar_make(
  names = tidyselect::ends_with("1_pre_aggregation_QC.my_GEM_well")
)
```

```text
<store>/plots/my_GEM_well/1_pre_aggregation_QC/
├── excluded_barcodes_by_type_upset.png
└── excluded_cellranger_only_barcodes_by_type_upset.png
```

**Revise:** choose cutoffs from the `per_aggregation_GEM_well_QC_comparisons/` distributions and enter them in `GEM_well_QC_exclude_list` for each GEM well; see [Set QC filters after the first run](reference_GEM_wells.md#qc-filters). After the rerun, these plots mark numeric cutoffs, and the UpSet and retention plots show the nuclei each filter removes. Remove a GEM well that fails QC from `aggregation_GEM_well_IDs`.

## Checkpoint 2: GEX dimension reduction {#checkpoint-2}

Normalizes the combined GEX data, selects variable genes and computes PCA, with Harmony correction when covariates are configured.

**Run pipeline:**

```{.r filename="R"}
targets::tar_make(
  names = tidyselect::ends_with("2_GEX_PCA_QC.my_aggregation")
)
```

**Review plots** ([examples](gallery.md#2-gex-pca-qc)):

```text
<store>/plots/my_aggregation/2_GEX_PCA_QC/
├── variable_feature_plot.png
├── VizDimLoadings_plots/
├── PCA_singular_values_elbow_plot.png
├── PCA_embedding_sdev_plot.png
└── PCA_metadata_association_barplots.png
```

**Revise:** choose the PCs used for clustering with [`aggregation_GEX_data_PCs`](parameters.html#aggregation_GEX_data_PCs) and for the UMAP with [`aggregation_UMAP_GEX_PCs`](parameters.html#aggregation_UMAP_GEX_PCs). When components track a batch variable such as `GEM_well_ID`, list it in [`aggregation_harmony_correction_metadata_col_names`](parameters.html#aggregation_harmony_correction_metadata_col_names); Harmony then corrects both modalities, and the embedding and association plots add Harmony panels.

## Checkpoint 3: GEX clusters and cell types {#checkpoint-3}

Clusters the GEX data, labels the clusters from your marker genes and removes GEX doublets.

**Run pipeline:**

```{.r filename="R"}
targets::tar_make(
  names = tidyselect::ends_with("3_GEX_QC.my_aggregation")
)
```

**Review plots** ([examples](gallery.md#3-gex-qc)):

```text
<store>/plots/my_aggregation/3_GEX_QC/
├── UMAPs/
│   ├── categorical/{harmony,non_harmony}/
│   ├── continuous/{harmony,non_harmony}/
│   └── cross/
├── markers_by_cluster_dot_plot.png
├── markers_by_cell_type_dot_plot.png
├── module_scores_by_cluster_dot_plot.png
├── module_scores_by_cell_type_dot_plot.png
├── cluster_UCell_advantage_plots/
├── cluster_marker_volcano_plots/
├── cell_type_marker_volcano_plots.png
├── categorical_by_cluster_bars_plots/
├── categorical_by_cell_type_bars_plots/
├── continuous_by_cluster_violin_plot/
├── continuous_by_cell_type_violin_plot/
└── cell_retention_flow_plot.png
```

**Revise:**

- [`aggregation_GEX_cluster_res`](parameters.html#aggregation_GEX_cluster_res) sets how finely the nuclei are clustered; judge it by marker consistency and cluster sizes. By default, the GEX clusters are also the peak-calling groups of checkpoint 4.
- [`aggregation_GEX_marker_genes`](parameters.html#aggregation_GEX_marker_genes) and [`aggregation_cluster_annotation_min_advantage`](parameters.html#aggregation_cluster_annotation_min_advantage) determine the cell-type labels; see [Cluster annotation](review_outputs.md#cluster-annotation).
- `UMAPs/cross/` compares layouts across PC and neighbor counts, for choosing [`aggregation_UMAP_GEX_PCs`](parameters.html#aggregation_UMAP_GEX_PCs) and [`aggregation_UMAP_nNNs`](parameters.html#aggregation_UMAP_nNNs).

## Checkpoint 4: Peak QC {#checkpoint-4}

Calls peaks within groups of nuclei, by default the GEX clusters, merges them into one consensus peak set and computes peak-based ATAC QC metrics.

**Run pipeline:**

```{.r filename="R"}
targets::tar_make(
  names = tidyselect::ends_with("4_peak_QC.my_aggregation")
)
```

**Review plots** ([examples](gallery.md#4-peak-qc)):

```text
<store>/plots/my_aggregation/4_peak_QC/
├── peaks_QC_violins_plot/
└── peaks_similarity_tiles_plot.png
```

**Revise:** choose ATAC filters from the `peaks_QC_violins_plot/` distributions and enter them in [`aggregation_QC_exclude_list_combined_object`](parameters.html#aggregation_QC_exclude_list_combined_object); see [QC filters after peak calling](reference_aggregations.md#qc-filters-after-peak-calling). Checkpoint 5 applies them.

## Checkpoint 5: ATAC filtering {#checkpoint-5}

Applies the ATAC filters from checkpoint 4 to the GEX-retained nuclei.

**Run pipeline:**

```{.r filename="R"}
targets::tar_make(
  names = tidyselect::ends_with("5_pre_LSI_QC.my_aggregation")
)
```

**Review plots** ([examples](gallery.md#5-pre-lsi-qc)):

```text
<store>/plots/my_aggregation/5_pre_LSI_QC/
├── QC_excluded_upset_plot.png
└── cell_retention_flow_plot.png
```

**Revise:** if the filters remove more nuclei than you can justify, or remove one GEM well disproportionately, revise [`aggregation_QC_exclude_list_combined_object`](parameters.html#aggregation_QC_exclude_list_combined_object).

## Checkpoint 6: ATAC dimension reduction {#checkpoint-6}

Computes LSI embeddings of the filtered peak matrix, with Harmony correction when covariates are configured.

**Run pipeline:**

```{.r filename="R"}
targets::tar_make(
  names = tidyselect::ends_with("6_ATAC_LSI_QC.my_aggregation")
)
```

**Review plots** ([examples](gallery.md#6-atac-lsi-qc)):

```text
<store>/plots/my_aggregation/6_ATAC_LSI_QC/
├── LSI_singular_values_elbow_plot.png
├── LSI_embedding_sdev_plot.png
├── VizDimLoadings_plots/
└── LSI_metadata_association_barplots.png
```

**Revise:** choose the LSI components used for clustering with [`aggregation_ATAC_data_PCs`](parameters.html#aggregation_ATAC_data_PCs) and for the UMAP with [`aggregation_UMAP_ATAC_PCs`](parameters.html#aggregation_UMAP_ATAC_PCs). Both start at component 2 by default; check the association plots for components that track sequencing depth or other technical metrics. Add ATAC-only Harmony covariates with [`aggregation_extra_harmony_covars_ATAC`](parameters.html#aggregation_extra_harmony_covars_ATAC).

## Checkpoint 7: ATAC clusters and motifs {#checkpoint-7}

Clusters the ATAC data, names the clusters from GEX marker scores, removes ATAC doublets, and computes gene activity and motif accessibility.

**Configure:** optionally, list expected transcription factors per cell type in [`aggregation_ATAC_marker_TFs`](parameters.html#aggregation_ATAC_marker_TFs) to organize the motif-accessibility plots.

**Run pipeline:**

```{.r filename="R"}
targets::tar_make(
  names = tidyselect::ends_with("7_ATAC_QC.my_aggregation")
)
```

**Review plots** ([examples](gallery.md#7-atac-qc)):

```text
<store>/plots/my_aggregation/7_ATAC_QC/
├── UMAPs/{categorical,continuous,cross}/
├── marker_gene_activity_dot_plot.png
├── motif_family_accessibility_by_ATAC_cluster_heatmap.png
├── motif_family_accessibility_by_GEX_cluster_heatmap.png
├── motif_family_accessibility_by_GEX_cell_type_heatmap.png
├── motif_family_accessibility_marker_volcano_plots.png
├── coverage_tracks_plots/
├── cluster_UCell_advantage_plots/
├── confusion_matrices_plots.png
├── categorical_by_cluster_bars_plots/
├── categorical_by_cell_type_bars_plots/
├── continuous_by_cluster_violin_plot/
├── continuous_by_cell_type_violin_plot/
└── cell_retention_flow_plot.png
```

**Revise:** choose [`aggregation_ATAC_cluster_res`](parameters.html#aggregation_ATAC_cluster_res) from the motif patterns and the agreement with GEX labels in `confusion_matrices_plots.png`.

## Checkpoint 8: WNN integration {#checkpoint-8}

Integrates the GEX and ATAC embeddings with weighted nearest neighbors (WNN), clusters and labels the integrated graph, and builds the final multimodal Seurat/Signac object.

**Run pipeline:**

```{.r filename="R"}
targets::tar_make(
  names = tidyselect::ends_with("8_multimodal_QC.my_aggregation")
)
```

**Review plots** ([examples](gallery.md#8-multimodal-qc)):

```text
<store>/plots/my_aggregation/8_multimodal_QC/
├── UMAPs/
│   ├── categorical/
│   ├── continuous/
│   ├── cross/
│   ├── cluster_named_dim_tri_plot.png
│   └── cluster_cell_type_dim_tri_plot.png
├── markers_by_cluster_dot_plot.png
├── module_scores_by_cluster_dot_plot.png
├── cluster_UCell_advantage_plots/
├── confusion_matrices_plots/
├── WNN_weight_metadata_associations_plot/
├── categorical_by_cluster_bars_plots/
├── categorical_by_cell_type_bars_plots/
├── continuous_by_cluster_violin_plot/
├── continuous_by_cell_type_violin_plot/
└── cell_retention_flow_plot.png
```

**Revise:** choose [`aggregation_WNN_cluster_res`](parameters.html#aggregation_WNN_cluster_res). The WNN clusters and cell types are the populations used in downstream summaries and comparisons.

The final object is `multimodal_Seurat_object.8_multimodal_QC.my_aggregation`; read it as in [Inspect the demo results](demo_outputs.md).

## Choose your next analysis

To build every remaining output, run `targets::tar_make()` without `names`. Then continue with an optional module:

- [Differential analyses](downstream_differential_analyses.md), if you have many donors and a condition of interest.
- [Genetic enrichment](downstream_genetic_enrichment.md), to relate fine-mapped GWAS variants for a human trait to accessibility in cell types.
- [Peak–gene correlation](downstream_peak_gene_correlation.md), to test associations between peak accessibility and gene expression.


## Part: Add an optional analysis


<!-- source: website/downstream_differential_analyses.md -->

# Differential analyses

## When to use

Use this module to test how cell-type proportions, gene expression or chromatin accessibility differ with a donor condition or phenotype. Donors are the biological replicates: proportions are modelled per donor, and molecular measurements are summed into **pseudobulks**, one per cell type and donor.

The module does not create replication. The donors, covariates, formula and contrasts must support the intended comparison. For this reason the public demo leaves the module off: one healthy PBMC donor and one lymphoma lymph-node donor cannot separate condition, donor and tissue effects.

## Prerequisites

Before enabling the module, confirm that:

- you have reviewed the aggregation through [checkpoint 8](main_running.md#checkpoint-8), including its cell-type annotations;
- the donor metadata has one row per `donor_id` and every variable used in a model;
- model variables describe donors, not individual nuclei; and
- each compared group has enough donors for the design and contrasts.

If the models need variables beyond the aggregation's donor metadata, supply a second table in [`differential_analyses_extended_donor_id_metadata_tsv`](parameters.html#differential_analyses_extended_donor_id_metadata_tsv) with the same unique `donor_id` key.

## Outputs

| Question | Output family |
|------------------------------------|------------------------------------|
| Do cell-type proportions differ? | `cell_type_composition` |
| Which genes change expression? | `gene_expression` |
| Which peaks change accessibility? | `chromatin_accessibility` |
| Which motif families change accessibility? | `motif_family_accessibility` (JASPAR) |
| Which regulators show altered expression-based activity? | `transcription_factor_activity` (CollecTRI) |

Motif-family accessibility summarizes ATAC data; transcription-factor activity is inferred from GEX data. The module also produces pseudobulk-depth and model diagnostics, comparisons across the feature families, and Hallmark and Reactome gene-set tests for gene expression.

## Configure

Add `differential_analyses` to the aggregation's [`modules`](parameters.html#modules) list, keeping its other settings:

```{.yaml filename="cfg_aggregations.yaml"}
my_aggregation:
  modules: [differential_analyses]
```

Then add an entry for `my_aggregation` to `cfg_module_differential_analyses.yaml` in the [selected configuration directory](main_overview.md#configuration-directory):

```{.yaml filename="cfg_module_differential_analyses.yaml"}
my_aggregation:
  differential_analyses_cell_type_composition_models:
    condition_abundance:
      formula: ~ condition
      contrast_specs_vec:
        treated_vs_control: conditiontreated
      plot_phenotype_vars: condition
      color_by: condition
  differential_analyses_pseudobulk_models:
    condition_model:
      formula: ~ 0 + cluster + condition
      contrast_specs_vec:
        treated_vs_control: conditiontreated
```

Each branch holds named models; the model and contrast names label the outputs.

- **Abundance models** ([`differential_analyses_cell_type_composition_models`](parameters.html#differential_analyses_cell_type_composition_models)) fit each cell type's share of a donor's nuclei. Use a one-sided, fixed-effects formula; the response is added for you. `plot_phenotype_vars` and `color_by` choose the plotted variables. Omit the field to skip this branch.
- **Feature models** ([`differential_analyses_pseudobulk_models`](parameters.html#differential_analyses_pseudobulk_models)) test the four feature families. The formula can use `cluster`, the cell type of each pseudobulk, and donor metadata variables. The example estimates one condition effect across cell types; to test within one cell type, add `cell_type_subset` and use `~ condition`.
- **Contrasts** in `contrast_specs_vec` are coefficient names of the model matrix, or linear combinations of them. A text variable's alphabetically first value is the reference, so `conditiontreated` is treated minus control. Replace the example variable, formula and contrast with your own.

Optional fields narrow a model:

- `donor_ids` restricts the donors in either branch;
- `cell_type_subset` restricts the cell types whose pseudobulks a feature model tests;
- `GEM_well_IDs` restricts the GEM wells that define an abundance model's population; and
- `cell_types_to_test` restricts the cell types an abundance model tests, while every retained nucleus still counts towards its donor's total.

Donors missing a model variable are excluded. The branches use different models, so their effect sizes are not directly comparable. The methods describe [abundance models](implementation/methods_differential_analyses.html#cell-type-composition) and the [feature-model routes](implementation/methods_differential_analyses.html#model-routes), including random effects, custom design and contrast functions and paired cell-type designs.

## Run

Preview the targets tagged for this module:

```{.r filename="R"}
targets::tar_manifest(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:differential_analyses")
  ) & tidyselect::ends_with(".my_aggregation"),
  callr_function = NULL
)[, c("name", "description")]
```

Then run the same selection:

```{.r filename="R"}
targets::tar_make(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:differential_analyses")
  ) & tidyselect::ends_with(".my_aggregation")
)
```

## Review

Review the plots ([examples](gallery.md#differential-analyses)); each plot's subtitle and caption say what to look for.

```text
<store>/plots/my_aggregation/differential_analyses/
├── pseudobulk_depth_distribution_plot.png
├── cell_type_composition/model_plots_condition_abundance/
├── gene_expression/
│   ├── volcano_plots/condition_model/
│   ├── gene_set_enrichment/{Hallmark,Reactome}/enrichment_plots/condition_model/
│   ├── significant_elements_plot.png
│   └── p_value_distribution_plot.png
├── chromatin_accessibility/
├── motif_family_accessibility/
├── transcription_factor_activity/
├── significant_elements_modality_distribution_plots/
└── CollecTRI_JASPAR/activity_accessibility_concordance_plots/
```

The other feature families follow the `gene_expression/` layout without gene sets.

Read feature-model results in R, for example for gene expression:

```{.r filename="R"}
targets::tar_read(
  results_tibble.gene_expression.differential_analyses.my_aggregation
)
```

The run also writes cohort tables recording each model's included and excluded donors, and the counts and results of each abundance model, under `<store>/files/my_aggregation/differential_analyses/`.

Use [Troubleshooting](troubleshooting.md) if a formula, contrast or metadata join fails.

## Parameter reference

[Open the searchable parameter browser](parameters.html#workflow=differential_analyses). The [differential analyses methods](implementation/methods_differential_analyses.html) list every fixed and configurable setting, and the [implementation graph](implementation/implementation_differential_analyses.html) shows the target structure.


<!-- source: website/downstream_genetic_enrichment.md -->

# Genetic enrichment



## When to use

Use this module to ask which cell types or nuclei have accessible chromatin overlapping genetic evidence for a human trait. It weights ATAC peaks by the fine-mapping probabilities of the GWAS variants they contain, scores enrichment per cell type and per nucleus, and uses the [`SCAVENGE`](https://github.com/sankaranlab/SCAVENGE) method to propagate nucleus scores through the WNN neighbor graph.

A **credible set** contains candidate causal variants at a GWAS locus, with probabilities from fine-mapping. Enrichment helps prioritize cellular contexts; it does not by itself identify a causal cell type, gene or mechanism.

The public demo leaves this module off. Its configuration already contains a study list for `immune_human_2x`, shown under [Parameter reference](#parameter-reference); add `genetic_enrichment` to that aggregation's `modules` to try it.

## Prerequisites

Before enabling the module, confirm that:

- the aggregation is human and you have reviewed it through [checkpoint 8](main_running.md#checkpoint-8), including its cell-type annotations;
- each configured study represents the intended trait and population;
- the machine can download and keep about 4 GB of Open Targets Parquet data, stored once under `<store>/files/OpenTargets/`; and
- you have decided whether to interpret cell-type summaries, individual nuclei or graph-propagated scores; these answer related but different questions.

## Outputs

| Question | Output |
|---|---|
| Which variants and peaks carry each study's evidence? | Peak-weight summaries and study and variant-to-peak tables |
| Which cell types are enriched? | Cell-type chromVAR deviation heatmaps |
| Which loci and variants drive a cell type's enrichment? | Locus-contribution heatmaps and bars, and variant detail plots |
| Which nuclei are trait-relevant? | SCAVENGE trait-relevance scores on the WNN graph, summarized by cell type and cluster |

## Configure

Add `genetic_enrichment` to the aggregation's [`modules`](parameters.html#modules) list, keeping its other settings:

```{.yaml filename="cfg_aggregations.yaml"}
my_aggregation:
  modules: [genetic_enrichment]
```

Then add an entry for `my_aggregation` to `cfg_module_genetic_enrichment.yaml` in the [selected configuration directory](main_overview.md#configuration-directory):

```{.yaml filename="cfg_module_genetic_enrichment.yaml"}
my_aggregation:
  genetic_enrichment_GWAS_studies:
    lymphocyte_count:
      Category: positive_control
      sourceId: GCST90002388
      finemappingMethod: auto
```

Each entry under [`genetic_enrichment_GWAS_studies`](parameters.html#genetic_enrichment_GWAS_studies) is one study. Its name labels the results, and `Category` groups studies in the plots.

`sourceId` values beginning with `GCST` are GWAS Catalog accessions, read from the Open Targets release `26.03` pinned by the workflow. `finemappingMethod: auto`, the default, selects the first available method in the order `SuSie`, `SuSiE-inf`, `PICS`. Name a method instead when the method is part of the analysis; the run then fails if that method is unavailable for the study.

Any other `sourceId` is the path of a local Parquet file, relative to the repository root. Omit `finemappingMethod` for local files: the method, like the study ID, genome build and provenance, is read from the file.

<details>
<summary>Required columns of a local fine-mapped GWAS file</summary>

The file holds one GWAS with one row per credible-set variant on GRCh38.

- **One value in every row:** `schemaVersion` (`1`), `studyId`, `studyType` (`gwas`), `finemappingMethod`, `confidence`, `credibleSetProbability` (greater than 0, at most 1), `genomeBuild` (`GRCh38`), `sampleSize`, `sourceUrl`, `sourceSha256` (lowercase SHA-256), `sourcePublication` and `sourceRelease`.
- **Per variant:** `studyLocusId`, `credibleSetIndex`, `variantId`, `chromosome` (without `chr`), `position`, `variantRepresentation`, `posteriorProbability` (0 to 1), `logBF`, `pValueMantissa`, `pValueExponent`, `beta`, `standardError`, `r2Overall`, `is95CredibleSet`, `is99CredibleSet`, `locusStart` and `locusEnd`.

Each `studyId`, `studyLocusId` and `variantId` combination must be unique, and each `position` must lie within its locus bounds. A file with `credibleSetProbability` 0.95 or 0.99 may contain only variants in that set.

</details>

## Run

Preview the targets tagged for this module:

```{.r filename="R"}
targets::tar_manifest(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:genetic_enrichment")
  ) & tidyselect::ends_with(".my_aggregation"),
  callr_function = NULL
)[, c("name", "description")]
```

Then run the same selection:

```{.r filename="R"}
targets::tar_make(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:genetic_enrichment")
  ) & tidyselect::ends_with(".my_aggregation")
)
```

Runtime and disk use grow with the number of studies, nuclei, SCAVENGE permutations and attributed loci.

## Review

Review the plots ([examples](gallery.md#genetic-enrichment)); each plot's subtitle and caption say what to look for.

```text
<store>/plots/my_aggregation/genetic_enrichment/
├── GWAS_peak_weights_barplot.png
├── cell_type_pseudobulk/
│   ├── chromVAR_deviation_heatmaps/
│   ├── chromVAR_locus_contribution_per_GWAS_heatmaps/
│   ├── chromVAR_locus_contribution_per_GWAS_faceted_bars_plots/
│   ├── chromVAR_absolute_effect_locus_contribution_per_GWAS_faceted_bars_plots/
│   └── chromVAR_variant_contribution_detail_plots/
└── single_nucleus/SCAVENGE/WNN_harmony_SNN/
    ├── TRS_heatmap/
    ├── TRS_UMAPs/
    ├── TRS_cluster_summary_plot/
    └── sig_prop_bars/
```

Read the tables behind the plots in R, for example the cell-type enrichment scores:

```{.r filename="R"}
targets::tar_read(
  chromVAR_deviation_tibble.cell_type_pseudobulk.genetic_enrichment.my_aggregation
)
```

## Parameter reference

[Open the searchable parameter browser](parameters.html#workflow=genetic_enrichment). The [genetic enrichment methods](implementation/methods_genetic_enrichment.html) list every fixed and configurable setting, and the [implementation graph](implementation/implementation_genetic_enrichment.html) shows the target structure.

<details>
<summary>Show the public <code>immune_human_2x</code> example</summary>

[Generated Quarto chunk omitted: `emit_yaml_entry(module_config_file, "immune_human_2x")`]

</details>


<!-- source: website/downstream_peak_gene_correlation.md -->

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

[Open the searchable parameter browser](parameters.html#workflow=peak_gene_correlation). The [peak–gene correlation methods](implementation/methods_peak_gene_correlation.html) list every fixed and configurable setting, including the support-filter thresholds, and the [implementation graph](implementation/implementation_peak_gene_correlation.html) shows the target structure.


## Part: Operation and scaling


<!-- source: website/performance_distributed_computing.md -->

# Choose where the analysis runs

A **worker** is an R process that builds targets. A **controller** starts and manages workers, either on your machine or through a cluster scheduler. multiomeR defines its controllers with `crew` in `crew_controllers.R` in the repository root; the [targets distributed-computing guide](https://books.ropensci.org/targets/crew.html) explains the general setup.

Independent GEM wells, modalities, and analysis branches run in parallel, so more workers shorten a run until the longest chain of dependent targets limits it. Use local workers on a workstation that meets the [system requirements](demo_installation.md#system-requirements). For large datasets, use a cluster scheduler and ask your support team which account and resource limits to use.

## What to expect

The table shows two recorded runs to `multimodal_Seurat_object.8_multimodal_QC.my_aggregation`, estimated from the recorded runtime of each target:

| Aggregation | GEM wells | Nuclei called by Cell Ranger ARC | Estimated critical path | Sum if targets ran one at a time |
|---|---:|---:|---:|---:|
| `immune_human_2x` | 2 | 17,277 | 17.3 minutes | 31.4 minutes |
| `PBMC_human_6x` | 6 | 51,291 | 23.8 minutes | 52.4 minutes |

`immune_human_2x` is the demo aggregation. `PBMC_human_6x` combines the six PBMC GEM wells in the public `cfg_GEM_wells.tsv`; the public `cfg_aggregations.yaml` does not define it.

The critical path is the longest chain of dependent targets: the shortest possible run time when enough workers are available. Fewer workers and scheduler queue time make real runs longer, as do more nuclei, peaks, plots, or modules. Treat the numbers as examples, not predictions for another machine or configuration, and time one representative aggregation before sizing a large run.

## Local execution

The committed `crew_controllers.R` is sized for a 16-CPU, 256-GB workstation: four `local-light` workers (1 core, 16 GB each) and two `local-heavy` workers (6 cores, 60 GB each). On a machine near the 60-GB minimum, lower the two `workers` values so that only one heavy target runs at a time:

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

The `RAM_GB` values in `controller_resources_tibble` route each target to a controller with enough declared memory; they do not limit memory use. Choose `workers` values so that the targets that can run at once fit in physical memory, and raise them only when you know the memory headroom.

After editing the file, restart R or reload the runtime, which also checks the file:

```{.r filename="R"}
load_project_runtime(force = TRUE)
```

## Scheduler execution

For SLURM, PBS, SGE, or LSF, replace the local controllers with the matching `crew.cluster` controllers. The commented SLURM example in `crew_controllers.R` defines light, heavy, and GPU tiers. For each tier:

1. Request, in the controller's scheduler options, the CPUs, memory, and GPUs that its `controller_resources_tibble` row declares.
2. Set the queue, account, wall time, modules, and worker start-up commands your cluster requires.
3. Give GPUs their own tier: only targets that request GPUs are routed to it.

Test a small target selection before increasing worker counts. For start-up and routing errors, see [Troubleshooting](troubleshooting.md#controller-and-scheduler-failures).

## Rules for `crew_controllers.R` {#controller-rules}

The last expression in the file must be a list with `controller_list` and `controller_resources_tibble`:

- Each `controller_name` in the table names one controller in `controller_list`. Names are unique.
- The table has exactly the columns `controller_name`, `cores`, `RAM_GB`, and `gpus`, in that order, with numeric, non-missing resources.
- Each target runs on the smallest tier that meets its CPU, RAM, and GPU request; targets without a request run on the smallest tier. Tiers are compared by GPUs, then cores, then RAM, so row order does not matter.

[Implementation conventions](implementation/implementation_conventions.html#runtime-bootstrap) describe how the runtime loads this file and how targets request resources.


<!-- source: website/troubleshooting.md -->

# Troubleshooting

Find the target name and first error message in the run output. Fix that cause, then [rerun the same selection](#rerun-safely).

| What happened? | Start here |
|---|---|
| Setup fails before any data processing | [The manifest does not build](#the-manifest-does-not-build) |
| A target reports an error | [A run reports errored targets](#a-run-reports-errored-targets) |
| You need to check a result that a failed target used | [Inspect one upstream value](#inspect-one-upstream-value) |
| An error names an input file, GEM well, donor, or metadata column | [Input and metadata failures](#input-and-metadata-failures) |
| Workers do not start or are killed | [Controller and scheduler failures](#controller-and-scheduler-failures) |
| A completed result needs rebuilding | [A target is unexpectedly outdated](#a-target-is-unexpectedly-outdated) |

Run the commands below in the repository-root R session. For general techniques, see the [targets debugging guide](https://books.ropensci.org/targets/debugging.html).

## The manifest does not build

Reload the runtime and build the manifest:

```{.r filename="R"}
load_project_runtime(force = TRUE)
targets::tar_manifest(callr_function = NULL)
```

Failures at this stage usually mean:

- `configuration.local` names a missing directory, or a `cfg_*` file is missing from the [selected configuration directory](main_overview.md#configuration-directory);
- a YAML parameter is unknown, misspelled, of the wrong type, or missing;
- an aggregation selects a GEM well that is absent from `cfg_GEM_wells.tsv`, or an active aggregation selects an inactive GEM well;
- an aggregation enables an unknown module, or the module's configuration has no entry for it;
- `crew_controllers.R` breaks the [controller rules](performance_distributed_computing.md#controller-rules); or
- an R package fails to load.

Correct the field or file named in the error, then repeat the manifest check.

## A run reports errored targets

List the errors, one target per distinct error message:

```{.r filename="R"}
list_distinct_errored_targets()
```

To add commands and stored tracebacks, filter errored targets by a regular expression that matches a target, aggregation, or module:

```{.r filename="R"}
list_distinct_errored_targets_w_tracebacks(
  target_name_pattern = "my_aggregation"
)
```

Some targets run separately for multiple groups; these runs are called **branches**, and their names end in a hash suffix. Copy the full target or branch name from the listing to inspect its saved workspace:

```{.r filename="R"}
inspect_target_workspace("categorical.UMAPs.8_multimodal_QC.my_aggregation")
```

The inspection reports the target command, stored error and traceback, and short summaries of its dependencies. Empty tables, zero-dimensional matrices, missing model columns, or unexpected labels usually point to an upstream data or configuration problem.

## Inspect one upstream value

Read a dependency by its full name and print only a small part of it:

```{.r filename="R"}
value <- targets::tar_read(metadata_w_cell_types_tibble.WNN.my_aggregation)
dim(value)
head(value)
```

## Input and metadata failures

Metadata files are checked when the targets that read them run. Compare the error with the [GEM well table](reference_GEM_wells.md), [donor metadata table](reference_donor_metadata.md), and [aggregation configuration](reference_aggregations.md) references, and check that:

- each `GEM_well_cellranger_arc_count_dir` contains the required `outs/` files, plus `atac_possorted_bam.bam` for genotype demultiplexing;
- `GEM_well_ID` values in `cfg_GEM_wells.tsv` and `donor_id` values in the donor metadata table are unique;
- every donor ID from the GEM well table or from genotype demultiplexing appears in the donor metadata table with identical spelling; and
- apart from these two keys, no column name appears in both tables.

## Controller and scheduler failures

If workers do not start, are killed, or no controller can run a target:

1. Check `crew_controllers.R` against the [controller rules](performance_distributed_computing.md#controller-rules), then reload it with `load_project_runtime(force = TRUE)`.
2. An error starting with `No controller found` means that no tier offers the requested cores, RAM, or GPUs; add or enlarge a tier.
3. For scheduler controllers, read the scheduler's output and error logs, and check the queue, account, wall time, memory, CPU, module, and file-system settings.
4. If local workers are killed for lack of memory, lower their `workers` values as described in [Local execution](performance_distributed_computing.md#local-execution).

## A target is unexpectedly outdated

Ask `targets` which part of a narrow selection needs rebuilding:

```{.r filename="R"}
targets::tar_outdated(
  names = tidyselect::ends_with("8_multimodal_QC.my_aggregation"),
  callr_function = NULL
)
```

Changes to code, configuration, input files, or any upstream target make downstream targets outdated. Look for the most upstream name in the result before assuming that the final target itself is the cause.

## Rerun safely

After fixing the cause, rerun the same selection; up-to-date results are reused:

```{.r filename="R"}
targets::tar_make(
  names = tidyselect::ends_with("8_multimodal_QC.my_aggregation")
)
```

Run `targets::tar_make()` without `names` only when you want every active aggregation and enabled module built. Do not delete the store to debug: it holds the error records, and deleting it forces unrelated recomputation.


## Part: Reference


<!-- source: website/gallery.md -->

# Output gallery



This page shows one saved plot from each plot target of the `mixed_human_31x` aggregation: 31 public GEM wells from 10x Genomics and ENCODE covering heart, blood, pancreas, liver, colon, lung and cerebellum, with all three optional modules enabled. The public default configuration includes this aggregation as an inactive example; its raw data are public, but the reprocessed `cellranger-arc count` outputs are not supplied.

Each card shows that plot, the target's description and its name without the `.mixed_human_31x` suffix. Select a preview to enlarge it. The [running guide](main_running.md#steps) explains when to review each checkpoint.

[Generated Quarto chunk omitted: `render_output_gallery()`]


<!-- source: website/review_outputs.md -->

# Output reference

This page collects facts about the saved outputs. For a guided first look, see [the demo results](demo_outputs.md); for what to review at each checkpoint, see the [running guide](main_running.md#steps). Plot subtitles and captions explain how to read each plot, and the [output gallery](gallery.md) shows an example of each.

## Output folders and metric selection {#output-folders}

Plots and files are saved under `<store>/plots/` and `<store>/files/`, in folders given by the target name read from right to left. For example, `cross.UMAPs.3_GEX_QC.my_aggregation` saves to `<store>/plots/my_aggregation/3_GEX_QC/UMAPs/cross/`. The first folder is the GEM well or aggregation. Existing files do not show which results are current: after a configuration change, run `targets::tar_outdated()` before reviewing them.

`QC_metric_manifest.tsv` in the repository root selects the plotted QC metrics, their display labels and plotting quantiles. `do_plot = FALSE` hides a metric. Plotting quantiles change the displayed range, not the nuclei retained; filters are set in `cfg_GEM_wells.tsv` and `cfg_aggregations.yaml`.

`UMAPs/cross/` at checkpoints 3, 7 and 8 redraws the cell-type UMAP for several numbers of dimensions and neighbours (neighbours only at checkpoint 8). Clusters and labels stay fixed, so these plots show whether the layout depends on the UMAP settings.

## Seurat objects {#seurat-objects}

Two targets export Seurat objects for work outside the pipeline:

- `GEX_Seurat_object.3_GEX_QC.my_aggregation`: GEX data only, available after checkpoint 3, before peak calling.
- `multimodal_Seurat_object.8_multimodal_QC.my_aggregation`: the final Seurat/Signac object with GEX, ATAC and WNN results.

## Cell retention tables {#cell-retention}

Append `.my_aggregation` to these target names when reading them with `targets::tar_read()`. The `cell_retention_flow_plot.png` of each checkpoint draws its table.

| Target | Checkpoint |
|---|---|
| `cell_retention_tibble.GEX_input` | 1 |
| `cell_retention_tibble.GEX` | 3 |
| `cell_retention_tibble.ATAC_input` | 5 |
| `cell_retention_tibble.ATAC` | 7 |
| `cell_retention_tibble.WNN` | 8 |

Each row is one GEM well and one exclusion action (`discard_action`, in `action_order` within its `stage`), with `input_cells`, `excluded_cells`, `retained_cells` and `retained_fraction`. A nucleus matching several actions counts under the first. Each table includes the stages of the tables before it, and GEM wells with no remaining nuclei keep their rows.

## Cluster numbering {#cluster-numbering}

Leiden clusters are numbered by size, starting with 1 for the largest, before clusters smaller than [`aggregation_cluster_min_barcodes`](parameters.html#aggregation_cluster_min_barcodes) are removed. The remaining clusters keep their numbers, so IDs can have gaps.

## Cluster annotation {#cluster-annotation}

GEX, ATAC and WNN clusters are all labelled from GEX data. Each cluster is scored against every marker set in [`aggregation_GEX_marker_genes`](parameters.html#aggregation_GEX_marker_genes), relative to random control genes matched for expression level and detection rate. A cluster is `Assigned` its best-scoring label when that label leads both the control background and the next-best label by at least [`aggregation_cluster_annotation_min_advantage`](parameters.html#aggregation_cluster_annotation_min_advantage); otherwise it is `Unassigned`. Raising the threshold can only withdraw assignments, and unassigned clusters keep their nuclei. The scores are cached separately, so changing the threshold does not recompute them. The scoring rule, control construction and fixed constants are described in [Cell-type annotation and motif accessibility](implementation/methods_annotation_and_motifs.html).

In the cell metadata, the cluster IDs are in `PCA_harmony_SNN_cluster` (GEX), `LSI_harmony_SNN_cluster` (ATAC) and `WNN_harmony_SNN_cluster` (WNN). Each has companion columns with the suffixes `_cell_type` (the label, or `Unassigned`), `_named` (cluster ID and label) and `_annotation_status` (`Assigned` or `Unassigned`).

The `cluster_UCell_diagnostics` targets write these files to `<store>/files/my_aggregation/3_GEX_QC/cluster_UCell_diagnostics/` and to the matching folders under `7_ATAC_QC/` and `8_multimodal_QC/`:

| File | Contents |
|---|---|
| `clusters.tsv` | One row per cluster: status, label, leading candidate, runner-up, the reason for an `Unassigned` status, and stability and GEM well agreement diagnostics. |
| `marker_evidence.tsv` | Each label's score, control background and adjusted score in each cluster. |
| `control_gene_matching.tsv` | Expression level and detection rate of each marker gene and its controls. |
| `GEM_well_agreement.tsv` | The same decision made separately within each GEM well of a cluster; written only when GEM wells contribute enough nuclei. |
| `settings.rds`, `method.txt` | The settings used and a short description of the rule. |

Stability and GEM well agreement are diagnostics; they never change an assignment. `<store>/files/my_aggregation/3_GEX_QC/marker_set_UCell_summary/marker_sets.tsv` summarizes each marker set across the GEX clusters before doublet filtering: how many clusters it leads or is assigned, and its closest competing set.

## Further methods

See [algorithm validation](implementation/algorithm_validation.html) for reference comparisons and deviations, and the [peak–gene module](downstream_peak_gene_correlation.md) for its model, support filters and output paths.


<!-- source: website/reference_GEM_wells.md -->

# GEM well table



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

[Generated Quarto chunk omitted: `emit_GEM_well_demo_table( GEM_well_config_file = "website/data/demo_GEM_wells.tsv", dictionary_file = "website/data/G...`]


<!-- source: website/reference_donor_metadata.md -->

# Donor metadata table

The donor metadata table is a TSV with one row per `donor_id`. Each aggregation points to one such file through [`aggregation_donor_id_metadata_tsv`](parameters.html#aggregation_donor_id_metadata_tsv) in the [aggregation configuration](reference_aggregations.md). Prepare it before the first run, as listed in [Plan your analysis](main_overview.md); the nuclei-per-donor plot at [checkpoint 1](main_running.md#checkpoint-1) shows the `donor_id` values the nuclei received.

## Minimal table

``` {.text filename="donor_metadata.tsv"}
donor_id	condition
donor_1	control
```

## Matching donors to nuclei

Every nucleus receives a `donor_id` from its GEM well: the configured `GEM_well_donor_id` for a non-multiplexed well, or a genotype-based assignment for a well with a configured VCF; see the [GEM well table](reference_GEM_wells.md#fill-in-a-row). Give each of those IDs exactly one row. Duplicated IDs stop the pipeline; nuclei whose donor is missing from the table get `NA` donor variables.

## Which variables belong here

Put donor-specific phenotypes and covariates in this table, for example condition, age, or sex. Put library-, run-, or batch-specific variables in the GEM well table with a `GEM_well_` prefix. Apart from their key columns, the two tables must not share column names. These rules are checked when the targets that read the table run, not when the pipeline starts.

## Extended table for differential analyses

The [differential analyses](downstream_differential_analyses.md) module can read additional donor-level model variables from a second table given in [`differential_analyses_extended_donor_id_metadata_tsv`](parameters.html#differential_analyses_extended_donor_id_metadata_tsv). It must keep the same unique `donor_id` key. If it is not set, the module uses the aggregation's donor table.


<!-- source: website/reference_aggregations.md -->

# Aggregation configuration



`cfg_aggregations.yaml` in the [selected configuration directory](main_overview.md#configuration-directory) has one top-level entry per aggregation: a joint GEX, ATAC and WNN analysis of one or more GEM wells. This page describes the entry structure; the [running guide](main_running.md#steps) explains when to set each parameter and how to review its effect.

The public configuration contains these entries; inactive ones can stay as examples:

- `template_aggregation` (inactive): a starting point for your own entry.
- `immune_human_2x` (active): the public demo, combining two human GEM wells with optional modules disabled.
- `brain_mouse` and `ENCODE_heart_LV_6x` (inactive): mouse and ENCODE validation examples.
- `mixed_human_31x` (inactive): the aggregation behind the [output gallery](gallery.md), with all optional modules enabled.

## Minimal entry

``` {.yaml filename="cfg_aggregations.yaml"}
my_aggregation:
  aggregation_GEM_well_IDs: [my_GEM_well]
  aggregation_donor_id_metadata_tsv: /path/to/donor_metadata.tsv
  aggregation_GEX_marker_genes:
    Cell_type_A: [GENE1, GENE2]
    Cell_type_B: [GENE3, GENE4]
  is_active: true
```

Every other parameter is optional or has a default, listed in the [parameter reference](#parameter-reference) below. Add a parameter to the entry only to change its default.

## Required keys

- [`aggregation_GEM_well_IDs`](parameters.html#aggregation_GEM_well_IDs): the `GEM_well_ID` values to combine. Each must be an active row of the [GEM well table](reference_GEM_wells.md), and all must use the same Cell Ranger ARC reference.
- [`aggregation_donor_id_metadata_tsv`](parameters.html#aggregation_donor_id_metadata_tsv): the [donor metadata table](reference_donor_metadata.md) for these GEM wells.
- [`aggregation_GEX_marker_genes`](parameters.html#aggregation_GEX_marker_genes): marker genes per expected cell type, used for cluster annotation and marker plots.

[`is_active`](parameters.html#is_active) (default `true`) controls whether targets are constructed for the aggregation. Set it to `false` for aggregations you are not ready to run before an unqualified `targets::tar_make()`; this does not delete existing results.

## Marker genes and transcription factors

Replace the placeholder genes with symbols appropriate for the tissue. List at least two cell types, and use gene names from the Cell Ranger ARC reference. A gene without a suffix or with a `+` suffix is a positive marker; a `-` suffix marks a gene that should be absent. Each cell type needs at least one positive marker. The optional [`aggregation_ATAC_marker_TFs`](parameters.html#aggregation_ATAC_marker_TFs) names transcription factors per cell type for the motif-accessibility plots at checkpoint 7. [Cluster annotation](review_outputs.md#cluster-annotation) describes how the marker lists are used.

## QC filters after peak calling

[`aggregation_QC_exclude_list_combined_object`](parameters.html#aggregation_QC_exclude_list_combined_object) lists filter expressions over the peak-based ATAC metrics. Nuclei for which an expression is `TRUE` are removed, for example:

``` {.yaml filename="cfg_aggregations.yaml"}
aggregation_QC_exclude_list_combined_object:
  - nCount_ATAC < 1000
  - atac_peak_counts_frac < 0.1
  - atac_peak_counts_blacklist_frac > 0.01
```

Omit it until the peak QC plots at [checkpoint 4](main_running.md#checkpoint-4) have shown the distributions; checkpoint 5 then shows which nuclei the filters remove.

## Optional modules

Omit [`modules`](parameters.html#modules) for the first run. After reviewing the primary module's results, list the optional modules to run:

``` {.yaml filename="cfg_aggregations.yaml"}
my_aggregation:
  modules: [differential_analyses]
```

Each listed module also needs an entry named after the aggregation in its own configuration file, in the same directory. The module pages describe these entries:

| Module | Configuration file | Page |
|---|---|---|
| `differential_analyses` | `cfg_module_differential_analyses.yaml` | [Differential analyses](downstream_differential_analyses.md) |
| `genetic_enrichment` | `cfg_module_genetic_enrichment.yaml` | [Genetic enrichment](downstream_genetic_enrichment.md) |
| `peak_gene_correlation` | `cfg_module_peak_gene_correlation.yaml` | [Peak–gene correlation](downstream_peak_gene_correlation.md) |

## Parameter reference {#parameter-reference}

The [parameter browser](parameters.html) lists every parameter of the primary module and the optional modules with its default, type and an example. Choose a workflow, then search by name or purpose. The [committed example](https://github.com/koefoeden/multiomeR/blob/main/configuration/cfg_aggregations.yaml) shows complete entries.

<details>

<summary>Show the public <code>immune_human_2x</code> example</summary>

[Generated Quarto chunk omitted: `emit_yaml_entry(aggregations_config_file, "immune_human_2x")`]

</details>


# Book: multiomeR Implementation


<!-- source: website/implementation/index.md -->

# Introduction

Use this book to trace a result back to its code or change how multiomeR works. For installation, configuration, execution, and output inspection, start with the [user manual](../). You do not need to read this book to run the demo.

## Design

multiomeR keeps the analysis steps in an editable repository. Configuration covers common choices such as inputs, markers, dimensions, and models; R helpers and target definitions are available when a study needs a change beyond those settings. This flexibility also means that users must review which methods and assumptions fit their data.

- **Reuse completed work.** [`targets`](https://books.ropensci.org/targets/) records dependencies between results, so a change rebuilds only the affected parts of an analysis, and independent tasks such as GEM wells run concurrently when workers are available.
- **Keep large matrices on disk.** [BPCells](https://bnprks.github.io/BPCells/) provides disk-backed matrices and streaming operations; some steps still need substantial RAM. [Choose where the analysis runs](../performance_distributed_computing.html#what-to-expect) gives multiomeR examples.
- **Keep the analysis inspectable.** Separate targets make intermediate tables, matrices, and files available for inspection, and Seurat/Signac exports allow exploration outside the pipeline.

## Where to start

1. Read [Reading the graph views](graph_methodology.md) and follow its configuration-to-target trace.
2. Open the [primary module](implementation_main.md) graph for the modality or checkpoint you plan to change. The [differential analyses](implementation_differential_analyses.md), [genetic enrichment](implementation_genetic_enrichment.md) and [peak–gene correlation](implementation_peak_gene_correlation.md) chapters cover the optional modules.
3. Use [Implementation conventions](implementation_conventions.md) for the manifest, mapping, symbol, tag, and runtime contracts.

The **Methods and parameters** chapters, from [Preprocessing and nucleus QC](methods_preprocessing_and_QC.md) to [Genetic enrichment](methods_genetic_enrichment.md), describe each stage and tabulate the settings that determine its results. Read them when you need the behaviour behind a result, or when deciding whether a change is a configuration edit or a code edit. [Algorithmic implementations, deviations and validation](algorithm_validation.md) records how the reimplemented reference algorithms differ from their references and how they are tested.

## Common entry points

| Change | Start with |
|---|---|
| Add or revise a YAML parameter | `cfg_pipeline_parameters.tsv`, then the owning config reader or target. |
| Find whether a threshold is configurable or fixed | The stage's chapter under **Methods and parameters**. |
| Change GEM well preprocessing | `_targets.R` mapping plus `extra_targets/per_GEM_well_targets.R`. |
| Change aggregation GEX, ATAC, or WNN processing | The corresponding graph section and `extra_targets/*_targets.R` file. |
| Inspect existing review selections | `[checkpoint:<name>]` description tags and the steps in [Run your own analysis](../main_running.html#steps). |
| Add a graph-visible target | Existing `[part_of_graph:<graph_id>]` tags and graph-pruning rules. |
| Change resource routing | `crew_controllers.R`, `packages/multiomeRCore/R/resource_helpers.R`, and the runtime bootstrap convention. |


## Part: Orientation


<!-- source: website/implementation/graph_methodology.md -->

# Reading the graph views



The graph chapters collect simplified views of the real `{targets}` dependency graph. They are meant to make the workflow easier to reason about before reading the target code directly.

The diagrams are generated from tagged target metadata and the real dependency graph, then simplified by pruning or bypassing lower-level nodes that would make each view harder to read. They keep real target names and preserve the dependency structure where practical, while staying compact enough to build intuition about the main control points.

The following chapters cover the primary module and the differential analyses, genetic enrichment and peak–gene correlation modules.

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

This chapter explains how configuration becomes target definitions: the parameter manifest supplies defaults and validation rules, mapping tables define repeated analyses, and target symbols connect their dependencies. Description tags support result selection and graph views. The last two sections cover project startup and the layout of the methods chapters. Preserve these contracts unless a change is meant to replace one of them.

## Target metadata tags

multiomeR stores lightweight target metadata in the `description` argument of `targets::tar_target()` and `tarchetypes::tar_file()` calls. The descriptions should remain readable prose, with bracketed tags appended when a target needs to be discoverable from the manifest.

``` r
targets::tar_target(
  name = harmony_embeddings_matrix.GEX,
  description = "Harmony-corrected SCTransform GEX PCA embeddings [part_of_graph:GEX] [part_of_graph:WNN]",
  command = ...
)
```

The tag families are:

``` text
[checkpoint:<name>]             review or execution checkpoint
[part_of_graph:<graph_id>]      curated membership in an implementation graph
[resource_observation:<note>]   compact empirical resource note
```

`[checkpoint:<name>]` marks targets selectable with `targets::tar_described_as()`. The numbered primary-module groups are listed in `QC_checkpoint_manifest.tsv`; optional module groups remain unnumbered. Selection matches description substrings; include the closing `]` to match a complete checkpoint tag. Dependencies still come from the target commands. [Run your own analysis](../main_running.html#steps) explains each boundary.

Numbered checkpoint plot targets end in the checkpoint name with hyphens replaced by underscores, before the mapped dataset or aggregation suffix. For example, `VizDimLoadings_plots.2_GEX_PCA_QC.my_aggregation` writes beneath `<store>/plots/my_aggregation/2_GEX_PCA_QC/`. Only plot targets use this naming convention; computational and metadata targets retain their modality suffixes.

`[part_of_graph:<graph_id>]` marks targets that should stay visible in a named implementation graph after graph-pruning helpers remove less informative intermediate nodes. This is the strictest tag family: `graph_id` must contain only letters, numbers, and underscores, and helper code parses these tags directly from target descriptions. A target may belong to several graph views.

``` r
description = paste(
  "Build the lightweight BPCells-backed chromVAR RSE",
  "[part_of_graph:ATAC]",
  "[part_of_graph:seurat_export]",
  "[part_of_graph:genetic_enrichment_single_nucleus]"
)
```

`[resource_observation:<note>]` keeps a short, dated runtime or memory observation next to the target that produced it. It is documentation, not a resource-estimation system.

Use tags only when they create a durable handle for readers, graph helpers, or checkpoint commands. Ordinary internal dependencies can stay untagged.

## Parameter manifest

`cfg_pipeline_parameters.tsv` is the schema for YAML-backed pipeline configuration. Each row defines one parameter for one scope: `aggregation`, or the name of an optional module. Its columns record the type, cardinality, default, missing-value rule, allowed values and description of the parameter. The YAML files then only need to specify values that differ from the manifest defaults, plus values that are required because their resolved value may not be missing.

Configuration readers resolve file names with `configuration_path()`, within `configuration/` or the directory named by the ignored `configuration.local`. The selection does not change data-path interpretation or the targets store. Aggregation and enabled-module settings are resolved during graph construction; disabled modules do not read their configuration.

At read time, the pipeline loads the manifest for a scope and parses each `default_value` as YAML. This allows defaults to be literal scalars, `NULL`, YAML lists, or evaluated YAML expressions such as `!expr 1:30`. Manifest defaults seed every config row before inheritance and row-specific overrides are applied.

The resolution order is:

1.  Start with manifest defaults for the requested scope.
2.  Resolve each parent listed in `inherits`.
3.  Overlay parent values onto the defaults.
4.  Overlay the child row onto the inherited values.
5.  Validate the fully resolved row.

``` yaml
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

``` text
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

1.  `GEM_well_tibble_all` reads only the pre-aggregation processing columns from every row in the canonical `cfg_GEM_wells.tsv`.
2.  `aggregation_tibble_all_from_yaml` is read from `cfg_aggregations.yaml`.
3.  `aggregation_tibble` keeps active aggregations, validates their GEM well references against the complete view, and adds upstream target-symbol columns.
4.  `GEM_well_tibble` keeps GEM wells whose `GEM_well_is_active` value is true.
5.  `_targets.R` expands active GEM wells and aggregations with `tar_map()`, then appends module target files. Cross-GEM-well QC summaries use the aggregation's selected wells.

Aggregation targets read GEM-well annotations through keyed projection targets that expose only the columns a consumer needs. These projections are cache boundaries: editing an unrelated column of `cfg_GEM_wells.tsv` does not invalidate expensive consumers.

``` r
tarchetypes::tar_map(
  values = GEM_well_tibble,
  names = GEM_well_ID,
  delimiter = ".",
  source("extra_targets/per_GEM_well_targets.R")$value
)
```

With `GEM_well_ID = "healthy_PBMC_human"`, a target named `cellranger_summary_file` becomes `cellranger_summary_file.healthy_PBMC_human`. The same dot-delimited suffix convention is used for aggregations, module targets, and nested module maps.

``` text
active cfg_GEM_wells.tsv row -> GEM_well_tibble row    -> per GEM well targets
cfg_aggregations.yaml key   -> aggregation_tibble row -> per-aggregation targets
```

Aggregation rows may opt into optional modules through `modules`. `_targets.R` validates module names against the known module list, and module target files then filter `aggregation_tibble` to the active aggregations that requested that module. Each opted-in aggregation must have a matching module config row.

The naming convention is therefore compositional:

``` text
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

``` r
aggregation_tibble |>
  add_target_sym_cols(
    aggregation_GEX_counts_BPCells_matrix_syms =
      target_sym_col("GEX_counts_BPCells_matrix", "aggregation_GEM_well_IDs")
  )
```

For an aggregation whose `aggregation_GEM_well_IDs` are `c("rx1", "rx2")`, this creates a row value equivalent to:

``` r
rlang::syms(c(
  "GEX_counts_BPCells_matrix.rx1",
  "GEX_counts_BPCells_matrix.rx2"
))
```

The aggregation target can then consume the row-local symbol list directly:

``` r
combined_counts_matrix <- purrr::reduce(
  aggregation_GEX_counts_BPCells_matrix_syms,
  cbind
)
```

Column names should describe the downstream scope, the upstream target, and the fact that the value is a symbol list. The `aggregation_*_syms` columns, such as `aggregation_GEX_counts_BPCells_matrix_syms`, splice per GEM well targets into aggregation-level targets.

Module target files also need aggregation-specific references to primary-module targets. For this, `add_aggregation_target_syms()` creates one symbol per row, suffixed by the aggregation name. These columns are named like the target they replace rather than with `*_syms`, because each cell is a single symbol rather than a list.

``` r
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

1.  loading core workflow packages and conflict preferences,
2.  sourcing generally reusable helpers from `packages/multiomeRCore/R`,
3.  sourcing pipeline-specific helpers from the root `R/` directory,
4.  applying global plotting and `{targets}` options,
5.  sourcing `crew_controllers.R` and installing controller resources.

The nested `multiomeRCore` directory is both ordinary editable pipeline source and an installable package boundary for standalone repositories. multiomeR does not install or attach that package itself: `targets::tar_source()` loads the same implementation files before the root helpers. Keep domain-specific code under `R/`, but do not duplicate the generally reusable implementations there.

For commands that intentionally bypass startup side effects, source the bootstrap helper directly and then load the runtime:

``` r
source("R/bootstrap_helpers.R")
load_project_runtime(force = TRUE)
targets::tar_manifest(callr_function = NULL)
```

Bootstrap state is cached in `bootstrap_state_env`. This avoids reloading packages, re-sourcing helpers, reapplying target options, reassigning patches, and reloading controllers on every call. Use `force = TRUE` when the current R session may be stale, such as after changing helper files, switching checkout roots, editing `crew_controllers.R`, or reusing a long-lived interactive session.

Project-root detection walks upward from the current working directory until it finds `pixi.toml`. Bootstrap commands should therefore be run from inside the multiomeR checkout.

Controller loading is part of the runtime contract, not a later execution detail. `crew_controllers.R` must return a named list with `controller_resources_tibble` and `controller_list`. The bootstrap validates that shape, installs a grouped `crew` controller into `{targets}`, and stores the resource table for `get_tar_resources()`.

``` r
targets::tar_target(
  example_target,
  example_function(),
  resources = get_tar_resources(cores_req = 6, RAM_GB_req = 60)
)
```

If `get_tar_resources()` is called before controller resources are loaded, it fails deliberately with an instruction to call `load_project_runtime()` first. Scheduler-specific examples belong in the main manual's [Choose where the analysis runs](../performance_distributed_computing.html) page; the implementation contract is that target code can request resources declaratively once the runtime has been loaded.

## Methods and parameter tables {#methods-and-parameter-tables}

The chapters in the **Methods and parameters** part pair a description of each analysis stage with a table of the settings that determine its results.

The descriptions state what each step does and why, without numerical values or links. They are kept in heading-less fragments under `_shared_methods/` and included both by these chapters and by the Supplementary Methods of the multiomeR manuscript, so the two texts cannot diverge while the manuscript is prepared. The submitted supplement cites an archived software release, which freezes the matching version of this book; afterwards the book continues to follow the code.

The tables record the values. Each uses one layout:

| Column | Content |
|---|---|
| Step | The analysis step, in the order the targets run. |
| Setting | The quantity or method choice. |
| Status | `Configurable` or `Fixed`. |
| Value | Shown only for fixed settings. |
| Source | Where the setting is controlled. |

A setting is **configurable** when a row of `cfg_pipeline_parameters.tsv` controls it, directly or as a field of a nested parameter such as a model specification, or when a column of `cfg_GEM_wells.tsv` controls it. Configurable rows link to the [parameter browser](../parameters.html), which renders the current default from the public manifest snapshot, and never repeat the default.

A setting is **fixed** when changing it requires a code edit. The Source cell names the symbol that fixes the value: the project function `f()` that contains it, the target whose command passes it, a vendored file, or `pkg::f()` when the value is a default of that package function. Package versions are locked by `pixi.lock`.

The tables list settings a methods section would report or a user might want to change. Parallelism, plotting style and input validation are left to the code. `pixi run --use-environment-activation-cache -e dev check-methods-parameters` checks that every manifest link resolves, every manifest parameter is linked and every cited symbol exists; it does not compare values with the code.


<!-- source: website/implementation/algorithm_validation.md -->

# Algorithmic implementations, deviations and validation

multiomeR reimplements a small number of reference algorithms so they can operate on the workflow's native matrices and graph state. The algorithms themselves are described in the **Methods and parameters** chapters. This page records why each was reimplemented, where it deliberately differs from its reference, and what the executable validation establishes.

The evidence labels are intentionally narrow:

- **Reference-parity tested** means the repository and named reference implementation run on the same deterministic fixture and their returned values are compared directly.
- **Reference-similarity tested** means exact equality is not an appropriate contract, so predefined similarity thresholds are checked against the named reference implementation.
- **Algorithmically derived** means the implementation is checked against an independent mathematical result, not against another software implementation.

Passing these fixtures does not validate every dataset, parameter regime, approximate-neighbour realization, biological interpretation, or downstream target. The test suite contains only such reference comparisons. Each test asserts the reference version it was written against; versions are locked by `pixi.lock`. The [CI workflow](https://github.com/koefoeden/multiomeR/blob/main/.github/workflows/algorithm-validation.yaml) runs the complete suite when tests, relevant helpers, or the Pixi environment change.

Run the complete suite with `pixi run --use-environment-activation-cache test`. The narrower `pixi run --use-environment-activation-cache test-algorithm-validation` task runs only the slow UCell, AMULET, WNN, and SCAVENGE tests.

| Implementation | Evidence status | Reference | Validation contract |
|---|---|---|---|
| BPCells-native UCell | Reference-parity tested | UCell | Identical values, dimensions, and dimnames |
| BPCells-native AMULET | Reference-parity tested | scDblFinder | Identical metrics and loci, including order |
| Native WNN | Reference-similarity tested | Seurat | Modality-weight Spearman and neighbour-overlap thresholds per fixture |
| Sparse SCAVENGE propagation | Algorithmically derived and reference-parity tested | SCAVENGE source at `8ee8b173d965` | Closed-form propagation; identical seed samples, exceedance counts and significant-cell calls |
| Peak–gene donor-slope REML and Kenward–Roger kernels | Reference-parity tested | lme4 and pbkrtest | Coefficients, df and P-values within 1e-6; identical fit statuses |
| Peak–gene HC3 statistics and compact BH breakpoints | Reference-parity tested | sandwich; `stats::p.adjust()` | HC3 coefficients, errors and P-values within 1e-10; identical FDR per chromosome slice |
| Sampled voom correlation | Reference-parity tested | edgeR and limma | Unsampled fits equal `voomLmFit()` within 1e-12; the sampled consensus correlation equals `duplicateCorrelation()` on the same features |

The peak–gene rows belong to the analyses in [Peak–gene correlation](methods_peak_gene_correlation.md); their tests are `test-peak-gene-hierarchical-parity.R` and `test-peak-gene-correlation-parity.R`. `fit_pseudobulk_voom()` can estimate the voom block correlation from a reproducible feature sample for downstream sensitivity analyses; the pipeline's own models use all features. Its test is `test-pseudobulk-correlation-sampling.R`.

## BPCells-native UCell scoring

**Reference algorithm.** [`UCell::ScoreSignatures_UCell()`](https://bioconductor.org/packages/release/bioc/html/UCell.html) calculates per-cell signature scores from descending feature ranks, caps ranks at `maxRank`, combines positive and negative signatures, and clips negative combined scores to zero.

**Reason for reimplementation.** The workflow keeps gene-by-cell counts in BPCells-backed matrices. Materializing the complete matrix or building a Seurat object solely for marker scoring would discard that storage contract, so multiomeR ranks bounded cell chunks and returns metadata-ready scores directly.

**Deliberate deviations and consequences.** Only one cell chunk is materialized at a time, and optional fork workers operate across chunks; this changes memory and execution behaviour but not the tested values. The helper returns a data frame instead of mutating a Seurat object. The target-level marker validator rejects configured genes missing from the reference, whereas the lower-level helper still exposes UCell's impute and skip modes. The production annotation reuses the chunked ranking helper and scores signed signatures per cell before aggregation, which costs more computation than a positive-only rank-summary shortcut.

**Implementation.** `calculate_BPCells_UCell_scores_from_matrix()` and `rank_UCell_count_chunk()` in `R/processing_GEX_helpers.R`; the cluster annotation built on them is in `R/cluster_annotation_helpers.R` and described in [Cell-type annotation and motif accessibility](methods_annotation_and_motifs.md#matched-control-cluster-annotation).

**Validation.** `tests/testthat/test-scoring-parity.R` compares signed signatures, with imputed and skipped missing genes, on a deterministic BPCells fixture and requires `identical()` values, dimensions, and dimnames. It also runs the production annotation path on unsigned, signed and negative-only signatures and compares per-cell scores, cluster means, matched-control summaries and marker-deletion effects with reference scores averaged within clusters, within 1e-12, and checks that chunking and fork workers leave them unchanged. The production cell-cycle scorer is compared with `Seurat::CellCycleScoring()`: phases are identical and scores agree within 1e-6, because BPCells normalizes counts at lower floating-point precision.

```bash
pixi run --use-environment-activation-cache test-scoring-parity
```

## BPCells-native AMULET

**Reference algorithm.** [`scDblFinder::amulet()`](https://bioconductor.org/packages/release/bioc/html/scDblFinder.html) detects likely scATAC-seq doublets from the number of genomic loci covered by more than two fragments, after filtering fragment sizes and excluded regions and removing loci recurrently covered across many cells, and derives Poisson P-values with Benjamini–Hochberg correction.

**Reason for reimplementation.** The per-GEM-well workflow already stores Cell Ranger ATAC fragments as compressed BPCells directories. Passing the fragment file to scDblFinder materializes chromosome-scale `GRanges` objects; the local implementation streams the BPCells fragments and retains only one chromosome's selected fragments while calculating coverage runs.

**Deliberate deviations and consequences.** Only unique-fragment operation is supported: BPCells fragments do not retain Cell Ranger's PCR-duplicate counts, so requesting non-unique expansion fails explicitly. Cell Ranger's inclusive end coordinate is shifted back by one base to match scDblFinder's BED import. BPCells does not export its fragment iterator header, so the native helper mirrors that private interface, verifies the pinned BPCells commit before use, and compiles a small shared library in each worker. A BPCells upgrade must revalidate this interface and the parity fixture before updating the pin.

**Implementation.** `src/amulet_bpcells.cpp` iterates fragments and calculates coverage runs; `calculate_amulet_metrics_BPCells()` in `R/amulet_BPCells_helpers.R` holds the interface check, high-overlap filtering and statistics; the `amulet_metrics_tibble` target consumes the prefixed BPCells fragments.

**Validation.** `tests/testthat/test-amulet-parity.R` requires `identical()` results against scDblFinder for the bundled fragment-file metrics, the production call with prefixed Cell Ranger barcodes, a deterministic multi-chromosome loci fixture, and the corresponding full AMULET metrics, and requires an explicit error for PCR-duplicate expansion.

```bash
pixi run --use-environment-activation-cache test-amulet-parity
```

## BPCells-backed ATAC scDblFinder feature aggregation

**Reference algorithm.** With `aggregateFeatures = TRUE`, [`scDblFinder::scDblFinder()`](https://bioconductor.org/packages/release/bioc/html/scDblFinder.html) performs its own TF-IDF-based feature clustering and sums peaks into feature groups before artificial-doublet classification.

**Reason for adaptation.** Materializing and transforming every peak within each GEM well can exhaust worker memory before classification. The pipeline instead derives feature groups once from the aggregation's LSI loadings, sums the disk-backed peak matrix with BPCells, and passes the compact matrix to scDblFinder with `aggregateFeatures = FALSE`.

**Deliberate deviations and consequences.** This is not a reimplementation of the classifier, which runs unchanged. Global LSI-derived groups replace scDblFinder's per-GEM-well feature groups, so the aggregated matrix, and therefore scores and calls, can differ; exact parity is not expected. The GEX path calls `scDblFinder::scDblFinder()` directly through the same memory-bounding per-GEM-well wrapper.

**Implementation.** `get_feature_groups_from_LSI_loadings()` and `aggregate_BPCells_rows_by_group()` in `R/processing_GEX_helpers.R`, called from `extra_targets/ATAC_targets.R`.

**Validation scope.** The reference-parity fixtures do not establish equivalence for this path. Its fixed settings are listed in [Preprocessing and nucleus QC](methods_preprocessing_and_QC.md#atac-qc-and-doublets).

## Native weighted nearest neighbors

**Reference algorithm.** [`Seurat::FindMultiModalNeighbors()`](https://satijalab.org/seurat/reference/findmultimodalneighbors) constructs cell-specific modality weights from within- and cross-modality neighbourhood prediction, collects candidate neighbours across modalities, and selects a weighted multimodal neighbour set.

**Reason for reimplementation.** The pipeline already holds aligned RNA and ATAC embeddings and needs reusable neighbour indices, distances, and modality weights without creating a Seurat object. Native graph state also feeds UMAP, Leiden clustering, SCAVENGE, and the Seurat/Signac export.

**Deliberate deviations and consequences.** BPCells HNSW replaces Seurat's Annoy search, so candidate sets need not be identical. The helper does not expose Seurat's optional smoothing or cross-constant list, and BPCells builds the downstream SNN graph rather than Seurat `Neighbor` and `Graph` objects. The helper's `seed` argument is not consulted by the HNSW calls and does not control neighbour-search randomness. These choices can change weights, selected neighbours, SNN edges, clusters, and UMAP coordinates, so correlation and overlap, not exact equality, are the validation contract.

**Implementation.** `weighted_nearest_neighbors_BPCells()` in `R/processing_multimodal_helpers.R` with the small-SNN bandwidth kernel in `src/wnn_snn_bandwidth.cpp`; `extra_targets/WNN_targets.R` aligns the embeddings and wires the result into clustering, UMAP and metadata targets.

**Validation.** `tests/testthat/test-wnn-parity.R` compares deterministic RNA/ATAC fixtures with Seurat at production search settings: a 400-cell fixture with candidate range 200 and `k` of 20 and 50, and a 160-cell stress fixture with `k = 15` and candidate range 50. Each case must exceed its thresholds for modality-weight Spearman correlation and mean and lower-quartile neighbour-set overlap, and one- and two-thread results must be identical. This does not assert equality of selected neighbours, SNN weights, clustering, or UMAP.

```bash
pixi run --use-environment-activation-cache test-algorithm-validation
```

## Sparse SCAVENGE propagation and significance

**Reference algorithm.** [SCAVENGE at commit `8ee8b173d965`](https://github.com/sankaranlab/SCAVENGE/tree/8ee8b173d965009a696b2a590d5b17b28b7cf851) selects high chromVAR z-score seed cells, constructs a binary mutual-nearest-neighbour graph, performs a column-normalized random walk with restart, caps and rescales the propagation score into a trait relevance score (TRS), and uses degree-matched seed permutations to identify significant cells.

**Reason for reimplementation.** The reference package's dependency stack predates the pipeline's R/Bioconductor environment. multiomeR needs sparse propagation over its native SNN graphs and must avoid materializing a cell-by-permutation score matrix for large cell sets.

**Deliberate deviations and consequences.** The reference builds a mutual-kNN graph, whereas multiomeR uses the binary support of its BPCells-derived SNN graph; edge weights are discarded, but topology can differ. Only per-cell exceedance counts and the cluster medians needed downstream are retained from the permutations, and random walks rather than random-number generation are parallelized, so the sampled null does not depend on the core count. Seed and scale-factor helpers guarantee at least one selected cell for small inputs, the degree sampler handles one-cell strata explicitly, and the random walk has a maximum-iteration guard. Cluster-level permutation medians, add-one P-values, and Benjamini–Hochberg adjustment within each grouping column are pipeline extensions.

**Implementation.** `R/SCAVENGE_helpers.R` and `src/scavenge_random_walk.cpp`; `module_genetic_enrichment/SCAVENGE_graph_targets.R` builds the graph and result records, and `SCAVENGE_group_targets.R` combines summaries and plots.

**Validation.** `tests/testthat/test-scavenge-parity.R` uses a deterministic fixture with heterogeneous-degree graph blocks and nonuniform edge weights, so it also tests conversion to binary adjacency. The iterative random walk must match the closed-form solution

\[ s = r\left(I - (1-r)P\right)^{-1}p_0 \]

within 1e-10. Compact local reference functions reproduce the relevant SCAVENGE code at the pinned commit without installing its dependency stack; propagation and transformed scores must agree within 1e-12, and the degree-matched seed samples, streamed exceedance counts, empirical P-values, significant-cell calls and one- versus two-core results must be identical. This does not establish parity of graph construction, chromVAR inputs, or biological interpretation.

```bash
pixi run --use-environment-activation-cache test-algorithm-validation
```


## Part: Target graph views


<!-- source: website/implementation/implementation_main.md -->

# Primary module



The root `_targets.R` creates GEM-well and aggregation mapping rows, then maps target fragments from `extra_targets/`. Use the diagrams to find the relevant stage, then inspect the corresponding source file for the complete command and resource declaration.

| Stage | Primary source |
|------------------------------------|------------------------------------|
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


<!-- source: website/implementation/implementation_differential_analyses.md -->

# Differential analyses



`module_differential_analyses/targets.R` filters the aggregations that enabled the module, joins their module configuration, attaches symbols for the accepted WNN metadata and pseudobulk inputs, and maps the target files in `module_differential_analyses/`. One generic pseudobulk model family is instantiated for each tested feature matrix.

The graph below is an orientation view; the target files hold the complete model and plotting commands. The methods, model routes and settings are in [Differential analyses methods](methods_differential_analyses.md), and the prerequisites and module selector in [Differential analyses](../downstream_differential_analyses.html).

[Mermaid graph omitted; source: `website/figures/human_curated/differential_analyses_v2.mmd`]


<!-- source: website/implementation/implementation_genetic_enrichment.md -->

# Genetic enrichment



`module_genetic_enrichment/targets.R` selects the aggregations whose `modules` include `genetic_enrichment`, resolves one configured Open Targets study set per aggregation, attaches symbols for the primary-module inputs it consumes, and maps the target files in `module_genetic_enrichment/`. It does not check the species; the GRCh38 GWAS inputs assume a human aggregation.

The methods and settings are in [Genetic enrichment methods](methods_genetic_enrichment.md), and the release, method-selection, and interpretation contracts in [Genetic enrichment](../downstream_genetic_enrichment.html).

## Single-nucleus and graph-based enrichment

This view covers the configured GWAS inputs, single-nucleus enrichment state, graph propagation, and downstream trait summaries. Additional cell-type contribution and locus-attribution branches may be pruned from this compact orientation view; use the manifest for the complete graph.

The sparse SCAVENGE reimplementation, deliberate graph and permutation differences, and validation evidence are recorded in [Algorithmic implementations, deviations and validation](algorithm_validation.md#sparse-scavenge-propagation-and-significance).

[Mermaid graph omitted; source: `website/figures/human_curated/genetic_enrichment_single_nucleus_v2.mmd`]

## Cell-type pseudobulk enrichment

This view covers the annotation-class pseudobulk deviations that are calculated separately from the nucleus-level results.

[Mermaid graph omitted; source: `website/figures/human_curated/genetic_enrichment_cell_type_absolute_effect_v2.mmd`]

## Cell-type contributions and locus attribution

This view covers the per-cell-type contribution and locus-attribution branches that are pruned from the single-nucleus view above.

[Mermaid graph omitted; source: `website/figures/human_curated/genetic_enrichment_cell_type_contributions_v2.mmd`]


<!-- source: website/implementation/implementation_peak_gene_correlation.md -->

# Peak–gene correlation



`module_peak_gene_correlation/targets.R` maps only opted-in aggregations, binds the primary-module inputs they consume, and maps the target files in `module_peak_gene_correlation/`. Parameters use the `peak_gene_correlation` manifest scope, targets end in `.peak_gene_correlation.<aggregation>`, and plot checkpoint tags use `peak_gene_correlation`, outside the numbered QC selections.

The methods and settings are in [Peak–gene correlation methods](methods_peak_gene_correlation.md), and the prerequisites and module selector in [Peak–gene correlation](../downstream_peak_gene_correlation.html).

## Candidate pairs, pseudobulks and tests

This view covers the TSS table, candidate peak–gene pairs, the broad WNN cell groups, donor–state pseudobulk construction, measurement-support filtering, the per-chromosome analysis branches and the diagnostic and top-link outputs.

[Mermaid graph omitted; source: `website/figures/human_curated/peak_gene_correlation_v2.mmd`]


## Part: Methods and parameters


<!-- source: website/implementation/methods_preprocessing_and_QC.md -->

# Preprocessing and nucleus QC

This chapter covers per-GEM-well preprocessing and the successive nucleus filters up to the final WNN cell set. The target structure is shown in the [primary-module graph](implementation_main.md), and GEM-well settings are columns of `cfg_GEM_wells.tsv`, described in [GEM well table](../reference_GEM_wells.html). The AMULET and ATAC scDblFinder adaptations are compared with their references in [Algorithmic implementations](algorithm_validation.md). The tables follow the layout in [Methods and parameter tables](implementation_conventions.md#methods-and-parameter-tables).


<!-- begin include: website/implementation/_shared_methods/quality_control.md -->

Each GEM well supplies a Cell Ranger ARC count directory. Gene-expression counts are imported into BPCells from the filtered feature-barcode matrix, or from a CellBender output when configured, and ATAC fragments from the Cell Ranger fragment file [@parks2025_bpcells; @fleming2023_cellbender]. Both carry a GEM-well-prefixed barcode, so nuclei remain distinct across wells. Per-well quality-control metrics are calculated for the barcodes that Cell Ranger called as cells: ATAC metrics such as TSS enrichment and nucleosome signal with BPCells, and gene-expression metrics such as the mitochondrial fraction from the imported matrix. Optional genotype demultiplexing with cellsnp-lite and Vireo replaces the GEM-well donor identifier with the assigned donor and records doublet and unassigned calls as metadata [@huang2021_cellsnp_lite; @huang2019_vireo]. A BPCells-native implementation of AMULET streams the stored fragments while reproducing the overlap metrics and q-values of `scDblFinder::amulet()` [@thibodeau2021_amulet; @germain2022_scdblfinder]; it supports unique fragments only and applies no automatic filter.

Configured GEM-well exclusions are filter expressions evaluated on the called cells of each well, so demultiplexing and AMULET results can serve as exclusion criteria before aggregation-level analysis. The remaining nuclei of the selected wells are combined without further filtering, and only barcodes represented in the aligned count data continue into dimensional reduction.

After GEX graph construction and Leiden clustering, clusters below the configured minimum size are removed; the same rule is later applied to the ATAC and WNN clusters. `scDblFinder` is then applied per GEM well to the raw counts of the retained nuclei [@germain2022_scdblfinder], using annotation-derived cluster labels: the cell-type label of an assigned cluster, otherwise its annotation status and identifier, optionally collapsed by a configured map. Doublet removal is configurable at the nucleus level and at the cluster level, where a Leiden cluster is removed when its doublet fraction exceeds a threshold.

The ATAC branch starts from the GEX-retained nuclei: peak-calling groups, the consensus peak matrix and peak-level quality-control metrics are built on that cell set, and the metrics are evaluated against configured aggregation-level exclusions. After ATAC clusters below the minimum size are removed, a separate ATAC-based `scDblFinder` analysis follows. Instead of scDblFinder's internal per-well feature aggregation, peaks are grouped once by *k*-means clustering of the aggregation-wide LSI loadings and summed with BPCells before the unchanged classifier is called. This adaptation bounds memory use; its scores are not expected to equal those obtained with the internal aggregation. ATAC doublets are removed at the nucleus and cluster level as for GEX.

WNN analysis uses the nuclei retained by the ATAC branch that are present in both embeddings and applies the minimum-cluster-size filter to the joint clusters; when nuclei are removed, the joint graph and embedding are recomputed on the retained nuclei while the cluster labels are kept. The GEX, ATAC and WNN metadata therefore represent successive cell universes rather than interchangeable versions of the original Cell Ranger calls.

<!-- end include: website/implementation/_shared_methods/quality_control.md -->


## Aggregation inputs and operational settings

These manifest parameters select inputs, plot variables and execution behaviour rather than algorithm settings.

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Inputs | GEM wells combined | Configurable |  | [`aggregation_GEM_well_IDs`](../parameters.html#aggregation_GEM_well_IDs) |
| Inputs | Donor metadata table | Configurable |  | [`aggregation_donor_id_metadata_tsv`](../parameters.html#aggregation_donor_id_metadata_tsv) |
| Inputs | Aggregation active | Configurable |  | [`is_active`](../parameters.html#is_active) |
| Inputs | Optional modules enabled | Configurable |  | [`modules`](../parameters.html#modules) |
| Plots | Categorical and continuous metadata plotted | Configurable |  | [`aggregation_categorical_vars`](../parameters.html#aggregation_categorical_vars), [`aggregation_continuous_vars`](../parameters.html#aggregation_continuous_vars) |
| Plots | Additional genes plotted | Configurable |  | [`aggregation_other_interesting_genes`](../parameters.html#aggregation_other_interesting_genes) |
| Tracks | Roadmap epigenome tracks | Configurable |  | [`aggregation_roadmap_EDACC_names`](../parameters.html#aggregation_roadmap_EDACC_names) |

## Per-GEM-well inputs and metrics

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| GEX input | CellBender versus Cell Ranger matrix | Configurable |  | `GEM_well_add_cellbender`, `GEM_well_cellbender_h5_file` |
| GEX input | Feature type kept | Fixed | `Gene Expression` | `GEX_counts_BPcells_matrix_dir` |
| Barcodes | Called-cell universe | Fixed | barcodes with `is_cell == 1` in the Cell Ranger `per_barcode_metrics.csv` | `per_barcode_metrics_tibble` |
| ATAC metrics | TSS enrichment | Fixed | `BPCells::qc_scATAC()` on Ensembl gene TSSs with an empty blacklist; package default: 101 bp centre window, 100 bp flanks 1.9–2 kb up- and downstream, flank signal floored at 0.1 | `ATAC_qc_metrics_tibble`, `BPCells::qc_scATAC()` |
| ATAC metrics | Nucleosome signal | Fixed | mono-nucleosomal / sub-nucleosomal fragment counts | `ATAC_qc_metrics_tibble` |
| GEX metrics | Mitochondrial genes | Fixed | gene names matching `(?i)^MT-` | `GEX_basic_metadata_tibble` |
| Demultiplexing | VCF, donor count and donor label | Configurable |  | `GEM_well_donors_VCF_file`, `GEM_well_n_donors`, `GEM_well_donor_id` |
| Demultiplexing | cellsnp-lite | Fixed | ATAC BAM, Cell Ranger-called barcodes, `--minMAF 0.1`, `--minCOUNT 20`, `--UMItag None` | `cellsnp_dir`, `get_cellsnp_dir()` |
| Demultiplexing | Vireo | Fixed | donor genotypes from the VCF (`-t GT`); no genotype learning | `get_vireo_donor_ids_tibble()` |
| AMULET | Nuclei scored | Fixed | Cell Ranger-called barcodes; no fragment minimum | `amulet_metrics_tibble` |
| AMULET | Fragment filters | Fixed | fragments ≤ 1,000 bp; chrM, chrX, chrY and their aliases excluded; Cell Ranger end coordinate shifted by −1 | `calculate_amulet_metrics_BPCells()`, `get_amulet_fragment_overlaps_BPCells()` |
| AMULET | High-overlap-site removal | Fixed | on; loci with Poisson P \< 0.01 across nuclei removed | `calculate_amulet_metrics_BPCells()`, `remove_high_overlap_amulet_loci()` |
| AMULET | Per-nucleus test | Fixed | upper-tail Poisson on the number of loci covered by more than two fragments; BH q-values | `calculate_amulet_metrics_BPCells()` |

## GEM-well exclusions and the GEX cell universe

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Exclusions | Per-well exclusion expressions | Configurable |  | `GEM_well_QC_exclude_list` |
| Cell set | Nuclei entering GEX PCA | Fixed | Cell Ranger-called nuclei passing the per-well exclusions, intersected with the combined matrix barcodes | `GEX_cellranger_kept_metadata_tibble`, `run_GEX_PCA_BPCells()` |

## Cluster-size filter and GEX doublets

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Cluster filter | Minimum cluster size, applied to GEX, ATAC, WNN and subgroups | Configurable |  | [`aggregation_cluster_min_barcodes`](../parameters.html#aggregation_cluster_min_barcodes) |
| Cluster filter | Comparison and disabling rule | Fixed | keep clusters with size ≥ threshold; disabled when NULL, NA or ≤ 1 | `filter_clusters_by_min_barcodes()` |
| GEX scDblFinder | Nuclei and counts | Fixed | GEX nuclei after the cluster-size filter and before doublet removal; raw counts; one run per GEM well | `scDblFinder_GEM_well_tibble.GEX`, `scDblFinder_results_by_GEM_well_tibble.GEX` |
| GEX scDblFinder | Label collapse map, reused for ATAC | Configurable |  | [`aggregation_scDblFinder_GEX_cell_type_collapse_list`](../parameters.html#aggregation_scDblFinder_GEX_cell_type_collapse_list) |
| GEX scDblFinder | Classifier arguments | Fixed | `dbr.sd = 1.0`; package defaults otherwise | `scDblFinder_results_by_GEM_well_tibble.GEX`, `scDblFinder::scDblFinder()` |
| GEX doublets | Nucleus-level removal | Configurable |  | [`aggregation_scDblFinder_GEX_remove_called_doublets`](../parameters.html#aggregation_scDblFinder_GEX_remove_called_doublets) |
| GEX doublets | Cluster-level removal | Configurable |  | [`aggregation_scDblFinder_GEX_max_doublet_fraction_per_cluster`](../parameters.html#aggregation_scDblFinder_GEX_max_doublet_fraction_per_cluster) |
| GEX doublets | Cluster fraction rule, also used for ATAC | Fixed | called-doublet fraction of each Leiden cluster, computed before nucleus-level removal; cluster removed when the fraction is strictly greater than the threshold | `filter_metadata_by_scDblFinder()` |

## ATAC QC and doublets

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Cell set | Nuclei entering the ATAC branch | Fixed | GEX nuclei after doublet removal | `BCs_per_peak_cluster_list.ATAC`, `consensus_peak_BPCells_matrix_dir.ATAC`, `metadata_w_QC_tibble.ATAC` |
| Peak QC | Metrics | Fixed | `nCount_ATAC` (counts in consensus peaks); blacklist counts / `nCount_ATAC`; `nCount_ATAC` / Cell Ranger ATAC fragments; that fraction divided by the genome fraction covered by peaks | `get_ATAC_QC_metadata_from_BPCells()` |
| Peak QC | Exclusion expressions | Configurable |  | [`aggregation_QC_exclude_list_combined_object`](../parameters.html#aggregation_QC_exclude_list_combined_object) |
| ATAC scDblFinder | Nuclei scored | Fixed | ATAC nuclei after the cluster-size filter that pass peak QC | `scDblFinder_GEM_well_tibble.ATAC` |
| ATAC scDblFinder | Feature groups | Fixed | 50 groups from `stats::kmeans(iter.max = 50, nstart = 1)`, seed 1, on the peak loadings of LSI dimensions `intersect(2:20, aggregation_ATAC_data_PCs)` | `scDblFinder_feature_groups.ATAC`, `get_feature_groups_from_LSI_loadings()` |
| ATAC scDblFinder | Classifier arguments | Fixed | `dbr.sd = 1.0`, `aggregateFeatures = FALSE`, `nfeatures = 50`, `processing = "normFeatures"`; package defaults otherwise | `scDblFinder_results_by_GEM_well_tibble.ATAC`, `scDblFinder::scDblFinder()` |
| ATAC scDblFinder | Cluster labels | Fixed | ATAC annotation-derived scDblFinder groups, collapsed with the GEX map | `scDblFinder_GEM_well_tibble.ATAC` |
| ATAC doublets | Nucleus-level removal | Configurable |  | [`aggregation_scDblFinder_ATAC_remove_called_doublets`](../parameters.html#aggregation_scDblFinder_ATAC_remove_called_doublets) |
| ATAC doublets | Cluster-level removal | Configurable |  | [`aggregation_scDblFinder_ATAC_max_doublet_fraction_per_cluster`](../parameters.html#aggregation_scDblFinder_ATAC_max_doublet_fraction_per_cluster) |

## WNN cell set

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Cell set | Nuclei entering WNN | Fixed | ATAC nuclei after doublet removal that have rows in both corrected embeddings | `embedding_matrices.WNN`, `get_WNN_embedding_matrices()` |
| Cluster filter | Threshold | Configurable |  | [`aggregation_cluster_min_barcodes`](../parameters.html#aggregation_cluster_min_barcodes) |


<!-- source: website/implementation/methods_GEX_ATAC_and_WNN.md -->

# GEX, ATAC, batch correction and WNN

This chapter covers normalization and dimensional reduction of both modalities, peak definition, batch correction, graph construction and clustering, and weighted nearest-neighbour (WNN) integration. The target structure is shown in the [primary-module graph](implementation_main.md), and the WNN implementation is compared with Seurat in [Algorithmic implementations](algorithm_validation.md#native-weighted-nearest-neighbors). The tables follow the layout in [Methods and parameter tables](implementation_conventions.md#methods-and-parameter-tables).

## GEX normalization and PCA


<!-- begin include: website/implementation/_shared_methods/GEX_normalization.md -->

For each aggregation, the quality-controlled gene-expression matrices are combined with their GEM-well-prefixed barcodes, and genes with low total counts are excluded before dimensional reduction. Normalization and dimensional reduction use either the disk-backed BPCells Pearson-residual workflow or Seurat `SCTransform` [@parks2025_bpcells; @hafemeister2019_sctransform]. In the BPCells branch, Pearson residuals are computed with a per-gene method-of-moments overdispersion estimate and clipped to a fixed range, configured cell-level covariates can be regressed out, the genes with the highest residual variance are retained and a truncated singular-value decomposition yields the principal components. The Seurat branch applies SCTransform to the same cells and genes, retains the same number of variable features and computes the principal components from the residual matrix. Regression of a cell-cycle difference score is the default; the S and G2M scores are computed from log-normalized counts with the Seurat cell-cycle gene sets.

<!-- end include: website/implementation/_shared_methods/GEX_normalization.md -->


| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Gene filter | Minimum total count per gene, both backends | Fixed | more than 50 counts over the PCA nuclei (`min_feature_count = 50`) | `run_GEX_PCA_BPCells()` |
| Backend | Normalization and PCA backend | Configurable |  | [`aggregation_GEX_PCA_backend`](../parameters.html#aggregation_GEX_PCA_backend) |
| Regression | Cell-level covariates regressed | Configurable |  | [`aggregation_SCT_regress_vars`](../parameters.html#aggregation_SCT_regress_vars) |
| Components | Number of PCs computed (last element of the list) | Configurable |  | [`aggregation_GEX_data_PCs`](../parameters.html#aggregation_GEX_data_PCs) |
| Variable genes | Genes retained by residual variance, both backends | Fixed | 3,000 (`n_variable_features = 3000`) | `run_GEX_PCA_BPCells()` |
| BPCells branch | Pearson residuals | Fixed | method-of-moments theta clamped to \[1e-6, 1e6\]; `clip_range = c(-10, 10)`; `min_var = 0` | `run_BPCells_native_GEX_PCA()` |
| BPCells branch | Regression and SVD | Fixed | `BPCells::regress_out(prediction_axis = "row")` on the residuals; `BPCells::svds()` without centring, embeddings = right singular vectors × singular values | `run_BPCells_native_GEX_PCA()` |
| Seurat branch | SCTransform | Fixed | `conserve.memory = TRUE`, `do.correct.umi = FALSE`; package defaults otherwise, including `vst.flavor = "v2"` and `ncells = 5000` | `run_Seurat_SCT_for_PCA()`, `Seurat::SCTransform()` |
| Seurat branch | PCA | Fixed | exact eigendecomposition of the gene × gene residual Gram matrix | `run_dense_feature_gram_PCA()` |
| Cell cycle | Gene sets and scoring | Fixed | `Seurat::cc.genes.updated.2019` for both organisms; log-normalization with scale factor 10,000; binned control-gene scores with seed 1; `CC.Difference = S − G2M` | `cell_cycle_gene_sets()`, `add_cell_cycle_scores_to_cell_attr()`, `calculate_BPCells_cell_cycle_scores_from_matrix()` |

## ATAC peak definition and dimensional reduction


<!-- begin include: website/implementation/_shared_methods/ATAC_peaks_and_LSI.md -->

Fragments on the standard chromosomes are used for peak calling. Peak-calling groups are defined by a configurable metadata column of the GEX-retained nuclei, and large groups are downsampled for peak discovery. Peaks are called per group either with MACS3, on exported group fragments with ATAC-style shift and extension and summit calling, or with the BPCells tile-based caller on the stored fragments [@zhang2008_macs; @parks2025_bpcells]. Consensus peaks follow the ArchR strategy of iterative fixed-width peak construction [@granja2021_archr]: summits are extended to fixed-width peaks, peaks overlapping the reference-genome blacklist are removed, and the group-level peak sets are combined into one non-overlapping consensus set, ranking overlapping peaks by summit significance within a group and by fold enrichment across groups. BPCells counts fragment insertions or overlaps, as configured, within these peaks for the GEX-retained nuclei.

ATAC dimensional reduction applies Signac's TF-IDF method 1 to the consensus peak matrix, followed by a BPCells singular-value decomposition to obtain latent semantic indexing (LSI) embeddings [@stuart2021_signac; @parks2025_bpcells]. No component is removed by the code; the default dimension settings omit the first component.

<!-- end include: website/implementation/_shared_methods/ATAC_peaks_and_LSI.md -->


### ATAC peak calling and consensus peaks

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Fragments | Chromosomes retained | Fixed | autosomes, X and Y | `combined_BPCells_fragment_obj.ATAC` |
| Groups | Peak-calling grouping column | Configurable |  | [`aggregation_call_peaks_by_cluster_col`](../parameters.html#aggregation_call_peaks_by_cluster_col) |
| Groups | Discovery cap per group | Fixed | 50,000 nuclei; larger groups randomly downsampled with seed 1 + group index | `build_peak_calling_cluster_discovery_tibble()` |
| Caller | Peak-calling method | Configurable |  | [`aggregation_ATAC_peak_calling_method`](../parameters.html#aggregation_ATAC_peak_calling_method) |
| MACS3 | Arguments | Fixed | `-f BED --nomodel --shift -75 --extsize 150 --call-summits --keep-dup all`; MACS3 default q-value cutoff 0.05 | `call_peaks_w_MACS3()` |
| Genome | Effective genome size, used by both callers and the peak-count enrichment | Fixed | GRCh38 2.913e9; mm10 and GRCm39 2.65e9 | `get_effective_genome_size()` |
| Tile caller | BPCells settings | Fixed | `peak_width = 500`, `peak_tiling = 3`, `fdr_cutoff = 0.01`, `merge_peaks = "none"` | `call_peaks_w_BPCells_tile()` |
| Peak shape | Summit extension | Fixed | 250 bp either side of the summit, giving 500 bp peaks | `get_peak_GRanges_w_fixed_width()` |
| Blacklist | Source and rule | Fixed | GRCh38 `hg38.Kundaje.GRCh38_unified_Excludable` and mm10 `mm10.Boyle.mm10-Excludable.v2` from AnnotationHub excluderanges; GRCm39 `resources/mm39.excluderanges.bed`; any overlap removes a peak | `get_blacklist_GRanges()`, `get_peak_GRanges_w_fixed_width()` |
| Consensus | Overlap ranking | Fixed | within groups `neg_log10pvalue_summit`, across groups `fold_change`, both decreasing | `within_clusters_collapsed_peaks_per_cluster_GRanges.ATAC`, `consensus_peak_GRanges.ATAC` |
| Peak matrix | Counting mode, also used for blacklist QC counts | Configurable |  | [`aggregation_ATAC_peak_matrix_mode`](../parameters.html#aggregation_ATAC_peak_matrix_mode) |
| Peak matrix | Nuclei | Fixed | GEX nuclei after doublet removal | `consensus_peak_BPCells_matrix_dir.ATAC` |

### ATAC TF-IDF and LSI

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| TF-IDF | Transform | Fixed | `log1p(10000 × TF × IDF)`, with TF = count / nucleus total and IDF = number of nuclei / peak total | `run_ATAC_LSI_BPCells()` |
| SVD | Components computed (last element of the list) | Configurable |  | [`aggregation_ATAC_data_PCs`](../parameters.html#aggregation_ATAC_data_PCs) |
| SVD | Decomposition | Fixed | `BPCells::svds()` without centring; embeddings = right singular vectors × singular values | `run_ATAC_LSI_BPCells()` |
| Dimensions | Components used downstream | Configurable |  | [`aggregation_ATAC_data_PCs`](../parameters.html#aggregation_ATAC_data_PCs), [`aggregation_UMAP_ATAC_PCs`](../parameters.html#aggregation_UMAP_ATAC_PCs) |

## Batch correction, clustering and weighted nearest neighbours


<!-- begin include: website/implementation/_shared_methods/batch_correction_clustering_WNN.md -->

Harmony correction can be applied separately to the selected GEX principal components and ATAC LSI dimensions [@patikas2026_harmony2]. Several configured metadata columns are combined into one interaction batch factor, nuclei with missing covariate values are removed with a warning, and the embeddings pass through unchanged when no columns are configured.

For each modality, approximate nearest neighbours are found on the selected dimensions with BPCells HNSW search, converted to a shared-nearest-neighbour (SNN) graph with Jaccard weights and clustered with the Leiden algorithm [@parks2025_bpcells; @traag2019_leiden]. The neighbour count includes the query cell itself, and clusters are renumbered by size. UMAP embeddings are computed with uwot on the same dimensions [@melville2026_uwot], and quality-control sweeps additionally render UMAPs over grids of dimensions and neighbour counts.

RNA and ATAC representations are combined with a BPCells-native implementation of the Seurat weighted-nearest-neighbour (WNN) strategy [@hao2021_multimodal; @parks2025_bpcells]. The GEX and ATAC embeddings, after optional Harmony correction, are aligned by barcode and L2-normalized, and candidate neighbours are found separately for each modality with HNSW search. A native routine applies Seurat's small-SNN bandwidth strategy to set a per-cell kernel width from the configured neighbour count, per-cell modality weights are derived from within- versus cross-modality prediction kernels, and the union of candidates is ranked by the weighted kernel score to select the final neighbours. BPCells builds the joint SNN graph, which is clustered with Leiden and embedded with UMAP using the precomputed neighbours. The implementation does not call `Seurat::FindMultiModalNeighbors()`; because HNSW replaces Seurat's Annoy search, exact equality is not expected, and its validation is a similarity contract on modality weights and neighbour overlap.

<!-- end include: website/implementation/_shared_methods/batch_correction_clustering_WNN.md -->


### Harmony batch correction

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Covariates | Shared correction columns | Configurable |  | [`aggregation_harmony_correction_metadata_col_names`](../parameters.html#aggregation_harmony_correction_metadata_col_names) |
| Covariates | Additional ATAC covariates | Configurable |  | [`aggregation_extra_harmony_covars_ATAC`](../parameters.html#aggregation_extra_harmony_covars_ATAC) |
| Dimensions | Corrected dimensions | Configurable |  | [`aggregation_GEX_data_PCs`](../parameters.html#aggregation_GEX_data_PCs), [`aggregation_ATAC_data_PCs`](../parameters.html#aggregation_ATAC_data_PCs) |
| Harmony | Arguments | Fixed | `max_iter = 25`, `lambda = 1`; package defaults otherwise | `run_harmony_on_embedding_matrix()`, `harmony::RunHarmony()` |

### Graph construction, Leiden clustering and UMAP

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| kNN | Neighbour count | Configurable |  | [`aggregation_data_nNNs`](../parameters.html#aggregation_data_nNNs) |
| kNN | Search | Fixed | `BPCells::knn_hnsw()` with cosine metric and `ef = 500` | `cluster_embedding_matrix_BPCells()` |
| SNN | Construction and pruning | Fixed | package default: Jaccard weights, `min_val = 1/15`, no self loops | `cluster_knn_snn_leiden()`, `BPCells::knn_to_snn_graph()` |
| Leiden | Resolution | Configurable |  | [`aggregation_GEX_cluster_res`](../parameters.html#aggregation_GEX_cluster_res), [`aggregation_ATAC_cluster_res`](../parameters.html#aggregation_ATAC_cluster_res), [`aggregation_WNN_cluster_res`](../parameters.html#aggregation_WNN_cluster_res) |
| Leiden | Objective, iterations and seed | Fixed | modularity with SNN edge weights, seed 1; package defaults: 2 iterations, `beta = 0.01` | `cluster_knn_snn_leiden()`, `igraph::cluster_leiden()` |
| UMAP | Dimensions | Configurable |  | [`aggregation_UMAP_GEX_PCs`](../parameters.html#aggregation_UMAP_GEX_PCs), [`aggregation_UMAP_ATAC_PCs`](../parameters.html#aggregation_UMAP_ATAC_PCs) |
| UMAP | Neighbours and minimum distance | Configurable |  | [`aggregation_UMAP_nNNs`](../parameters.html#aggregation_UMAP_nNNs), [`aggregation_UMAP_min_dist`](../parameters.html#aggregation_UMAP_min_dist) |
| UMAP | Fixed arguments | Fixed | cosine metric, 2 components, seed 1; package defaults otherwise, including spectral initialisation | `run_UMAP_from_embedding_matrix()`, `uwot::umap()` |
| UMAP sweeps | Grids | Fixed | 3 dimension counts from 5 to the number of data dimensions; 3 neighbour counts from 10 to the configured UMAP neighbour count | `UMAP_n_dims_seq.GEX`, `UMAP_n_dims_seq.ATAC`, `UMAP_neighbors_seq` |

### Weighted nearest neighbours

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Inputs | Embeddings after optional Harmony, and their dimensions | Configurable |  | [`aggregation_GEX_data_PCs`](../parameters.html#aggregation_GEX_data_PCs), [`aggregation_ATAC_data_PCs`](../parameters.html#aggregation_ATAC_data_PCs) |
| Candidates | Candidates per modality | Fixed | `candidate_k = 200`, from `BPCells::knn_hnsw()` on the L2-normalized embeddings with Euclidean metric, `k = candidate_k + 1` and `ef = 500` | `WNN_results_raw`, `WNN_results`, `weighted_nearest_neighbors_BPCells()` |
| Final neighbours | Neighbour count, also the bandwidth and imputation neighbourhood | Configurable |  | [`aggregation_data_nNNs`](../parameters.html#aggregation_data_nNNs) |
| Bandwidth | Small-SNN kernel width | Fixed | mean distance, minus the nearest non-self distance, to the k cells with the fewest shared neighbours among those sharing at least one; `sd_scale = 1`; floored at machine epsilon | `weighted_nearest_neighbors_BPCells()`, `calculate_small_SNN_bandwidth()` |
| Weights | Modality weight kernel | Fixed | `exp(−d/σ)`; ratio `within / (cross + 1e-4)` clipped to \[0, 200\]; softmax across modalities | `weighted_nearest_neighbors_BPCells()` |
| Selection | Weighted score and distance | Fixed | `Σ w_m · exp(−d_m/σ_m)` (`kernel_power = 1`); `nn_dist = sqrt((1 − score)/2)` | `weighted_nearest_neighbors_BPCells()` |
| SNN and Leiden | Graph and clustering | Fixed | as for GEX and ATAC: package-default SNN, Leiden modularity, seed 1 | `cluster_WNN_graph()`, `cluster_knn_snn_leiden()` |
| UMAP | Neighbours (capped at the final neighbour count) and minimum distance | Configurable |  | [`aggregation_UMAP_nNNs`](../parameters.html#aggregation_UMAP_nNNs), [`aggregation_UMAP_min_dist`](../parameters.html#aggregation_UMAP_min_dist) |
| UMAP | Fixed arguments | Fixed | precomputed WNN neighbours, 2 components, seed 1; package defaults otherwise | `run_WNN_UMAP()`, `uwot::umap()` |


<!-- source: website/implementation/methods_annotation_and_motifs.md -->

# Cell-type annotation and motif accessibility

This chapter covers marker-signature cluster annotation and motif-family accessibility. The BPCells-native UCell scorer is compared with UCell in [Algorithmic implementations](algorithm_validation.md#bpcells-native-ucell-scoring). The tables follow the layout in [Methods and parameter tables](implementation_conventions.md#methods-and-parameter-tables).

## Cell-type annotation


<!-- begin include: website/implementation/_shared_methods/cell_type_annotation.md -->

Cell-type annotation is driven by user-supplied GEX marker signatures. Signatures accept unsigned genes, positive markers and genes expected to be absent; genes missing from the count matrix fail validation before scoring. Scores are computed with a BPCells-native implementation of UCell scoring [@andreatta2021_ucell; @parks2025_bpcells], which reproduces per-cell descending ranks, rank truncation, positive and negative signatures and lower-bound clipping on bounded chunks of the raw GEX counts; the production workflow does not call the UCell package, and its validation requires values identical to the UCell reference. Signed signatures are clipped at zero per cell before averaging. The same scoring, with the same GEX control reference, annotates the GEX, ATAC and WNN cluster partitions.

Each label's cluster-level mean score is compared with random control signatures matched on gene abundance and detection. The control reference samples cells per GEM well from the GEX metadata before doublet filtering. For each control replicate, markers are visited in random order and each is replaced by a gene drawn from its nearest eligible candidates not yet used in that replicate; candidates exclude all marker genes and undetected genes. A label's adjusted score in a cluster is its observed mean minus an upper quantile of its matched controls. The label with the highest adjusted score is the candidate, and its advantage is the smaller of its lead over zero and its lead over the runner-up. The candidate is assigned when its advantage is positive, untied and at least the configured minimum; otherwise the cluster remains unassigned, with the candidate and reason retained, so raising the minimum can only withdraw assignments. Because annotation depends on the supplied signatures and their level of detail, labels can represent either cell types or broader source classes.

Three diagnostics accompany each decision without vetoing it. Marker stability is the fraction of leave-one-marker-out variants, each removing one candidate marker and its matched control, in which the candidate keeps a positive advantage of at least the configured minimum. Cell stability is the fraction of leave-one-block-out replicates, with nuclei split into blocks stratified by cluster and GEM well, that assign the same candidate. GEM-well agreement is the fraction of sufficiently large per-well subgroups that assign it. Detection counts additionally report the positive markers detected in a minimum fraction of nuclei.

<!-- end include: website/implementation/_shared_methods/cell_type_annotation.md -->


### Signature scoring

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Signatures | Marker genes per label | Configurable |  | [`aggregation_GEX_marker_genes`](../parameters.html#aggregation_GEX_marker_genes) |
| Ranks | Rank cap | Fixed | `min(1500, number of genes)` | `prepare_cluster_UCell_controls()` |
| Ranks | Direction and ties | Fixed | descending raw counts per nucleus; ties averaged | `rank_UCell_count_chunk()` |
| Scores | Signed signatures | Fixed | positive-gene score minus negative-gene score, clipped with `pmax(0, score)` per nucleus before averaging | `score_signed_UCell_cells()` |
| Scores | Partitions annotated | Fixed | GEX `PCA_harmony_SNN_cluster` before doublet removal; ATAC `LSI_harmony_SNN_cluster` before ATAC doublet removal; WNN `WNN_harmony_SNN_cluster` | `cluster_UCell_evidence.GEX`, `cluster_UCell_evidence.ATAC`, `cluster_UCell_evidence.WNN` |

### Matched-control cluster annotation

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Reference | Nuclei and counts | Fixed | up to 50 nuclei per GEM well, sampled from the pre-doublet-filter GEX nuclei; full aggregated GEX counts | `cluster_UCell_controls.GEX`, `prepare_cluster_UCell_controls()` |
| Controls | Random mappings | Fixed | 999 | `build_UCell_controls()` |
| Controls | Matching coordinates | Fixed | `log1p(abundance × 1e4)` and `asin(sqrt(detection))`, each standardized; Euclidean distance | `build_UCell_controls()` |
| Controls | Candidate pool and draw | Fixed | `neighbours = 50`; pool of the max(4 × neighbours, 2 × marker genes) nearest eligible genes, which is 200 for panels of up to 100 marker genes; each marker drawn uniformly from its `neighbours` nearest candidates not yet used in the replicate | `build_UCell_controls()` |
| Controls | Random seed | Fixed | 20260910 | `prepare_cluster_UCell_controls()`, `build_UCell_controls()` |
| Adjusted score | Background | Fixed | 0.95 quantile of the label's matched-control scores | `score_UCell_group_evidence()` |
| Assignment | Advantage rule | Fixed | `min(best, best − second)`; requires best \> 0 and no exact tie | `assign_UCell_cluster_evidence()` |
| Assignment | Minimum advantage | Configurable |  | [`aggregation_cluster_annotation_min_advantage`](../parameters.html#aggregation_cluster_annotation_min_advantage) |
| Diagnostics | Marker deletion | Fixed | each marker of the leading label, with its matched control, omitted in turn (labels with at least 2 markers); stability = fraction of deletions keeping a positive advantage at least the minimum | `score_cluster_UCell_summaries()`, `evaluate_cluster_UCell_evidence()` |
| Diagnostics | Cell-deletion blocks | Fixed | 10 blocks stratified by cluster and GEM well, each left out once; stability reported when at least 2 blocks are assessable | `make_annotation_blocks()`, `summarize_cluster_UCell_counts()`, `evaluate_cluster_UCell_evidence()` |
| Diagnostics | Marker detection | Fixed | positive markers detected in at least 10% of a cluster's nuclei | `score_UCell_group_evidence()` |
| Diagnostics | GEM-well agreement | Fixed | GEM wells with at least 25 nuclei in the cluster; agreement reported when at least 2 wells are assessed | `score_cluster_UCell_summaries()`, `evaluate_cluster_UCell_evidence()` |
| Plots | Marker-set order | Fixed | `hclust(method = "ward.D2")` on Euclidean distances between adjusted-score profiles | `plot_UCell_annotation_dot()` |

## Motif families and motif accessibility


<!-- begin include: website/implementation/_shared_methods/motif_accessibility.md -->

Transcription-factor motifs are represented by the sequence-similarity families of the JASPAR 2026 CORE vertebrate collection [@ovekbaydar2026_jaspar]. The pipeline vendors the familial root motifs and the family membership table and matches the root motifs to the consensus peaks with `motifmatchr` [@schep2025_motifmatchr]. Configured transcription factors are resolved to families by name, or by name and motif identifier when a symbol belongs to more than one family. betterChromVAR applies the chromVAR workflow with GC-bias correction to obtain analytic motif-family deviations and z-scores for individual nuclei [@schep2017_chromvar; @germain2026_betterchromvar]. All families are summarized per ATAC annotation class with a Wilcoxon marker test and as cell-weighted means per cluster and cell type; the configured transcription factors select only the families shown on embeddings. These values measure accessibility associated with a motif family and are not direct measurements of transcription-factor activity.

<!-- end include: website/implementation/_shared_methods/motif_accessibility.md -->


| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Families | JASPAR 2026 CORE vertebrate files | Fixed | 233 familial root motifs; family membership of 1,019 motifs | `resources/JASPAR2026_vertebrate_familial_root_motifs.tf`, `resources/JASPAR2026_vertebrate_motif_families.tsv` |
| Families | Configured transcription factors | Configurable |  | [`aggregation_ATAC_marker_TFs`](../parameters.html#aggregation_ATAC_marker_TFs) |
| Families | Symbol resolution | Fixed | case-insensitive TF name, or `TF__motifID` for symbols in several families | `resolve_marker_motif_families()` |
| Scanning | Motif matrices | Fixed | PFMs with a uniform 0.25 background on the `+` strand; all-zero columns dropped | `read_JASPAR_familial_root_PFMatrixList()` |
| Scanning | Match thresholds | Fixed | package defaults: `p.cutoff = 5e-5`, `bg = "subject"`, `w = 7` | `get_motif_matrix_from_peak_ranges()`, `motifmatchr::matchMotifs()` |
| Scanning | Genome | Fixed | BSgenome UCSC hg38, mm10 or mm39, matching the reference | `get_chromVAR_genome_obj()` |
| chromVAR | Nuclei | Fixed | ATAC nuclei after doublet removal | `chromVAR_obj.ATAC` |
| chromVAR | Peaks, bias and expectation | Fixed | zero-count peaks removed; `betterChromVAR::addGCBias()` with missing values set to 0; expectation = mean count per peak over nuclei | `get_chromVAR_obj_from_peak_matrix()`, `get_chromVAR_peak_expectation()` |
| chromVAR | Background | Fixed | package defaults: default background bins, no shrinkage (`shrinkage = "none"`) | `betterChromVAR::getBackgroundBins()`, `betterChromVAR::computeBackgrounds()` |
| chromVAR | Deviations | Fixed | analytic deviations and z-scores; package default `denominator = "global"` | `compute_chromVAR_annotation_chunk_result()`, `betterChromVAR::computeDeviationsAnalytic()` |
| Summaries | Per-cell-type test | Fixed | `BPCells::marker_features(method = "wilcoxon")` on z-scores by ATAC cell type; BH across all family-by-group tests; effect = difference in mean z-score | `get_marker_motif_family_accessibility_from_chromVAR_BPCells_z_scores()` |


<!-- source: website/implementation/methods_peak_gene_correlation.md -->

# Peak–gene correlation

This chapter covers the optional `peak_gene_correlation` module. The target structure is shown in the [peak–gene correlation graph](implementation_peak_gene_correlation.md), and the configuration in [Peak–gene correlation](../downstream_peak_gene_correlation.html). The compiled kernels are compared with their references in [Algorithmic implementations](algorithm_validation.md). The tables follow the layout in [Methods and parameter tables](implementation_conventions.md#methods-and-parameter-tables).


<!-- begin include: website/implementation/_shared_methods/peak_gene_correlation.md -->

The optional peak–gene correlation module runs on the accepted WNN cell set with its own configuration, review checkpoint and output directory. It operates separately within sufficiently large WNN annotation classes and requires donor identifiers and per-nucleus GEX and ATAC depth metadata; classes lacking them are skipped with a diagnostic. Consensus ATAC peaks are paired with gene transcription start sites (TSSs) on the same chromosome within a fixed window measured from the peak centre, and each pair is classified as self-promoter, gene-body, proximal or distal. Within each class, nuclei from eligible donors are partitioned into mutually exclusive ATAC-state bins by *k*-means clustering of scaled ATAC LSI or Harmony dimensions, with the number of bins adapting to the median number of nuclei per donor. Counts are summed within each donor–state combination to create non-overlapping pseudobulks; pseudobulks below a minimum size are dropped, and donors and states are pruned iteratively until every donor contributes a minimum number of states.

GEX and ATAC pseudobulk counts are separately scaled to counts per million, using the metadata-derived pseudobulk depth, and log1p-transformed. Genes and peaks must be detected in a minimum fraction of pseudobulks, and a chromosome branch requires minimum numbers of pseudobulks and residual degrees of freedom. A configurable measurement-support filter then removes hypotheses whose gene and peak are not both supported in enough shared donors, with count thresholds scaled by each pseudobulk's depth relative to the class median. Excluded hypotheses are never tested and do not enter the multiple-testing family.

GEX and ATAC values are residualized against donor and scaled log RNA and ATAC depth, retaining variation between ATAC states; the Pearson correlation of the residuals and a heteroskedasticity-robust (HC3) regression test give conditional association summaries. The hierarchical analysis fits donor fixed intercepts, depth covariates and a donor-varying slope for the within-donor-centred peak value, using project-owned compiled kernels for profiled restricted maximum likelihood and Kenward–Roger inference for the average slope; `lme4` and `pbkrtest` serve only as test references. It requires within-donor peak variation in at least two donors, and because every support preset requires at least two shared donors, single-donor classes yield diagnostics but no tests. Fits failing numerical diagnostics retain their estimates but no inferential P-value. Hierarchical P-values are Benjamini–Hochberg-corrected within each annotation class over the complete eligible pair family, counting unreliable tests; conditional P-values are corrected within each class over the non-missing values. SuSiE fine-mapping prioritizes peaks for genes with at least one conditional link, using the donor- and depth-residualized values [@wang2020_susie]. These model-based associations do not establish causal enhancer–gene regulation.

Top-link plots rank estimable positive hierarchical associations by nominal P-value, without a significance cutoff, excluding self-promoter peaks; gene-body peaks remain eligible. Each plot shows the gene context, the focal class's insertion coverage and the donor-residual scatter, and donor-level direction and covariance diagnostics help identify associations dominated by one donor.

<!-- end include: website/implementation/_shared_methods/peak_gene_correlation.md -->


## Cell groups, candidate pairs and donor–state pseudobulks

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Cell groups | Annotation column and minimum class size | Fixed | `WNN_harmony_SNN_cluster_cell_type`; at least 200 nuclei | `make_peak_gene_correlation_cell_groups()` |
| Pairs | TSS and search window | Fixed | TSS at the gene start, or the gene end on the minus strand; peak centre within 250 kb of the TSS | `make_peak_gene_correlation_gene_TSS_tibble()`, `make_peak_gene_correlation_candidate_pairs()` |
| Pairs | Self-promoter window and link classes | Fixed | strand-aware −1,500 to +500 bp around the TSS; self-promoter, gene body, proximal (≤ 10 kb), distal, in that precedence | `make_peak_gene_correlation_candidate_pairs()` |
| Donors | Donor and depth columns | Fixed | `donor_id`; first available of `nCount_RNA`, `gex_umis_count` and of `nCount_ATAC`, `atac_fragments` | `make_peak_gene_correlation_donor_state_record()` |
| States | Embedding | Configurable |  | [`aggregation_ATAC_data_PCs`](../parameters.html#aggregation_ATAC_data_PCs), [`aggregation_harmony_correction_metadata_col_names`](../parameters.html#aggregation_harmony_correction_metadata_col_names), [`aggregation_extra_harmony_covars_ATAC`](../parameters.html#aggregation_extra_harmony_covars_ATAC) |
| States | Dimensions and scaling | Fixed | LSI dimensions 2–20 present in the embedding; zero-variance dimensions dropped and the rest standardized within the class | `make_peak_gene_correlation_donor_state_record()` |
| States | Number of bins | Fixed | `min(20, floor(median nuclei per eligible donor / 20))`; classes with fewer than 2 bins are skipped | `make_peak_gene_correlation_donor_state_record()` |
| States | k-means | Fixed | Lloyd algorithm, `iter.max = 1000`, `nstart = 1`, seed 1 | `make_peak_gene_correlation_donor_state_record()` |
| Pseudobulks | Minimum sizes | Fixed | donors need at least 40 nuclei before binning; at least 20 nuclei per donor–state, 2 states per donor, and 1 donor per class and per state | `make_peak_gene_correlation_donor_state_record()` |

## Normalization, eligibility and measurement-support filtering

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Normalization | Scaling | Fixed | counts per million (`scale_factor = 1e6`) of the metadata-derived pseudobulk depth, then `log1p` | `normalize_peak_gene_correlation_aggregate_matrices()` |
| Eligibility | Detection | Fixed | genes detected and peaks accessible in at least 5 % of pseudobulks | `prepare_peak_gene_correlation_branch()` |
| Eligibility | Branch size | Fixed | at least 10 pseudobulks; pseudobulks − design rank − 1 ≥ 5 | `prepare_peak_gene_correlation_branch()` |
| Nuisance design | Covariates | Fixed | donor fixed effects; standardized `log1p` GEX and ATAC depth when non-constant; QR-pruned to full rank | `make_peak_gene_correlation_design_matrix()` |
| Support filter | Preset | Configurable |  | [`peak_gene_correlation_filter`](../parameters.html#peak_gene_correlation_filter) |
| Support filter | Preset thresholds | Fixed | RNA count, ATAC count, supporting aggregates, shared donors, supporting aggregates per donor. Lenient: 5, 3, `max(6, 10%)`, 2, 2; moderate: 10, 5, `max(6, 10%)`, 2, 2; strict: 10, 5, `max(10, 20%)`, 3, 3 | `peak_gene_filter_settings()` |
| Support filter | Depth scaling of count thresholds | Fixed | `max(2, count × aggregate depth / median depth)` | `filter_peak_gene_candidate_pairs()` |

## Conditional and hierarchical tests

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| HC3 | Test distribution | Fixed | two-sided t with pseudobulks − design rank − 1 degrees of freedom | `score_peak_gene_correlations_for_cell_group()` |
| Multiplicity | BH family | Fixed | within each class; HC3 over the non-missing P-values; hierarchical over all eligible pairs, including unreliable fits | `finalize_peak_gene_correlation_results()`, `finalize_peak_gene_hierarchical_results()` |
| Links | Conditional link | Fixed | correlation ≥ 0.15, FDR \< 0.05, not self-promoter | `make_peak_gene_correlation_links()` |
| Hierarchical | Donor requirement | Fixed | at least 2 donors with within-donor peak variation | `score_peak_gene_hierarchical_associations()` |
| Hierarchical | Unreliable fits | Fixed | kernel diagnostic raised or Kenward–Roger df \< 1: estimate kept, P-value missing | `score_peak_gene_hierarchical_associations()` |

## Prioritization, top links and plots

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| SuSiE | Genes and candidate peaks | Fixed | genes with a conditional link; up to 50 genes per branch by best FDR and 500 peaks per gene by absolute correlation; at least 2 variable peaks | `finemap_peak_gene_correlations_for_branch()` |
| SuSiE | Model settings | Fixed | `L = min(10, n_peaks)`, `intercept = FALSE`, `standardize = TRUE`, `estimate_residual_variance = TRUE`, `max_iter = 100`; credible-set coverage 0.95 | `finemap_peak_gene_correlations_for_branch()` |
| SuSiE | Records retained | Fixed | PIP ≥ 0.01 or credible-set member, otherwise the top peak | `finemap_peak_gene_correlations_for_branch()` |
| Top links | Links per cell group | Configurable |  | [`peak_gene_correlation_top_links_per_cell_group`](../parameters.html#peak_gene_correlation_top_links_per_cell_group) |


<!-- source: website/implementation/methods_differential_analyses.md -->

# Differential analyses

This chapter covers the optional `differential_analyses` module. The target structure is shown in the [differential analyses graph](implementation_differential_analyses.md), and the prerequisites and configuration in [Differential analyses](../downstream_differential_analyses.html). The tables follow the layout in [Methods and parameter tables](implementation_conventions.md#methods-and-parameter-tables).

## Cell-type composition


<!-- begin include: website/implementation/_shared_methods/cell_type_composition.md -->

The differential-analysis module treats donors, rather than nuclei, as the biological replicates and takes annotation classes from the WNN cell-type label of the final WNN metadata. For cell-type composition, nuclei are counted by donor and annotation class, completing absent combinations with zero counts, so no class is dropped for being observed in few donors; every observed class is tested unless the model restricts the response classes. Each class is modelled separately, with the numbers of nuclei in that class and in all other classes as a two-column response, in a fixed-effects beta-binomial model with logit link fitted by `glmmTMB` using the configured formula [@brooks2017_glmmtmb]; random effects and custom design functions are not supported in this branch. Named linear contrasts of the fixed effects are tested with Wald statistics, and fits that fail, do not converge or lack a positive-definite Hessian are reported as non-estimable. The pipeline does not impose a universal replication threshold, so the configured covariates must be supported by the donor count and study design. Benjamini–Hochberg correction is applied across the tested annotation classes within each model and contrast; there is no adjustment across contrasts or models.

<!-- end include: website/implementation/_shared_methods/cell_type_composition.md -->


| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Input | Annotation column counted per donor | Fixed | `WNN_harmony_SNN_cluster_cell_type` | `model_data.cell_type_composition` |
| Population | GEM wells defining the population | Configurable |  | `GEM_well_IDs` inside [`differential_analyses_cell_type_composition_models`](../parameters.html#differential_analyses_cell_type_composition_models) |
| Population | Donor restriction | Configurable |  | `donor_ids` inside [`differential_analyses_cell_type_composition_models`](../parameters.html#differential_analyses_cell_type_composition_models) |
| Response | Classes tested; denominators always use all retained nuclei | Configurable |  | `cell_types_to_test` inside [`differential_analyses_cell_type_composition_models`](../parameters.html#differential_analyses_cell_type_composition_models) |
| Model | Fixed-effects formula | Configurable |  | `formula` inside [`differential_analyses_cell_type_composition_models`](../parameters.html#differential_analyses_cell_type_composition_models) |
| Model | Family, link and response | Fixed | beta-binomial with logit link on `cbind(n_nuclei, n_other_nuclei)`; no dispersion or zero-inflation formula, so package defaults apply | `fit_cell_type_composition_model()`, `glmmTMB::glmmTMB()` |
| Contrasts | Named linear contrasts | Configurable |  | `contrast_specs_vec` inside [`differential_analyses_cell_type_composition_models`](../parameters.html#differential_analyses_cell_type_composition_models) |
| Plot | Phenotype panels and colour variable | Configurable |  | `plot_phenotype_vars`, `color_by` inside [`differential_analyses_cell_type_composition_models`](../parameters.html#differential_analyses_cell_type_composition_models) |

## Molecular pseudobulk analyses


<!-- begin include: website/implementation/_shared_methods/pseudobulk_differential_analyses.md -->

For molecular analyses, GEX and ATAC counts are summed within each donor–annotation-class combination to form pseudobulk samples, the same pseudobulks used for the Seurat/Signac export. Four feature matrices are tested:

- **Gene expression (DGE):** gene pseudobulk counts.
- **Chromatin accessibility (DCA):** consensus-peak pseudobulk counts after peak-level quality control.
- **Motif-family accessibility (DTFA):** betterChromVAR deviations of the JASPAR motif families calculated from the pseudobulk ATAC counts, column-centred and quantile-normalized across samples, so this matrix is continuous [@germain2026_betterchromvar].
- **Transcription-factor activity (DCTA):** signed CollecTRI regulator activities inferred with the `decoupleR` univariate linear model from filtered, library-size-normalized log-CPM GEX pseudobulks, requiring a minimum number of measured targets per regulator [@muller_dott2023_collectri; @badia2022_decoupler]. The CollecTRI network is accepted only when it matches a pinned checksum, and complexes such as AP1 and NF-κB remain intact as regulons.

For each model, donors missing a model variable or outside an optional donor list and samples outside an optional annotation-class subset are removed; the DTFA branch additionally removes samples below a configured minimum ATAC depth. For count matrices, zero-depth samples are removed, features are filtered with the design-aware `edgeR::filterByExpr()`, and library sizes are normalized with `edgeR::normLibSizes()` [@robinson2009_edger].

Model matrices are generated from a configured formula or custom design function, and contrasts from named linear expressions or custom contrast functions. A configured random effect selects voom for counts followed by `limma::duplicateCorrelation()` with the random effect as block [@ritchie2015_limma]. Otherwise, count designs with more than one coefficient use edgeR quasi-likelihood dispersion estimation with robust fitting and testing, and continuous matrices and single-coefficient designs use limma; a single-coefficient count design is fitted on the raw counts without the filtering and normalization above. An optional paired-cell-type route fits each annotation class separately, permits at most one pseudobulk per pairing unit and class, estimates residual correlations between classes from shared donors, and incorporates those correlations when testing cross-class contrasts.

For each model and contrast, Benjamini–Hochberg FDR is calculated across the tested features, and a fixed FDR threshold defines significant results. DGE statistics are also tested against the MSigDB Hallmark and Reactome collections using competitive `cameraPR` tests with a fixed inter-gene correlation and a minimum number of tested genes per set [@dolgalev2026_msigdbr]. The test statistic is the moderated *t*-statistic where the route provides one and otherwise a signed normal quantile of the nominal P-value; FDR is calculated within each contrast and collection. A cross-modality comparison maps CollecTRI regulators to JASPAR families and compares the model *t*-statistics of transcription-factor activity, motif-family accessibility and TF expression rather than their raw scales; family-level summaries use the median regulator *t*-statistic.

Diagnostic outputs report pseudobulk depth, retained sample and donor counts, paired-donor support, P-value distributions, effect directions and significant-feature counts. When a trait identifier is configured, the top DGE features are annotated with Open Targets evidence. A branch fails rather than returning a model when no samples or testable features remain, residual degrees of freedom are below one, a correlation block or contrast is invalid, or too few shared donors remain for a requested paired-cell-type comparison.

<!-- end include: website/implementation/_shared_methods/pseudobulk_differential_analyses.md -->


### Pseudobulk construction

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Summation | Grouping and inputs | Fixed | counts summed (`method = "sum"`) per `WNN_harmony_SNN_cluster_cell_type` and donor; the ATAC input is the peak-QC-filtered consensus-peak matrix | `get_BPCells_pseudobulk_matrix()`, `pseudobulk_counts_BPCells_matrix_dir.GEX`, `pseudobulk_counts_BPCells_matrix_dir.ATAC` |
| DTFA | Motif families | Fixed | 233 JASPAR 2026 CORE vertebrate familial root motifs; membership map of 1,019 motifs | `resources/JASPAR2026_vertebrate_familial_root_motifs.tf`, `resources/JASPAR2026_vertebrate_motif_families.tsv` |
| DTFA | Scores | Fixed | betterChromVAR analytic deviations with `compute = c("deviations", "z")` on peaks with non-zero pseudobulk counts; z-scores column-centred and quantile-normalised with `limma::normalizeBetweenArrays(method = "quantile")`, switched on by the project helper argument `normalize = TRUE` | `get_pseudobulk_chromVAR_background_record()`, `compute_pseudobulk_chromVAR_deviation_SE()`, `get_pseudobulk_chromVAR_accessibility_matrix()`, `get_pseudobulk_motif_family_accessibility_matrix()` |
| DCTA | Regulon network | Fixed | CollecTRI from `https://rescued.omnipathdb.org/CollecTRI.csv`, SHA-256 `86c90b30f2cc75c189da1f0a8c353d1547287cd656a9fac1c678634285bcb4e0`; 43,536 signed interactions from 1,189 regulators | `CollecTRI_human_network_csv`, `read_CollecTRI_human_network()` |
| DCTA | Expression preprocessing | Fixed | `edgeR::filterByExpr(group = cluster)`; TMM normalisation with `normLibSizes()`; `cpm(log = TRUE, prior.count = 2)` | `get_pseudobulk_CollecTRI_TF_activity_matrix()` |
| DCTA | Inference | Fixed | `decoupleR::run_ulm()` with at least 5 measured targets per regulator (`min_targets = 5`) | `get_pseudobulk_CollecTRI_TF_activity_matrix()` |

### Sample and feature filtering

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Donors | Variables that must be non-missing | Configurable |  | `formula` inside [`differential_analyses_pseudobulk_models`](../parameters.html#differential_analyses_pseudobulk_models) |
| Donors | Donor restriction | Configurable |  | `donor_ids` inside [`differential_analyses_pseudobulk_models`](../parameters.html#differential_analyses_pseudobulk_models) |
| Donors | Extended donor metadata table | Configurable |  | [`differential_analyses_extended_donor_id_metadata_tsv`](../parameters.html#differential_analyses_extended_donor_id_metadata_tsv) |
| Samples | Annotation-class subset | Configurable |  | `cell_type_subset` inside [`differential_analyses_pseudobulk_models`](../parameters.html#differential_analyses_pseudobulk_models) |
| Samples | Minimum ATAC depth, DTFA only | Configurable |  | [`differential_analyses_motif_family_accessibility_min_ATAC_counts`](../parameters.html#differential_analyses_motif_family_accessibility_min_ATAC_counts) |
| Count routes | Sample and feature filtering | Fixed | samples with zero counts removed; `filterByExpr(design = design_matrix)` and TMM `normLibSizes()`, otherwise package defaults | `fit_pseudobulk_feature_matrix_model()`, `fit_pseudobulk_cell_type_matrix()`, `edgeR::filterByExpr()`, `edgeR::normLibSizes()` |

### Model routes

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Design | Formula or custom design function | Configurable |  | `formula`, `design_matrix_func_name` inside [`differential_analyses_pseudobulk_models`](../parameters.html#differential_analyses_pseudobulk_models) |
| Design | Contrasts | Configurable |  | `contrast_specs_vec` and custom contrast functions inside [`differential_analyses_pseudobulk_models`](../parameters.html#differential_analyses_pseudobulk_models) |
| Route | Count detection | Fixed | integer check on up to 10 × 10 randomly sampled entries | `is_count_matrix()` |
| Route | Correlation route trigger | Configurable |  | `random_effect` inside [`differential_analyses_pseudobulk_models`](../parameters.html#differential_analyses_pseudobulk_models) |
| edgeR route | Dispersion and test | Fixed | `estimateDisp()` with package defaults; `glmQLFit(robust = TRUE)`; `glmQLFTest()` | `fit_pseudobulk_feature_matrix_model()`, `get_pseudobulk_feature_model_results()`, `edgeR::estimateDisp.DGEList()` |
| limma routes | Moderation | Fixed | `eBayes()` with package defaults, without `robust` or `trend` | `get_pseudobulk_feature_model_results()`, `get_pseudobulk_cell_type_contrast_statistics()`, `limma::eBayes()` |
| Paired route | Trigger and pairing | Configurable |  | `cell_type_formula`, `pairing_variable`, `correlation_block` inside [`differential_analyses_pseudobulk_models`](../parameters.html#differential_analyses_pseudobulk_models) |
| Paired route | Per-class count fit | Fixed | `edgeR::voomLmFit()` with `normalize.method = "none"` after the filtering and TMM normalisation above | `fit_pseudobulk_cell_type_matrix()` |
| Paired route | Residual correlation between classes | Fixed | shared donors, more than coefficients + 2 required; up to 2,000 evenly spaced common features; Fisher-z mean trimmed at 0.15 | `estimate_pseudobulk_cell_type_residual_correlations()` |
| Paired route | Cross-class test | Fixed | t from the two class estimates and their covariance; df is the smaller per-class total df | `get_pseudobulk_paired_cell_type_contrast_statistics()` |

### Multiplicity, significance and gene sets

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Significance | FDR and significant features | Fixed | BH across the tested features of each model and contrast; significant when FDR \< 0.05 and `logFC != 0` | `get_pseudobulk_feature_model_results()`, `get_pseudobulk_cell_type_model_results()`, `get_pseudobulk_differential_significant_elements_tibble()` |
| Top features | Features labelled and queried per contrast | Fixed | 40 with the smallest nominal P | `top_features_tibble` |
| Gene sets | Collections | Fixed | MSigDB Hallmark (`H`) and Reactome (`C2`, `CP:REACTOME`); human gene sets for `Homo_sapiens`, mouse for `Mus_musculus` | `gene_sets`, `get_msigdb_gene_sets()` |
| Gene sets | Test | Fixed | `limma::cameraPR(inter.gene.cor = 0.01)`; at least 10 tested genes per set; BH within contrast and collection | `get_gene_set_enrichment_results()` |
| Open Targets | Trait identifier; empty skips the query | Configurable |  | [`differential_analyses_pseudobulk_OT_GWAS_efo_id`](../parameters.html#differential_analyses_pseudobulk_OT_GWAS_efo_id) |
| Open Targets | Query | Fixed | Open Targets Platform GraphQL API (`https://api.platform.opentargets.org/api/v4/graphql`), queried at run time for the top gene-expression features only | `top_feature_open_targets_evidence_tibble`, `get_OT_GWAS_gene_evidence_tibble()` |

### Cross-modality comparison

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Crosswalk | Complex members | Fixed | AP1: FOS, FOSB, FOSL1, FOSL2, JUN, JUNB, JUND; NFKB: NFKB1, NFKB2, REL, RELA, RELB | `get_CollecTRI_JASPAR_family_map()` |
| Summary | Family statistics and concordance | Fixed | median CollecTRI and TF-expression t per family; Spearman correlation with the motif-family t when at least 3 families are mapped in a contrast | `get_CollecTRI_JASPAR_family_comparison_tibble()`, `get_CollecTRI_JASPAR_concordance_tibble()` |


<!-- source: website/implementation/methods_genetic_enrichment.md -->

# Genetic enrichment

This chapter covers the optional `genetic_enrichment` module. The target structure is shown in the [genetic enrichment graph](implementation_genetic_enrichment.md), and the configuration in [Genetic enrichment](../downstream_genetic_enrichment.html). The sparse SCAVENGE implementation is compared with its reference in [Algorithmic implementations](algorithm_validation.md#sparse-scavenge-propagation-and-significance). The tables follow the layout in [Methods and parameter tables](implementation_conventions.md#methods-and-parameter-tables).


<!-- begin include: website/implementation/_shared_methods/genetic_enrichment.md -->

The implementation pins one Open Targets Platform release and downloads its study, credible-set, credible-set evidence and target datasets [@buniello2025_open_targets_platform]; local fine-mapping results can be supplied as Parquet files satisfying a fixed schema. For each configured GWAS, multiomeR retrieves the requested study and fine-mapping result, taking the first available fine-mapping method in a fixed priority order when selection is automatic. Variants of all credible sets, without a 95% credible-set restriction, are filtered by posterior probability and mapped to consensus ATAC peaks. Variant weights that map to the same peak are summed and capped at one. These peak weights form the trait annotation from which analytic nucleus-level chromVAR z-scores are calculated, reusing the GC-bias background of the motif analysis [@schep2017_chromvar; @ulirsch2019_gchromvar; @germain2026_betterchromvar].

Annotation-class pseudobulk deviations are calculated separately rather than by averaging the nucleus-level results. ATAC counts are summed by the GEX-derived annotation carried into the final WNN metadata, and a betterChromVAR background model is fitted to the resulting peak-by-class matrix. The raw deviation is the observed-minus-background accessibility of the weighted peaks relative to their expected accessibility, and a relative deviation additionally standardizes the raw deviations across annotation classes within each trait. The analytic z-score uses the background variance of the weighted peaks; one-sided P-values and Benjamini–Hochberg-adjusted values are retained, and plot labels mark fixed unadjusted z-score thresholds. Thus, the heatmap fill, within-trait standardization and support statistic are separate quantities. When effect sizes are available, an absolute-effect branch weights variants by posterior probability times effect size. Locus-level attribution decomposes each class deviation into contributing loci and variants, reconciled against the class totals and labelled with the highest-scoring Open Targets locus-to-gene genes. Detail plots are drawn for classes with a positive deviation and a z-score at or above a configurable screen, showing the top loci by absolute contribution and by combined contribution and effect-size rank.

SCAVENGE-style trait-relevance scores are calculated from the nucleus-level z-scores on the WNN graph with a local sparse-matrix implementation of the SCAVENGE propagation strategy rather than the reference package [@yu2022_scavenge]. Nuclei whose one-sided normal-tail probability is at or below a fixed cutoff are seeds, subject to the configured maximum seed fraction; when no nucleus qualifies, all scores are zero. The nonzero support of the WNN graph is converted to binary adjacency, degree-zero nuclei are excluded, and seed signal is propagated by a random walk with the configured restart probability until convergence. Degree-matched seed permutations are sampled sequentially as in the reference implementation, while parallel native random walks stream per-cell exceedance counts without materializing the cell-by-permutation score matrix. Cell-level empirical P-values are the exceedance fraction, with a fixed cutoff defining significant cells. Scores are capped at an upper quantile, min–max scaled and multiplied by the mean z-score of the top-scoring nuclei. As pipeline extensions, cluster-level permutation medians receive add-one P-values with Benjamini–Hochberg adjustment within each metadata grouping, and each summarized grouping reports the number of nuclei, the number and proportion of significant nuclei, and the median, mean, interquartile range and range of the scores.

<!-- end include: website/implementation/_shared_methods/genetic_enrichment.md -->


## GWAS inputs and peak weights

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Studies | Study label, category, source identifier, fine-mapping method | Configurable |  | fields inside [`genetic_enrichment_GWAS_studies`](../parameters.html#genetic_enrichment_GWAS_studies) |
| Open Targets | Platform release and datasets | Fixed | release `26.03`: `credible_set`, `study`, `evidence_gwas_credible_sets`, `target` | `open_targets_credible_set_dataset_path`, `open_targets_study_dataset_path`, `open_targets_gwas_credible_sets_evidence_dataset_path`, `open_targets_target_dataset_path` |
| Open Targets | Study selection | Fixed | identifiers matching `^GCST[0-9]+$`, anything else is a local file; `studyType == "gwas"` | `classify_GWAS_source()`, `resolve_open_targets_GWAS_input_tibble()` |
| Fine-mapping | Automatic priority | Fixed | SuSie, SuSiE-inf, PICS | `resolve_open_targets_GWAS_input_tibble()` |
| Credible sets | Open Targets variants | Fixed | every variant in each credible-set `locus`, without an `is95CredibleSet` filter, on autosomes, X, Y or MT | `get_open_targets_credible_set_variants_tibble()` |
| Variants | Posterior probability cutoff | Configurable |  | [`genetic_enrichment_posterior_probability_cutoff`](../parameters.html#genetic_enrichment_posterior_probability_cutoff) |
| Variants | Cutoff comparison | Fixed | posterior probability strictly greater than the cutoff | `filter_credible_set_variants()` |
| Peak weights | Peaks and combination | Fixed | peaks of the ATAC chromVAR object; weights of variants in the same peak summed and capped at 1 (`weight_transform = "cap_1"`) | `genetic_enrichment_peak_ranges`, `get_GWAS_chromVAR_peak_weight_record()` |

## Nucleus-level deviations

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Deviations | Nucleus-level statistic | Fixed | analytic z-scores only (`compute = "z"`) | `get_GWAS_chromVAR_z_score_chunk_record()` |

## Annotation-class pseudobulk deviations

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Pseudobulks | Grouping | Fixed | ATAC counts summed (`method = "sum"`) per `PCA_harmony_SNN_cluster_cell_type` from the WNN metadata | `cell_type_pseudobulk_counts_BPCells_matrix_dir.ATAC`, `get_BPCells_group_pseudobulk_matrix()` |
| Pseudobulks | Peak filter | Fixed | peaks with zero pseudobulk counts removed before the background | `get_pseudobulk_chromVAR_background_record()` |
| Plots | Support labels | Fixed | `**` for z ≥ 2.326, `*` for z ≥ 1.645, unadjusted | `chromVAR_Z_support_labels()` |
| Plots | Compartment grouping | Configurable |  | [`genetic_enrichment_compartment_patterns`](../parameters.html#genetic_enrichment_compartment_patterns) |
| Absolute effect | Eligibility | Fixed | variant-level effects when every variant has one, otherwise locus-level effects when every locus has one; otherwise skipped | `infer_GWAS_absolute_effect_weighting()` |

## Locus attribution and detail plots

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| L2G | Gene labels | Fixed | L2G score ≥ 0.05; top 3 genes per locus | `get_open_targets_GWAS_locus_to_gene_tibble()`, `get_GWAS_chromVAR_locus_contribution_tibble()` |
| Detail plots | Minimum z | Configurable |  | [`genetic_enrichment_variant_detail_min_z`](../parameters.html#genetic_enrichment_variant_detail_min_z) |
| Detail plots | Loci | Fixed | classes with a positive deviation passing the z screen; per class the top 3 loci by absolute contribution plus the top 3 by combined contribution and effect-size percentile; 25 kb flank | `prepare_GWAS_variant_contribution_detail_records()`, `select_GWAS_detail_loci()` |
| Attribution plots | Loci shown individually | Fixed | 15 in heatmaps; 5 in bar plots | `chromVAR_locus_contribution_per_GWAS_heatmaps.cell_type_pseudobulk`, `chromVAR_locus_contribution_per_GWAS_faceted_bars_plots.cell_type_pseudobulk` |

## SCAVENGE trait-relevance propagation

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Graph | Representation | Fixed | WNN SNN graph (`WNN_harmony_SNN`) only | `graph_matrix` |
| Seeds | Selection | Fixed | one-sided normal P ≤ 0.05; nuclei with non-finite z or z \> 1000 excluded first | `get_SCAVENGE_seed_index()`, `get_SCAVENGE_result_from_chromVAR_z_score_record()` |
| Seeds | Maximum seed fraction | Configurable |  | [`genetic_enrichment_SCAVENGE_seed_percent`](../parameters.html#genetic_enrichment_SCAVENGE_seed_percent) |
| Walk | Restart probability | Configurable |  | [`genetic_enrichment_SCAVENGE_restart_prob`](../parameters.html#genetic_enrichment_SCAVENGE_restart_prob) |
| Walk | Convergence | Fixed | L1 change ≤ 1e-5; at most 10,000 iterations | `run_sparse_random_walk_with_restart()`, `run_SCAVENGE_permutation_statistics()` |
| Permutations | Count | Configurable |  | [`genetic_enrichment_SCAVENGE_permutation_times`](../parameters.html#genetic_enrichment_SCAVENGE_permutation_times) |
| Significance | Cell and cluster P-values | Fixed | cell: strictly greater permuted scores / permutations, significant at P ≤ 0.05; cluster: median score, (exceedances + 1) / (permutations + 1) counting ties, BH within each grouping | `get_SCAVENGE_result_from_chromVAR_z_score_record()`, `summarize_SCAVENGE_cluster_permutations()` |
| Score | Cap and scale | Fixed | capped at the 0.95 quantile, min–max scaled, multiplied by the mean z of the top 1 % of nuclei | `get_SCAVENGE_result_from_chromVAR_z_score_record()`, `get_SCAVENGE_scale_factor()` |
| Summaries | Groupings | Fixed | `WNN_harmony_SNN_cluster_named`, `WNN_harmony_SNN_cluster_cell_type` | `summarize_SCAVENGE_TRS_by_groups()`, `get_SCAVENGE_cluster_index_record()` |
| Plots | Heatmap stars | Fixed | on BH-adjusted cluster P: `***` ≤ 0.001, `**` ≤ 0.01, `*` ≤ 0.05 | `add_SCAVENGE_heatmap_significance()` |


# Orphaned Markdown Pages

Tracked Markdown files not reached from the Quarto book graph or include graph.

- `website/data/README.md`
- `website/figures/human_curated/README.md`
