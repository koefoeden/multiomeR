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

This chapter explains how configuration becomes target definitions: the parameter manifest supplies defaults and validation rules, mapping tables define repeated analyses, and target symbols connect their dependencies. Description tags support result selection and graph views. The final section covers project startup and resource configuration.

## Target metadata tags

multiomeR stores lightweight target metadata in the `description` argument of `targets::tar_target()` and `tarchetypes::tar_file()` calls. The descriptions should remain readable prose, with bracketed tags appended when a target needs to be discoverable from the manifest.

``` r
targets::tar_target(
  name = harmony_embeddings_matrix.GEX,
  description = "Harmony-corrected SCTransform GEX PCA embeddings [part_of_graph:GEX] [part_of_graph:WNN]",
  command = ...
)
```

The currently meaningful tag families are:

``` text
[checkpoint:<name>]             review or execution checkpoint
[part_of_graph:<graph_id>]      curated membership in an implementation graph
[resource_observation:<note>]   compact empirical resource note
```

`[checkpoint:<name>]` marks targets selectable with `targets::tar_described_as()`. The eight numbered primary-module groups are listed in `QC_checkpoint_manifest.tsv`; optional module groups remain unnumbered. UMAP parameter sweeps belong to their modality's numbered checkpoint, and compatibility objects belong to GEX checkpoint 3 or multimodal checkpoint 8. Selection matches description substrings; include the closing `]` to match a complete checkpoint tag. Dependencies still come from the target commands. [Run your own analysis](../main_running.html#steps) explains each boundary; acceptance criteria depend on the study.

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

`[resource_observation:<note>]` is currently best treated as provisional documentation. It is useful when a target has a compact empirical runtime or memory observation worth keeping near the target definition, but it is not yet a structured resource-estimation system. Keep these notes short, dated when relevant, and self-explanatory.

Use tags only when they create a durable handle for readers, graph helpers, or checkpoint commands. Ordinary internal dependencies can stay untagged.

## Parameter manifest

`cfg_pipeline_parameters.tsv` is the schema for YAML-backed pipeline configuration. Each row defines one parameter for one scope:

``` text
aggregation
differential_analyses
genetic_enrichment
peak_gene_correlation
```

For each parameter, the manifest records its name, type, cardinality, default value, missing-value rule, allowed values, example values, topic, graph/module ownership, and human description. The YAML files then only need to specify values that differ from the manifest defaults, plus values that are required because their resolved value may not be missing.

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

Configuration readers use `configuration_path()` to resolve a basename within `configuration/` or the directory named by the ignored `configuration.local`. Relative selections are anchored at the repository root. This selection does not change data-path interpretation or `_targets.yaml`. The selected GEM-well path is a graph global consumed by the file target, so changing directories also changes its dependency. Aggregation and enabled-module settings are resolved during graph construction. Disabled modules do not read their configuration.

1.  `GEM_well_tibble_all` reads only the pre-aggregation processing columns from every row in the canonical `cfg_GEM_wells.tsv`.
2.  `aggregation_tibble_all_from_yaml` is read from `cfg_aggregations.yaml`.
3.  `aggregation_tibble` keeps active aggregations, validates their GEM well references against the complete view, and adds upstream target-symbol columns.
4.  `GEM_well_tibble` keeps GEM wells whose `GEM_well_is_active` value is true.
5.  `_targets.R` expands active GEM wells and aggregations with `tar_map()`, then appends module target files. Cross-GEM-well QC summaries use the aggregation's selected wells.

Within each aggregation, `GEM_well_metadata_tibble` reads the same canonical file, subsets it to `aggregation_GEM_well_IDs`, and preserves that order. Cheap keyed projection targets then expose only the columns requested for SCT, Harmony or configured analyses. Complete non-processing annotations are joined only for explicit export objects. These projection targets are cache boundaries: a newly added or edited online column can update the canonical table without changing expensive consumers whose selected view is identical.

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

The chapters in the **Methods and parameters** part describe each analysis stage and list every setting that determines its result. They are the single place where exact values are recorded; the manuscript supplement describes the same algorithms without values. Each chapter uses one table layout:

| Column | Content |
|---|---|
| Step | The analysis step, in the order the targets run. |
| Setting | The quantity or method choice. |
| Status | `Configurable` or `Hardcoded`, as defined below. |
| Value | Shown only for hardcoded settings. |
| Source | Where the setting lives. |

A setting is **configurable** when a row in `cfg_pipeline_parameters.tsv` controls it, directly or as a field inside a nested manifest parameter such as a model specification, or when a column of `cfg_GEM_wells.tsv` controls it. Configurable rows link to the [parameter browser](../parameters.html), which renders the current default from the public manifest snapshot. They never repeat the default in prose, so a default change needs no edit outside the manifest.

A setting is **hardcoded** when no manifest row or GEM-well column controls it. Hardcoded rows show the value and the source, which is one of:

``` text
helper default   default argument of an R helper that the calling target does not override
target literal   a literal passed explicitly in a target command
inline literal   a constant inside a function body
```

Helper defaults are the easiest to expose as parameters later; inline literals require a code change. Changing any hardcoded value invalidates the affected targets on the next run, like any other code change.

The demonstration settings in the manuscript supplement are resolved values for one aggregation and are not repeated here.

## How to read the rest of the implementation book

These conventions are the connective tissue behind the graph chapters. The parameter manifest explains why config rows can be compact. Mapping tibbles explain why target names have stable suffixes. Target-symbol columns explain how mapped targets pass sets of upstream targets across graph levels. Target metadata tags explain why some nodes remain visible in curated graph views. The bootstrap contract explains why helper functions, controller resources, and target options are available before `_targets.R` is evaluated.

When modifying the implementation, preserve these contracts unless the change is explicitly meant to replace one of them.


<!-- source: website/implementation/algorithm_validation.md -->

# Algorithmic implementations, deviations and validation

multiomeR reimplements a small number of reference algorithms so they can operate on the workflow's native matrices and graph state. This page is the maintained record of what those implementations preserve, where they deliberately differ, and what the executable validation establishes.

The evidence labels used below are intentionally narrow:

- **Reference-parity tested** means the repository and named reference implementation run on the same deterministic fixture and their returned values are compared directly.
- **Reference-similarity tested** means exact equality is not an appropriate contract, so predefined similarity thresholds are checked against the named reference implementation.
- **Algorithmically derived** means the implementation is checked against an independent mathematical result, not against another software implementation.

Passing these fixtures does not validate every dataset, parameter regime, approximate-neighbor realization, biological interpretation, or downstream target. The test suite contains only such reference comparisons of the repository's reimplementations. The [CI workflow](https://github.com/koefoeden/multiomeR/blob/main/.github/workflows/algorithm-validation.yaml) runs the complete suite when tests, relevant helpers, or the Pixi environment change.

Run the complete suite with `pixi run --use-environment-activation-cache test`. The narrower `pixi run --use-environment-activation-cache test-algorithm-validation` task runs only the slow UCell, AMULET, WNN, and SCAVENGE parity and acceptance tests.

| Implementation | Evidence status | Maintained reference | Current fixed-fixture result |
|---|---|---|---|
| BPCells-native UCell | Reference-parity tested | UCell 2.14.0 | Exact values, dimensions, and dimnames |
| BPCells-native AMULET | Reference-parity tested | scDblFinder 1.24.0 | Exact metrics and multi-chromosome loci, including order |
| Native WNN | Reference-similarity tested | Seurat 5.5.0 | Production settings, `k` 20/50: weight Spearman 0.983/0.991; mean neighbor overlap 0.992/0.998 |
| Sparse SCAVENGE propagation | Algorithmically derived and reference-parity tested | SCAVENGE 1.0.2 at `8ee8b173d965` | Closed-form delta 8.61e-13; pinned-reference propagation delta 1.11e-16; exact streamed exceedance counts and significant-cell calls |
| Peak-gene donor-slope REML and Kenward-Roger kernels | Reference-parity tested | lme4 2.0.1 and pbkrtest 0.5.5 | Production scan coefficients, df and p-values within 1e-6; identical fit statuses |
| Sampled voom correlation | Reference-parity tested | edgeR 4.8.2 and limma 3.66.0 | Unsampled fits equal `voomLmFit()`; sampled consensus correlation equals `duplicateCorrelation()` on the same features |
| Peak-gene HC3 statistics and compact BH breakpoints | Reference-parity tested | sandwich 3.1.1; `stats::p.adjust()` | HC3 coefficients, errors and p-values within 1e-10 for one to six donors; identical FDR per chromosome slice |

The peak-gene and voom rows are described with their analyses in [Peak-gene correlation](methods_peak_gene_correlation.md) and [Differential analyses](methods_differential_analyses.md); their tests are `test-peak-gene-hierarchical-parity.R`, `test-peak-gene-correlation-parity.R` and `test-pseudobulk-correlation-sampling.R`.

## BPCells-native UCell scoring

**Reference algorithm.** [`UCell::ScoreSignatures_UCell()`](https://bioconductor.org/packages/release/bioc/html/UCell.html) calculates per-cell signature scores from descending feature ranks, caps ranks at `maxRank`, combines positive and negative signatures, and clips negative combined scores to zero. The maintained comparison also covers `UCell::AddModuleScore_UCell()`.

**Reason for reimplementation.** The workflow keeps gene-by-cell counts in BPCells-backed matrices. Materializing the complete matrix in memory or building a Seurat object solely for marker scoring would discard that storage contract, so multiomeR ranks bounded cell chunks and returns metadata-ready scores directly.

**Behavior preserved.** The implementation preserves UCell signature syntax (`+` and `-` suffixes), descending per-cell ranks, configurable tie handling, `maxRank` capping, impute/skip behavior for missing genes, negative-signature weighting, lower-bound clipping, signature names, and cell order.

**Deliberate deviations and consequences.** Matrix materialization is limited to one cell chunk at a time, and optional fork workers operate across chunks. This changes memory and execution behavior but not the tested score values. The helper returns a data frame instead of mutating a Seurat object. The target-level marker validator rejects configured genes missing from the Cell Ranger reference before normal pipeline scoring, whereas the lower-level helper still exposes UCell's impute/skip modes for explicit use.

**Implementation and wiring.** The scorer is [`calculate_BPCells_UCell_scores_from_matrix()` in `R/processing_GEX_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/processing_GEX_helpers.R). [`extra_targets/general_aggregation_targets.R`](https://github.com/koefoeden/multiomeR/blob/main/extra_targets/general_aggregation_targets.R) validates `UCell_GEX_marker_genes_list`. The production annotation method in `R/cluster_annotation_helpers.R` does not call that scorer; it reuses its chunked ranking helper and ranks bounded GEX count chunks once per modality's cluster partition, retaining exact sufficient statistics for positive signatures and per-cell scores for GEX metadata. For signed signatures, it scores and clips each cell before aggregation, preserving positive/negative membership in observed signatures, matched controls and marker-deletion variants. This costs more computation than the positive-only rank-summary shortcut; cell chunks bound temporary memory.

**Matched-control cluster annotation.** The annotation method built on these scores, its fixed constants and its single configurable margin are documented in [Cell-type annotation and motif accessibility](methods_annotation_and_motifs.md#matched-control-cluster-annotation). The scoring-parity test runs this production path on unsigned, signed and negative-only signatures and compares per-cell scores, observed cluster means, matched-control means and 95th percentiles, and marker-deletion excess with UCell 2.14.0 scores averaged within clusters, within 1e-12. Chunking and fork workers leave every statistic unchanged.

**Validation.** [`tests/testthat/test-scoring-parity.R`](https://github.com/koefoeden/multiomeR/blob/main/tests/testthat/test-scoring-parity.R) creates a deterministic 500-gene by 37-cell matrix, writes the project input as BPCells, and compares signed signatures with imputed and skipped missing genes. It requires `identical()` values, dimensions, and dimnames against UCell 2.14.0. It also runs the production cell-cycle scorer on the Seurat 2019 cell-cycle genes and compares it with `Seurat::CellCycleScoring()`: phases are identical, and S and G2M scores agree within 1e-6 because BPCells normalizes counts at lower floating-point precision than `NormalizeData()`.

The UCell 2.14.0 reference call with imputed missing genes can emit non-fatal R stack-imbalance warnings under the repository's R 4.5 environment. A narrow diagnostic reproduced them in the UCell reference call but not in the repository scorer. The validation therefore runs all reference calls in one disposable `callr` process and compares its returned matrix in the clean parent session; this isolates the package warning without weakening the equality assertion or changing the runtime.

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

**Validation.** [`tests/testthat/test-amulet-parity.R`](https://github.com/koefoeden/multiomeR/blob/main/tests/testthat/test-amulet-parity.R) compares against scDblFinder 1.24.0 in disposable `callr` processes. It requires `identical()` results for the bundled fragment-file metrics, the production call with prefixed Cell Ranger barcodes, a deterministic 12,000-fragment multi-chromosome loci fixture, and the corresponding full AMULET metrics. It also requires an explicit error for unsupported PCR-duplicate expansion. The disposable reference processes isolate stack-imbalance warnings emitted by the current scDblFinder reference under R 4.5 without weakening the returned-object comparison.

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

**Validation scope.** The public reference-parity fixtures do not establish equivalence for this alternative feature-aggregation path. The fixed classifier and grouping settings are listed in [Preprocessing and nucleus QC](methods_preprocessing_and_QC.md#atac-qc-and-doublets).

## Native weighted nearest neighbors

**Reference algorithm.** [`Seurat::FindMultiModalNeighbors()`](https://satijalab.org/seurat/reference/findmultimodalneighbors) constructs cell-specific modality weights from within- and cross-modality neighborhood prediction, collects candidate neighbors across modalities, and selects a weighted multimodal neighbor set.

**Reason for reimplementation.** The pipeline already has aligned RNA PCA/Harmony and ATAC LSI/Harmony matrices and needs reusable neighbor indices, distances, and modality weights without creating a Seurat object. Native graph state also feeds UMAP, Leiden clustering, SCAVENGE, and the optional Seurat/Signac export.

**Behavior preserved.** [`weighted_nearest_neighbors_BPCells()`](https://github.com/koefoeden/multiomeR/blob/main/R/processing_multimodal_helpers.R) preserves row-wise L2 normalization, per-modality nearest neighbors, nearest-distance correction, Seurat's small-SNN far-neighbour kernel bandwidth, within/cross prediction kernels, capped modality affinity ratios, normalized cell-specific modality weights, candidate-set union, weighted neighbor ranking, and Seurat's transformation from weighted affinity to neighbor distance.

**Deliberate deviations and consequences.** BPCells HNSW replaces Seurat's Annoy search, so approximate candidate sets need not be identical. The helper does not expose Seurat's optional smoothing or cross-constant list, and BPCells builds downstream SNN state rather than storing Seurat `Neighbor` and `Graph` objects. Its current `seed` argument is not consulted by the HNSW calls, so it must not be interpreted as controlling neighbor-search randomness. These choices can change weights, selected neighbors, SNN edges, clusters, and UMAP coordinates; correlation and overlap are therefore the validation contract, not exact equality.

**Implementation and wiring.** The implementation and graph consumers are in [`R/processing_multimodal_helpers.R`](https://github.com/koefoeden/multiomeR/blob/main/R/processing_multimodal_helpers.R); the tracked project-owned small-SNN kernel is in [`src/wnn_snn_bandwidth.cpp`](https://github.com/koefoeden/multiomeR/blob/main/src/wnn_snn_bandwidth.cpp). [`extra_targets/WNN_targets.R`](https://github.com/koefoeden/multiomeR/blob/main/extra_targets/WNN_targets.R) aligns modality embeddings, creates `WNN_results_raw`, filters small clusters, optionally recomputes `WNN_results`, and wires that state into WNN UMAP, clustering, metadata, and cell-type targets.

**Validation.** [`tests/testthat/test-wnn-parity.R`](https://github.com/koefoeden/multiomeR/blob/main/tests/testthat/test-wnn-parity.R) compares two deterministic RNA/ATAC fixtures with Seurat 5.5.0 using the production search settings: the helper's default HNSW `ef`, two threads and, for the 400-cell, 12-dimension fixture, candidate range 200 with `k = 20` (the manifest default) and `k = 50` (the most common configured value). Modality-weight Spearman correlations are 0.983 and 0.991 and mean neighbour-set overlaps 0.992 and 0.998; a 160-cell stress fixture with `k = 15` and candidate range 50 gives 0.934 and 0.990. One- and two-thread results are identical. Each acceptance threshold sits about 0.01 below the observed value. Before the small-SNN migration, the 26,667-cell production object gave RNA/ATAC weight Spearman 0.988/0.988 and mean neighbour overlap 0.979. The current locked BPCells 0.3.1 build is pinned at `28759cdd5125`.

**Status and rerun.** Reference-similarity tested at production settings on the fixtures above. This does not assert exact equality of selected neighbours, SNN weights, clustering, or UMAP.

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

\[ s = r\left(I - (1-r)P\right)^{-1}p_0, \]

with a maximum absolute tolerance of 1e-10; the current delta is 8.61e-13. Second, compact local reference functions reproduce the relevant SCAVENGE 1.0.2 code at the pinned commit without installing its historical dependency stack. The random-walk delta against that reference is 1.11e-16, the transformed-score delta is 3.33e-16, and all 199 fixed-RNG degree-matched seed samples, streamed per-cell exceedance counts, empirical P-values, and significant-cell calls are identical. One- and two-core native results are also identical.

**Status and rerun.** Random-walk propagation is algorithmically derived against the closed form. Seed selection, binary propagation, transformed scores, sequential permutation sampling, streamed exceedance counts, empirical P-values, and significant-cell calls are reference-parity tested against the pinned source calculation; cluster summaries are documented pipeline extensions. This does not establish parity of mutual-kNN versus pipeline graph construction, chromVAR inputs, or biological interpretation.

```bash
pixi run --use-environment-activation-cache test-algorithm-validation
```

Run the complete parity suite with:

```bash
pixi run --use-environment-activation-cache test
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

This chapter describes per-GEM-well preprocessing and the successive nucleus filters up to the final WNN cell set, and lists every setting that determines them, using the layout defined in [Methods and parameter tables](implementation_conventions.md#methods-and-parameter-tables). The target structure is shown in the [primary-module graph](implementation_main.md). GEM-well-level settings are columns of `cfg_GEM_wells.tsv`, described in [GEM well table](../reference_GEM_wells.html); aggregation-level settings are manifest parameters.

## Aggregation inputs and operational settings

These manifest parameters select inputs, plot variables and execution behaviour rather than algorithm settings. They are listed here so that every manifest row has a home in this book.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| Inputs | GEM wells combined | Configurable |  | [`aggregation_GEM_well_IDs`](../parameters.html#aggregation_GEM_well_IDs) |
| Inputs | Donor metadata table | Configurable |  | [`aggregation_donor_id_metadata_tsv`](../parameters.html#aggregation_donor_id_metadata_tsv) |
| Inputs | Aggregation active | Configurable |  | [`is_active`](../parameters.html#is_active) |
| Inputs | Optional modules enabled | Configurable |  | [`modules`](../parameters.html#modules) |
| Plots | Categorical and continuous metadata plotted | Configurable |  | [`aggregation_categorical_vars`](../parameters.html#aggregation_categorical_vars), [`aggregation_continuous_vars`](../parameters.html#aggregation_continuous_vars) |
| Plots | Additional genes plotted | Configurable |  | [`aggregation_other_interesting_genes`](../parameters.html#aggregation_other_interesting_genes) |
| Tracks | Roadmap epigenome tracks | Configurable |  | [`aggregation_roadmap_EDACC_names`](../parameters.html#aggregation_roadmap_EDACC_names) |

## Per-GEM-well inputs and metrics

Each GEM well supplies a Cell Ranger ARC count directory. The GEX matrix is imported into BPCells from the filtered feature-barcode matrix, or from a CellBender output when configured, keeping only gene-expression features. ATAC fragments are imported from the Cell Ranger fragment file. Both carry a GEM-well-prefixed barcode so nuclei stay distinct across wells. The called-cell universe is the set of barcodes that Cell Ranger flagged as cells in its per-barcode metrics. ATAC QC metrics come from BPCells, GEX metrics from the imported matrix. Optional genotype demultiplexing runs cellsnp-lite on the ATAC BAM and Vireo against the configured VCF; the resulting donor label replaces the GEM-well donor identifier, and doublet or unassigned calls are recorded as metadata only. The BPCells-native AMULET implementation, described in [Algorithmic implementations](algorithm_validation.md#bpcells-native-amulet), calculates overlap metrics and q-values on the called cells; no automatic AMULET filter is applied.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| GEX input | CellBender versus Cell Ranger matrix | Configurable |  | `GEM_well_add_cellbender`, `GEM_well_cellbender_h5_file` in `cfg_GEM_wells.tsv` |
| GEX input | Feature type kept | Hardcoded: target literal | `Gene Expression` | `GEX_counts_BPcells_matrix_dir` in `extra_targets/per_GEM_well_targets.R` |
| Barcodes | Prefix on GEX and fragment barcodes | Hardcoded: target literal | `<GEM_well_ID>_` | `extra_targets/per_GEM_well_targets.R` |
| Barcodes | Called-cell universe | Hardcoded: target literal | rows of `per_barcode_metrics.csv` with `is_cell == 1` | `cellranger_kept_metadata_tibble` in `extra_targets/per_GEM_well_targets.R` |
| ATAC metrics | QC function and blacklist | Hardcoded: target literal | `BPCells::qc_scATAC()` with Ensembl genes and an empty blacklist | `ATAC_qc_metrics_tibble` in `extra_targets/per_GEM_well_targets.R` |
| ATAC metrics | TSS enrichment window | Hardcoded: environment pin | BPCells 0.3.1: 101 bp centre window, 100 bp flanks at ±1.9–2 kb, denominator floor 0.1 | BPCells `qc_scATAC()` |
| ATAC metrics | Nucleosome signal | Hardcoded: target literal | mono-nucleosomal / sub-nucleosomal fragment counts | `extra_targets/per_GEM_well_targets.R` |
| GEX metrics | Mitochondrial gene pattern | Hardcoded: target literal | `(?i)^MT-` | `extra_targets/per_GEM_well_targets.R` |
| Demultiplexing | VCF, donor count and donor label | Configurable |  | `GEM_well_donors_VCF_file`, `GEM_well_n_donors`, `GEM_well_donor_id` in `cfg_GEM_wells.tsv` |
| Demultiplexing | cellsnp-lite settings | Hardcoded: helper default | `--minMAF 0.1`, `--minCOUNT 20`, `--UMItag None`, ATAC BAM | `R/parallel_GEM_well_preprocessing_helpers.R` |
| Demultiplexing | Vireo settings | Hardcoded: helper default | genotypes from the VCF (`-t GT`), no genotype learning | `R/parallel_GEM_well_preprocessing_helpers.R` |
| Demultiplexing | Cores | Hardcoded: target literal | cellsnp 6, Vireo 4 | `extra_targets/per_GEM_well_targets.R` |
| Demultiplexing | Use of doublet and unassigned calls | Hardcoded: inline literal | metadata only; removable through the GEM-well exclusion list | `extra_targets/per_GEM_well_targets.R` |
| AMULET | Barcodes scored | Hardcoded: target literal | Cell Ranger called cells | `amulet_metrics_tibble` in `extra_targets/per_GEM_well_targets.R` |
| AMULET | Maximum fragment size | Hardcoded: helper default | 1,000 bp | `calculate_amulet_metrics_BPCells()` in `R/amulet_BPCells_helpers.R` |
| AMULET | Excluded regions | Hardcoded: helper default | chrM, chrX, chrY and their aliases | `R/amulet_BPCells_helpers.R` |
| AMULET | Cell Ranger end-inclusive shift | Hardcoded: target literal | end − 1 | `extra_targets/per_GEM_well_targets.R` |
| AMULET | High-overlap-site removal | Hardcoded: helper default | on; Poisson P \< 0.01 | `remove_high_overlap_amulet_loci()` in `R/amulet_BPCells_helpers.R` |
| AMULET | Per-cell test | Hardcoded: inline literal | upper-tail Poisson on loci covered by more than two fragments; BH q-values | `R/amulet_BPCells_helpers.R` |
| AMULET | Nuclei scored | Hardcoded: target literal | Cell Ranger-called barcodes; no fragment minimum | `extra_targets/per_GEM_well_targets.R` |
| AMULET | Automatic filter | Hardcoded: inline literal | none; q-values usable in the GEM-well exclusion list | `R/processing_and_aggregation_constants.R` |

## GEM-well exclusions and the GEX cell universe

Configured GEM-well exclusions are dplyr filter expressions evaluated on the called cells of each well; matching nuclei are removed from the union of called cells before any aggregation-level analysis. The per-well GEX matrices are then column-bound without further filtering. At PCA, the metadata are intersected with the matrix barcodes, so only barcodes present in the aligned count data continue.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| Exclusions | Per-well exclusion expressions | Configurable |  | `GEM_well_QC_exclude_list` in `cfg_GEM_wells.tsv` |
| Exclusions | Where applied | Hardcoded: target literal | union of called cells, before GEX combination | `GEX_cellranger_kept_metadata_tibble` in `extra_targets/general_aggregation_targets.R` |
| Combination | Matrix merge | Hardcoded: target literal | column bind of every per-well matrix, no gene or cell filter | `extra_targets/GEX_merge_and_dim_reduc_targets.R` |
| Retention | Barcodes carried into PCA | Hardcoded: inline literal | intersection of metadata and matrix barcodes | `run_GEX_PCA_BPCells()` in `R/processing_GEX_helpers.R` |

## Cluster-size filter and GEX doublets

After GEX graph construction and Leiden clustering, clusters below the configured minimum size are removed; the same threshold is applied to the ATAC and WNN clusters and to the optional subgroup clusterings. The filter keeps clusters at or above the threshold and is disabled when the value is missing or at most one. scDblFinder is then run per GEM well on the raw GEX counts of the retained nuclei. Its cluster labels are the annotation-derived scDblFinder groups: the cell-type label when the cluster was assigned, otherwise the status and cluster identifier, optionally collapsed with the configured map. Doublet removal is controlled separately at the nucleus level and at the cluster level; the cluster-level rule drops a Leiden cluster whose doublet fraction exceeds the threshold.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| Cluster filter | Minimum cluster size, applied to GEX, ATAC, WNN and subgroups | Configurable |  | [`aggregation_cluster_min_barcodes`](../parameters.html#aggregation_cluster_min_barcodes) |
| Cluster filter | Comparison and disabling rule | Hardcoded: inline literal | keep size ≥ threshold; disabled when NULL, NA or ≤ 1 | `cluster_embedding_matrix_BPCells()` in `R/processing_ATAC_helpers.R` |
| GEX scDblFinder | Input | Hardcoded: target literal | per-well slice of the raw aggregated GEX counts | `extra_targets/GEX_graph_and_cluster_targets.R` |
| GEX scDblFinder | Cluster labels | Hardcoded: inline literal | annotation-derived scDblFinder groups | `R/cluster_annotation_helpers.R`, `R/processing_GEX_helpers.R` |
| GEX scDblFinder | Label collapse map, reused for ATAC | Configurable |  | [`aggregation_scDblFinder_GEX_cell_type_collapse_list`](../parameters.html#aggregation_scDblFinder_GEX_cell_type_collapse_list) |
| GEX scDblFinder | `dbr.sd` | Hardcoded: target literal | 1.0 | `extra_targets/GEX_graph_and_cluster_targets.R` |
| GEX scDblFinder | Return type and parallelism | Hardcoded: inline literal | scores; `BiocParallel::SerialParam()` | `R/processing_GEX_helpers.R` |
| GEX scDblFinder | Other arguments | Hardcoded: environment pin | scDblFinder 1.24.0 defaults | `extra_targets/GEX_graph_and_cluster_targets.R` |
| GEX doublets | Nucleus-level removal | Configurable |  | [`aggregation_scDblFinder_GEX_remove_called_doublets`](../parameters.html#aggregation_scDblFinder_GEX_remove_called_doublets) |
| GEX doublets | Cluster-level removal | Configurable |  | [`aggregation_scDblFinder_GEX_max_doublet_fraction_per_cluster`](../parameters.html#aggregation_scDblFinder_GEX_max_doublet_fraction_per_cluster) |
| GEX doublets | Cluster fraction rule | Hardcoded: inline literal | fraction over all calls on the raw Leiden cluster; strictly greater than the threshold | `R/processing_GEX_helpers.R` |

## ATAC QC and doublets

The ATAC branch starts from the GEX-retained nuclei: peak-calling groups, the consensus peak matrix and peak-level QC metrics are all built on that cell set. Peak-level metrics are evaluated against the configured aggregation-level exclusions, ATAC clusters below the minimum size are removed, and a separate scDblFinder run uses the feature-aggregation adaptation described in [Algorithmic implementations](algorithm_validation.md#bpcells-backed-atac-scdblfinder-feature-aggregation): peaks are grouped by k-means on the aggregation-wide LSI loadings, the groups are summed with BPCells, and the compact matrix is passed to the unchanged classifier.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| Cell set | Nuclei entering the ATAC branch | Hardcoded: target literal | post-doublet-filter GEX metadata | `extra_targets/ATAC_targets.R` |
| Peak QC | Metrics | Hardcoded: inline literal | ATAC counts in peaks, blacklist counts and fraction, peak-count fraction of fragments, peak-count enrichment | `R/processing_ATAC_helpers.R` |
| Peak QC | Exclusion expressions | Configurable |  | [`aggregation_QC_exclude_list_combined_object`](../parameters.html#aggregation_QC_exclude_list_combined_object) |
| ATAC scDblFinder | Nuclei scored | Hardcoded: target literal | post-cluster-filter ATAC cells passing peak QC | `extra_targets/ATAC_targets.R` |
| ATAC scDblFinder | LSI dimensions for feature groups | Hardcoded: target literal | `intersect(2:20, aggregation_ATAC_data_PCs)` | `extra_targets/ATAC_targets.R` |
| ATAC scDblFinder | Number of feature groups and seed | Hardcoded: target literal | 50 groups; seed 1 | `extra_targets/ATAC_targets.R` |
| ATAC scDblFinder | Grouping method | Hardcoded: inline literal | `stats::kmeans(iter.max = 50, nstart = 1)` on loadings; sums via `BPCells::pseudobulk_matrix(method = "sum")` | `get_feature_groups_from_LSI_loadings()`, `aggregate_BPCells_rows_by_group()` in `R/processing_GEX_helpers.R` |
| ATAC scDblFinder | Classifier arguments | Hardcoded: target literal | `dbr.sd = 1.0`, `aggregateFeatures = FALSE`, `nfeatures = 50`, `processing = "normFeatures"` | `extra_targets/ATAC_targets.R` |
| ATAC scDblFinder | Cluster labels | Hardcoded: inline literal | ATAC annotation-derived groups collapsed with the GEX map | `extra_targets/ATAC_targets.R` |
| ATAC doublets | Nucleus-level removal | Configurable |  | [`aggregation_scDblFinder_ATAC_remove_called_doublets`](../parameters.html#aggregation_scDblFinder_ATAC_remove_called_doublets) |
| ATAC doublets | Cluster-level removal | Configurable |  | [`aggregation_scDblFinder_ATAC_max_doublet_fraction_per_cluster`](../parameters.html#aggregation_scDblFinder_ATAC_max_doublet_fraction_per_cluster) |

## WNN cell set

WNN integration uses the nuclei retained by the ATAC branch that have rows in both corrected embeddings. The joint clusters are filtered with the minimum-cluster-size rule; when nuclei are dropped, the WNN graph and UMAP are recomputed on the retained nuclei while the cluster labels from the first run are kept. The GEX, ATAC and WNN metadata therefore describe successive cell universes.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| Cell set | Nuclei entering WNN | Hardcoded: target literal | post-doublet-filter ATAC metadata ∩ rows of both embeddings | `extra_targets/WNN_targets.R`, `R/processing_multimodal_helpers.R` |
| Cluster filter | Threshold and recomputation | Configurable |  | [`aggregation_cluster_min_barcodes`](../parameters.html#aggregation_cluster_min_barcodes) |
| Cluster filter | Recompute rule | Hardcoded: target literal | graph and UMAP recomputed when any nucleus is dropped; labels kept | `extra_targets/WNN_targets.R` |


<!-- source: website/implementation/methods_GEX_ATAC_and_WNN.md -->

# GEX, ATAC, batch correction and WNN

This chapter describes normalization and dimensional reduction of both modalities, peak definition, batch correction, weighted nearest-neighbour (WNN) integration, and the shared graph, clustering and UMAP steps. It lists every setting that determines them, using the layout defined in [Methods and parameter tables](implementation_conventions.md#methods-and-parameter-tables). The target structure is shown in the [primary-module graph](implementation_main.md). Library versions are pinned by the Pixi environment: BPCells 0.3.1, igraph 2.3.0, harmony 2.0.2, uwot 0.2.4 and Seurat 5.5.0 at the time of writing.

## GEX normalization and PCA

The combined GEX matrix keeps every gene. At PCA, genes whose total count over the Cell Ranger-kept cells is at or below the minimum are excluded, for both backends. The backend is configurable. The BPCells-native branch computes Pearson residuals with a per-gene method-of-moments theta, clips them, optionally regresses configured cell-level covariates, keeps the genes with the highest residual variance and takes a truncated SVD. The Seurat branch runs SCTransform v2 on the same cells and genes, keeps the same number of variable features and computes the PCA from the dense residual matrix. Regression of the cell-cycle difference is the manifest default; when a cell-cycle column is regressed, S and G2M scores are computed from log-normalized counts with the Seurat 2019 gene sets.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| Gene filter | Minimum total count per gene, both backends | Hardcoded: helper default | `min_feature_count = 50`, genes with more than 50 counts kept | `run_GEX_PCA_BPCells()` in `R/processing_GEX_helpers.R` |
| Backend | Normalization and PCA backend | Configurable |  | [`aggregation_GEX_PCA_backend`](../parameters.html#aggregation_GEX_PCA_backend) |
| Regression | Cell-level covariates regressed | Configurable |  | [`aggregation_SCT_regress_vars`](../parameters.html#aggregation_SCT_regress_vars) |
| Components | Number of PCs computed | Configurable |  | last element of [`aggregation_GEX_data_PCs`](../parameters.html#aggregation_GEX_data_PCs) |
| Variable genes | Genes retained by residual variance | Hardcoded: helper default | `n_variable_features = 3000` | `R/processing_GEX_helpers.R` |
| BPCells branch | Residual clip range | Hardcoded: helper default | `clip_range = c(-10, 10)` | `R/processing_GEX_helpers.R` |
| BPCells branch | Minimum variance passed to `sctransform_pearson` | Hardcoded: helper default | `min_var = 0` | `R/processing_GEX_helpers.R` |
| BPCells branch | Theta estimate | Hardcoded: inline literal | method of moments, clamped to \[1e-6, 1e6\] | `R/processing_GEX_helpers.R` |
| BPCells branch | Regression | Hardcoded: inline literal | `BPCells::regress_out(prediction_axis = "row")` on residuals | `R/processing_GEX_helpers.R` |
| BPCells branch | SVD | Hardcoded: inline literal | `BPCells::svds()`, embeddings = right singular vectors × singular values, no centring | `R/processing_GEX_helpers.R` |
| Seurat branch | SCTransform arguments | Hardcoded: inline literal | `conserve.memory = TRUE`, `do.correct.umi = FALSE`, otherwise Seurat 5.5.0 defaults (v2, 5,000 model cells) | `run_Seurat_SCT_for_PCA()` in `R/processing_GEX_helpers.R` |
| Seurat branch | PCA | Hardcoded: inline literal | eigendecomposition of the residual gram matrix | `R/processing_GEX_helpers.R` |
| Cell cycle | Gene sets and scoring | Hardcoded: inline literal | `Seurat::cc.genes.updated.2019`, both organisms; log-normalization with scale factor 10,000; `CC.Difference = S − G2M`; seed 1 | `R/processing_GEX_helpers.R` |

## ATAC peak calling and consensus peaks

Fragments are restricted to the standard chromosomes. Peak-calling groups are defined by a configurable metadata column of the GEX-retained nuclei; groups above the discovery cap are randomly downsampled and fragments are exported per group. MACS3 is run with ATAC-style shift and extension and summit calling, or the BPCells tile caller is used. Summits are extended to fixed-width peaks, peaks overlapping the reference blacklist are removed, and an ArchR-style iterative overlap removal produces one non-overlapping consensus set: within a group, overlapping peaks are ranked by summit significance; across groups, by MACS3 fold enrichment. The peak matrix counts insertions or fragment overlaps, as configured, for the GEX-retained nuclei.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| Fragments | Chromosomes retained | Hardcoded: target literal | autosomes, X and Y | `extra_targets/ATAC_targets.R` |
| Groups | Peak-calling grouping column | Configurable |  | [`aggregation_call_peaks_by_cluster_col`](../parameters.html#aggregation_call_peaks_by_cluster_col) |
| Groups | Discovery cap per group | Hardcoded: helper default | `max_cells_per_cluster = 50000`, seed 1 + group index | `R/processing_ATAC_helpers.R` |
| Caller | Peak-calling method | Configurable |  | [`aggregation_ATAC_peak_calling_method`](../parameters.html#aggregation_ATAC_peak_calling_method) |
| MACS3 | Command-line arguments | Hardcoded: inline literal | `-f BED --nomodel --shift -75 --extsize 150 --call-summits --keep-dup all`; default q = 0.05 | `R/processing_ATAC_helpers.R` |
| MACS3 | Effective genome size | Hardcoded: inline literal | GRCh38 2.913e9; mm10 and GRCm39 2.65e9 | `R/processing_ATAC_helpers.R` |
| Tile caller | BPCells settings | Hardcoded: helper default | `peak_width = 500`, `peak_tiling = 3`, `fdr_cutoff = 0.01`, `merge_peaks = "none"` | `R/processing_ATAC_helpers.R` |
| Peak shape | Summit extension | Hardcoded: helper default | `extend_summits = 250`, giving 500 bp peaks | `R/processing_ATAC_helpers.R` |
| Blacklist | Source per genome | Hardcoded: inline literal | GRCh38 Kundaje unified; mm10 Boyle v2; `resources/mm39.excluderanges.bed` | `R/processing_ATAC_helpers.R` |
| Blacklist | Rule | Hardcoded: inline literal | any overlap removes the peak | `R/processing_ATAC_helpers.R` |
| Consensus | Within-group ranking | Hardcoded: target literal | `neg_log10pvalue_summit`, decreasing | `extra_targets/ATAC_targets.R` |
| Consensus | Cross-group ranking | Hardcoded: target literal | `fold_change`, decreasing | `extra_targets/ATAC_targets.R` |
| Consensus | Overlap removal | Hardcoded: inline literal | iterative: reduce, keep best per cluster, drop overlaps, repeat | `R/processing_ATAC_helpers.R` |
| Peak matrix | Counting mode, also used for blacklist QC counts | Configurable |  | [`aggregation_ATAC_peak_matrix_mode`](../parameters.html#aggregation_ATAC_peak_matrix_mode) |
| Peak matrix | Cells | Hardcoded: target literal | post-doublet-filter GEX metadata barcodes | `extra_targets/ATAC_targets.R` |

## ATAC TF-IDF and LSI

The QC-filtered peak matrix is transformed with Signac's TF-IDF method 1 and decomposed with BPCells SVD to obtain latent semantic indexing (LSI) embeddings. No code rule drops the first component; the manifest defaults for the data and UMAP dimensions start at the second component.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| TF-IDF | Method and scale factor | Hardcoded: helper default | term frequency × inverse document frequency, `log1p(scale_factor × TF·IDF)`, `scale_factor = 10000` | `R/processing_ATAC_helpers.R` |
| SVD | Components computed | Configurable |  | last element of [`aggregation_ATAC_data_PCs`](../parameters.html#aggregation_ATAC_data_PCs) |
| SVD | Function | Hardcoded: inline literal | `BPCells::svds()`, embeddings = right singular vectors × singular values | `R/processing_ATAC_helpers.R` |
| Dimensions | Components used downstream | Configurable |  | [`aggregation_ATAC_data_PCs`](../parameters.html#aggregation_ATAC_data_PCs), [`aggregation_UMAP_ATAC_PCs`](../parameters.html#aggregation_UMAP_ATAC_PCs) |

## Harmony batch correction

Harmony can be applied separately to the selected GEX PCs and ATAC LSI dimensions. When no covariates are configured, the uncorrected embedding is returned unchanged. Several covariates are collapsed into one interaction batch factor, and nuclei with missing covariate values are dropped with a warning.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| Covariates | Shared correction columns | Configurable |  | [`aggregation_harmony_correction_metadata_col_names`](../parameters.html#aggregation_harmony_correction_metadata_col_names) |
| Covariates | Additional ATAC covariates | Configurable |  | [`aggregation_extra_harmony_covars_ATAC`](../parameters.html#aggregation_extra_harmony_covars_ATAC) |
| Dimensions | Corrected dimensions | Configurable |  | [`aggregation_GEX_data_PCs`](../parameters.html#aggregation_GEX_data_PCs), [`aggregation_ATAC_data_PCs`](../parameters.html#aggregation_ATAC_data_PCs) |
| Function | Call | Hardcoded: inline literal | `harmony::RunHarmony(max_iter = 25, lambda = 1)`, other arguments harmony 2.0.2 defaults | `R/processing_ATAC_helpers.R` |
| Covariates | Combination rule | Hardcoded: inline literal | interaction of all covariates as one batch factor | `R/processing_ATAC_helpers.R` |
| Covariates | Missing values | Hardcoded: inline literal | nuclei removed with a warning | `R/processing_ATAC_helpers.R` |
| Resources | Cores | Hardcoded: target literal | 6 | `extra_targets/GEX_merge_and_dim_reduc_targets.R`, `extra_targets/ATAC_targets.R` |

## Graph construction, Leiden clustering and UMAP

For GEX and ATAC, approximate nearest neighbours are found with BPCells HNSW on the selected dimensions, converted to a shared-nearest-neighbour (SNN) graph with Jaccard weights and clustered with Leiden. Because the query cell is its own first neighbour, a neighbour count of k yields k − 1 non-self neighbours. Clusters are renumbered by size. UMAP is computed with uwot on the selected dimensions. QC parameter sweeps additionally render UMAPs over grids of dimensions and neighbour counts.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| kNN | Neighbour count | Configurable |  | [`aggregation_data_nNNs`](../parameters.html#aggregation_data_nNNs) |
| kNN | Function, metric and search effort | Hardcoded: inline literal | `BPCells::knn_hnsw(metric = "cosine", ef = 500)` | `R/processing_ATAC_helpers.R` |
| SNN | Construction and pruning | Hardcoded: environment pin | `BPCells::knn_to_snn_graph()` defaults: Jaccard weights, `min_val = 1/15`, no self loops | `R/processing_ATAC_helpers.R` |
| Leiden | Resolution | Configurable |  | [`aggregation_GEX_cluster_res`](../parameters.html#aggregation_GEX_cluster_res), [`aggregation_ATAC_cluster_res`](../parameters.html#aggregation_ATAC_cluster_res), [`aggregation_WNN_cluster_res`](../parameters.html#aggregation_WNN_cluster_res) |
| Leiden | Objective, weights and seed | Hardcoded: inline literal | `igraph::cluster_leiden(objective_function = "modularity")` with SNN weights; seed 1 | `R/processing_ATAC_helpers.R` |
| Leiden | Iterations and randomness | Hardcoded: environment pin | igraph 2.3.0 defaults: 2 iterations, beta 0.01 | `R/processing_ATAC_helpers.R` |
| Leiden | Relabelling | Hardcoded: inline literal | renumbered by size, largest first | `R/processing_ATAC_helpers.R` |
| UMAP | Dimensions | Configurable |  | [`aggregation_UMAP_GEX_PCs`](../parameters.html#aggregation_UMAP_GEX_PCs), [`aggregation_UMAP_ATAC_PCs`](../parameters.html#aggregation_UMAP_ATAC_PCs) |
| UMAP | Neighbours and minimum distance | Configurable |  | [`aggregation_UMAP_nNNs`](../parameters.html#aggregation_UMAP_nNNs), [`aggregation_UMAP_min_dist`](../parameters.html#aggregation_UMAP_min_dist) |
| UMAP | Function and fixed arguments | Hardcoded: inline literal | `uwot::umap(metric = "cosine", n_components = 2, n_sgd_threads = 0)`, seed 1 | `R/processing_ATAC_helpers.R` |
| UMAP | Other arguments | Hardcoded: environment pin | uwot 0.2.4 defaults, spectral initialisation | `R/processing_ATAC_helpers.R` |
| UMAP sweeps | Grids | Hardcoded: target literal | dimensions from 5 to the configured count in 3 steps; neighbours from 10 to the configured UMAP count in 3 steps | `extra_targets/GEX_graph_and_cluster_targets.R`, `extra_targets/ATAC_targets.R` |
| Resources | Threads | Hardcoded: target literal | 6 | `extra_targets/GEX_graph_and_cluster_targets.R`, `extra_targets/ATAC_targets.R` |

## Weighted nearest neighbours

The BPCells-native WNN implementation is described, with its deviations from Seurat and its validation, in [Algorithmic implementations](algorithm_validation.md#native-weighted-nearest-neighbors). The selected GEX and ATAC dimensions after optional Harmony correction are aligned by barcode and L2-normalized. HNSW finds a large candidate set per modality; Seurat's small-SNN bandwidth strategy sets the kernel width per cell from the configured neighbour count; per-cell modality weights come from within- versus cross-modality prediction kernels; and the union of candidates is ranked by the weighted kernel score to select the final neighbours. BPCells builds the SNN graph, which is clustered with Leiden and embedded with UMAP using the precomputed neighbours.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| Inputs | Embeddings and dimensions | Configurable |  | [`aggregation_GEX_data_PCs`](../parameters.html#aggregation_GEX_data_PCs), [`aggregation_ATAC_data_PCs`](../parameters.html#aggregation_ATAC_data_PCs) after Harmony |
| Inputs | Normalization | Hardcoded: helper default | row-wise L2 | `weighted_nearest_neighbors_BPCells()` in `R/processing_multimodal_helpers.R` |
| Candidates | Candidates per modality | Hardcoded: target literal | `candidate_k = 200` | `extra_targets/WNN_targets.R` |
| Candidates | Per-modality search | Hardcoded: inline literal | `BPCells::knn_hnsw(k = candidate_k + 1, metric = "euclidean", ef = 500)` | `R/processing_multimodal_helpers.R` |
| Final neighbours | Neighbour count, also the bandwidth and imputation neighbourhood | Configurable |  | [`aggregation_data_nNNs`](../parameters.html#aggregation_data_nNNs) |
| Bandwidth | Small-SNN kernel width | Hardcoded: inline literal | mean distance to the k lowest-shared-neighbour cells after subtracting the nearest non-self distance; `sd_scale = 1`; floor at machine epsilon | `src/wnn_snn_bandwidth.cpp`, `R/processing_multimodal_helpers.R` |
| Weights | Modality weight kernel | Hardcoded: inline literal | `exp(−d/σ)`; ratio `within / (cross + 1e-4)` clipped to \[0, 200\]; softmax-normalized | `R/processing_multimodal_helpers.R` |
| Selection | Weighted score and distance | Hardcoded: inline literal | `Σ w_m · exp(−(d_m/σ_m))`, `kernel_power = 1`; `nn_dist = sqrt((1 − score)/2)` | `R/processing_multimodal_helpers.R` |
| SNN and Leiden | Graph and clustering | Hardcoded: environment pin | `knn_to_snn_graph(min_val = 1/15)`; Leiden modularity, seed 1 | `R/processing_multimodal_helpers.R` |
| UMAP | Neighbours and minimum distance | Configurable |  | [`aggregation_UMAP_nNNs`](../parameters.html#aggregation_UMAP_nNNs) capped at the neighbour count; [`aggregation_UMAP_min_dist`](../parameters.html#aggregation_UMAP_min_dist) |
| UMAP | Fixed arguments | Hardcoded: inline literal | precomputed neighbours, 2 components, seed 1 | `R/processing_multimodal_helpers.R` |
| Resources | Threads | Hardcoded: target literal | 6 | `extra_targets/WNN_targets.R` |


<!-- source: website/implementation/methods_annotation_and_motifs.md -->

# Cell-type annotation and motif accessibility

This chapter describes marker-signature cluster annotation and motif-family accessibility, and lists every setting that determines them, using the layout defined in [Methods and parameter tables](implementation_conventions.md#methods-and-parameter-tables). The BPCells-native UCell scorer and its validation are described in [Algorithmic implementations](algorithm_validation.md#bpcells-native-ucell-scoring).

## Signature scoring

Annotation is driven by the configured GEX marker signatures. Signatures accept unsigned genes, `+` suffixes for positive markers and `-` suffixes for genes expected to be absent. Genes missing from the reference fail validation before any scoring; the production path never imputes missing genes. Scores are UCell-style capped-rank statistics computed on bounded chunks of the raw GEX counts, with signed signatures clipped at zero per cell before averaging. The same scoring evidence is computed for the GEX, ATAC and WNN cluster partitions, reusing the GEX control reference.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| Signatures | Marker genes per label | Configurable |  | [`aggregation_GEX_marker_genes`](../parameters.html#aggregation_GEX_marker_genes) |
| Signatures | Panel constraints | Hardcoded: inline literal | at least 2 labels; no signature longer than the rank cap | `R/cluster_annotation_helpers.R` |
| Signatures | Missing genes | Hardcoded: target literal | validation error before scoring | `UCell_GEX_marker_genes_list` in `extra_targets/general_aggregation_targets.R` |
| Ranks | Rank cap | Hardcoded: inline literal | `min(1500, n_genes)` | `prepare_cluster_UCell_controls()` in `R/cluster_annotation_helpers.R` |
| Ranks | Direction and ties | Hardcoded: helper default | descending counts, `ties.method = "average"` | `rank_UCell_count_chunk()` in `R/processing_GEX_helpers.R` |
| Ranks | Cell chunk size | Hardcoded: helper default | 250 cells | `summarize_cluster_UCell_counts()` in `R/cluster_annotation_helpers.R` |
| Ranks | Fork workers | Hardcoded: target literal | 2 | `extra_targets/GEX_graph_and_cluster_targets.R`, `extra_targets/ATAC_targets.R`, `extra_targets/WNN_targets.R` |
| Scores | Lower-bound clipping of signed scores | Hardcoded: inline literal | `pmax(0, score)` per cell | `R/cluster_annotation_helpers.R` |
| Scores | Cluster columns annotated | Hardcoded: target literal | `PCA_harmony_SNN_cluster`, `LSI_harmony_SNN_cluster`, `WNN_harmony_SNN_cluster` | the three target files above |
| Scores | Per-cell scores retained | Hardcoded: target literal | GEX only | `extra_targets/GEX_graph_and_cluster_targets.R` |

## Matched-control cluster annotation

Each label is compared with random control signatures matched on gene abundance and detection. The control reference samples cells per GEM well from the pre-doublet-filter GEX metadata. For each replicate, markers are visited in random order and each is replaced by one gene drawn from its nearest eligible candidates not yet used in that replicate; candidates exclude all marker genes and undetected genes. The adjusted score of a label in a cluster is its observed mean score minus the upper quantile of its matched controls. The highest adjusted score nominates the candidate; its advantage is the smaller of its lead over zero and its lead over the runner-up. Assignment requires a positive best score, no exact tie and an advantage at least the configured margin; otherwise the cluster is `Unassigned` with the candidate and reason retained. Raising the margin can only withdraw assignments.

Diagnostics never veto an assignment: marker-deletion blocks remove one marker and its control from each signature and record whether the same candidate would still be assigned; detection counts report positive markers detected in a minimum fraction of cells; GEM-well agreement compares subgroups of sufficient size. GEX-module dot plots use the pre-doublet-filter evidence, order marker sets by Ward clustering of the adjusted profiles and colour by the cached adjusted scores.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| Reference | Cells per GEM well | Hardcoded: helper default | 50 | `R/cluster_annotation_helpers.R` |
| Reference | Cell set and matrix | Hardcoded: target literal | pre-doublet `metadata_w_clusters_tibble.GEX`; full aggregated GEX counts | `cluster_UCell_controls.GEX` in `extra_targets/GEX_graph_and_cluster_targets.R` |
| Controls | Random mappings | Hardcoded: helper default | 999 | `build_UCell_controls()` in `R/cluster_annotation_helpers.R` |
| Controls | Matching coordinates | Hardcoded: inline literal | `log1p(abundance × 1e4)` and `asin(sqrt(detection))`, standardized, Euclidean distance | `R/cluster_annotation_helpers.R` |
| Controls | Candidate pool and draw | Hardcoded: helper default | pool of 200 nearest; draw from the 50 nearest unused | `R/cluster_annotation_helpers.R` |
| Controls | Exclusions | Hardcoded: inline literal | all marker genes; genes with zero detection | `R/cluster_annotation_helpers.R` |
| Controls | Random seed | Hardcoded: helper default | 20260910 | `R/cluster_annotation_helpers.R` |
| Adjusted score | Background quantile | Hardcoded: inline literal | 0.95 | `R/cluster_annotation_helpers.R` |
| Assignment | Advantage rule | Hardcoded: inline literal | `min(best, best − second)`; requires best \> 0 and no tie | `R/cluster_annotation_helpers.R` |
| Assignment | Minimum advantage | Configurable |  | [`aggregation_cluster_annotation_min_advantage`](../parameters.html#aggregation_cluster_annotation_min_advantage) |
| Diagnostics | Marker-deletion blocks | Hardcoded: helper default | 10, stratified by cluster and GEM well | `R/cluster_annotation_helpers.R` |
| Diagnostics | Minimum assessable deletions | Hardcoded: inline literal | 2 | `R/cluster_annotation_helpers.R` |
| Diagnostics | Marker detection fraction | Hardcoded: inline literal | 0.10 | `R/cluster_annotation_helpers.R` |
| Diagnostics | GEM-well subgroup size and count | Hardcoded: inline literal | at least 25 cells; at least 2 wells | `R/cluster_annotation_helpers.R` |
| Plots | Dot-plot cell set | Hardcoded: target literal | GEX: pre-doublet-filter; ATAC and WNN: post-filter | `extra_targets/GEX_graph_and_cluster_targets.R`, `extra_targets/ATAC_targets.R`, `extra_targets/WNN_targets.R` |
| Plots | Marker-set ordering | Hardcoded: inline literal | Euclidean distance of adjusted profiles, `hclust(method = "ward.D2")` | `R/cluster_annotation_helpers.R` |

## Motif families and motif accessibility

Motif families are the sequence-similarity clusters of the JASPAR 2026 CORE vertebrate collection. The pipeline vendors the 233 familial root motifs and the family membership table under `resources/` and scans the root motifs directly against the consensus peaks with `motifmatchr`. Configured transcription factors are resolved to families by name, or by name and motif identifier when a symbol belongs to more than one family. betterChromVAR computes analytic deviations and z-scores per nucleus with GC-bias correction; no fragment-length bias term is used in the pinned version. Per-cell-type summaries test all families with a Wilcoxon marker test and report mean differences; cell-weighted mean heatmaps cover all families by ATAC cluster, GEX cluster and GEX cell type. The configured transcription factors select only the UMAP feature colourings. These values measure accessibility associated with a motif family, not transcription-factor activity.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| Families | JASPAR source | Hardcoded: target literal | `resources/JASPAR2026_vertebrate_familial_root_motifs.tf` (233 roots); `resources/JASPAR2026_vertebrate_motif_families.tsv` (1,019 motifs) | `extra_targets/setup_targets.R` |
| Families | Integrity check | Hardcoded: inline literal | counts only: 233 roots in order, 1,019 rows, 233 families | `extra_targets/setup_targets.R` |
| Families | Configured transcription factors | Configurable |  | [`aggregation_ATAC_marker_TFs`](../parameters.html#aggregation_ATAC_marker_TFs) |
| Families | Symbol resolution | Hardcoded: inline literal | upper-case `TF` or `TF__motifID`; ambiguous symbols error | `R/celltype_labeling_helpers.R` |
| Scanning | PFM construction | Hardcoded: inline literal | uniform 0.25 background, `+` strand, zero-sum columns dropped | `R/celltype_labeling_helpers.R` |
| Scanning | `motifmatchr::matchMotifs()` | Hardcoded: environment pin | package defaults: `p.cutoff = 5e-5`, `bg = "subject"`, `w = 7` | `R/celltype_labeling_helpers.R` |
| Scanning | Genome | Hardcoded: inline literal | BSgenome hg38, mm10 or mm39 by reference | `R/celltype_labeling_helpers.R` |
| chromVAR | Package | Hardcoded: environment pin | betterChromVAR 0.99.41 at commit `82ae1e4` | `scripts/github_packages.R` |
| chromVAR | Bias | Hardcoded: inline literal | `addGCBias()`; missing bias set to 0 | `R/celltype_labeling_helpers.R` |
| chromVAR | Expectation and peak filter | Hardcoded: inline literal | row mean over cells; zero-count peaks removed | `R/celltype_labeling_helpers.R` |
| chromVAR | Background bins and shrinkage | Hardcoded: environment pin | `getBackgroundBins()` defaults; `computeBackgrounds(shrinkage = "none")` | `R/celltype_labeling_helpers.R` |
| chromVAR | Deviations | Hardcoded: inline literal | `computeDeviationsAnalytic(denominator = "global")`, deviations and z | `R/celltype_labeling_helpers.R`, `extra_targets/ATAC_targets.R` |
| chromVAR | Cell set and chunking | Hardcoded: target literal | post-doublet-filter ATAC metadata; `chunk_nonzero_limit = 2^27` | `extra_targets/ATAC_targets.R` |
| Summaries | Per-cell-type test | Hardcoded: inline literal | `BPCells::marker_features(method = "wilcoxon")`, BH, on the ATAC cell-type column, all families | `R/celltype_labeling_helpers.R` |
| Summaries | Heatmaps | Hardcoded: target literal | cell-weighted means of all families by ATAC cluster, GEX cluster and GEX cell type | `extra_targets/ATAC_targets.R` |
| Summaries | Use of configured families | Hardcoded: target literal | UMAP feature colourings only | `extra_targets/ATAC_targets.R` |


<!-- source: website/implementation/methods_peak_gene_correlation.md -->

# Peak–gene correlation

This chapter describes the optional `peak_gene_correlation` module and lists every setting that determines its results, using the layout defined in [Methods and parameter tables](implementation_conventions.md#methods-and-parameter-tables). The target structure is shown in the [peak–gene correlation graph](implementation_peak_gene_correlation.md); the user-facing configuration is in [Peak–gene correlation](../downstream_peak_gene_correlation.html). Only two manifest parameters belong to this module; everything else is fixed. The module also depends indirectly on the ATAC embedding settings of the aggregation.

## Cell groups, candidate pairs and donor–state pseudobulks

The analysis runs separately within WNN annotation classes with enough nuclei. Consensus peaks are paired with gene transcription start sites (TSSs) on the same chromosome within a fixed window, measured from the peak centre; each pair is classified as self-promoter, gene-body, proximal or distal. Within a class, nuclei from eligible donors are partitioned into mutually exclusive ATAC-state bins by k-means on the scaled ATAC LSI or Harmony dimensions; the number of bins adapts to the median number of nuclei per eligible donor. Counts are summed within each donor–state combination. Pseudobulks below the minimum size are dropped, and donors and states are pruned iteratively until every donor contributes the minimum number of states. Per-nucleus GEX and ATAC depth columns are required; groups lacking donor or depth metadata are skipped with a diagnostic.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| Cell groups | Annotation column | Hardcoded: helper default | `WNN_harmony_SNN_cluster_cell_type` | `make_peak_gene_correlation_cell_groups()` in `R/peak_gene_correlation_helpers.R` |
| Cell groups | Minimum nuclei per class | Hardcoded: helper default | 200 | `R/peak_gene_correlation_helpers.R` |
| Pairs | TSS definition | Hardcoded: inline literal | gene start, or end on the minus strand | `R/peak_gene_correlation_helpers.R` |
| Pairs | Maximum peak-centre to TSS distance | Hardcoded: helper default | `max_distance = 250000` | `R/peak_gene_correlation_helpers.R` |
| Pairs | Self-promoter window | Hardcoded: inline literal | strand-aware −1,500 to +500 bp around the TSS | `R/peak_gene_correlation_helpers.R` |
| Pairs | Link classes | Hardcoded: inline literal | self-promoter, gene body, proximal (≤ 10 kb), distal, in that precedence | `R/peak_gene_correlation_helpers.R` |
| Donors | Donor column | Hardcoded: helper default | `donor_id` | `make_peak_gene_correlation_donor_state_record()` in `R/peak_gene_correlation_helpers.R` |
| Donors | Depth columns | Hardcoded: inline literal | first of `nCount_RNA`, `gex_umis_count`; first of `nCount_ATAC`, `atac_fragments` | `R/peak_gene_correlation_helpers.R` |
| States | Embedding source | Configurable |  | ATAC LSI or Harmony dimensions from [`aggregation_ATAC_data_PCs`](../parameters.html#aggregation_ATAC_data_PCs) and the Harmony covariates |
| States | Dimensions used | Hardcoded: helper default | `2:20`, intersected with available dimensions | `R/peak_gene_correlation_helpers.R` |
| States | Preprocessing | Hardcoded: inline literal | zero-variance dimensions dropped; columns standardized within the group | `R/peak_gene_correlation_helpers.R` |
| States | Maximum bins and bin rule | Hardcoded: helper default | 20; `min(20, floor(median nuclei per eligible donor / 20))`; fewer than 2 skips the group | `R/peak_gene_correlation_helpers.R` |
| States | k-means | Hardcoded: inline literal | Lloyd, `iter.max = 1000`, `nstart = 1`, seed 1; non-convergence is an error | `R/peak_gene_correlation_helpers.R` |
| Pseudobulks | Minimum nuclei per donor–state | Hardcoded: helper default | 20 | `R/peak_gene_correlation_helpers.R` |
| Pseudobulks | Minimum states per donor | Hardcoded: helper default | 2 | `R/peak_gene_correlation_helpers.R` |
| Pseudobulks | Donor eligibility before binning | Hardcoded: inline literal | at least 40 nuclei | `R/peak_gene_correlation_helpers.R` |
| Pseudobulks | Minimum donors per group and per state | Hardcoded: helper default | 1 and 1 | `R/peak_gene_correlation_helpers.R` |

## Normalization, eligibility and measurement-support filtering

GEX and ATAC pseudobulk counts are separately scaled to counts per million using the metadata-derived pseudobulk depth and log1p-transformed. A chromosome branch requires a minimum number of pseudobulks and residual degrees of freedom after the nuisance design, and genes and peaks must be detected in a minimum fraction of pseudobulks. The configurable measurement-support filter then removes hypotheses whose gene and peak are not both supported in enough shared donors; count thresholds scale with each aggregate's depth relative to the class median. Excluded hypotheses never enter the BH family.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| Normalization | Scale factor and transform | Hardcoded: helper default | `scale_factor = 1e6`, `log1p` | `R/peak_gene_correlation_helpers.R` |
| Eligibility | Gene detection fraction | Hardcoded: helper default | 0.05 | `prepare_peak_gene_correlation_branch()` in `R/peak_gene_correlation_helpers.R` |
| Eligibility | Peak accessibility fraction | Hardcoded: helper default | 0.05 | `R/peak_gene_correlation_helpers.R` |
| Eligibility | Minimum pseudobulks per branch | Hardcoded: helper default | 10 | `R/peak_gene_correlation_helpers.R` |
| Eligibility | Minimum residual df | Hardcoded: inline literal | pseudobulks − design rank − 1 ≥ 5 | `R/peak_gene_correlation_helpers.R` |
| Nuisance design | Covariates | Hardcoded: inline literal | donor fixed effects; standardized `log1p` GEX and ATAC depth when non-constant; QR-pruned to full rank | `make_peak_gene_correlation_design_matrix()` in `R/peak_gene_correlation_helpers.R` |
| Support filter | Preset | Configurable |  | [`peak_gene_correlation_filter`](../parameters.html#peak_gene_correlation_filter) |
| Support filter | Preset thresholds | Hardcoded: inline literal | lenient: RNA 5, ATAC 3 counts, `max(6, 10%)` aggregates, 2 shared donors, 2 aggregates per donor; moderate: 10, 5, `max(6, 10%)`, 2, 2; strict: 10, 5, `max(10, 20%)`, 3, 3 | `peak_gene_filter_settings()` in `R/peak_gene_filter_helpers.R` |
| Support filter | Depth scaling of count thresholds | Hardcoded: inline literal | `max(2, count × depth / median depth)` | `R/peak_gene_filter_helpers.R` |

## Conditional and hierarchical tests

The conditional analysis residualizes GEX and ATAC values against the nuisance design and reports the Pearson correlation of the residuals with an HC3 heteroskedasticity-robust regression test. The hierarchical analysis fits donor fixed intercepts, the depth covariates and a donor-varying slope for the within-donor-centred peak value, using project-owned compiled kernels for profiled restricted maximum likelihood and Kenward–Roger inference; the lme4 and pbkrtest route exists only as a test reference. It requires within-donor peak variation in at least two donors. Fits that fail the kernel diagnostics keep their estimates but no inferential P-value. BH correction runs within each annotation class over the complete eligible pair family, counting unreliable tests. Because every support preset requires at least two shared donors, single-donor classes yield diagnostics but no tests.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| HC3 | Correlation and test | Hardcoded: inline literal | Pearson on residuals; OLS slope with HC3 sandwich SE; two-sided t on residual df | `R/peak_gene_correlation_helpers.R` |
| HC3 | Numerical guards | Hardcoded: inline literal | leverage ≥ 1 − 1e-8 gives NA; residual-variation tolerance 1e-12 | `R/peak_gene_correlation_helpers.R` |
| HC3 | Multiplicity | Hardcoded: inline literal | BH by class over non-missing P-values | `R/peak_gene_correlation_helpers.R` |
| HC3 | Link definition for summary outputs | Hardcoded: inline literal | correlation ≥ 0.15, FDR \< 0.05, not self-promoter | `make_peak_gene_correlation_links()` in `R/peak_gene_correlation_helpers.R` |
| Hierarchical | Model | Hardcoded: inline literal | `y ~ 0 + design + x + (0 + x | donor)`, REML, x within-donor centred | `score_peak_gene_hierarchical_associations()` in `R/peak_gene_hierarchical_helpers.R` |
| Hierarchical | Engine | Hardcoded: inline literal | `src/peak_gene_REML.cpp`, `src/peak_gene_KR.cpp` (pbkrtest 0.5.5 equations) | `R/peak_gene_hierarchical_helpers.R` |
| Hierarchical | Minimum donors and within-donor variation | Hardcoded: inline literal | ≥ 2 donors; ≥ 2 donors with within-donor peak variation | `R/peak_gene_hierarchical_helpers.R` |
| Hierarchical | Variance-ratio search | Hardcoded: inline literal | grid over `expm1(0..32)` then golden section, tolerance 1e-9; zero allowed | `src/peak_gene_REML.cpp` |
| Hierarchical | Unreliable-fit rules | Hardcoded: inline literal | kernel diagnostics (no residual variation, ratio out of range, conditioning \< 1e-12, invalid KR covariance) or KR df \< 1 | `R/peak_gene_hierarchical_helpers.R`, `src/peak_gene_KR.cpp` |
| Hierarchical | P-value | Hardcoded: inline literal | F test with 1 and KR degrees of freedom | `src/peak_gene_KR.cpp` |
| Hierarchical | Multiplicity | Hardcoded: inline literal | BH by class, family size including unreliable tests | `R/peak_gene_hierarchical_helpers.R` |

## Prioritization, top links and plots

SuSiE fine-mapping prioritizes peaks for genes that have at least one conditional link, using the donor- and depth-residualized values. Top-link figures rank estimable positive hierarchical slopes by nominal P-value, excluding self-promoter peaks; gene-body peaks remain eligible, and no significance cutoff is applied. Each figure shows the gene context, the focal cell type's insertion coverage and the donor-residual scatter.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| SuSiE | Input and gene screen | Hardcoded: inline literal | HC3 results; genes with a link at correlation ≥ 0.15, FDR \< 0.05, not self-promoter | `R/peak_gene_finemapping_helpers.R` |
| SuSiE | Genes per branch and peaks per gene | Hardcoded: helper default | 50 genes by best FDR; 500 peaks by absolute correlation; at least 2 variable peaks | `R/peak_gene_finemapping_helpers.R` |
| SuSiE | Model settings | Hardcoded: inline literal | `L = min(10, n_peaks)`, `intercept = FALSE`, `standardize = TRUE`, `estimate_residual_variance = TRUE`, `max_iter = 100`, credible-set coverage 0.95 | `R/peak_gene_finemapping_helpers.R` |
| SuSiE | Records retained | Hardcoded: inline literal | PIP ≥ 0.01 or credible-set members, else the top peak | `R/peak_gene_finemapping_helpers.R` |
| Top links | Links per cell group | Configurable |  | [`peak_gene_correlation_top_links_per_cell_group`](../parameters.html#peak_gene_correlation_top_links_per_cell_group) |
| Top links | Selection | Hardcoded: inline literal | estimable, coefficient \> 0, not self-promoter; ordered by hierarchical P, gene, peak | `R/peak_gene_hierarchical_helpers.R` |
| Plots | Window padding and coverage | Hardcoded: inline literal | 2 kb padding; 500 bins; focal WNN cell type; fragments as read counts | `R/peak_gene_plot_helpers.R` |
| Plots | Coverage clip | Hardcoded: helper default | `clip_quantile = 0.999` | `R/ATAC_tracks_helpers.R` |
| Plots | Histogram bin width and distance bins | Hardcoded: target literal | 0.025; 5 kb bins capped at 245 kb | `module_peak_gene_correlation/correlation_targets.R`, `R/peak_gene_correlation_helpers.R` |


<!-- source: website/implementation/methods_differential_analyses.md -->

# Differential analyses

This chapter describes the optional `differential_analyses` module and lists every setting that determines its results, using the layout defined in [Methods and parameter tables](implementation_conventions.md#methods-and-parameter-tables). The target structure is shown in the [differential analyses graph](implementation_differential_analyses.md); the user-facing prerequisites and configuration are in [Differential analyses](../downstream_differential_analyses.html).

Donors, not nuclei, are the biological replicates. Both branches take their annotation classes from the WNN cell-type label carried in the final WNN metadata. The module retains its own metadata file and full-tibble targets, then projects a canonical analysis view containing only donors in the aggregation and the columns required by the configured models and composition plots. Rows are ordered by `donor_id` and non-key columns by name, so changes to unused columns or out-of-aggregation donors stop at this projection boundary.

## Cell-type composition

Nuclei are counted by donor and annotation class. Every observed class is tested unless `cell_types_to_test` restricts the response classes; missing donor–class combinations are completed with zero counts, so no class is dropped for being observed in few donors. Each class is fitted separately as a beta-binomial model with logit link on the two-column response of nuclei in the class versus nuclei in all other classes, using the configured one-sided formula. Random effects, custom design functions and two-sided formulas are rejected in this branch. Contrasts are named linear combinations of the fixed effects and are tested with normal Wald statistics. Fits that error, fail to converge or lack a positive-definite Hessian are marked non-estimable. Benjamini–Hochberg (BH) correction is applied across the tested classes within each model and contrast; there is no adjustment across contrasts or models.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| Input | Annotation column counted per donor | Hardcoded: target literal | `WNN_harmony_SNN_cluster_cell_type` | `model_data.cell_type_composition` in `module_differential_analyses/setup_and_cell_type_composition_targets.R` |
| Population | GEM wells defining the population | Configurable |  | `GEM_well_IDs` inside [`differential_analyses_cell_type_composition_models`](../parameters.html#differential_analyses_cell_type_composition_models) |
| Population | Donor restriction | Configurable |  | `donor_ids` inside [`differential_analyses_cell_type_composition_models`](../parameters.html#differential_analyses_cell_type_composition_models) |
| Response | Classes tested; denominators always use all retained nuclei | Configurable |  | `cell_types_to_test` inside [`differential_analyses_cell_type_composition_models`](../parameters.html#differential_analyses_cell_type_composition_models) |
| Response | Two-column response prepended to the formula | Hardcoded: inline literal | `cbind(n_nuclei, n_other_nuclei)` | `fit_cell_type_composition_model()` in `R/differential_analysis_helpers.R` |
| Model | Fixed-effects formula | Configurable |  | `formula` inside [`differential_analyses_cell_type_composition_models`](../parameters.html#differential_analyses_cell_type_composition_models) |
| Model | Family and link; no dispersion or zero-inflation formula | Hardcoded: inline literal | `glmmTMB::betabinomial(link = "logit")`, package defaults | `fit_cell_type_composition_model()` in `R/differential_analysis_helpers.R` |
| Model | Random effects, design functions, contrast functions, paired fields | Hardcoded: inline literal | rejected with an error | `fit_cell_type_composition_model()` in `R/differential_analysis_helpers.R` |
| Model | Design validity | Hardcoded: inline literal | full rank, finite, more donors than coefficients | `validate_differential_design()` in `R/differential_analysis_helpers.R` |
| Contrasts | Named linear contrasts | Configurable |  | `contrast_specs_vec` inside [`differential_analyses_cell_type_composition_models`](../parameters.html#differential_analyses_cell_type_composition_models) |
| Test | Wald test from the conditional covariance | Hardcoded: inline literal | two-sided normal | `fit_cell_type_composition_model()` in `R/differential_analysis_helpers.R` |
| Test | Non-estimable rule | Hardcoded: inline literal | fit error, non-zero convergence code, or non-positive-definite Hessian | `fit_cell_type_composition_model()` in `R/differential_analysis_helpers.R` |
| Multiplicity | BH across tested classes within model and contrast | Hardcoded: inline literal | `p.adjust(method = "BH")` | `R/differential_analysis_helpers.R` |
| Plot | Interval and significance colour | Hardcoded: inline literal | ±1.96 SE; FDR \< 0.05 | `R/differential_analysis_helpers.R` |
| Plot | Phenotype panels and colour variable | Configurable |  | `plot_phenotype_vars`, `color_by` inside [`differential_analyses_cell_type_composition_models`](../parameters.html#differential_analyses_cell_type_composition_models) |

## Pseudobulk construction

GEX and ATAC counts are summed within each donor and annotation class in the primary module, so this module reuses the same pseudobulk targets as the compatibility export. Sample identifiers combine the class and the donor. Four feature matrices are tested:

- **Gene expression (DGE)**: the GEX pseudobulk count matrix.
- **Chromatin accessibility (DCA)**: the consensus-peak pseudobulk count matrix after peak-level QC.
- **Motif-family accessibility (DTFA)**: betterChromVAR deviations of the JASPAR 2026 familial root motifs computed on the pseudobulk ATAC counts. Peaks with zero pseudobulk counts are dropped before the background model. The z-scores are column-centred and quantile-normalised across samples, so this matrix is continuous.
- **Transcription-factor activity (DCTA)**: signed CollecTRI regulator activities inferred from the GEX pseudobulks with the `decoupleR` univariate linear model. Genes are filtered by expression across classes, library sizes are normalised and log-CPM values are computed before inference. Regulators need a minimum number of measured targets. The CollecTRI network is downloaded from the OmniPath rescue archive and accepted only when it matches the pinned checksum.

Motif families come from the official JASPAR 2026 CORE vertebrate clustering. Each of the 233 families is represented by its published root motif, which is scanned directly against the consensus peaks; individual member motifs are used only as family metadata. The same family-level matrix supports marker plots and the Seurat export.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| Summation | Aggregation method | Hardcoded: inline literal | `BPCells::pseudobulk_matrix(method = "sum")` | `get_BPCells_pseudobulk_matrix()` in `R/pseudobulk_helpers.R` |
| Summation | Grouping column for GEX and ATAC | Hardcoded: target literal | `WNN_harmony_SNN_cluster_cell_type` | `pseudobulk_counts_matrix.GEX` and `.ATAC` in `extra_targets/general_aggregation_targets.R` |
| Summation | ATAC input matrix | Hardcoded: target literal | `peak_QC_filtered_BPCells_matrix.ATAC` | `extra_targets/general_aggregation_targets.R` |
| DTFA | Deviation computation | Hardcoded: inline literal | betterChromVAR `compute = c("deviations", "z")`, `normalize = TRUE` | `R/pseudobulk_helpers.R` |
| DTFA | Post-processing of z-scores | Hardcoded: inline literal | column centring, `limma::normalizeBetweenArrays(method = "quantile")` | `R/pseudobulk_helpers.R` |
| DTFA | Peaks dropped before background | Hardcoded: inline literal | zero pseudobulk row sum | `R/pseudobulk_helpers.R` |
| DTFA | Motif universe | Hardcoded: target literal | `resources/JASPAR2026_vertebrate_familial_root_motifs.tf`, 233 families checked | `extra_targets/setup_targets.R`, `extra_targets/ATAC_targets.R` |
| DTFA | Family membership table | Hardcoded: target literal | `resources/JASPAR2026_vertebrate_motif_families.tsv`, 1,019 motifs | `extra_targets/setup_targets.R` |
| DCTA | Network source and checksum | Hardcoded: target literal | `https://rescued.omnipathdb.org/CollecTRI.csv`, SHA-256 `86c90b30f2cc75c189da1f0a8c353d1547287cd656a9fac1c678634285bcb4e0` | `module_differential_analyses/shared_targets.R` |
| DCTA | Network validation | Hardcoded: inline literal | 43,536 interactions, 1,189 sources, `mor` in {−1, 1} | `R/TF_activity_helpers.R` |
| DCTA | Gene filter before inference | Hardcoded: inline literal | `edgeR::filterByExpr(group = cluster)` | `R/TF_activity_helpers.R` |
| DCTA | Normalisation and log-CPM | Hardcoded: inline literal | `normLibSizes()` (TMM), `cpm(log = TRUE, prior.count = 2)` | `R/TF_activity_helpers.R` |
| DCTA | Minimum measured targets per regulator | Hardcoded: helper default | `min_targets = 5` | `R/TF_activity_helpers.R`, not overridden in `module_differential_analyses/targets.R` |

## Sample and feature filtering

For each configured pseudobulk model, donors missing any variable in the formula, donors outside an optional `donor_ids` list, and samples outside an optional annotation-class subset are removed. The DTFA branch additionally removes samples below the configured minimum ATAC depth. The remaining steps depend on the model route, described next. For the routes that model counts, zero-depth samples are removed, features are filtered with the design-aware `edgeR::filterByExpr()`, and library sizes are normalised with `edgeR::normLibSizes()`.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| Donors | Variables that must be non-missing | Configurable |  | `formula` inside [`differential_analyses_pseudobulk_models`](../parameters.html#differential_analyses_pseudobulk_models) |
| Donors | Donor restriction | Configurable |  | `donor_ids` inside [`differential_analyses_pseudobulk_models`](../parameters.html#differential_analyses_pseudobulk_models) |
| Donors | Extended donor metadata table | Configurable |  | [`differential_analyses_extended_donor_id_metadata_tsv`](../parameters.html#differential_analyses_extended_donor_id_metadata_tsv) |
| Samples | Annotation-class subset | Configurable |  | `cell_type_subset` inside [`differential_analyses_pseudobulk_models`](../parameters.html#differential_analyses_pseudobulk_models) |
| Samples | Minimum ATAC depth, DTFA only | Configurable |  | [`differential_analyses_motif_family_accessibility_min_ATAC_counts`](../parameters.html#differential_analyses_motif_family_accessibility_min_ATAC_counts) |
| Samples | Zero-depth removal, count routes | Hardcoded: inline literal | column sum \> 0 | `R/pseudobulk_helpers.R` |
| Features | Expression filter, count routes | Hardcoded: inline literal | `edgeR::filterByExpr(design = design_matrix)`, package defaults | `R/pseudobulk_helpers.R` |
| Features | Library-size normalisation, count routes | Hardcoded: inline literal | `edgeR::normLibSizes()`, TMM | `R/pseudobulk_helpers.R` |
| Residual df | Minimum residual degrees of freedom, checked before and after filtering | Hardcoded: inline literal | samples − coefficients \> 0 | `R/pseudobulk_helpers.R` |

## Model routes

Model matrices come from the configured formula or a custom design function, and contrasts from named linear expressions or custom contrast functions. Whether a matrix holds counts is decided by checking a random sample of entries for integers. The route is then selected as follows:

1.  A configured `random_effect` selects the correlation route: `limma::voom()` for counts, then one pass of `limma::duplicateCorrelation()` with the random effect as block, then `lmFit()` with the consensus correlation.
2.  Otherwise, a count matrix with more than one coefficient selects the edgeR route: `estimateDisp()` with package defaults, `glmQLFit(robust = TRUE)` and `glmQLFTest()`.
3.  Otherwise, the basic limma route: `lmFit()`, `contrasts.fit()` and `eBayes()` with package defaults. Note that a count matrix with a single coefficient takes this route on the raw counts, without the filtering and normalisation described above.

An optional paired-cell-type route is selected by a `cell_type_formula`. It fits each annotation class separately, permits at most one pseudobulk per pairing unit and class, and requires at least two classes. Count matrices are fitted per class with `edgeR::voomLmFit()` blocked on the correlation block; continuous matrices use `duplicateCorrelation()` and `lmFit()`. Residual correlations between classes are estimated on shared donors from the common features and combined into cross-class contrast tests. A class pair needs more shared donors than coefficients plus two.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| Design | Formula or custom design function | Configurable |  | `formula`, `design_matrix_func_name` inside [`differential_analyses_pseudobulk_models`](../parameters.html#differential_analyses_pseudobulk_models) |
| Design | Contrasts | Configurable |  | `contrast_specs_vec` and custom contrast functions inside [`differential_analyses_pseudobulk_models`](../parameters.html#differential_analyses_pseudobulk_models) |
| Design | Coefficient-name sanitising | Hardcoded: inline literal | `:` to `.`, `-` to `_` | `R/pseudobulk_helpers.R` |
| Route | Count detection | Hardcoded: inline literal | integer check on up to 10 × 10 sampled entries | `is_count_matrix()` in `R/data_transformations.R` |
| Route | Correlation route trigger | Configurable |  | `random_effect` inside [`differential_analyses_pseudobulk_models`](../parameters.html#differential_analyses_pseudobulk_models) |
| Correlation route | Fit | Hardcoded: inline literal | `limma::voom()` for counts, one `duplicateCorrelation()` pass, `lmFit(correlation, block)` | `R/pseudobulk_helpers.R` |
| edgeR route | Dispersion and fit | Hardcoded: inline literal | `estimateDisp()` defaults, `glmQLFit(robust = TRUE)`, `glmQLFTest()` | `R/pseudobulk_helpers.R` |
| limma routes | Moderation | Hardcoded: inline literal | `eBayes()` defaults, no `robust` or `trend` | `R/pseudobulk_helpers.R` |
| Paired route | Trigger and pairing | Configurable |  | `cell_type_formula`, `pairing_variable`, `correlation_block` inside [`differential_analyses_pseudobulk_models`](../parameters.html#differential_analyses_pseudobulk_models) |
| Paired route | One pseudobulk per pairing unit and class; at least two classes | Hardcoded: inline literal | error otherwise | `R/pseudobulk_helpers.R` |
| Paired route | Per-class count fit | Hardcoded: inline literal | `edgeR::voomLmFit(block, normalize.method = "none")` | `R/pseudobulk_helpers.R` |
| Paired route | Minimum shared donors per class pair | Hardcoded: inline literal | more than coefficients + 2 | `R/pseudobulk_helpers.R` |
| Paired route | Residual-correlation estimate | Hardcoded: helper default | up to 2,000 evenly spaced common features; Fisher-z trimmed mean, `trim = 0.15` | `R/pseudobulk_helpers.R` |
| Paired route | Cross-class statistic | Hardcoded: inline literal | t from the two class estimates and their covariance; df is the minimum per-class total df | `R/pseudobulk_helpers.R` |
| Paired route | Parallel class fits | Hardcoded: inline literal | up to 6 forks | `R/pseudobulk_helpers.R`; `cores_req = 6` in `module_differential_analyses/pseudobulk_differential_targets.R` |
| Failure | Conditions that fail the branch | Hardcoded: inline literal | no samples or features, residual df below one, invalid block or contrast, too few shared donors | `R/pseudobulk_helpers.R` |

## Multiplicity, significance and gene sets

For each model and contrast, BH FDR is calculated across the tested features. A feature is counted as significant when its FDR is below the threshold and its effect is non-zero. Differential-expression statistics are tested against the MSigDB Hallmark and Reactome collections with the competitive `cameraPR` test. The statistic is the moderated t where the route provides one and otherwise a signed normal quantile derived from the nominal P-value, which is the case for the edgeR route. A set must contain the minimum number of tested genes, and FDR is calculated within each contrast and collection. The optional Open Targets annotation queries the platform for the top features of the gene-expression branch when an EFO identifier is configured.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| Multiplicity | FDR across features per contrast | Hardcoded: inline literal | `p.adjust(method = "BH")` | `R/pseudobulk_helpers.R` |
| Significance | FDR threshold for significant counts | Hardcoded: helper default | `FDR_threshold = 0.05`, plus `logFC != 0` | `R/pseudobulk_helpers.R`, not overridden in `module_differential_analyses/pseudobulk_differential_targets.R` |
| Significance | Threshold used in plots and comparisons | Hardcoded: inline literal | 0.05 | `R/pseudobulk_helpers.R`, `R/TF_activity_helpers.R` |
| Top features | Features labelled and queried per contrast | Hardcoded: target literal | `n = 40` by nominal P | `module_differential_analyses/pseudobulk_differential_targets.R` |
| Gene sets | Collections | Hardcoded: target literal | MSigDB `H`; `C2` with `CP:REACTOME` | `module_differential_analyses/targets.R` |
| Gene sets | Species mapping | Hardcoded: inline literal | `Homo_sapiens` to human, `Mus_musculus` to mouse | `R/pseudobulk_helpers.R` |
| Gene sets | msigdbr version | Hardcoded: environment pin | `r-msigdbr >=26.1.0,<27` | `pixi.toml` |
| Gene sets | Test | Hardcoded: inline literal | `limma::cameraPR(inter.gene.cor = 0.01)` | `R/pseudobulk_helpers.R` |
| Gene sets | Minimum tested genes per set | Hardcoded: helper default | `min_genes_per_set = 10` | `R/pseudobulk_helpers.R`, not overridden in `module_differential_analyses/gene_set_enrichment_targets.R` |
| Gene sets | Multiplicity | Hardcoded: inline literal | BH within contrast and collection | `R/pseudobulk_helpers.R` |
| Gene sets | Terms shown per plot | Hardcoded: inline literal | 25 smallest P | `R/pseudobulk_helpers.R` |
| Open Targets | Trait identifier; empty skips the query | Configurable |  | [`differential_analyses_pseudobulk_OT_GWAS_efo_id`](../parameters.html#differential_analyses_pseudobulk_OT_GWAS_efo_id) |
| Open Targets | Scope and endpoint | Hardcoded: inline literal | gene-expression branch only; `https://api.platform.opentargets.org/api/v4/graphql` | `module_differential_analyses/pseudobulk_differential_targets.R`, `R/pseudobulk_helpers.R` |

## Cross-modality comparison

The cross-modality fragment maps CollecTRI regulators to JASPAR families and compares model t-statistics, not raw activity scales. AP1 and NFKB stay intact as complex regulons during activity inference; their canonical members are used only to associate the complexes with motif families. Source-level tables retain TF expression as a third reference. Family-level summaries use the median CollecTRI regulator t-statistic and report whether any mapped source is FDR-significant.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| Crosswalk | Complex members | Hardcoded: inline literal | AP1: FOS, FOSB, FOSL1, FOSL2, JUN, JUNB, JUND; NFKB: NFKB1, NFKB2, REL, RELA, RELB | `R/TF_activity_helpers.R` |
| Summary | Family-level statistics | Hardcoded: inline literal | median CollecTRI t, median TF-expression t, first motif-family t and FDR | `R/TF_activity_helpers.R` |
| Summary | Concordance | Hardcoded: inline literal | Spearman correlation, at least 3 mapped families per contrast | `R/TF_activity_helpers.R` |

## Diagnostics and resources

Diagnostic outputs report pseudobulk depth with the ATAC threshold marked, cohort tables per model, contrast support (samples, donors, paired donors, smallest group), P-value distributions, signed significant-feature counts and volcano plots.

| Step | Setting | Status | Value | Source |
|---------------|---------------|---------------|---------------|---------------|
| Resources | Filtering, fitting and DCTA allocations | Hardcoded: target literal | filter 60 GB; fit 6 cores, 60 GB; DCTA 32 GB | `module_differential_analyses/pseudobulk_differential_targets.R`, `module_differential_analyses/targets.R` |
| Environment | Package pins | Hardcoded: environment pin | `bioconductor-edger >=4.8.2`, `bioconductor-limma >=3.66.0`, `r-glmmtmb >=1.1.14`, `bioconductor-decoupler >=2.16.0` | `pixi.toml` |


<!-- source: website/implementation/methods_genetic_enrichment.md -->

# Genetic enrichment

This chapter describes the optional `genetic_enrichment` module and lists every setting that determines its results, using the layout defined in [Methods and parameter tables](implementation_conventions.md#methods-and-parameter-tables). The target structure is shown in the [genetic enrichment graph](implementation_genetic_enrichment.md); the user-facing configuration is in [Genetic enrichment](../downstream_genetic_enrichment.html). The sparse SCAVENGE reimplementation and its validation are described in [Algorithmic implementations](algorithm_validation.md#sparse-scavenge-propagation-and-significance).

## GWAS inputs and peak weights

Study identifiers of the Open Targets form are resolved against the pinned platform release, whose study, credible-set, credible-set evidence and target datasets are downloaded once. Any other identifier is treated as a local Parquet file that must satisfy the local schema. With automatic fine-mapping selection, the first available method in a fixed priority order is used. Variants with posterior probability above the configured cutoff are retained and mapped to the consensus peaks that carry chromVAR state. Variant weights mapping to the same peak are summed and capped at one.

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Studies | Study label, category, source identifier, fine-mapping method | Configurable | | fields inside [`genetic_enrichment_GWAS_studies`](../parameters.html#genetic_enrichment_GWAS_studies) |
| Open Targets | Platform release | Hardcoded: target literal | `26.03` | `module_genetic_enrichment/shared_targets.R` |
| Open Targets | Datasets | Hardcoded: target literal | `credible_set`, `study`, `evidence_gwas_credible_sets`, `target` | `module_genetic_enrichment/shared_targets.R` |
| Open Targets | Source classification | Hardcoded: inline literal | `^GCST[0-9]+$` is Open Targets; anything else is a local path | `R/GWAS_chromVAR_input_helpers.R` |
| Open Targets | Study type | Hardcoded: inline literal | `studyType == "gwas"` | `R/GWAS_chromVAR_input_helpers.R` |
| Fine-mapping | Automatic priority | Hardcoded: inline literal | SuSie, SuSiE-inf, PICS | `resolve_open_targets_GWAS_input_tibble()` in `R/GWAS_chromVAR_input_helpers.R` |
| Credible sets | Probability | Hardcoded: inline literal | 0.95 | `R/GWAS_chromVAR_input_helpers.R` |
| Credible sets | Chromosomes | Hardcoded: inline literal | autosomes, X, Y, MT | `R/GWAS_chromVAR_input_helpers.R` |
| Local files | Validator | Hardcoded: inline literal | 29 required columns, schema version 1, GRCh38, SHA-256 field | `validate_local_finemapped_GWAS_tibble()` in `R/GWAS_chromVAR_input_helpers.R` |
| Variants | Posterior probability cutoff | Configurable | | [`genetic_enrichment_posterior_probability_cutoff`](../parameters.html#genetic_enrichment_posterior_probability_cutoff) |
| Variants | Comparison | Hardcoded: inline literal | strictly greater than the cutoff | `R/GWAS_chromVAR_helpers.R` |
| Peak weights | Peak set | Hardcoded: target literal | rows of the ATAC chromVAR object | `module_genetic_enrichment/gchromVAR_targets.R` |
| Peak weights | Combination | Hardcoded: helper default | `weight_transform = "cap_1"`: sum, capped at 1 | `R/GWAS_chromVAR_helpers.R` |

## Nucleus-level deviations

The trait peak weights form a chromVAR annotation. Deviations and z-scores per nucleus are computed analytically with betterChromVAR on the ATAC chromVAR object reused from the primary module, with the GC-bias background described in [Cell-type annotation and motif accessibility](methods_annotation_and_motifs.md#motif-families-and-motif-accessibility).

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Object | chromVAR object and chunk context | Hardcoded: target literal | `chromVAR_obj.ATAC`, `chromVAR_chunk_context_records.ATAC` | `module_genetic_enrichment/targets.R` |
| Deviations | Method | Hardcoded: inline literal | `computeDeviationsAnalytic(compute = "z")` | `R/GWAS_chromVAR_helpers.R` |
| Background | Bias and bins | Hardcoded: environment pin | `addGCBias()`; `getBackgroundBins()` defaults | `R/celltype_labeling_helpers.R` |

## Annotation-class pseudobulk deviations

Annotation-class deviations are calculated separately rather than by averaging nucleus-level results. ATAC counts are summed by the GEX-derived cell-type label carried into the final WNN metadata, a betterChromVAR background is fitted to the peak-by-class matrix, and the raw deviation is the observed-minus-background accessibility of the weighted peaks relative to their expected accessibility. A relative deviation standardizes the raw deviations across classes within each trait. The analytic z-score uses the background variance; one-sided P-values and BH-adjusted values over the whole table are retained. Plot labels use unadjusted z thresholds. An absolute-effect branch weights variants by posterior probability times effect size when effect sizes are available.

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Pseudobulks | Grouping column | Hardcoded: target literal | `PCA_harmony_SNN_cluster_cell_type` from the WNN metadata | `module_genetic_enrichment/GWAS_chromVAR_cell_type_targets.R` |
| Pseudobulks | Aggregation | Hardcoded: inline literal | `BPCells::pseudobulk_matrix(method = "sum")`, 6 threads | `R/pseudobulk_helpers.R` |
| Pseudobulks | Peak filter | Hardcoded: inline literal | zero-count peaks removed before the background | `R/pseudobulk_helpers.R` |
| Deviations | Raw, relative and z | Hardcoded: inline literal | weighted observed minus background over expected; standardized within trait; z from background variance | `R/GWAS_chromVAR_contribution_helpers.R` |
| Deviations | P-values and adjustment | Hardcoded: inline literal | one-sided normal; BH over all trait × class rows | `R/GWAS_chromVAR_contribution_helpers.R` |
| Plots | Support labels | Hardcoded: inline literal | `**` for z ≥ 2.326, `*` for z ≥ 1.645, unadjusted | `R/GWAS_chromVAR_contribution_helpers.R` |
| Plots | Compartment grouping | Configurable | | [`genetic_enrichment_compartment_patterns`](../parameters.html#genetic_enrichment_compartment_patterns) |
| Plots | Unmatched classes and ordering | Hardcoded: inline literal | `Other`; hierarchical clustering within compartments of more than two classes | `R/GWAS_plot_helpers.R` |
| Absolute effect | Eligibility | Hardcoded: inline literal | variant-level effects when all variants have them, else locus-level, else skipped | `R/GWAS_chromVAR_absolute_effect_helpers.R` |

## Locus attribution and detail plots

Per-class contributions are attributed to loci and variants and reconciled against the class totals. Loci are labelled with the highest-scoring locus-to-gene (L2G) genes from the Open Targets evidence. Detail plots are drawn for the top loci of classes passing a configurable z screen.

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| L2G | Score filter and label | Hardcoded: inline literal | score ≥ 0.05; top 3 genes | `R/GWAS_chromVAR_contribution_helpers.R` |
| Reconciliation | Tolerances | Hardcoded: inline literal | 1e-10 peak–variant; 1e-8 level sums | `R/GWAS_chromVAR_contribution_helpers.R` |
| Detail plots | Minimum z | Configurable | | [`genetic_enrichment_variant_detail_min_z`](../parameters.html#genetic_enrichment_variant_detail_min_z) |
| Detail plots | Loci, flank, coverage | Hardcoded: helper default | 3 loci per class; 25 kb flank; 500 bins; 0.999 coverage cap | `R/GWAS_chromVAR_contribution_helpers.R` |
| Attribution plots | Loci shown | Hardcoded: target literal | heatmaps 15; bar plots 5 | `module_genetic_enrichment/GWAS_chromVAR_contribution_targets.R` |

## SCAVENGE trait-relevance propagation

Trait-relevance scores are computed from the nucleus-level z-scores on the WNN SNN graph only. Nuclei whose one-sided normal-tail probability is at or below the seed cutoff are seeds, capped at the configured fraction; when no nucleus qualifies, all scores are zero. The graph's nonzero support becomes binary adjacency, degree-zero nuclei are removed, and seed mass is propagated by a random walk with restart until the L1 change is below the tolerance. Degree-matched seed permutations are sampled sequentially in R and evaluated by a native worker that streams per-nucleus exceedance counts. Cell-level empirical P-values are the exceedance fraction; cluster-level tests compare medians with add-one P-values and BH within each grouping. Scores are capped at an upper quantile, min–max scaled and multiplied by the mean z of the top nuclei.

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Graph | Representation | Hardcoded: target literal | `WNN_harmony_SNN` only | `module_genetic_enrichment/targets.R` |
| Graph | Binarization | Hardcoded: inline literal | nonzero support set to 1; column-normalized transition matrix | `R/SCAVENGE_helpers.R` |
| Seeds | P cutoff | Hardcoded: helper default | `p_value_cutoff = 0.05` | `R/SCAVENGE_helpers.R` |
| Seeds | Maximum seed fraction | Configurable | | [`genetic_enrichment_SCAVENGE_seed_percent`](../parameters.html#genetic_enrichment_SCAVENGE_seed_percent) |
| Seeds | z pre-filter | Hardcoded: helper default | non-finite dropped; `max_z_score = 1000` | `R/SCAVENGE_helpers.R` |
| Walk | Restart probability | Configurable | | [`genetic_enrichment_SCAVENGE_restart_prob`](../parameters.html#genetic_enrichment_SCAVENGE_restart_prob) |
| Walk | Convergence | Hardcoded: helper default | L1 change ≤ 1e-5; at most 10,000 iterations | `R/SCAVENGE_helpers.R`, `src/scavenge_random_walk.cpp` |
| Permutations | Count | Configurable | | [`genetic_enrichment_SCAVENGE_permutation_times`](../parameters.html#genetic_enrichment_SCAVENGE_permutation_times) |
| Permutations | Cores and chunks | Hardcoded: target literal | 15 cores; 4 chunks per core | `module_genetic_enrichment/SCAVENGE_graph_targets.R` |
| Cell P | Definition and call | Hardcoded: inline literal | exceedances (strictly greater) / permutations; significant at P ≤ 0.05 | `R/SCAVENGE_helpers.R` |
| Cluster P | Definition | Hardcoded: inline literal | median score; (exceedances + 1)/(permutations + 1) with ≥; BH within grouping | `R/SCAVENGE_helpers.R` |
| Score | Cap, scaling and scale factor | Hardcoded: helper default | cap at the 0.95 quantile; min–max; × mean z of the top 1 % | `R/SCAVENGE_helpers.R` |
| Summaries | Groupings | Hardcoded: helper default | `WNN_harmony_SNN_cluster_named`, `WNN_harmony_SNN_cluster_cell_type` | `R/SCAVENGE_helpers.R` |
| Summaries | Statistics | Hardcoded: inline literal | count, significant count and proportion, median, mean, quartiles, range, null median and 0.95 quantile | `R/SCAVENGE_helpers.R` |
| Plots | Heatmap labels | Hardcoded: inline literal | on BH-adjusted cluster P: `***` ≤ 0.001, `**` ≤ 0.01, `*` ≤ 0.05 | `R/GWAS_plot_helpers.R` |


# Orphaned Markdown Pages

Tracked Markdown files not reached from the Quarto book graph or include graph.

- `website/data/README.md`
- `website/figures/human_curated/README.md`
