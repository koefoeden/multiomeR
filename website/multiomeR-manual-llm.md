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
multiomeR is in beta and may introduce breaking changes between releases. The steps in [Run your own analysis](main_running.md#steps) show what to inspect before each stage, but acceptance criteria still depend on the tissue, study design, and intended use. Report problems or questions through [GitHub issues](https://github.com/koefoeden/multiomeR/issues).
:::

## Your first analysis

Start with the public demo: two human GEM wells with supplied configuration. You will install the software, run one joint analysis, and read its cell metadata and multimodal Seurat object. This gives you a working example before you choose settings for your own study.

You need basic R skills, a Linux terminal, and a machine with sufficient [memory and disk space](demo_installation.md#system-requirements). You do not need to know how to write a `targets` pipeline. Commands labeled **Bash** run in the terminal; commands labeled **R** run in the R session opened during installation. Run both from the repository folder unless stated otherwise.

## Terms used in this manual

In this manual, a **GEM well** is one configured 10x library and output directory, an **aggregation** is a joint analysis of one or more GEM wells, and a **donor** is the individual identified by `donor_id`. One GEM well may contain multiple donors.

A **target** is a named result, such as a metadata table, matrix directory, or plot. Its **dependencies** are the inputs and earlier results needed to build it. You request the result you want; `targets` works out the order and reuses results that are up to date. The **store** is the folder where it keeps results and the records needed for reruns. For a small worked introduction, see the [targets walkthrough](https://books.ropensci.org/targets/walkthrough.html).

## Workflow at a glance

The **main pipeline** processes each GEM well, aggregates selected GEM wells, and builds multimodal RNA/ATAC outputs for clustering, cell typing, and WNN integration. Three optional modules extend completed aggregations with differential analyses, genetic enrichment or peak–gene correlation. The [output gallery](gallery.md) previews every plot they and the main pipeline save.

[Image omitted; source: `figures/multiomeR_overview_simplified.drawio.svg`; alt: multiomeR workflow from Cell Ranger ARC GEM well outputs through per GEM well processing, aggregation-level GEX and A...]

Continue to [Install and prepare the demo](demo_installation.md).


## Part: Try the public demo


<!-- source: website/demo_installation.md -->

# Install and prepare the demo



## System requirements

- Linux system with `git` and `curl`
- At least 60 GB of RAM. This is enough for one heavy target at a time; machines near the minimum should reduce concurrent workers in `crew_controllers.R`.
- Multiple CPU cores are strongly recommended. The timing quoted in the next chapter was measured with 16 logical threads.

::: {.callout-tip title="Machines with less than 256 GB of RAM"}
The committed `crew_controllers.R` is sized for a 16-CPU, 256-GB workstation: four light workers and two heavy workers that may each use 60 GB. On a machine near the 60-GB minimum, edit the two `workers` values in `crew_controllers.R` before running the demo so that only one heavy target runs at a time:

``` {.r filename="crew_controllers.R"}
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

The `RAM_GB` values in the same file describe routing capacity, not enforced limits, so the workers that can run at once must fit in physical memory. To run on a SLURM or other scheduler instead, see [Choose where the analysis runs](performance_distributed_computing.md).
:::

## Set up the demo

Run this block from the directory where you want to clone multiomeR to. The single `pixi run` setup command installs the locked environment before its `setup-demo` task downloads the two configured public inputs and installs the pinned GitHub-only R packages.

``` {.bash filename="Bash"}
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

The download task is restart-safe: non-empty files already present under `example_data` are skipped.

Continue to [Run the demo](demo_running.md) from the R prompt.


<!-- source: website/demo_running.md -->

# Run the demo



In the R session opened during installation, run the command below to process the `immune_human_2x` aggregation. It combines the two demo GEM wells, produces a Seurat/Signac object containing the multimodal results, and draws the integrated WNN UMAPs colored by cluster, cell type, and the other categorical metadata.

`names` selects these two targets by their exact names using `all_of()`. `tar_make()` also builds the dependencies needed for them, but does not build every plot in the [output gallery](gallery.md).

``` {.r filename="R"}
demo_targets <- c(
  "multimodal_Seurat_object.8_multimodal_QC.immune_human_2x",
  "categorical.UMAPs.8_multimodal_QC.immune_human_2x"
)

targets::tar_make(names = tidyselect::all_of(demo_targets))
```

Keep the R session open until the command finishes. Progress messages report targets being dispatched, completed, or skipped because they are already up to date. Using 16 threads, this should take \~ 30 minutes, writing about 6 GB to disk.

## Confirm success

After the run, the following command should return `character(0)`, meaning the requested results and their dependencies are up to date:

``` {.r filename="R"}
targets::tar_outdated(
  names = tidyselect::all_of(demo_targets),
  callr_function = NULL
)
```

If names are returned, those results still need building. If the run failed, follow [Troubleshooting](troubleshooting.md), fix the reported cause, and run the same `tar_make()` command again. Completed results can be reused.

Continue to [Inspect the demo results](demo_outputs.md) to read the object and find the associated files.


<!-- source: website/demo_outputs.md -->

# Inspect the demo results



The pipeline saves four main kinds of output in its targets store, normally `outputs/`. The `store` setting in `_targets.yaml` selects this folder.

- Serialized R objects in `objects/`, managed by targets.
- Data files in `files/`, grouped by target and analysis.
- Plot images in `plots/`, normally in PNG format.
- Editable plot objects in `plot_objects/`, saved as RDS files alongside the corresponding image hierarchy.

## Objects

Most intermediate and final result objects are saved automatically by `targets` and can be loaded into any repository-root R-session using `targets::tar_read()`:

``` {.r filename="R"}
cell_metadata <- targets::tar_read(
  metadata_w_cell_types_tibble.WNN.immune_human_2x
)
dim(cell_metadata)
head(cell_metadata)

demo_object <- targets::tar_read(multimodal_Seurat_object.8_multimodal_QC.immune_human_2x)
demo_object
```

The metadata table describes the retained nuclei and their annotations. WNN means *weighted nearest neighbors*: the integrated representation uses information from both RNA and ATAC. The Seurat/Signac object is a convenient export for further exploration; the pipeline also retains its matrices in BPCells format on disk.

## Files

File targets also load with `targets::tar_read()`, but their value is a path rather than an in-memory result, so you will have to load them yourself using the appropriate tool.

``` {.r filename="R"}
targets::tar_read(cellranger_barcodes_tsv.healthy_PBMC_human)
targets::tar_read(aggregated_GEX_BPCells_matrix_dir.GEX.immune_human_2x)
targets::tar_read(consensus_peak_BPCells_matrix_dir.ATAC.immune_human_2x)
```

As you can see from the output above, pipeline-generated files generally follow a folder hierarchy derived from their target names; `tar_read()` gives their actual paths.

## Plots

The demo command also built `categorical.UMAPs.8_multimodal_QC.immune_human_2x`. Plots use the same scope-based layout as files, but under `outputs/plots/` instead. Note that a target might produce multiple files, as seen in the example below:

``` {.r filename="R"}
targets::tar_read(categorical.UMAPs.8_multimodal_QC.immune_human_2x)
```

Open `WNN_harmony_SNN_cluster_cell_type.png` there to see the integrated clusters and cell-type labels from your own run. It should resemble this documentation snapshot:

[Image omitted; source: `figures/demo_WNN_cell_type_UMAP.png`; alt: WNN UMAP of the two demo GEM wells, colored by cluster and cell type]

## Plot objects

Each plot saved by the standard plotting helper also has an `.rds` copy under `plot_objects/`, unless plot-object saving was disabled. This lets you reopen a plot in R without repeating the analysis. For example:

``` {.r filename="R"}
plot_file <- file.path(
  targets::tar_config_get("store"),
  "plot_objects/immune_human_2x/8_multimodal_QC/UMAPs/categorical",
  "WNN_harmony_SNN_cluster_cell_type.rds"
)
p <- readRDS(plot_file)
p
```

For a ggplot object, edit it with the usual ggplot2 functions and save a separate copy:

``` {.r filename="R"}
p <- p + ggplot2::labs(title = "My integrated cell types")
ggplot2::ggsave("my_cell_types.png", p, width = 10, height = 8)
```

Some outputs are composite plots rather than ordinary ggplot objects and need their own editing methods. Keep custom exports separate from pipeline outputs, which can be overwritten on a rerun.

## Possible next steps

- To continue with the demo-aggregation, and explore other outputs, run `targets::tar_make()` without `names`.
- To get started with your own data, please continue at [Plan your analysis](main_overview.md).
- To diagnose a failed or unexpectedly stale target, use [Troubleshooting](troubleshooting.md).


## Part: Analyze your own data


<!-- source: website/main_overview.md -->

# Plan your analysis

Before configuring anything, please check that you have the following inputs available:

**cellranger-arc count directories:** The output directories produced by `cellranger-arc count` containing the `outs/`-folder with `summary.csv`, `filtered_feature_bc_matrix.h5`, `atac_fragments.tsv.gz`, its `.tbi` index, and `per_barcode_metrics.csv` inside. All count directories must have been generated using the same reference if you want to combine later in the pipeline.

**donor metadata TSV** with one row per donor and the phenotypes or covariates you will use; see the [donor metadata table](reference_donor_metadata.md).

**VCF files for demultiplexing by genotype (optional):** If you have multiplexed several donors on one or more GEM-wells, and you wish to demultiplex them using reference genotypes, prepare a VCF file containing the donors in each pool, following the [Vireo genotype-input documentation](https://vireosnp.readthedocs.io/en/stable/manual.html). This functionality also requires `atac_possorted_bam.bam` in the cellranger-arc count directories.

**CellBender H5 files (optional)**: If you wish the pipeline to use gene-expression data that has been filtered for ambient RNA, run [CellBender remove-background](https://cellbender.readthedocs.io/en/latest/usage/) beforehand and supply the resulting H5 file in the GEM well table as described later.

**Compute setup:** If your dataset is large, it is highly recommended to run the pipeline on a compute cluster with a job-scheduler available. See [Choose where the analysis runs](performance_distributed_computing.md) for more info.

When these things are in order, please continue to [Run your own analysis](main_running.md#steps).


<!-- source: website/main_running.md -->

# Run your own analysis

After completing the demo, replace its configuration and work through **Configure → Run pipeline → Review plots**, repeating each step as needed. The plots explain what to inspect and which settings to revise.

Run commands from the repository-root R session. Replace `my_GEM_well` with your configured identifiers. `<store>` means the folder selected in `_targets.yaml`, which is normally `outputs` in the root of the repository.

## 1. Pre-process the cellranger-arc count dirs (GEM wells) {#steps}

**Initial configuration:**

- Add entries for each GEM-well (cellranger-arc count dir) you want to process inside `cfg_GEM_wells.tsv`. See [GEM well table](reference_GEM_wells.md) for more information

- Locate the `template_aggregation` inside `cfg_aggregations.yaml` and replace the GEM_well_ID-placeholders with the GEM-well IDs that you just added - these GEM-wells are now officially part of this aggregation. Finish by renaming the aggregation-entry using a short, descriptive name. Throughout the rest of this page, replace the \<my_aggregation\>-placeholders with your custom name.

**Run pipeline:**

``` {.r filename="R"}
targets::tar_make(
  names = tidyselect::ends_with("1_pre_aggregation_QC.<my_aggregation>")
)
```

**Review plots** (examples: [checkpoint 1](gallery.md#1-pre-aggregation-qc)):

``` text
<store>/plots/my_aggregation/1_pre_aggregation_QC/
├── per_aggregation_GEM_well_QC_comparisons/
├── aggregation_excluded_cellranger_only_barcodes_by_type_upset.png
├── aggregation_excluded_barcodes_by_type_upset.png
├── nuclei_per_donor_id_bars.png
└── cell_retention_flow_plot.png
```

For exclusion overlaps within an individual GEM well, you can also request its plots:

``` {.r filename="R"}
targets::tar_make(
  names = tidyselect::ends_with("1_pre_aggregation_QC.<my_GEM_well_ID>")
)
```

``` text
<store>/plots/my_GEM_well/1_pre_aggregation_QC/
├── excluded_barcodes_by_type_upset.png
└── excluded_cellranger_only_barcodes_by_type_upset.png
```

## 2. Normalize, reduce dimensions, and batch-correct the gene-expression data

**Initial configuration:**

- Start with default settings and revise as described in the plots if necessary.

**Run pipeline:**

``` {.r filename="R"}
targets::tar_make(
  names = tidyselect::ends_with("2_GEX_PCA_QC.my_aggregation")
)
```

**Review plots** (examples: [checkpoint 2](gallery.md#2-gex-pca-qc)):

``` text
<store>/plots/my_aggregation/2_GEX_PCA_QC/
├── variable_feature_plot.png
├── VizDimLoadings_plots/
├── PCA_singular_values_elbow_plot.png
├── PCA_embedding_sdev_plot.png
└── PCA_metadata_association_barplots/
```

## 3. Cluster and label cell types with GEX

**Initial configuration:**

Set [`aggregation_GEX_marker_genes`](parameters.html#aggregation_GEX_marker_genes), [`aggregation_categorical_vars`](parameters.html#aggregation_categorical_vars) & [`aggregation_continuous_vars`](parameters.html#aggregation_continuous_vars)

``` {.r filename="R"}
targets::tar_make(
  names = tidyselect::ends_with("3_GEX_QC.my_aggregation")
)
```

**Review plots** (examples: [checkpoint 3](gallery.md#3-gex-qc)):

``` text
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

## 4. Call peaks, inspect ATAC quality, and filter nuclei

**Initial configuration:**

- Start with default settings and revise as described in the plots if necessary.

**Run pipeline:**

``` {.r filename="R"}
targets::tar_make(
  names = tidyselect::ends_with(c(
    "4_peak_QC.my_aggregation",
    "5_pre_LSI_QC.my_aggregation"
  ))
)
```

**Review plots** (examples: [checkpoint 4](gallery.md#4-peak-qc) and [checkpoint 5](gallery.md#5-pre-lsi-qc)):

``` text
<store>/plots/my_aggregation/
├── 4_peak_QC/
│   ├── peaks_QC_violins_plot/
│   └── peaks_similarity_tiles_plot.png
└── 5_pre_LSI_QC/
    ├── QC_excluded_upset_plot.png
    └── cell_retention_flow_plot.png
```

## 5. Normalize, reduce dimensions and batch-correct the ATAC-data

Computes LSI on the retained nuclei with optional Harmony correction.

**Initial configuration:**

- Start with default settings and revise as described in the plots if necessary.

**Run pipeline:**

``` {.r filename="R"}
targets::tar_make(
  names = tidyselect::ends_with("6_ATAC_LSI_QC.my_aggregation")
)
```

**Review plots** (examples: [checkpoint 6](gallery.md#6-atac-lsi-qc)):

``` text
<store>/plots/my_aggregation/6_ATAC_LSI_QC/
├── LSI_singular_values_elbow_plot.png
├── LSI_embedding_sdev_plot.png
├── VizDimLoadings_plots/
└── LSI_metadata_association_barplots/
```

## 6. Generate ATAC-clusters and motif accessibility

**Initial configuration:**

- `Set aggregation_ATAC_marker_TFs`

**Run pipeline:**

``` {.r filename="R"}
targets::tar_make(
  names = tidyselect::ends_with("7_ATAC_QC.my_aggregation")
)
```

**Review plots** (examples: [checkpoint 7](gallery.md#7-atac-qc)):

``` text
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

## 7. Integrate GEX and ATAC

**Initial configuration:**

- Start with default settings and revise as described in the plots if necessary.

**Run pipeline:**

``` {.r filename="R"}
targets::tar_make(
  names = tidyselect::ends_with("8_multimodal_QC.my_aggregation")
)
```

**Review plots** (examples: [checkpoint 8](gallery.md#8-multimodal-qc)):

``` text
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

The final object is `multimodal_Seurat_object.8_multimodal_QC.my_aggregation`.

## 8. Choose your next analysis!

Continue to one of the three modules:

- [Differential analyses](downstream_differential_analyses.md), if you have many donors and a condition of interest

- [Genetic enrichment](downstream_genetic_enrichment.md), if you are interested in pinpointing cell-type-level genetic enrichment

- [Peak–gene correlation](downstream_peak_gene_correlation.md), to test associations between peak accessibility and gene expression.


## Part: Add an optional analysis


<!-- source: website/downstream_differential_analyses.md -->

# Differential analyses



## When to use this module

Use this module to ask how cell-type proportions, gene expression, or chromatin accessibility differ with a condition or donor phenotype. It requires biological replication across donors or samples. Molecular measurements are combined into **pseudobulks**: counts summarized for each cell type within a donor or sample.

The module does not create biological replication. The donor structure, covariates, design formula, and contrasts must be defensible for the intended analysis before the workflow is run.

See the [example plots](gallery.md#differential-analyses) for representative diagnostics and the [implementation graph](implementation/implementation_differential_analyses.html) for target structure.

## Prerequisites

Before enabling the module, confirm that:

- you have reviewed the main aggregation and its cell-type annotations;
- WNN cell-type metadata and GEX and ATAC pseudobulk matrices are available;
- the donor metadata contains one unique row per `donor_id` and every variable used in a model;
- model variables are donor- or pseudobulk-sample-level variables, not duplicated cell-level measurements; and
- the number and distribution of donors support the specified design and contrasts.

Use [`differential_analyses_extended_donor_id_metadata_tsv`](parameters.html#differential_analyses_extended_donor_id_metadata_tsv) when the modelling table needs variables beyond the aggregation's normal donor metadata. It must retain the same unique `donor_id` key.

## Outputs

Choose the output that matches your question:

| Question | Output family |
|------------------------------------|------------------------------------|
| Do cell-type proportions differ? | `cell_type_composition` |
| Which genes change expression? | `gene_expression` |
| Which peaks change accessibility? | `chromatin_accessibility` |
| Which motif families change accessibility? | `motif_family_accessibility` (JASPAR) |
| Which regulators show altered expression-based activity? | `transcription_factor_activity` (CollecTRI) |

The module also produces model diagnostics, comparisons across modalities, and gene-set tests for Hallmark and Reactome pathways. Motif-family accessibility summarizes ATAC evidence; transcription-factor activity is inferred from gene expression using CollecTRI. Interpret each in the context of its measurement.

Plot directories use these descriptive family names below `plots/<aggregation>/differential_analyses/`. Gene-set plots appear under `gene_expression/gene_set_enrichment/Hallmark/enrichment_plots/<model>/` or the corresponding `Reactome` directory. Volcano outputs use `<family>/volcano_plots/<model>/`; saved plot targets omit redundant `_file` and `_files` suffixes. Renaming targets creates new cache entries and output paths on the next run; existing output directories are not migrated.

See the [differential analyses methods](implementation/methods_differential_analyses.html) for activity inference, motif-family definitions, gene-set testing and the fixed and configurable settings.

## Configure

Add [`modules`](parameters.html#modules) to the existing aggregation entry, keeping its input and marker settings:

``` {.yaml filename="cfg_aggregations.yaml"}
your_aggregation:
  modules: [differential_analyses]
```

Then create a matching row directly in `configuration/cfg_module_differential_analyses.yaml`.

``` {.yaml filename="configuration/cfg_module_differential_analyses.yaml"}
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

Both branches use named models, donor eligibility checks and named contrasts. Abundance models use [`differential_analyses_cell_type_composition_models`](parameters.html#differential_analyses_cell_type_composition_models); feature models use [`differential_analyses_pseudobulk_models`](parameters.html#differential_analyses_pseudobulk_models). Omitting abundance models disables that branch. The former aggregation-wide cell-type composition formula, phenotype and colour settings have been replaced by fields inside each named model.

For mixed tissues, set `GEM_well_IDs` inside an abundance model to define its population, for example the six left-ventricle wells. Optional `donor_ids` can further restrict donors in either branch. Donors with missing model metadata or no selected samples are excluded and recorded in model-specific cohort TSVs. Feature cohorts also report retained pseudobulk sample counts and depth-filter exclusions.

By default, abundance models test all cell-type labels observed in their eligible population, including unassigned labels. Optional `cell_types_to_test` restricts the response cell types **without changing the denominator**: every retained nucleus in the selected wells contributes to its donor's total. Zero donor–cell-type counts remain in the analysis. In contrast, feature-model `cell_type_subset` selects the cells represented by the pseudobulks. Feature matrices already pool wells within donors and cell types, so they cannot support a late `GEM_well_IDs` filter; the module rejects that field for feature models.

Abundance models fit a separate fixed-effects beta-binomial logit model per cell type. Use a one-sided predictor formula, `formula: ~ ...`, and named `contrast_specs_vec`; two-sided formulas, random effects and custom design or contrast functions are rejected in this branch. Contrast tables report log-odds effects, Wald uncertainty, donor counts and BH FDR across tested cell types within each model and contrast, and failed fits are marked non-estimable. The response construction, test and every fixed setting are documented in [Differential analyses methods](implementation/methods_differential_analyses.html#cell-type-composition).

The module selection below requests both configured abundance and pseudobulk outputs. Formula terms and contrast coefficients must match columns produced by the model matrix. The two branches retain their distinct response construction and fitting methods; sharing configuration does not make their effect estimates interchangeable.

The model example assumes `condition` distinguishes treated and control donors. Check which group is the reference and what each model coefficient represents before using `conditiontreated` as a contrast. Replace the example formula and contrast to match your study.

## Run

Preview the selected module outputs before running them:

``` {.r filename="R"}
targets::tar_manifest(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:differential_analyses")
  ) & tidyselect::ends_with(".your_aggregation"),
  callr_function = NULL
)[, c("name", "description")]
```

Then run the same selection:

``` {.r filename="R"}
targets::tar_make(
  names = targets::tar_described_as(
    tidyselect::contains("checkpoint:differential_analyses")
  ) & tidyselect::ends_with(".your_aggregation")
)
```

## Review

Open the configured model outputs listed above; the [output gallery](gallery.md#differential-analyses) shows one example per plot. Interpretation and method details are included in the plot subtitles and captions.

Runtime depends on donors, cell types, models, contrasts, and gene-set analyses. Use [Troubleshooting](troubleshooting.md) if a formula, contrast, or metadata join fails.

## Parameter reference

The OLINK and bulk-RNA path fields are reserved optional integration inputs and are not consumed by the current public differential-analysis selection. Leave them `NULL` unless the corresponding integration is implemented in your downstream workflow.

[Open the searchable parameter browser](parameters.html#workflow=differential_analyses).

The public demos leave this module disabled. Comparing one healthy PBMC donor with one lymphoma lymph-node donor cannot separate condition, donor, and tissue effects. Configure differential analyses for a design with biological replication.


<!-- source: website/downstream_genetic_enrichment.md -->

# Genetic enrichment



## When to use this module

Use this module to ask which cell types or nuclei have accessible regions overlapping genetic evidence for a human trait. It connects fine-mapped GWAS variants to ATAC peaks, calculates accessibility-based enrichment, and uses [`SCAVENGE`](https://github.com/sankaranlab/SCAVENGE) to summarize trait relevance across related nuclei.

A **credible set** contains candidate causal variants at a GWAS locus, with probabilities from fine-mapping. Enrichment helps prioritize cellular contexts; it does not by itself identify a causal cell type, gene, or mechanism.

See the [example plots](gallery.md#genetic-enrichment) for representative results and the [implementation graph](implementation/implementation_genetic_enrichment.html) for upstream ATAC and WNN dependencies.

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

Add [`modules`](parameters.html#modules) to the existing human aggregation entry, keeping its other settings:

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

`sourceId` values beginning with `GCST` use the pinned Open Targets datasets. Every other value is a local Parquet filename, resolved from the project root and tracked as a file target. Local files must satisfy the schema enforced by `validate_local_finemapped_GWAS_tibble()`; their study ID, fine-mapping method, build, credible-set probability, and provenance are read from the file rather than repeated in YAML.

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

Open the study-selection summaries, chromVAR summaries, SCAVENGE heatmaps and locus-contribution plots produced for your configured studies. The [output gallery](gallery.md#genetic-enrichment) shows one example per plot; interpretation belongs to each plot.

Runtime and disk use grow with studies, cells, graph representations, permutations, and attributed loci.

## Parameter reference

[Open the searchable parameter browser](parameters.html#workflow=genetic_enrichment).

<details>
<summary>Show the public <code>immune_human_2x</code> example</summary>

[Generated Quarto chunk omitted: `emit_yaml_entry(module_config_file, "immune_human_2x")`]

</details>


<!-- source: website/downstream_peak_gene_correlation.md -->

# Peak–gene correlation

Run this optional module after accepting the final WNN cell set. It relates ATAC accessibility to RNA expression within broad GEX-derived cell types, using the retained WNN nuclei; it is no longer part of checkpoint 8.

## Configure

Add `peak_gene_correlation` to the aggregation's existing [`modules`](parameters.html#modules) list in `cfg_aggregations.yaml`, and add a matching row in `cfg_module_peak_gene_correlation.yaml` in the selected configuration directory:

```yaml
my_aggregation:
  peak_gene_correlation_top_links_per_cell_group: 3
```

The top-link count controls the number of detail figures per cell type. [`peak_gene_correlation_filter`](parameters.html#peak_gene_correlation_filter) selects `lenient` (default), `moderate`, or `strict` measurement-support filtering. Disabled aggregations contribute no module targets.

## Run

```r
targets::tar_make(
  names = tidyselect::ends_with(".peak_gene_correlation.my_aggregation")
)
```

## Review

The [output gallery](gallery.md#peak-gene-correlation) shows one example per plot. All module targets have description tag `[checkpoint:peak_gene_correlation]`. Their paths are `<store>/plots/my_aggregation/peak_gene_correlation/`; file exports use the corresponding `files` directory, with any modality suffixes as deeper subdirectories. For example, read selected links with:

```r
targets::tar_read(
  peak_gene_correlation_links_tibble.WNN.peak_gene_correlation.my_aggregation
)
```

## Model and scope

Peak–gene links are candidate regulatory relationships. Cells are aggregated within donor and ATAC-defined state using WNN cell-type annotations, without reusing a cell across aggregates. Measurement-support filters select hypotheses before fitting; the default requires shared support from at least two donors, and the strict preset requires three. Review `filter_retention_plot` first.

The hierarchical analysis fits a mean peak effect with donor-specific slope variation, donor intercepts and RNA/ATAC depth adjustment. It reports Kenward–Roger p-values, BH FDR and numerical reliability diagnostics. Top-link figures rank positive, estimable nonpromoter slopes by p-value without a significance cutoff, so appearing in a figure is not evidence of significance. They combine focal-cell-type coverage, gene context, hierarchical evidence and adjusted aggregate scatterplots.

The existing HC3 correlation summaries remain a separate conditional analysis. Neither analysis establishes causal regulation; numerical reference parity does not establish statistical calibration across datasets. See the [peak–gene correlation methods](implementation/methods_peak_gene_correlation.html) for inference limits, support thresholds and the fixed and configurable settings, and the [implementation graph](implementation/implementation_peak_gene_correlation.html) for target structure.


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



A **worker** is an R process that runs an analysis task. A **controller** starts and manages those workers, either on your machine or through a cluster scheduler. multiomeR uses `crew` for this. The committed configuration runs the demo unchanged on a machine that meets the [system requirements](demo_installation.md#system-requirements); adjust it as described below before running on a smaller machine or a scheduler.

Use local workers on a suitable workstation. On a shared cluster, ask your support team which scheduler, account, and resource limits to use. The [targets distributed-computing guide](https://books.ropensci.org/targets/crew.html) explains the general setup; this page covers multiomeR's configuration file.

## Local execution

A fresh clone includes a local `crew_controllers.R` sized for a 16-CPU, 256-GB workstation, with four light workers and two heavy workers. A machine near the 60-GB minimum should reduce concurrency to one heavy worker and should not run several memory-intensive targets simultaneously. Edit the worker counts and resource tiers directly in `crew_controllers.R`, keeping controller names identical between `controller_list` and `controller_resources_tibble`.

After changing the file, restart R or reload the project runtime explicitly:

```{.r filename="R"}
load_project_runtime(force = TRUE)
```

Rebuild a narrow manifest selection before starting the data run to validate the controller contract.

## Scheduler execution

For SLURM, PBS, SGE, or LSF, replace the local controllers with the corresponding `crew.cluster` controllers. The commented SLURM section in `crew_controllers.R` shows the expected shape.

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

Find the target name and first error message in the run output. Fix that cause, then rerun the same selection: `targets` can reuse completed work. Keep the store intact, since it also contains the records needed to diagnose the failure.

| What happened? | Start here |
|---|---|
| Setup fails before any data processing | [The manifest does not build](#the-manifest-does-not-build) |
| A target reports an error | [A run reports errored targets](#a-run-reports-errored-targets) |
| A completed result needs rebuilding | [A target is unexpectedly outdated](#a-target-is-unexpectedly-outdated) |
| Workers do not start or are killed | [Controller and scheduler failures](#controller-and-scheduler-failures) |

Run these commands in the repository's Pixi R session. For general debugging techniques beyond the project helpers below, see the [targets debugging guide](https://books.ropensci.org/targets/debugging.html).

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

Correct the field or file named in the error, then repeat the manifest check before running the analysis.

## A run reports errored targets

List the errors with the project's helper. It groups repeated failures so you can start with their common cause:

```{.r filename="R"}
list_distinct_errored_targets()
```

For commands and stored tracebacks matching a target, aggregation, or module:

```{.r filename="R"}
list_distinct_errored_targets_w_tracebacks(
  target_name_pattern = "your_target_or_aggregation"
)
```

Some targets run separately for multiple groups; these runs are called **branches**. Copy the full target or branch name from the error listing to inspect its saved workspace:

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

After fixing the cause, rerun the same target selection. Successful upstream results remain cached.

```{.r filename="R"}
targets::tar_make(
  names = tidyselect::matches("your_target.*your_aggregation")
)
```

Use an unqualified `targets::tar_make()` only when every active aggregation and enabled module is intentionally in scope. Avoid deleting the target store as a debugging step: it removes evidence and forces unrelated recomputation.


## Part: Reference


<!-- source: website/gallery.md -->

# Output gallery



Each plot target has one example here, taken from `mixed_human_31x`: 31 public GEM wells from 10x Genomics and ENCODE covering heart, blood, pancreas, liver, colon, lung and cerebellum, with all optional analyses enabled. Its configuration is included as an inactive example; the raw data are public, but the reprocessed Cell Ranger ARC inputs are not supplied. A card shows one of the files its target saves, the target's description and its name; select a preview to enlarge it. [Run your own analysis](main_running.md#steps) explains when to review each checkpoint.

[Generated Quarto chunk omitted: `render_output_gallery()`]


<!-- source: website/review_outputs.md -->

# Output files and metadata

The [running guide](main_running.md#steps) lists the files to open after each checkpoint. Read each plot's subtitle and caption for interpretation and method details.

## Output folders and metric selection {#output-folders}

Plots live under `<store>/plots/<scope>/<checkpoint>/`. The store comes from `_targets.yaml`; scope is a GEM well or aggregation. Existing files are not a record of which targets are current: use `targets::tar_outdated()` before reviewing results after a configuration change.

`QC_metric_manifest.tsv` selects metrics, display labels and plotting quantiles. `do_plot = FALSE` hides a metric; plotting quantiles change the displayed range, not the cells retained by the pipeline. Filtering is configured separately in `cfg_GEM_wells.tsv` and `cfg_aggregations.yaml`.

GEX, ATAC and WNN parameter sweeps belong to checkpoints 3, 7 and 8. The GEX compatibility object belongs to checkpoint 3 and the multimodal object to checkpoint 8. Peak–gene correlation is a separate optional module.

## Cell retention tables {#cell-retention}

Append `.my_aggregation` to these target names when reading them with `targets::tar_read()`:

| Target | Checkpoint |
|---|---|
| `cell_retention_tibble.GEX_input` | 1 |
| `cell_retention_tibble.GEX` | 3 |
| `cell_retention_tibble.ATAC_input` | 5 |
| `cell_retention_tibble.ATAC` | 7 |
| `cell_retention_tibble.WNN` | 8 |

Tables accumulate stages and retain wells with zero surviving cells. Rows identify an exclusion action and its order, input count, excluded count and retained count.

## Cluster numbering {#cluster-numbering}

Leiden clusters receive size-ranked IDs before minimum-size filtering. Filtering does not renumber survivors, so IDs can have gaps.

## Cluster annotation {#cluster-annotation}

Cluster labelling uses adjusted-score advantage over matched random controls. Marker lists accept unsigned genes, positive `+` suffixes and negative `-` suffixes. The same GEX control reference annotates GEX, ATAC and WNN clusters, and each cluster receives an `Assigned` or `Unassigned` status; an unassigned cluster retains its candidate label and the reason in the diagnostics, and no cells are removed. Set [`aggregation_cluster_annotation_min_advantage`](parameters.html#aggregation_cluster_annotation_min_advantage) to the required lead over both background and the next-best label; raising it only withdraws assignments. The scoring rule, control construction and all fixed constants are documented in [Cell-type annotation and motif accessibility](implementation/methods_annotation_and_motifs.html). The former module-score labeller and its `aggregation_allow_multiple_cell_types` and `aggregation_cluster_annotation_method` settings have been removed; delete them from older configuration files.

Annotation exports include `clusters.tsv`, `marker_evidence.tsv`, `control_gene_matching.tsv` and, where eligible, `GEM_well_agreement.tsv`. The `cluster_UCell_evidence` targets cache scoring separately from the threshold-dependent annotation. `marker_set_UCell_summary.3_GEX_QC` exports `marker_sets.tsv` and `method.txt` before GEX doublet filtering.

## Further methods

See [algorithm validation](implementation/algorithm_validation.html) for reference comparisons and deviations, and the [peak–gene module](downstream_peak_gene_correlation.md) for its model, support filters and output paths.


<!-- source: website/reference_GEM_wells.md -->

# GEM well table



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

[Generated Quarto chunk omitted: `emit_GEM_well_demo_table( GEM_well_config_file = "website/data/demo_GEM_wells.tsv", dictionary_file = "website/data/G...`]


<!-- source: website/reference_donor_metadata.md -->

# Donor metadata table

The donor metadata table is a TSV with one unique row per `donor_id`. Each aggregation points to one such file through [`aggregation_donor_id_metadata_tsv`](parameters.html#aggregation_donor_id_metadata_tsv) in the [aggregation configuration](reference_aggregations.md). It is created in step 2 of [Run your own analysis](main_running.md#steps).

## Minimal table

``` {.text filename="donor_metadata.tsv"}
donor_id	condition
donor_1	control
```

## Matching donors to nuclei

Every nucleus receives a `donor_id` from its GEM well: the configured `GEM_well_donor_id` for a non-multiplexed well, or a genotype-based assignment for a well with a configured VCF. Each of those IDs must appear exactly once in this table; see the [GEM well table](reference_GEM_wells.md#fill-in-a-row).

## Which variables belong here

Put donor-specific phenotypes and covariates in this table, for example condition, age, or sex. Put library-, run-, or batch-specific variables in the GEM well table with a `GEM_well_` prefix. Apart from their key columns, the two tables must not reuse column names. Metadata file contents are validated when their targets run, not when the manifest is built.

## Extended table for differential analyses

The [differential analyses](downstream_differential_analyses.md) module can read additional donor-level model variables from a second table given in [`differential_analyses_extended_donor_id_metadata_tsv`](parameters.html#differential_analyses_extended_donor_id_metadata_tsv). It must keep the same unique `donor_id` key. If it is not set, the module inherits the aggregation's donor table.


<!-- source: website/reference_aggregations.md -->

# Aggregation configuration



`cfg_aggregations.yaml` has one top-level entry per aggregation: a joint GEX, ATAC, and WNN analysis of one or more GEM wells. The committed file enables the two human GEM wells in `immune_human_2x`, with optional modules disabled. The mouse and ENCODE validation examples, and the `mixed_human_31x` aggregation behind the [output gallery](gallery.md), are inactive by default. Edit the file directly; the demo entries can stay as worked examples. This page describes the entry structure and lists every parameter. When to set each parameter, and how to review the effect, is given step by step in [Run your own analysis](main_running.md#steps).

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

Every other parameter has a default from `cfg_pipeline_parameters.tsv`, listed in the [parameter reference](#parameter-reference) below. Add a parameter to the entry only when you want to change its default.

## Required keys

- [`aggregation_GEM_well_IDs`](parameters.html#aggregation_GEM_well_IDs): the `GEM_well_ID` values to combine. Each must be an active row of the [GEM well table](reference_GEM_wells.md), and all must use the same Cell Ranger reference.
- [`aggregation_donor_id_metadata_tsv`](parameters.html#aggregation_donor_id_metadata_tsv): the [donor metadata table](reference_donor_metadata.md) for these GEM wells.
- [`aggregation_GEX_marker_genes`](parameters.html#aggregation_GEX_marker_genes): a named list of marker genes per expected cell type, used for cluster annotation and marker plots.
- [`is_active`](parameters.html#is_active): whether targets are constructed for the aggregation. Deactivate aggregations you are not ready to run before an unqualified `targets::tar_make()`.

## Marker genes and transcription factors

Replace the placeholder genes with symbols appropriate for the tissue and reference. A gene listed without a suffix or with a `+` suffix is a positive marker; a `-` suffix marks a gene that should be absent. The optional [`aggregation_ATAC_marker_TFs`](parameters.html#aggregation_ATAC_marker_TFs) list names transcription factors per cell type for the motif-activity plots. How the annotation uses these lists is described in [Output files and metadata](review_outputs.md#cluster-annotation).

## QC filters after peak calling

[`aggregation_QC_exclude_list_combined_object`](parameters.html#aggregation_QC_exclude_list_combined_object) lists dplyr filter expressions applied to the peak-based ATAC metrics of the combined object, for example:

``` {.yaml filename="cfg_aggregations.yaml"}
aggregation_QC_exclude_list_combined_object:
  - nCount_ATAC < 1000
  - atac_peak_counts_frac < 0.1
  - atac_peak_counts_blacklist_frac > 0.01
```

Omit it or set it to `null` until step 4 of [Run your own analysis](main_running.md#steps) has shown the distributions.

## Optional modules

Omit [`modules`](parameters.html#modules) for the first run. After reviewing the main results, enable an optional analysis by listing its name and adding a matching entry for the aggregation in the module's own configuration file:

``` {.yaml filename="cfg_aggregations.yaml"}
your_aggregation:
  modules: [differential_analyses]
```

See [Differential analyses](downstream_differential_analyses.md) and [Genetic enrichment](downstream_genetic_enrichment.md) for the module entries and their parameters.

## Parameter reference {#parameter-reference}

The [standalone parameter browser](parameters.html) is generated from `cfg_pipeline_parameters.tsv`, using a shared snapshot of the public runtime defaults and validation schema. Choose the main workflow or an optional module, then search by name or purpose. Cards are grouped by whether a value is required, defaulted, or optional. Defaults are visible beside each parameter; open a row for its type and example. See the [committed example](https://github.com/koefoeden/multiomeR/blob/main/configuration/cfg_aggregations.yaml) for a complete configuration.

<details>

<summary>Show the public <code>immune_human_2x</code> example</summary>

[Generated Quarto chunk omitted: `emit_yaml_entry(aggregations_config_file, "immune_human_2x")`]

</details>


# Book: multiomeR Implementation


<!-- source: website/implementation/index.md -->

# Introduction

Use this book to trace a result back to its code or change how multiomeR works. For installation, configuration, execution, and output inspection, start with the [user manual](../). You do not need to read this book to run the demo.

Use this book when you need to trace a configuration value into mapped targets, understand how the simplified graph views relate to the real `{targets}` graph, or decide where an implementation change belongs.

## Where to start

For a first implementation pass:

1. Read [Reading the graph views](graph_methodology.md) and follow its configuration-to-target trace.
2. Open the [main pipeline](implementation_main.md) graph for the modality or checkpoint you plan to change.
3. Use [Implementation conventions](implementation_conventions.md) to understand the relevant manifest, mapping, symbol, tag, and runtime contracts.
4. Read [Background and design philosophy](background_philosophy.md) when you need the rationale for the editable-workflow design.

The [differential analyses](implementation_differential_analyses.md), [genetic enrichment](implementation_genetic_enrichment.md) and [peak–gene correlation](implementation_peak_gene_correlation.md) chapters cover the optional module graphs.

The **Methods and parameters** chapters, from [Preprocessing and nucleus QC](methods_preprocessing_and_QC.md) to [Genetic enrichment](methods_genetic_enrichment.md), describe each stage and list every fixed value and every configurable setting, using the table layout defined in [Implementation conventions](implementation_conventions.md#methods-and-parameter-tables). Read them when you need the exact behaviour behind a result, or when deciding whether a change is a configuration edit or a code edit.

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

`[checkpoint:<name>]` marks targets selectable with `targets::tar_described_as()`. The eight numbered main-pipeline groups are listed in `QC_checkpoint_manifest.tsv`; optional module groups remain unnumbered. UMAP parameter sweeps belong to their modality's numbered checkpoint, and compatibility objects belong to GEX checkpoint 3 or multimodal checkpoint 8. Selection matches description substrings; include the closing `]` to match a complete checkpoint tag. Dependencies still come from the target commands. [Run your own analysis](../main_running.html#steps) explains each boundary; acceptance criteria depend on the study.

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

Module target files also need aggregation-specific references to main-pipeline targets. For this, `add_aggregation_target_syms()` creates one symbol per row, suffixed by the aggregation name. These columns are named like the target they replace rather than with `*_syms`, because each cell is a single symbol rather than a list.

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


## Part: Background


<!-- source: website/implementation/background_philosophy.md -->

# Why an editable workflow?

multiomeR keeps the analysis steps in an editable repository. Configuration covers common choices such as inputs, markers, dimensions, and models; R helpers and target definitions are available when a study needs a change beyond those settings. This flexibility also means that users must review which methods and assumptions fit their data.

## Reuse completed work

`targets` records dependencies between results so that a change can rebuild the affected parts of an analysis. Independent tasks can run concurrently when worker capacity permits. This is useful when processing several GEM wells or repeating analyses with revised settings. The [targets manual](https://books.ropensci.org/targets/) explains the execution model and its limits.

## Keep large matrices on disk

BPCells provides disk-backed matrices and streaming operations that can reduce the need to hold full matrices in memory. Some analysis steps still need substantial RAM, and performance depends on the data, storage, and available workers. See the [BPCells documentation](https://bnprks.github.io/BPCells/) for its matrix operations and [Performance and scaling](../performance_overview.html) for multiomeR examples.

## Keep the analysis inspectable

Separate targets make intermediate tables, matrices, and files available for inspection. Seurat/Signac exports provide another way to explore completed results. The [implementation conventions](implementation_conventions.md) explain where to change parameters, helpers, and target definitions; the [user manual](../) covers running an existing configuration.


## Part: Target graph views


<!-- source: website/implementation/implementation_main.md -->

# Main pipeline



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



`module_differential_analyses/targets.R` filters aggregations that enabled the module, joins their module config, attaches symbols for accepted WNN metadata and pseudobulk inputs, and maps the composition, pseudobulk, gene-set enrichment, and cross-modality target fragments. The generic pseudobulk model family is instantiated for gene expression, chromatin accessibility, motif-family accessibility, and expression-derived CollecTRI activity (transcription-factor activity). transcription-factor activity first converts filtered, normalized GEX pseudobulks to signed ULM scores and then reuses the same model and contrast machinery.

The cross-modality fragment creates a CollecTRI-to-JASPAR family crosswalk, a detailed regulator-level table containing transcription-factor activity, motif-family accessibility, and TF-expression results, a family-level comparison table, a contrast-level concordance summary, and its plot. CollecTRI complexes remain intact in transcription-factor activity; complex-member mappings are introduced only by the comparison crosswalk.

The graph below is an orientation view. Inspect `setup_and_cell_type_composition_targets.R`, `pseudobulk_differential_targets.R`, `gene_set_enrichment_targets.R`, and `cross_modality_targets.R` for the complete model and plotting commands. The methods, model routes and every fixed or configurable setting are listed in [Differential analyses methods](methods_differential_analyses.md). The user-facing prerequisites and module selector are documented in [Differential analyses](../downstream_differential_analyses.html).

[Mermaid graph omitted; source: `website/figures/human_curated/differential_analyses_v2.mmd`]


<!-- source: website/implementation/implementation_genetic_enrichment.md -->

# Genetic enrichment



`module_genetic_enrichment/targets.R` filters enabled human aggregations, resolves one configured Open Targets study set per aggregation, and attaches symbols for WNN metadata, graphs, embeddings, consensus peaks, chromVAR state, and ATAC fragments.

The main target fragments live in `setup_targets.R`, `gchromVAR_targets.R`, `SCAVENGE_graph_targets.R`, `SCAVENGE_group_targets.R`, `GWAS_chromVAR_cell_type_targets.R`, and `GWAS_chromVAR_contribution_targets.R`. The user-facing release, method-selection, and interpretation contracts are documented in [Genetic enrichment](../downstream_genetic_enrichment.html).

## Single-nucleus and graph-based enrichment

This view covers the configured GWAS inputs, single-nucleus enrichment state, graph propagation, and downstream trait summaries. Additional cell-type contribution and locus-attribution branches may be pruned from this compact orientation view; use the manifest for the complete graph.

The sparse SCAVENGE reimplementation, deliberate graph and permutation differences, and validation evidence are recorded in [Algorithmic implementations, deviations and validation](algorithm_validation.md#sparse-scavenge-propagation-and-significance).

[Mermaid graph omitted; source: `website/figures/human_curated/genetic_enrichment_single_nucleus_v2.mmd`]

The fixed thresholds and the configurable settings of this module are listed in [Genetic enrichment methods](methods_genetic_enrichment.md).

## Cell-type pseudobulk enrichment

This view covers the annotation-class pseudobulk deviations that are calculated separately from the nucleus-level results.

[Mermaid graph omitted; source: `website/figures/human_curated/genetic_enrichment_cell_type_absolute_effect_v2.mmd`]

## Cell-type contributions and locus attribution

This view covers the per-cell-type contribution and locus-attribution branches that are pruned from the single-nucleus view above.

[Mermaid graph omitted; source: `website/figures/human_curated/genetic_enrichment_cell_type_contributions_v2.mmd`]


<!-- source: website/implementation/implementation_peak_gene_correlation.md -->

# Peak–gene correlation



`module_peak_gene_correlation/targets.R` maps only opted-in aggregations and binds their existing WNN metadata, GEX and ATAC matrices, ATAC embeddings, fragments and reference annotations. `correlation_targets.R` owns the candidate-pair construction, donor–state pseudobulking, filtering, the conditional and hierarchical analyses, SuSiE prioritization, exports and plots.

Parameters use the `peak_gene_correlation` manifest scope and matching module YAML rows. Targets end in `.peak_gene_correlation.<aggregation>`, with `.WNN` before that suffix for intermediate results. Plot checkpoint tags use `peak_gene_correlation`, keeping this analysis outside the numbered QC selections.

The fixed thresholds and the configurable settings of this module are listed in [Peak–gene correlation methods](methods_peak_gene_correlation.md). The user-facing prerequisites and module selector are documented in [Peak–gene correlation](../downstream_peak_gene_correlation.html).

## Candidate pairs, pseudobulks and tests

This view covers the TSS table, candidate peak–gene pairs, the broad WNN cell groups, donor–state pseudobulk construction, measurement-support filtering, the per-chromosome analysis branches and the diagnostic and top-link outputs.

[Mermaid graph omitted; source: `website/figures/human_curated/peak_gene_correlation_v2.mmd`]


## Part: Methods and parameters


<!-- source: website/implementation/methods_preprocessing_and_QC.md -->

# Preprocessing and nucleus QC

This chapter describes per-GEM-well preprocessing and the successive nucleus filters up to the final WNN cell set, and lists every setting that determines them, using the layout defined in [Methods and parameter tables](implementation_conventions.md#methods-and-parameter-tables). The target structure is shown in the [main pipeline graph](implementation_main.md). GEM-well-level settings are columns of `cfg_GEM_wells.tsv`, described in [GEM well table](../reference_GEM_wells.html); aggregation-level settings are manifest parameters.

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

This chapter describes normalization and dimensional reduction of both modalities, peak definition, batch correction, weighted nearest-neighbour (WNN) integration, and the shared graph, clustering and UMAP steps. It lists every setting that determines them, using the layout defined in [Methods and parameter tables](implementation_conventions.md#methods-and-parameter-tables). The target structure is shown in the [main pipeline graph](implementation_main.md). Library versions are pinned by the Pixi environment: BPCells 0.3.1, igraph 2.3.0, harmony 2.0.2, uwot 0.2.4 and Seurat 5.5.0 at the time of writing.

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

GEX and ATAC counts are summed within each donor and annotation class in the main pipeline, so the module reuses the same pseudobulk targets as the compatibility export. Sample identifiers combine the class and the donor. Four feature matrices are tested:

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

The trait peak weights form a chromVAR annotation. Deviations and z-scores per nucleus are computed analytically with betterChromVAR on the ATAC chromVAR object reused from the main pipeline, with the GC-bias background described in [Cell-type annotation and motif accessibility](methods_annotation_and_motifs.md#motif-families-and-motif-accessibility).

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
