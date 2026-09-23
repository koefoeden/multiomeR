![](figures/multiomeR-logo.svg){fig-alt="multiomeR logo" width="100%"}

# Start here

## What is multiomeR?

multiomeR is a [targets](https://books.ropensci.org/targets/)-based workflow for single-nucleus 10x Genomics Multiome ATAC + Gene Expression data. It starts from `cellranger-arc count` outputs and is meant to be adapted to your own data, compute setup, and biological questions.

::: {.callout-warning title="Beta software"}
multiomeR is in beta and may introduce breaking changes between releases. The [running guide](main_running.md#steps) shows what to inspect at each checkpoint, but acceptance criteria still depend on the tissue, study design, and intended use. Report problems or questions through [GitHub issues](https://github.com/koefoeden/multiomeR/issues).
:::

## Workflow at a glance

The **main pipeline** processes each GEM well, combines selected GEM wells into an aggregation, clusters and labels cell types in its gene-expression (GEX) and ATAC data, and integrates both modalities with weighted nearest neighbors (WNN). Three optional modules extend a completed aggregation with differential analyses, genetic enrichment for human traits, or peak–gene correlation.

![](figures/multiomeR_overview_simplified.drawio.svg){fig-alt="multiomeR workflow from cellranger-arc count outputs through per-GEM-well processing, aggregation-level GEX and ATAC analysis, and WNN integration to optional downstream modules"}

## How this manual is organized

- **Try the public demo:** [install multiomeR](demo_installation.md), then run and inspect a small example analysis.
- **Analyze your own data:** [check your inputs](main_overview.md), then configure, run, and review the main pipeline one checkpoint at a time.
- **Add an optional analysis:** run [differential analyses](downstream_differential_analyses.md), [genetic enrichment](downstream_genetic_enrichment.md), or [peak–gene correlation](downstream_peak_gene_correlation.md) on a completed aggregation.
- **Operation and scaling:** [run locally or on a scheduler](performance_distributed_computing.md), and [troubleshoot](troubleshooting.md) failed or outdated targets.
- **Reference:** browse the [output gallery](gallery.md) for an example of each plot, and look up output files, configuration tables, and methods.

## Your first analysis

Start with the public demo, even if you plan to analyze your own data. It runs one aggregation of two human GEM wells with supplied configuration and shows how to read its cell metadata and multimodal Seurat object.

You need basic R skills, a Linux terminal, and a machine with sufficient [memory and disk space](demo_installation.md#system-requirements). You do not need to know how to write a `targets` pipeline. Run blocks labeled **Bash** in the terminal and blocks labeled **R** in the R session opened during installation, both from the repository folder unless stated otherwise. A block labeled with a file name, such as `cfg_aggregations.yaml`, shows content for that file.

## Terms used in this manual

- A **GEM well** is one 10x Chromium chip channel and its `cellranger-arc count` output directory. It may contain nuclei from several **donors**, the individuals identified by `donor_id`.
- An **aggregation** is a joint analysis of one or more GEM wells.
- A **checkpoint** is one of the eight stages of the main pipeline, each ending with plots to review before you continue. Its name, such as `8_multimodal_QC`, appears in target names and plot folders.
- A **target** is a named result, such as a metadata table, matrix directory, or plot. You request the targets you want; `targets` builds the earlier results they depend on, in order, and reuses those that are up to date. Most target names end with their GEM well or aggregation, as in `multimodal_Seurat_object.8_multimodal_QC.immune_human_2x`.
- The **store** is the folder where `targets` keeps results and the records needed for reruns. Paths in this manual write it as `<store>`.

For a small worked introduction to `targets`, see the [targets walkthrough](https://books.ropensci.org/targets/walkthrough.html).

Continue to [Install and prepare the demo](demo_installation.md).
