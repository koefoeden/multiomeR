![](figures/multiomeR-logo.svg){fig-alt="multiomeR logo" width="100%"}

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

The **main pipeline** processes each GEM well, aggregates selected GEM wells, and builds multimodal RNA/ATAC outputs for clustering, cell typing, and WNN integration. Two optional modules extend completed aggregations with differential analyses or genetic enrichment.

![multiomeR workflow from Cell Ranger ARC GEM well outputs through per GEM well processing, aggregation-level GEX and ATAC analysis, WNN integration, and optional downstream modules](figures/multiomeR_overview_simplified.drawio.svg)

Continue to [Install and prepare the demo](demo_installation.md).