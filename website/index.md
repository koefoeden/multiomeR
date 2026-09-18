![](figures/multiomeR-logo.svg){fig-alt="multiomeR logo" width="100%"}

# Start here

## What is multiomeR?

multiomeR is a [targets](https://books.ropensci.org/targets/)-based workflow for processing and analyzing single-nucleus 10x Genomics Multiome ATAC + Gene Expression datasets. It is meant to be adapted to your own data, compute setup, and biological questions.

The workflow starts from `cellranger-arc count` outputs, processes gene-expression (GEX) and ATAC data, builds multimodal aggregations, and supports optional downstream modules for differential analyses and genetic enrichment for human datasets.

::: {.callout-warning title="Beta software"}
multiomeR is in beta and may introduce breaking changes between releases. The steps in [Run your own analysis](main_running.md#steps) show what to inspect before each stage, but acceptance criteria still depend on the tissue, study design, and intended use. Report problems or questions through [GitHub issues](https://github.com/koefoeden/multiomeR/issues).
:::

## Your first analysis

Start with the public demo: two human GEM wells with supplied configuration. You will install the software, run one joint analysis, and read its cell metadata and multimodal object. This gives you a working example before you choose settings for your own study.

You need basic R skills, a Linux terminal, and a machine with sufficient [memory and disk space](demo_installation.md#system-requirements). You do not need to know how to write a `targets` pipeline. Commands labelled **Bash** run in the terminal; commands labelled **R** run in the R session opened during installation. Run both from the repository folder unless stated otherwise.

## Find what you need

| If you want to... | Start here |
|------------------------------------|------------------------------------|
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

In this manual, a **GEM well** is one configured 10x library and output directory, an **aggregation** is a joint analysis of one or more GEM wells, and a **donor** is the individual identified by `donor_id`. One GEM well may contain multiple donors.

A **target** is a named result, such as a metadata table, matrix directory, or plot. Its **dependencies** are the inputs and earlier results needed to build it. You request the result you want; `targets` works out the order and reuses results that are up to date. The **store** is the folder where it keeps results and the records needed for reruns. For a small worked introduction, see the [targets walkthrough](https://books.ropensci.org/targets/walkthrough.html).

## Workflow at a glance

The **main pipeline** processes each GEM well, aggregates selected GEM wells, and builds multimodal RNA/ATAC outputs for clustering, cell typing, and WNN integration. Two optional modules extend completed aggregations with differential analyses or genetic enrichment.

![multiomeR workflow from Cell Ranger ARC GEM well outputs through per GEM well processing, aggregation-level GEX and ATAC analysis, WNN integration, and optional downstream modules](figures/multiomeR_overview_simplified.drawio.svg)

Continue to [Install and prepare the demo](demo_installation.md).