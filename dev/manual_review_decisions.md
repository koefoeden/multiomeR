# Manual Review Decision Checklist

Temporary working checklist for reviewing subagent suggestions before editing the published manual.

Status key:

- `pending`: not reviewed yet
- `accept`: implement substantially as suggested
- `revise`: implement a modified version
- `reject`: do not implement
- `defer`: keep for later, outside the next docs pass

## Batch A: Correctness And Broken Links

| ID | Status | Page | Issue | Suggested direction |
|---|---|---|---|---|
| A1 | pending | `website/demo_installation.qmd` | Fresh quickstart may load the tracked Esrum/SLURM `crew_controllers.R` before users choose a local controller. | Add an explicit local-demo controller step before entering R, likely `cp crew_controllers_template.R crew_controllers.R`, with wording that scheduler users should adapt it instead. |
| A2 | pending | `website/demo_outputs.qmd` | `consensus_peak_BED.ATAC.immune_human_2x` does not exist in the current target manifest. | Replace with a real target, likely `consensus_peak_BPCells_matrix_dir.ATAC.immune_human_2x` or `consensus_peak_tibble.ATAC.immune_human_2x`, depending on what the example should teach. |
| A3 | pending | `website/demo_outputs.qmd` | Page implies the quickstart creates plot outputs, but the previous quickstart command only builds the Seurat export and dependencies. | State that review plots are created only when plot/checkpoint targets are requested; link to `main_running.qmd` and `gallery_main.qmd`. |
| A4 | pending | `website/main_inputs.qmd` | Required `outs/` file list omits `summary.csv` and treats `atac_possorted_bam.bam` as generally required. | Document default non-demultiplexed files separately from genotype-demultiplexing BAM requirement. |
| A5 | pending | `website/main_inputs.qmd` | Aggregations default to active, so template rows can accidentally enter unqualified runs. | Warn readers to set `is_active: false` or narrow target selection before unqualified `tar_make()`. |
| A6 | pending | `website/main_running.qmd` | `targets::tar_make()` can run enabled optional modules, not only the main pipeline. | Reword the "full main pipeline" section and possibly show a main-export selector separately from full root workflow execution. |
| A7 | pending | `website/_quarto.yml` | Rendered "View source" links point to root-level QMD paths instead of `website/*.qmd`. | Add `repo-subdir: website` under `book:` if source links should stay enabled. |
| A8 | pending | `website/performance_overview.qmd` | Benchmark image links to ignored `../outputs/benchmark/...`, likely broken on public docs. | Move/copy a finalized image into tracked website assets or remove the figure until a tracked asset exists. |
| A9 | pending | `website/downstream_differential_analyses.qmd` | Rendered source/GitHub links may be broken from public unauthenticated clients. | Likely covered by A7; verify after `repo-subdir` or decide whether source links should be disabled. |

## Batch B: Reader Routing And Prerequisites

| ID | Status | Page | Issue | Suggested direction |
|---|---|---|---|---|
| B1 | pending | `website/index.qmd` | Landing page gives no explicit next step after introducing the workflow. | Add a compact "Where to start" section linking to manual usage, quickstart, own-data setup, optional modules. |
| B2 | pending | `website/index.qmd` | Optional modules are introduced without prerequisites. | Add one sentence saying modules run from completed aggregations and need module-specific metadata/inputs. |
| B3 | pending | `website/index.qmd` | "10x Genomics Epi Multiome" wording is inconsistent with README/manual terminology. | Rename to "10x Genomics Multiome ATAC + Gene Expression" unless there is a deliberate brand reason. |
| B4 | pending | `website/index.qmd` | Overview figure lacks useful alt text. | Add descriptive alt text for the workflow figure. |
| B5 | pending | `website/intro_manual_usage.qmd` | Quickstart route links only to installation despite covering install, run, inspect. | Split quickstart into install/run/output links. |
| B6 | pending | `website/intro_manual_usage.qmd` | Own-data route links directly to inputs, bypassing overview/running pages. | Link to `main_overview.qmd`, `main_inputs.qmd`, and `main_running.qmd` as separate steps. |
| B7 | pending | `website/intro_manual_usage.qmd` | Manual routing page omits basic prerequisites. | Add a short "Before you start" sentence with Linux, R/targets familiarity, demo inputs or Cell Ranger ARC outputs. |
| B8 | pending | `website/intro_manual_usage.qmd` | Definition of "reaction" is too hardware-specific and may confuse sequencing lanes. | Define reaction as one configured 10x Multiome library/run row plus one `cellranger-arc count` output directory. |
| B9 | pending | `website/intro_manual_usage.qmd` | Manual map omits Performance and Output gallery sections. | Add low-detail bullets for these sections. |
| B10 | pending | `website/main_overview.qmd` | First own-data page lacks entry requirements and next action. | Add a short bridge to installation, configuration, and running workflow. |
| B11 | pending | `website/main_overview.qmd` | Role is ambiguous because intro page bypasses it. | Resolve with B6, or decide whether overview should be merged/shortened. |
| B12 | pending | `website/main_running.qmd` | Controller/resource prerequisite appears later in the manual. | Add a short reminder before first run to check `crew_controllers.R`, linking to distributed computing. |
| B13 | pending | `website/main_running.qmd` | "Inspect" after checkpoints is underspecified. | Add a bridge to output locations, `tar_read()`, `demo_outputs.qmd`, and gallery pages. |
| B14 | pending | `website/main_running.qmd` | Checkpoint selector mechanism is opaque. | Consider adding a `tar_manifest()` preview command before `tar_make()`. |

## Batch C: Configuration Semantics

| ID | Status | Page | Issue | Suggested direction |
|---|---|---|---|---|
| C1 | pending | `website/demo_installation.qmd` | "Download the two outputs" overstates what is downloaded; commands fetch a subset. | Say the demo downloads the required subset, while user data should keep complete `outs/` directories. |
| C2 | pending | `website/demo_installation.qmd` | Parallel data-download scripts/manifests may be stale relative to docs and config. | Decide source of truth: either update scripts/manifests or avoid pointing readers to parallel mechanisms. |
| C3 | pending | `website/demo_installation.qmd` | System requirements understate disk and command prerequisites. | Add `git`, `curl`, `tar`, HTTPS access, and a conservative disk estimate. |
| C4 | pending | `website/demo_installation.qmd` | Install timing claim is too crisp because GitHub R packages build from source. | Remove timing or scope it to pixi solve/install vs GitHub package install. |
| C5 | pending | `website/main_inputs.qmd` | Metadata requirements omit unique keys and non-overlapping non-key columns. | Document uniqueness and donor-level vs reaction-level variable placement. |
| C6 | pending | `website/main_inputs.qmd` | Demultiplexing semantics need a clearer donor-ID contract. | Explain dummy donor assignment when no VCF is configured and donor metadata matching for VCF path. |
| C7 | pending | `website/performance_distributed_computing.qmd` | Local template guidance is too optimistic for 60 GB machines. | Clarify template is sized for larger workstation and users may need fewer workers near minimum RAM. |
| C8 | pending | `website/performance_distributed_computing.qmd` | Page omits reload/restart requirement after editing `crew_controllers.R`. | Add restart or `load_project_runtime(force = TRUE)` instruction. |
| C9 | pending | `website/performance_distributed_computing.qmd` | Controller return contract is under-specified. | Add exact `controller_resources_tibble` column/order/type/name requirements. |
| C10 | pending | `website/performance_distributed_computing.qmd` | "Local default" block looks like Bash but is not a command. | Convert to prose or show actual copy/edit commands. |

## Batch D: Module-Specific Nuance

| ID | Status | Page | Issue | Suggested direction |
|---|---|---|---|---|
| D1 | pending | `website/downstream_differential_analyses.qmd` | Full checkpoint hides branch-specific config requirements for DCTC and pseudobulk. | Add pre-run note about required DCTC formula/phenotype fields and pseudobulk model entries, or narrow target selection. |
| D2 | pending | `website/downstream_differential_analyses.qmd` | Metadata prerequisites are underexplained. | Specify donor- or pseudobulk-sample-level variables joined by `donor_id`. |
| D3 | pending | `website/downstream_differential_analyses.qmd` | Parameter overview includes apparently unused OLINK/bulk RNA fields. | Remove from public reference, label as reserved, or decide to keep because they are near-term planned. |
| D4 | pending | `website/downstream_differential_analyses.qmd` | Module config page says edit symlinked/template config directly. | Mirror main config guidance: keep template as reference and maintain project-specific module rows. |
| D5 | pending | `website/downstream_genetic_enrichment.qmd` | Running prerequisite is too vague. | State dependency on GEX, ATAC, WNN checkpoints plus Open Targets study set and human aggregation. |
| D6 | pending | `website/downstream_genetic_enrichment.qmd` | Open Targets release pin and `finemappingMethod: auto` behavior are hidden. | Add reproducibility sentence under configuration. |
| D7 | pending | `website/downstream_genetic_enrichment.qmd` | `gchromVAR` wording reads like a direct external package claim. | Reword as GWAS-chromVAR-style scores using pipeline chromVAR/betterChromVAR machinery. |
| D8 | pending | `website/downstream_genetic_enrichment.qmd` | Optional gene-level branch is not tied to enabling parameter. | Note it runs only when `genetic_enrichment_psbulk_GWAS_gene_chromVAR_GWAS_IDs` lists GWAS IDs. |
| D9 | pending | `website/downstream_genetic_enrichment.qmd` | Missing implementation graph cross-link. | Add link to genetic enrichment implementation diagrams. |

## Batch E: Gallery And Output Presentation

| ID | Status | Page | Issue | Suggested direction |
|---|---|---|---|---|
| E1 | pending | `website/gallery_main.qmd` | Gallery assets are snapshots, but page may imply quickstart recreates all plots. | Add reproducibility warning and explain plot targets must be run explicitly. |
| E2 | pending | `website/gallery_main.qmd` | Dataset-level QC cards appear under aggregation demo without explaining dataset/aggregation boundary. | Add note that first QC cards are for `immune_human_dataset`, then aggregation-level cards follow. |
| E3 | pending | `website/gallery_main.qmd` | Missing setup/output cross-links. | Link to installation, demo running, demo outputs. |
| E4 | pending | `website/gallery_differential_analyses.qmd` | No context before generated cards. | Add intro saying outputs are representative module snapshots, curated subset, target names trace to workflow. |
| E5 | pending | `website/gallery_differential_analyses.qmd` | Gallery page shares title with module page. | Rename H1 to `Differential Analyses Gallery`. |
| E6 | pending | `website/output_gallery.yaml` or gallery renderer | Differential volcano cards expose dynamic branch hashes. | Prefer stable display target names plus "representative branch" wording. |
| E7 | pending | `website/gallery_differential_analyses.qmd` | Gallery under-represents module outputs promised elsewhere. | Either add more cards or explicitly call gallery curated subset. |
| E8 | pending | `website/gallery_genetic_enrichment.qmd` | Gallery silently switches from quickstart `immune_human_2x` to larger `PBMC_human_6x`. | Add intro explaining scope and that minimal quickstart does not produce all shown outputs. |
| E9 | pending | `website/gallery_genetic_enrichment.qmd` | Most rendered sections lack reader context. | Add one-sentence descriptions for all sections or move descriptions into `output_gallery.yaml`. |
| E10 | pending | `website/gallery_genetic_enrichment.qmd` | "Raw GWAS-specific trait relevance scores" is misleading for SCAVENGE UMAPs. | Reword to "per-cell graph-propagated trait relevance scores". |

## Batch F: Drift Reduction And Overview Polish

| ID | Status | Page | Issue | Suggested direction |
|---|---|---|---|---|
| F1 | pending | `website/main_overview.qmd` | Feature list is exhaustive and drift-prone. | Replace with stable stage-based summary. |
| F2 | pending | `website/main_overview.qmd` | Optional/config-dependent behavior is mixed with baseline pipeline behavior. | Add "when configured" language and cross-links to optional module pages. |
| F3 | pending | `website/main_overview.qmd` | Opening sentence blurs root workflow and main pipeline. | Say page covers the main branch of root `_targets.R`, not all targets. |
| F4 | pending | `website/demo_outputs.qmd` | Path-layout explanation duplicates implementation internals. | Keep one example and point to `tar_read()` / `tar_meta(fields = path)` for exact paths. |
| F5 | pending | `website/performance_overview.qmd` | Scaling claim is overconfident from two benchmark points. | Reword as current example benchmark, not full scaling result. |
| F6 | pending | `website/performance_overview.qmd` | Future benchmark promise reads like internal TODO. | Delete or replace with method note. |
| F7 | pending | `website/performance_overview.qmd` | Missing method/prerequisite context. | Add short "How to read this" paragraph and link to distributed computing. |
| F8 | pending | `website/performance_distributed_computing.qmd` | Cross-links are sparse. | Add opening links to installation, running workflow, and implementation conventions. |

## Review Log

Use this section to record decisions as we go.

| ID | Decision | Notes |
|---|---|---|
| A1 | pending | Start here. |
