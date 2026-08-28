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
| A1 | reject | `website/demo_installation.qmd` | Fresh quickstart may load an Esrum/SLURM `crew_controllers.R` before users choose a local controller. | The committed file now contains the local controller setup; documentation explains how to adjust it for another machine. |
| A2 | accept | `website/demo_outputs.qmd` | `consensus_peak_BED.ATAC.immune_human_2x` does not exist in the current target manifest. | Replaced with the verified file target `consensus_peak_BPCells_matrix_dir.ATAC.immune_human_2x`. |
| A3 | accept | `website/demo_outputs.qmd` | Page implies the quickstart creates plot outputs, but the previous quickstart command only builds the Seurat export and dependencies. | The quickstart/output contract now distinguishes endpoint dependencies from checkpoint plot targets. |
| A4 | accept | `website/main_inputs.qmd` | Required `outs/` file list omits `summary.csv` and treats `atac_possorted_bam.bam` as generally required. | Baseline and VCF-demultiplexing file requirements are now separate. |
| A5 | accept | `website/main_inputs.qmd` | Aggregations default to active, so example rows can accidentally enter unqualified runs. | Only the quickstart GEM wells and aggregation are active by default, with explicit scope guidance. |
| A6 | accept | `website/main_running.qmd` | `targets::tar_make()` can run enabled optional modules, not only the main pipeline. | Added a single-aggregation endpoint selector and clarified root-workflow scope. |
| A7 | accept | `website/_quarto.yml` | Rendered "View source" links point to root-level QMD paths instead of `website/*.qmd`. | Added and rendered with `repo-subdir: website`. |
| A8 | revise | `website/performance_overview.qmd` | Benchmark image links to ignored `../outputs/benchmark/...`, likely broken on public docs. | Removed the unpublished image and replaced it with an evidence-bounded benchmark table. |
| A9 | accept | `website/downstream_differential_analyses.qmd` | Rendered source/GitHub links may be broken from public unauthenticated clients. | Covered by the verified `repo-subdir` fix. |

## Batch B: Reader Routing And Prerequisites

| ID | Status | Page | Issue | Suggested direction |
|---|---|---|---|---|
| B1 | accept | `website/index.qmd` | Landing page gives no explicit next step after introducing the workflow. | Added a goal-based route table. |
| B2 | accept | `website/index.qmd` | Optional modules are introduced without prerequisites. | Added completed-aggregation and module-input prerequisites. |
| B3 | accept | `website/index.qmd` | "10x Genomics Epi Multiome" wording is inconsistent with README/manual terminology. | Standardized the product wording. |
| B4 | accept | `website/index.qmd` | Overview figure lacks useful alt text. | Added workflow-specific alt text. |
| B5 | accept | `website/intro_manual_usage.qmd` | Quickstart route links only to installation despite covering install, run, inspect. | Added the full three-step quickstart route. |
| B6 | accept | `website/intro_manual_usage.qmd` | Own-data route links directly to inputs, bypassing overview/running pages. | Added overview, configuration, and checkpoint links. |
| B7 | accept | `website/intro_manual_usage.qmd` | Manual routing page omits basic prerequisites. | Added a dedicated prerequisites section. |
| B8 | accept | `website/intro_manual_usage.qmd` | Definition of "GEM well" is too hardware-specific and may confuse sequencing lanes. | Replaced it with the configuration/output-directory definition. |
| B9 | accept | `website/intro_manual_usage.qmd` | Manual map omits Performance and Output gallery sections. | Added output, operation, scaling, and troubleshooting routes. |
| B10 | accept | `website/main_overview.qmd` | First own-data page lacks entry requirements and next action. | Added entry requirements and a recommended path. |
| B11 | revise | `website/main_overview.qmd` | Role is ambiguous because intro page bypasses it. | Kept the page but made it a stable stage-based orientation page. |
| B12 | accept | `website/main_running.qmd` | Controller/resource prerequisite appears later in the manual. | Added a pre-run controller check and link. |
| B13 | accept | `website/main_running.qmd` | "Inspect" after checkpoints is underspecified. | Added concrete review criteria, key outputs, and feedback decisions. |
| B14 | accept | `website/main_running.qmd` | Checkpoint selector mechanism is opaque. | Added a verified `tar_manifest()` preview. |

## Batch C: Configuration Semantics

| ID | Status | Page | Issue | Suggested direction |
|---|---|---|---|---|
| C1 | accept | `website/demo_installation.qmd` | "Download the two outputs" overstates what is downloaded. | The Pixi task now says it downloads the required subsets of the two public outputs. |
| C2 | accept | `website/demo_installation.qmd` | Inline download commands can drift from the repository manifests. | Replaced them with the manifest-backed `download-demo-data` Pixi task. |
| C3 | accept | `website/demo_installation.qmd` | System requirements understate disk and command prerequisites. | Added commands, network, disk, RAM, and CPU guidance. |
| C4 | accept | `website/demo_installation.qmd` | Install timing claim is too crisp because GitHub R packages build from source. | Removed the fixed install-time promise. |
| C5 | accept | `website/main_inputs.qmd` | Metadata requirements omit unique keys and non-overlapping non-key columns. | Added keyed examples and column-ownership rules. |
| C6 | accept | `website/main_inputs.qmd` | Demultiplexing semantics need a clearer donor-ID contract. | Documented VCF and non-VCF donor assignment. |
| C7 | accept | `website/performance_distributed_computing.qmd` | Local controller guidance is too optimistic for 60 GB machines. | Added default sizing and safe-concurrency guidance. |
| C8 | accept | `website/performance_distributed_computing.qmd` | Page omits reload/restart requirement after editing `crew_controllers.R`. | Added restart/reload instructions. |
| C9 | accept | `website/performance_distributed_computing.qmd` | Controller return contract is under-specified. | Documented the validated shape, order, types, names, and default-row behavior. |
| C10 | accept | `website/performance_distributed_computing.qmd` | "Local default" block looks like Bash but is not a command. | Replaced it with an actual copy/customization workflow. |

## Batch D: Module-Specific Nuance

| ID | Status | Page | Issue | Suggested direction |
|---|---|---|---|---|
| D1 | accept | `website/downstream_differential_analyses.qmd` | Full checkpoint hides branch-specific config requirements for DCTC and pseudobulk. | Added full-checkpoint configuration requirements and a minimum model. |
| D2 | accept | `website/downstream_differential_analyses.qmd` | Metadata prerequisites are underexplained. | Added donor/sample-level modelling and replication requirements. |
| D3 | revise | `website/downstream_differential_analyses.qmd` | Parameter overview includes apparently unused OLINK/bulk RNA fields. | Kept the manifest fields but labeled them reserved and unused by the public checkpoint. |
| D4 | accept | `website/downstream_differential_analyses.qmd` | Module config setup previously required a separate template workflow. | The module config is now a regular active file and the guidance says to edit it directly. |
| D5 | accept | `website/downstream_genetic_enrichment.qmd` | Running prerequisite is too vague. | Added human, checkpoint, target, study, download, and interpretation prerequisites. |
| D6 | accept | `website/downstream_genetic_enrichment.qmd` | Open Targets release pin and `finemappingMethod: auto` behavior are hidden. | Documented release 26.03 and the verified method priority. |
| D7 | accept | `website/downstream_genetic_enrichment.qmd` | `gchromVAR` wording reads like a direct external package claim. | Reframed as pipeline GWAS chromVAR-style deviations. |
| D9 | accept | `website/downstream_genetic_enrichment.qmd` | Missing implementation graph cross-link. | Added a rendered cross-book graph link. |

## Batch E: Gallery And Output Presentation

| ID | Status | Page | Issue | Suggested direction |
|---|---|---|---|---|
| E1 | accept | `website/gallery_main.qmd` | Gallery assets are snapshots, but page may imply quickstart recreates all plots. | Added snapshot and checkpoint-selection guidance. |
| E2 | accept | `website/gallery_main.qmd` | Dataset-level QC cards appear under aggregation demo without explaining dataset/aggregation boundary. | Explained dataset versus aggregation cards. |
| E3 | accept | `website/gallery_main.qmd` | Missing setup/output cross-links. | Added install, run, output, and checkpoint links. |
| E4 | accept | `website/gallery_differential_analyses.qmd` | No context before generated cards. | Added curated-snapshot and interpretation context. |
| E5 | accept | `website/gallery_differential_analyses.qmd` | Gallery page shares title with module page. | Renamed all gallery pages unambiguously. |
| E6 | accept | `website/output_gallery.yaml` or gallery renderer | Differential volcano cards expose dynamic branch hashes. | Replaced hashes with stable parent target names. |
| E7 | accept | `website/gallery_differential_analyses.qmd` | Gallery under-represents module outputs promised elsewhere. | Explicitly labeled the gallery as a curated subset. |
| E8 | accept | `website/gallery_genetic_enrichment.qmd` | Gallery silently switches from quickstart `immune_human_2x` to larger `PBMC_human_6x`. | Added an explicit scope note for six GEM wells. |
| E9 | accept | `website/gallery_genetic_enrichment.qmd` | Most rendered sections lack reader context. | Added page and subsection context and removed stale missing-asset cards. |
| E10 | accept | `website/gallery_genetic_enrichment.qmd` | "Raw GWAS-specific trait relevance scores" is misleading for SCAVENGE UMAPs. | Reworded as per-cell graph-propagated scores. |

## Batch F: Drift Reduction And Overview Polish

| ID | Status | Page | Issue | Suggested direction |
|---|---|---|---|---|
| F1 | accept | `website/main_overview.qmd` | Feature list is exhaustive and drift-prone. | Replaced with a stage-based overview. |
| F2 | accept | `website/main_overview.qmd` | Optional/config-dependent behavior is mixed with baseline pipeline behavior. | Separated baseline stages from configuration-dependent capabilities. |
| F3 | accept | `website/main_overview.qmd` | Opening sentence blurs root workflow and main pipeline. | Clarified the main branch versus root workflow. |
| F4 | accept | `website/demo_outputs.qmd` | Path-layout explanation duplicates implementation internals. | Reduced it to one path example and the user-relevant grouping. |
| F5 | accept | `website/performance_overview.qmd` | Scaling claim is overconfident from two benchmark points. | Replaced the claim with an explicitly limited two-run snapshot. |
| F6 | accept | `website/performance_overview.qmd` | Future benchmark promise reads like internal TODO. | Removed the internal promise. |
| F7 | accept | `website/performance_overview.qmd` | Missing method/prerequisite context. | Added estimation assumptions and usage guidance. |
| F8 | accept | `website/performance_distributed_computing.qmd` | Cross-links are sparse. | Added quickstart, troubleshooting, and implementation links. |

## Review Log

Use this section to record decisions as we go.

| ID | Decision | Notes |
|---|---|---|
| A1 | reject | Revalidated after materializing `crew_controllers.R` as the local default configuration. |
| Implementation | complete | Reader-journey revision implemented and rendered on `codex/revise-doc-reader-journey`. |
| Concision pass | complete | Removed the redundant manual-routing chapter, compressed repeated summaries, and moved parameter details behind searchable disclosure controls. |
