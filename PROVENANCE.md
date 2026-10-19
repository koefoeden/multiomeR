# FLINC freeze 2.24: code snapshot

This tag holds the code and configuration that produced freeze
`2.24_w_Seurat_SCT` of the FLINC (CoCoMASLD) liver 10x Multiome data. The
freeze was completed on 10 June 2026. The code is a development version of
multiomeR. The maintained pipeline is on `main` and in the
[releases](https://github.com/koefoeden/multiomeR/releases).

## What the snapshot contains

The snapshot starts from public commit `f1a41827` (10 June 2026). One commit
adds the changes from the private development repository that the freeze
used. That private commit (`fc7ab78`) is the one the freeze recorded:

- the code in `R/gchromvar_SCAVENGE_helpers.R`,
  `R/parallel_reaction_preprocessing_helpers.R`, `R/targets_helpers.R`,
  `extra_targets/per_reaction_targets.R`,
  `module_genetic_enrichment/gchromVAR_targets.R` and
  `module_genetic_enrichment/psbulk_GWAS_gene_chromVAR_targets.R`, plus the
  empty placeholder `NULL.txt` that this code uses;
- the FLINC configuration: `cfg_aggregations.yaml`, `cfg_datasets.yaml`,
  `cfg_reactions.tsv` and `module_genetic_enrichment/cfg.yaml`.

Left out, because none of it computes the frozen results:

- private helper functions and the Slurm controller settings of the CBMR
  cluster (copy `crew_controllers_template.R` to `crew_controllers.R` to run);
- one private export of genome-browser and scE2G inputs;
- the script that copied the outputs into the freeze folder.

## What differs from the run

- The configuration keeps only the `default` and `FLINC` entries. The reaction
  table keeps the seven columns of `cfg_reactions_template.tsv`, the only ones
  the code reads.
- Donor identifiers are removed. Two wells that each hold one donor were named
  after that donor's study number. They are called `liver_single_donor_1` and
  `liver_single_donor_2` here, and their donor IDs `single_donor_1` and
  `single_donor_2`. Their paths are renamed the same way.
- One donor was excluded in the per-well QC filter of `cfg_datasets.yaml`. Its
  study ID is replaced by the placeholder `withheld_donor_ID`.
- Two inactive wells, which were commented out, are removed.
- The Galaxy API key and history ID for genome-browser track uploads are
  removed.

The paths point to the cluster where the data are stored. The data are not
included.

## Software

- `pixi.toml` and `pixi.lock` give the conda environment: R 4.5.3, Seurat
  5.5.0, SeuratObject 5.4.0, Harmony 2.0.2, chromVAR 1.32.0, MACS3 3.0.4,
  cellsnp-lite 1.2.3, vireo 0.5.9, scDblFinder 1.24.0 and targets 1.12.0.
- BPCells, Signac (`v2` branch) and betterChromVAR were installed from GitHub
  outside the lock file. Their exact commits were not recorded.
- Before multiomeR, Cell Ranger ARC 2.1.0 aligned the reads to GRCh38-2024-A,
  and CellBender 0.3.0 removed ambient RNA from the gene expression counts.

## Build history

The targets store kept results whose code it judged unchanged, so not every
step was rerun at the recorded commit:

| Built | Steps |
|---|---|
| 17–31 May 2026 | per-well inputs, GEX aggregation, peak calling and the first QC steps |
| 3 June 2026 | `GEX_Seurat_object` |
| 4 June 2026 | `multimodal_Seurat_object` |
| 10 June 2026 | single-nucleus gchromVAR z-scores and SCAVENGE scores |

## Frozen outputs

The freeze records SHA-256 checksums of its copy manifests:

| Part | Files | Manifest SHA-256 |
|---|---|---|
| Objects | 6 | `ff2430e5abf2e3a52ff1e47d15aa27130824b940b12bf8113c450fd108003f49` |
| Matrix folders | 3 | `42199b9c350590820db5d33cc1ec57c647af293abcaeb61c3ff759ec90a608d1` |
| Plots | 651 | `d2c0d4befd99c76605629c54c4b7ed35e05130aff4a1a4892f98a49bbcd9aeaa` |
| Path updates | 6 | `141366b0fdc250a3a9c3800677ed8d3f02107a0af610474c506d0ef35cd20376` |

The objects are `GEX_Seurat_object`, `multimodal_Seurat_object`,
`GWAS_chromVAR_ZScores_tibble.single_cell` and the three SCAVENGE score tables
(LSI, PCA and WNN graphs). The matrix folders hold the aggregated GEX counts,
the consensus peak counts and the TF activities.

## Later changes

Later releases revised several methods; see the release notes of v1.0.0 and
v1.2.0. Two affect these outputs:

- v1.2.0 corrected the background variance of gchromVAR z-scores when peak
  weights are fine-mapping probabilities below one. The gchromVAR z-scores and
  SCAVENGE scores in this freeze predate that correction.
- v1.2.0 estimates WNN modality weights within subsets of at most 10,000
  nuclei. This freeze estimated them across all nuclei.
