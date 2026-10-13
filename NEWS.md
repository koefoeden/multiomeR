# multiomeR 1.2.1 (2026-10-12)

multiomeR 1.2.1 lowers the memory requirement from 60 GB to 32 GB of RAM, labels
the `brain_mouse` example with embryonic marker sets, and fixes the UMAP
parameter sweeps on large aggregations. The results of your aggregations do not
change: existing stores stay current apart from a few plots.

## What you need to do

1. **Update to this release.** Releases tag stable, validated commits, so update
   to the tag rather than pulling `main`, which can carry development work that
   is not yet validated. Commit or stash local changes first; resolve any rebase
   conflicts as usual. While no other R process uses the environment, run:
   ```bash
   cd multiomeR  # your clone
   git fetch origin tag v1.2.1
   git rebase v1.2.1
   pixi install --locked --run-post-link-scripts
   pixi run --use-environment-activation-cache install-r-github-packages
   ```
   Nothing else is required: the environment and the parameter schema are
   unchanged.

2. **Coming from 1.1.0?** Also follow steps 2 to 6 of the [1.2.0
   notes](https://github.com/koefoeden/multiomeR/releases/tag/v1.2.0). They now
   include the rename of AMULET's fragment-count column from `nFrags` to
   `amulet_nFrags`, which 1.2.0 made like the other AMULET columns but its notes
   first omitted.

## What happens to your aggregations

- **Existing stores stay current.** Resource requests, controllers and target
  descriptions do not invalidate targets.
- **Only plots rebuild.** The GEX and ATAC UMAP sweeps (`cross.UMAPs.3_GEX_QC`,
  `cross.UMAPs.7_ATAC_QC`) and `TRS_cluster_summary_plot` rerun with new plots.
  `UMAP_embeddings_tibble.GEX`, `.GEX_non_harmony` and `.ATAC` rerun and
  reproduce their values.
- **`brain_mouse`** rebuilds from its marker lists onward.
- **In `configuration_dev/`**, only targets of
  `pipeline_manuscript_figures/targets.R` rerun: Figure 1 panels A and C–E with
  its plot and data, and the Table S4 writer. `figure_1_traits`,
  `supplementary_table_S3` and `supplementary_figure_S3` are new.

## What changed

### Memory requirement

- Each target's `get_tar_resources()` request is about 1.5 times the peak memory
  it reached on aggregations of 213,000 and 713,000 nuclei, in 16- and 32-GB
  steps. With the default GEX backend, no target requests more than 32 GB,
  against 49 requests of 40 or 60 GB before; the largest targets need about 7 GB
  plus 2.5 GB per 100,000 nuclei. Unmeasured targets request at most 32 GB.
- GEX normalization and the two Seurat objects request 60 GB only when an
  aggregation of the configuration uses the optional `Seurat_SCT` GEX backend,
  whose GEX normalization needs about 0.25 GB per 1,000 nuclei.
- The default `local-heavy` controller is declared with 32 GB, and a new
  `local-large` controller (6 cores, 60 GB, one worker) takes those requests. A
  custom `crew_controllers.R` needs a tier of at least 60 GB only for
  aggregations that use `Seurat_SCT`.
- Built from an empty store in a 16-CPU Slurm job limited to 32 GiB, the demo's
  processes used at most 17 GiB.
- Target descriptions record the measurements behind the requests as
  `[resource_observation: ...]` tags. A target may carry several, and the manual
  suggests common keys.

### Examples

- The inactive `brain_mouse` example (10x Genomics' E18 mouse brain) has
  embryonic marker sets: radial glia, intermediate progenitors, upper- and
  deep-layer and hippocampal excitatory neurons, and inhibitory neurons. The
  adult glial sets left every WNN cluster unassigned; the new sets label these
  populations, leaving at most a small mixed cluster unassigned, and its ATAC
  marker TFs follow.

### Fixes

- The GEX and ATAC cross-parameter UMAP sweeps run over the first n configured
  UMAP components, so the ATAC sweep no longer includes the LSI depth component
  and its largest setting matches the main UMAP. They are initialized with
  spectral PCA (`init = "spca"`), because uwot's default spectral initialization
  could crash on embeddings of several hundred thousand nuclei.
- The SCAVENGE TRS interval plot uses at most four facet columns, fewer for more
  than 30 groups, so its width stays within the 50,000-pixel image limit.

### Manuscript results

- Figure 1C shows the cell-type enrichments at q < 0.05 that rest on many loci
  (top locus below half of the signal, at least 20 effective weighted peaks),
  one per study and per cell type and trait category, with stars at q < 0.05;
  panels D and E show LDL cholesterol in hepatocytes and colorectal cancer in
  colonocytes. Trait names are in sentence case, and panel A names its cell type
  as a plural noun.
- New Supplementary Figure S3 shows every tested study's cell-type deviations on
  Figure 1C's scale, and Supplementary Table S3 is generated from the stored
  GWAS inputs of `mixed_human_31x`.
- Supplementary Figure S2's reserved CPU and RAM are computed from the worker
  sizes of the measured runs, recorded in
  `seurat_signac_comparison_allocations.tsv`, so later request changes leave
  them unchanged.
- The Table S4 data leave out the compute node's name.
- The manuscript figure helpers wrap patchwork panels with
  `patchwork::wrap_elements()`, so their legends and labels are measured on the
  saving device and no longer overlap or clip.

## Development

- Agents start long pipeline drivers detached and follow their log, since an
  agent tool's time limit stops the driver and its running targets.

# multiomeR 1.2.0 (2026-10-12)

multiomeR 1.2.0 revises the statistical methods after an audit: GEX
normalization, WNN modality weights, clustering, consensus peaks, doublet
detection, and the differential, gene-set, composition, genetic-enrichment and
peak–gene tests. Every aggregation's results change, from QC and cell-type
labels to all module results, so plan a full rebuild after updating. The release
also adds an optional WNN-specific marker list, per-cluster composition tables,
lower memory use for CellBender input, targets stores that can move between
disks, and the public configuration, figures and results of the multiomeR
manuscript.

## What you need to do

1. **Update to this release.** Releases tag stable, validated commits, so update
   to the tag rather than pulling `main`, which can carry development work that
   is not yet validated. Commit or stash local changes first; resolve any rebase
   conflicts as usual. While no other R process uses the environment, run:
   ```bash
   cd multiomeR  # your clone
   git fetch origin tag v1.2.0
   git rebase v1.2.0
   pixi install --locked --run-post-link-scripts
   pixi run --use-environment-activation-cache install-r-github-packages
   ```
   This installs Seurat 5.6.0 from GitHub, with its build dependencies
   `r-pbapply` and `r-rcppthread` from the lock file.

2. **Give your configuration directory its own `crew_controllers.R`.** The root
   `crew_controllers.R` is gone: each configuration directory now carries
   controllers that suit its aggregations, and the default moved to
   `configuration/crew_controllers.R`. If you selected your own directory in
   `configuration.local`, copy `configuration/crew_controllers.R` into it and
   adapt it, moving in any changes you made to the former root file.

3. **Update moved paths** in your configuration and scripts:
   - `cfg_pipeline_parameters.tsv`, `QC_checkpoint_manifest.tsv`,
     `QC_metric_manifest.tsv` and the demo download lists are in `manifests/`.
   - `reference_metadata/` is `resources/reference_metadata/`; add new Cell
     Ranger ARC references there.
   - `module_<name>/` is `extra_targets/module_<name>/`.
   - The example donor metadata tables formerly in `example_data/` are in
     `configuration/` and `configuration_dev/`, and the download script is in
     `scripts/`. `example_data/` now only receives downloads; three unused
     tables were removed.
   - `mixed_human_31x` and the nested `comparison_1x` to `comparison_20x`
     aggregations moved to the new `configuration_dev/`, with their donor
     metadata and GWAS selection; select it to rebuild them.
     `manuscript_figures/` is replaced by `pipeline_manuscript_figures/`.

4. **Check your configuration against the changed defaults and checks:**
   - The cell-cycle score is no longer regressed out by default. Set
     `aggregation_SCT_regress_vars: [CC.Difference]` to keep it.
   - A differential-analysis formula that uses `cluster` without
     `cell_type_formula` or `random_effect` now stops with an error: it treated
     each donor's cell types as independent samples. Rewrite such formulas.
   - CellBender counts are used whenever a GEM well gives a CellBender H5 file;
     `GEM_well_add_cellbender` is no longer needed and a leftover column is
     ignored. Set `GEM_well_cellbender_h5_file` to `NA` to keep Cell Ranger
     counts.
   - Marker symbols that several genes share now stop marker validation, as
     missing symbols do, because they were scored on only one of the genes.
     Replace them.
   - The example configurations now exclude nuclei with more than 5% of their
     peak counts in blacklisted regions (was 1%, which removed cell types
     unevenly) and AMULET multiplets at the authors' q < 0.01 (was 0.05).
     Configurations copied from the old examples keep the old limits; change
     them to follow the examples.

5. **Update code that reads pipeline outputs:**
   - Each unassigned cluster is labeled `Unassigned_cluster_<id>` instead of one
     merged `Unassigned`.
   - AMULET's fragment-count column is `amulet_nFrags` instead of `nFrags`, like
     the other AMULET columns.

6. **Demo users:** run `pixi run --use-environment-activation-cache setup-demo`
   again to download `controller_PBMC_human`.

## What happens to your aggregations

- **Everything after the per-GEM-well inputs rebuilds.** Every GEM well reruns
  `amulet_metrics_tibble` and `full_metadata_tibble`, and every aggregation
  target rebuilds. Cell Ranger, CellBender, cellsnp-lite, Vireo and the other
  per-GEM-well inputs stay current. The first run downloads the AMULET
  excluded-region list.
- **Cell-type labels may change.** GEX normalization, WNN weights, clustering
  and doublet calls all change, and each unassigned cluster keeps its own label.
  Review the labels at checkpoints 3, 7 and 8 before the modules run; the new
  `aggregation_WNN_marker_genes` lets you tune the final labels without
  rebuilding the GEX and ATAC processing.
- **Peaks change.** Each group's peaks must be called in both random halves of
  its nuclei, so small groups lose their weakly supported peaks.
- **Module results change:** differential expression and accessibility, gene-set
  and composition tests, genetic enrichment, and peak–gene links, which are now
  calibrated and fewer.
- **Targets that rerun and keep their results.** Targets that open BPCells
  matrices or fragments now record project-relative paths, so they rerun with
  unchanged data and their consumers rerun with them:
  `GEX_counts_BPCells_matrix`, `fragments_w_prefix_bpcells`,
  `aggregated_counts_BPCells_matrix.GEX`, `consensus_peak_BPCells_matrix.ATAC`,
  `motif_family_accessibility_BPCells_matrix.ATAC`,
  `gene_score_archr_BPCells_matrix.ATAC` and the pseudobulk count matrices.
  `GEM_well_config_tsv`, `GEM_well_metadata_tibble`,
  `GEM_well_annotation_metadata_tibble`, `GEX_counts_BPcells_matrix_dir`,
  `chromHMMs_list_proj.ATAC`, the five downloaded reference-file targets, the
  file targets of the moved manifests, reference metadata and donor metadata,
  and each GEM well's Cell Ranger reference lookup
  (`cellranger_reference_json_file`) rerun and reproduce their results. The
  exclusion lists (`excluded_barcodes_by_type_list`,
  `excluded_cellranger_only_barcodes_by_type_list`) rerun for GEM wells whose QC
  list lacks `not_found_in_GEX_matrix`, which includes every GEM well of the
  public configurations; their changes stop at the checkpoint 1 plots and
  retention tables. The new composition tables build without rerunning anything
  else.
- **Demo:** its aggregation targets rebuild with the new `controller_PBMC_human`
  well.
- **In `configuration_dev/`**, the Figure 1 targets rebuild. Supplementary
  Figure S2 reads the recorded runtimes of the compared workflows from the
  store, so rebuild its data before rerunning targets of the `comparison_`
  aggregations.

## What changed

Every change under Processing, Quality control, Differential analyses, Genetic
enrichment and Peak–gene links alters results.

### Processing

- The BPCells backend's Pearson residuals use one overdispersion, theta = 100,
  for every gene and are clipped at ±√N, the analytic Pearson residuals of Lause
  et al. (2021). The per-gene method-of-moments theta counted differences
  between cell types as overdispersion and, with the ±10 clip, flattened marker
  variation: on a multi-tissue test set, 24 of 27 canonical markers now enter
  the 3,000 variable genes, against 10.
- `aggregation_SCT_regress_vars` is unset by default, so the cell-cycle score,
  mostly noise in post-mitotic nuclei, is no longer regressed out;
  `CC.Difference` remains an example value.
- WNN modality weights are estimated within seeded random subsets of at most
  10,000 nuclei, so they no longer drift towards 0 or 1 as an aggregation grows.
  Neighbours are still chosen among all nuclei, and smaller aggregations are
  unchanged. Each nucleus is now first in its own neighbour list.
- Leiden clustering runs for ten iterations, as Seurat does, instead of igraph's
  default two.
- Consensus peaks are built as ArchR does. Each group's discovery nuclei are
  split into two random halves, MACS3 calls peaks in each half, and only peaks
  called in both are kept; each group then keeps its strongest peaks up to 500
  per nucleus and 150,000 in all, and overlapping peaks across groups are ranked
  by the quantile of their significance within their group. Deep groups no
  longer dominate, and small groups lose their weakly supported tail (a
  1,200-nucleus group keeps 54% of its peaks). A half takes about half the time
  and memory of a full call, so MACS3 workers request 16 GB. Consensus peaks
  keep their group in `cluster`.
- The JASPAR familial root motifs lose flank columns with less than half of the
  largest column count; these low-support columns restricted the matches of 56
  of 233 families to a sub-family.
- Each unassigned cluster is labeled `Unassigned_cluster_<id>` instead of all
  being merged into one `Unassigned` group, which peak calling, pseudobulks,
  composition and genetic enrichment analysed as a real cell type.

### Quality control

- AMULET runs with chrM, chrX and chrY excluded in full and, for GRCh38, the
  authors' list of repeats, segmental duplications and blacklisted regions,
  downloaded once into the store and checked against its SHA-256 checksum. The
  example configurations exclude nuclei at the authors' q < 0.01 instead of
  0.05.
- In insertions mode, the peak fraction is divided by twice the fragment count,
  so `atac_peak_counts_frac` is a fraction and filters on it take effect.
- ATAC scDblFinder groups peaks on singular-value-scaled LSI loadings, as
  scDblFinder itself does.
- The example configurations exclude nuclei with more than 5% of their peak
  counts in blacklisted regions instead of 1%, which removed cell types
  unevenly.
- Marker symbols that several genes share are rejected; they were scored on one
  of the genes.

### Differential analyses

- GEX markers and ATAC differential peaks are tested on counts per 10,000
  instead of raw counts, with Seurat-style fold changes, so they no longer
  follow each group's sequencing depth.
- Motif-family pseudobulks are tested on bias-corrected chromVAR deviations
  without centring or quantile normalization, which collapsed cell-type-defining
  families and produced false positives.
- Contrasts between annotation classes in the paired route are fitted on each
  donor's combined response, weighted by its precision weights, instead of using
  one residual correlation between classes for all features, which made
  low-expression features anti-conservative.
- The Hallmark and Reactome tests run with camera, estimating each set's
  inter-gene correlation, instead of cameraPR's default of 0.01; in permutations
  of a donor covariate, the old test called 17% of sets at P < 0.05. Sets need
  at least ten tested genes.
- The random-effect route uses edgeR's `voomLmFit()`, so precision weights
  account for the correlation block.
- A joint formula that uses `cluster` without `cell_type_formula` or
  `random_effect` stops with an error. The example formula now tests one
  annotation class.
- Composition models test by default only assigned cell types present in at
  least half of the donors.

### Genetic enrichment

- The cell-type pseudobulks are summed by the final WNN labels, as SCAVENGE, the
  differential analyses and the peak–gene module do, instead of by the GEX
  labels carried into the WNN metadata, so labels set only at the WNN level
  reach the pseudobulk enrichment too.
- gchromVAR's background variance is scaled with squared peak weights, the
  variance of a weighted sum, in both the cell-type contributions and the
  single-nucleus z-scores. Summed weights deflated the z-scores by a
  trait-dependent factor whenever fine-mapping probabilities were below one.

### Peak–gene links

- The model includes up to ten principal components of each cell type's residual
  peak accessibility, one per ten residual degrees of freedom when fewer, from
  the new target `peak_gene_correlation_latent_factors.WNN`, so peaks and genes
  that vary along the same chromatin states no longer appear linked.
- Each candidate pair is calibrated against its own gene, following the idea of
  Signac's LinkPeaks: up to 200 of the branch's candidate peaks more than 10 Mb
  from the gene's TSS are scored with the same model, the pair's statistic is
  standardized by this background (`hierarchical_gene_z`) and rescaled within
  peak-detection quintiles, and its empirical P-value against the cell type's
  background is Benjamini–Hochberg-corrected. Links also need a positive gene z.
  Before, the model P-values were corrected directly, and random pairs on
  different chromosomes passed the FDR cutoff almost as often as cis pairs. On
  held-out background pairs, the new P-values fall below 0.05 and 0.01 within 5%
  of the nominal rates.
- The random-slope test requires within-donor peak variation in at least eight
  donors and leaves unassigned clusters out.
- SuSiE's purity filter (minimum absolute correlation 0.5) is applied, so
  credible sets no longer span most candidate peaks.

### Options and outputs

- The optional `aggregation_WNN_marker_genes` replaces, adds or, given as `~`,
  drops marker sets for the WNN labels only, with their own matched controls.
  The final cell types can then be tuned without rebuilding the GEX and ATAC
  processing, whose doublet detection and peak calling use the GEX labels.
- `cluster_composition_table.tsv` is written at checkpoints 3, 7 and 8: for each
  cluster and each categorical variable, `GEM_well_ID` and `donor_id`, the
  cluster size, its most frequent category and that category's share, and every
  category above 5%.
- The peak–gene results gain the gene-standardized statistic
  `hierarchical_gene_z`, the calibrated `hierarchical_background_pvalue` and
  each cell type's `n_background_pairs`.
- `crew_controllers.R` is read from the selected configuration directory, and
  CellBender counts follow the H5 file alone (see What you need to do).
- The new `configuration_dev/` builds the manuscript's results in its own
  targets store (`TAR_PROJECT=dev`, store `outputs_dev`), with Slurm workers on
  one node. `configuration/` keeps the demo and the small examples, and
  validation covers whichever validation aggregations the selected directory
  defines.

### Performance and fixes

- CellBender matrices are read in blocks of 2,000 barcodes, keeping only the
  gene features, so a 34,498-nucleus GEM well peaks at 5.2 GB instead of 17.1
  GB. The written matrices are byte-identical.
- Seurat 5.6.0, installed from GitHub, has threaded C++ kernels for
  normalization, PCA, neighbour search and UMAP; they get the cores of each
  Slurm worker through `Seurat.nthreads`.
- Opened BPCells matrices and fragments and downloaded reference files are
  recorded by unresolved, project-relative paths, so a targets store keeps
  working after it moves, for example from node-local disk to shared storage.
- Roadmap chromHMM annotations load per aggregation, so activating another
  aggregation or selecting another configuration directory no longer rebuilds
  them.
- Nuclei missing from the GEX count matrix are always excluded at checkpoint 1,
  so the checkpoint 1 plots and retention tables report them instead of the
  checkpoint 3 retention table counting them as small GEX clusters.
- All unassigned clusters are drawn as one track in the ATAC coverage plots.

### Demo

- `controller_PBMC_human`, 10x Genomics' 10k PBMC set from the Chromium
  Controller, replaces `unsorted_PBMC_human`. The unsorted 3k well came from the
  same donor as `healthy_PBMC_human`, so the demo's donors `pbmc1` and `pbmc2`
  are now two people, and its categorical UMAPs and composition tables show the
  donor.
- Harmony corrects for the GEM well, which now separates donors and runs, and
  cell types need a 0.02 margin over the marker background instead of 0.05,
  which the T-cell clusters, inflated by the larger well's T cells, did not
  reach. The NK markers drop KLRD1, which effector CD8 T cells share, for
  SH2D1B, and a MAIT set (SLC4A10, ZBTB16, KLRB1, NCR3) labels the MAIT cells.
  With the revised methods, three small WNN clusters (725 of 7,838 nuclei) match
  no marker set above background and stay unassigned.
- The download grows from about 1.3 to 4.6 GB, the stored results from about 3.3
  to 6.6 GB and the run from about 20 to 30 minutes with 16 threads.

### Manuscript

- The six 10x Genomics PBMC GEM wells of `mixed_human_31x` are labeled with
  their two donors: `pbmc1` for the four from one woman and `pbmc2` for the two
  from one man, as 10x Genomics' dataset descriptions, sex markers and
  mitochondrial genotypes agree. They were six donors before. The `comparison_`
  copies keep their labels, since those aggregations run no donor-level
  analyses.
- `pipeline_manuscript_figures/` builds Figure 1, Supplementary Figure S2 and
  Supplementary Table S4 of the manuscript from public data in
  `configuration_dev/`, with the results and their data under `results/`.
- Figure 1, Supplementary Figure S2 and Supplementary Table S4 are rebuilt with
  the revised methods. Figure 1A–B now shows a peak 72 kb from SGPL1 in NK
  cells, its trait panels show assigned cell types only, and `figure_1_data`
  writes the values the manuscript quotes. Panels C to E spell cell types alike,
  and panel D shows protein-coding genes with labels that no longer overlap.
- Supplementary Figure S2 compares multiomeR with a current Seurat v5/Signac v2
  workflow on BPCells-backed assays, counting only the analyses both workflows
  perform. Its peak quantification uses fragtk, which `pixi run
  --use-environment-activation-cache install-fragtk` installs.

### Documentation

- The In depth pages and the shared methods text describe the revised methods,
  and the checkpoint 8 guidance mentions the WNN marker list.
- New guidance covers aggregations that combine tissues, previewing the
  per-GEM-well targets of checkpoint 1, and leaving the costly Seurat exports
  out until the cell-type labels are accepted.
- The manual describes the per-cluster composition tables.
- The workflow overview figure no longer shows the dropped HC3 scan.
- Bibliography entries are complete and case-protected.
- `CITATION.cff` is added; GitHub shows it as "Cite this repository", and Zenodo
  uses it for the authors of archived releases. The README and the manual's home
  page ask users to cite multiomeR, targets and BPCells, and the methods behind
  the results they report; `CITATION.cff` lists targets and BPCells as
  references.

## Development

- The WNN and cell-cycle reference tests require the pinned Seurat 5.6.0, and
  the peak–gene parity test scores background pairs.
- `validation/` and `RELEASES.md` moved into `dev/`, release validation reports
  are kept in `dev/validation/reports/`, the scripts that make Cell Ranger ARC
  test subsets moved to `dev/`, and the unused `air.toml` was removed.

# multiomeR 1.1.0 (2026-09-29)

A documentation and demo release. The public demo aggregation is rebuilt with
new data and markers; other aggregations need no recomputation.

## Demo

- Combine two healthy PBMC GEM wells from the same 10x Genomics series in the
  demo: the granulocyte-sorted well as before and the unsorted 3k well, now
  `unsorted_PBMC_human`, which replaces the lymphoma lymph node. The download
  shrinks from about 4.1 to 1.3 GB and the stored results from about 6.4 to
  3.3 GB.
- Annotate the demo with marker sets whose genes single nuclei capture: B,
  Mono, NK, DC and T. The previous two-gene markers left most demo nuclei
  `Unassigned`; every GEX, ATAC and WNN cluster is now labeled.
- Plot only the sorting annotation in the demo's categorical UMAPs, because its
  other GEM-well annotations are now constant.

## Documentation

- Merge the implementation book into the manual. Its method descriptions form
  an In depth part with one primary-module page organized by checkpoint stage,
  one page per optional module and the algorithm-validation page; the
  implementation conventions become one Change the pipeline page. Each
  checkpoint of the running guide and each module page links to its methods.
- Mark the reading layers in the sidebar: the guides, a Reference part that now
  includes the parameter browser, and the In depth part with an advanced badge.
- Show the printed results of every command on the demo-results page, refreshed
  from the built demo with `pixi run refresh-demo-outputs`.
- Give the GEM well, donor and aggregation reference pages one layout, and
  document pooled GEM wells, the rules checked and when, and aggregation
  inheritance.
- Use US spelling throughout the manual.
- Deploy the documentation only from `main`.

## Validation

- Retry each live preflight check of `validate-local` up to four times, also
  after refused connections, so an intermittent outage of a remote source no
  longer fails the whole validation.

## Migration

- Run `pixi run --use-environment-activation-cache setup-demo` again to
  download the new demo well; `example_data/lymphoma_lymph_human` is no longer
  used and can be deleted.
- Target invalidation: every `immune_human_2x` target and the per-GEM-well
  targets of `unsorted_PBMC_human` build anew. Changed target descriptions
  invalidate nothing, so other aggregations stay current.

# multiomeR 1.0.0 (2026-09-29)

The first stable release. The manuscript describing multiomeR is under peer
review (link to come). From this release on, incompatible changes to
configuration, target names or output schemas require a new major version; see
[RELEASES.md](dev/RELEASES.md).

This release changes cell-type labels, peak–gene inference, genetic-enrichment
inputs and several target and parameter names. Expect a full rebuild of existing
aggregations and review the migration notes first.

## Configuration

- Keep each project's settings in its own configuration directory, selected by
  a single path in the ignored `configuration.local`; `configuration/` holds the
  public defaults and examples. Module settings use flat
  `cfg_module_<module>.yaml` files, and disabled modules need no file.
- Pass scalar component counts to GEX PCA and ATAC LSI, so changing only the
  selected dimensions, for example dropping the first LSI component, no longer
  recomputes either reduction.
- Stop before building the target graph, with the install command, when the
  Bioconda BSgenome data packages are missing because Pixi skipped their
  post-link scripts.
- Add the output gallery's `mixed_human_31x` aggregation of public 10x Genomics
  and ENCODE data as an inactive example with its module settings. Its GEM wells
  are processed locally with Cell Ranger ARC 2.1.0 and GRCh38-2024-A.

## QC and cell-type annotation

- Title the QC-exclusion UpSet plots with the number of retained barcodes and
  plot the nuclei per donor that reach the final WNN object.
- Centre the motif-family accessibility heatmaps on each family's mean across
  groups, so abundant groups do not set their zero point.
- Unify pre-filter QC plots across modalities, show doublet evidence beside GEX
  cluster markers, and report the selected wells, QC rules and configuration
  values in checkpoint captions.
- Diagnose WNN modality weights by their association with metadata within
  clusters, and show WNN markers and UCell evidence as cluster dot plots.

## Differential analyses

- Use WNN-derived cell-type labels throughout, with named abundance and feature
  models, predictor-only abundance formulas and descriptive target names.
- Test gene sets on the existing contrast results.
- Give donors absent from the extended donor metadata missing model variables,
  so models exclude them instead of the module stopping.
- Estimate intra-block correlations on at most 50,000 evenly spaced features and
  give each forked cell-type fit one BLAS thread, so a 590,000-peak
  chromatin-accessibility fit takes minutes instead of hours.

## Genetic enrichment

- Report SCAVENGE trait-relevance scores without permutation P-values, summarized
  by WNN cluster and cell type in four heatmaps.
- Decompose cell-type chromVAR deviations additively into peak, variant and
  locus contributions, with ordinary and automatic absolute-effect weighting,
  and screen locus detail plots with a configurable z-score threshold.
- Compare the credible sets of all configured GWAS by their shared normalized
  posterior-probability mass, to show which traits give non-independent
  enrichment results.
- Build each GWAS's peak weights once and derive the chromVAR annotation from
  the peak-to-variant allocation. This also fixes the ordinary posterior-
  probability path, whose normalized inputs had lost their `posteriorProbability`
  column, and reports PIPs rather than effect weights as absolute-effect locus
  leads.
- Sort the Open Targets credible-set records deterministically and allocate
  peak contributions to variants without per-variant summaries.
- Weight every cell type equally in the expected accessibility of the
  cell-type pseudobulk chromVAR background, so an abundant cell type is no
  longer compared mainly with itself and pinned near zero.

## Peak–gene correlation

- Make peak–gene correlation an optional module on WNN cell types.
- Replace the HC3 scan with a hierarchical donor-slope model fitted by native
  REML with Kenward–Roger inference over all eligible pairs; links require a
  positive, reliable hierarchical estimate at FDR < 0.05.
- Filter hypotheses by RNA, ATAC and shared-donor measurement support, and show
  genomic context, coverage and donor-adjusted scatterplots in top-link figures.
  The moderate support preset is the default.
- Fit the hierarchical model in donor space, without aggregate-by-aggregate
  matrices, so a branch projected at about a day takes minutes.

## Runtime and dependencies

- Export normalized GEX data in the Seurat objects: `RNA` gains a lazy
  log-normalized `data` layer, the BPCells backend adds the regressed Pearson
  residuals as `RNA` `scale.data`, and the Seurat backend adds the
  `SCTransform()` assay as `SCT`. The cell metadata includes the Seurat
  cell-cycle scores.
- Build the multimodal Seurat export without rehashing every fragment file and
  count each assay once, which takes about a fifth of the time for large
  aggregations.
- Save plots through one staged path that removes only obsolete outputs recorded
  in each target's inventory, and build each ggplot once.
- Save plots as images only. The mirrored `plot_objects/*.rds` copies
  re-serialized each plot's inputs, reaching hundreds of gigabytes for large
  aggregations, and took about half of all plot-saving time.
- Score the cluster UCell evidence with a compiled kernel, and compute the
  consensus peak matrix, ATAC blacklist counts and gene scores per GEM well;
  both return identical results, several times faster.
- Compute ATAC coverage tracks from an in-memory copy of the plotted regions,
  which takes seconds instead of minutes to hours per region set on large
  merged fragment objects and returns identical data.
- Make the project bootstrap stateless and share one checksum-keyed loader for
  the standalone native sources.
- Rerun the dependents of opened BPCells directories when the directory
  content changes; previously they could keep results computed from the old
  content, and a cheap CI check now rejects bare BPCells opens in targets.
- Update the locked environment within R 4.5 and Bioconductor 3.22 (among them
  Seurat 5.5.1, scDblFinder 1.24.10, arrow 25 and Python 3.13), move BPCells,
  Signac and betterChromVAR to current revisions, and drop 31 unused
  dependencies. multiomeRCore 0.2.0 narrows `get_tar_resources()`.
- Keep only reference-parity tests of the reimplemented algorithms, which pass
  with the updated reference packages.
- Run only cheap checks in cloud CI: documentation source links and, also
  weekly, the example-data download URLs. Run the reference-parity tests
  locally, and build the documentation site with the repository builder.

## Documentation

- Restructure the manual around Configure, Run and Review steps, with one page
  per checkpoint and module, a standalone parameter browser, an output gallery
  generated from saved plots and an implementation book with one methods page
  per stage.
- Share the method descriptions with the manuscript supplement and describe the
  validation contract of each reimplemented algorithm.

## Migration

- The Seurat objects no longer contain the counts-only `SCTregr` assay. Use
  `RNA` (BPCells backend) or `SCT` (Seurat backend); the WNN weight column is
  `RNA.weight` or `SCT.weight` instead of `SCTregr.weight`.
- Move existing settings into a configuration directory and select it in
  `configuration.local`.
- Remove these parameters, which configuration validation now rejects:
  `aggregation_GALAXY_track_upload_API_KEY`,
  `aggregation_GALAXY_track_upload_HISTORY_ID`,
  `aggregation_tar_make_skip_regex_patterns`,
  `differential_analyses_bulk_RNA_rds_file_path`,
  `differential_analyses_OLINK_parquet_file_path`,
  `genetic_enrichment_posterior_probability_weighting_function_name` and
  `genetic_enrichment_SCAVENGE_permutation_times`, and the per-study
  `variant_weighting_mode` GWAS field.
- Reinstall the environment with `pixi install --locked --run-post-link-scripts`
  and the pinned GitHub packages with `pixi run install-r-github-packages`.
- Update scripts that read HC3 peak–gene results, SCAVENGE P-values or renamed
  differential-analysis targets.
- Redraw plots from their data instead of reading `plot_objects/*.rds`, as the
  demo-output page shows.

# multiomeR 0.5.0 (2026-09-11)

Released through [PR #5](https://github.com/koefoeden/multiomeR/pull/5).

This release changes cluster annotations, cluster IDs, peak–gene inference and
several target/output names. Review the migration notes before rerunning an
existing analysis.

## QC and interpretation

- Organize review outputs into eight numbered QC checkpoints, separating GEX
  PCA and ATAC LSI review from clustering. Add cumulative, per-action nuclei
  retention flows and validate matching gene definitions before aggregation.
- Select QC violin metrics from the shared manifest. Correct UpSet percentages,
  align RNA/ATAC/WNN confusion matrices, improve marker-set dot plots, and show
  all computed PCA/LSI dimensions in review plots.
- Add within-cell-type cluster-marker comparisons and volcano plots. Number
  GEX, ATAC and WNN Leiden clusters by decreasing size without changing their
  membership; break size ties by the original cluster ID.
- Add readable motif-family labels, palettes that accommodate many cell groups,
  and aligned genetic-enrichment panels and shared legends.

## Cluster annotation

- Replace the former standardized module-score labeller with matched-control
  UCell evidence across GEX, ATAC and WNN. All three use a shared GEX control
  reference. There is one annotation method, with no legacy selector.
- Assign the leading label only when its adjusted score exceeds both matched
  background and competing labels by `aggregation_cluster_annotation_min_advantage`
  (default 0.05). Otherwise retain `Unassigned`, the candidate and a reason.
  Stability and GEM-well agreement are diagnostics, not assignment gates.
- Preserve positive and negative UCell marker signatures. Signed scores are
  clipped per cell before averaging for observed and control signatures,
  including marker-deletion diagnostics. Positive-only panels retain the
  faster exact rank-summary calculation.
- Export cluster evidence, marker-set support and competition, and clearer
  pre-doublet-filtering review plots. Unassigned clusters remain separate for
  doublet detection. Adjusted scores are not identity probabilities or validated
  biological error rates.

## Peak–gene analysis

- Allow exploratory discovery with few donors, including one donor, when
  aggregate counts and residual variation permit estimation. Adjust for donor
  and library depth; report conditional HC3 association statistics and separate
  donor-support diagnostics. State is no longer a nuisance-design term.
- Retain nonpromoter gene-body candidates and plot the donor/depth residuals
  used in scoring. These associations do not establish population-level
  replication or causal enhancer–gene relationships.
- Add SuSiE conditional peak prioritization and compact global-FDR reconstruction
  from chromosome branch inputs. Handle empty candidate sets with explicit
  plot outputs. SuSiE PIPs remain exploratory conditional prioritization.

## Migration

- Remove `aggregation_cluster_annotation_method` and
  `aggregation_allow_multiple_cell_types` from configuration. Signed marker
  panels remain supported; check the new evidence before interpreting labels.
- Recheck external references to numeric cluster IDs. Annotation changes can
  propagate through doublet grouping, retained cells and downstream analyses.
- Find plots and compatibility exports in numbered checkpoint folders, including
  `GEX_Seurat_object.3_GEX_QC` and `multimodal_Seurat_object.8_multimodal_QC`.
  Update scripts selecting the former target names. Old output files are not
  automatically deleted; consult current target metadata when reviewing results.
- Subgroup reprocessing is removed from the public workflow. Remove its six
  `aggregation_subgroups_*`/contamination-filter settings from public configs.
- Peak–gene results replace `cluster_robust_SE` with `association_SE` and add
  donor-support and inference-status fields. Recompute affected targets before
  comparing results; this is a method change, not just a speed improvement.

Principal recomputation boundaries (destination manifests provide full names):

cascading_target_breaking: cell_retention_tibble.GEX_input
cascading_target_breaking: QC_metric_manifest_tsv
cascading_target_breaking: PCA_clusters.GEX
cascading_target_breaking: LSI_clusters.ATAC
cascading_target_breaking: clusters_tibble_raw.WNN
cascading_target_breaking: metadata_w_cell_types_unfiltered_tibble.GEX
cascading_target_breaking: metadata_w_cell_types_unfiltered_tibble.ATAC
cascading_target_breaking: metadata_w_cell_types_tibble.WNN
cascading_target_breaking: peak_gene_correlation_donor_state_records.peak_gene_correlation.WNN
contained_target_breaking: GEX_Seurat_object.3_GEX_QC
contained_target_breaking: multimodal_Seurat_object.8_multimodal_QC

# Earlier untagged development (0.4.1.9000)

- Redesigned the searchable parameter reference with topic filters, visible
  defaults, readable expanded details, and a responsive layout.
- Integrated the public demo gallery refresh and seven numbered QC reviews,
  retaining a separate GEX PCA review before clustering.

- Added ordered QC checkpoints, per-GEM-well cutoff comparisons, and clearer
  PCA review plots, with an informational QC metric manifest.

- Fixed peak retention when no blacklist intervals overlap. Limited
  differential-analysis metadata dependencies to the configured donors and variables.

- Split chromVAR detail preparation from rendering and cache details per GWAS;
  added optional deterministic sampling for pseudobulk correlation estimation.

- Moved validation into a root testthat suite and simplified helper ownership,
  implementation diagrams, and reader documentation.

- Made fresh clones directly runnable after downloading the public demo data:
  configuration and local controller files are now regular files, and only the
  two `immune_human_2x` GEM wells and that aggregation are active by default.

- Added a restart-safe `setup-demo` Pixi task that installs the environment,
  downloads the two configured public inputs, and installs GitHub-only R
  dependencies behind one quickstart command.

- Replaced the public reaction-based API and configuration vocabulary with 10x
  Genomics GEM well terminology. Existing configurations must migrate to
  `cfg_GEM_wells.tsv`, `GEM_well_ID`, and the corresponding `GEM_well_*`
  parameter names documented in the migration table.

- Made manuscript benchmark wall-time results portable across targets stores
  and removed machine-specific figure paths ([#4](https://github.com/koefoeden/multiomeR/pull/4)).

- Replaced Seurat-backed cell-cycle module scoring with BPCells-native control-binned scoring, restored BPCells-native UCell marker scoring for signed marker sets, and added synthetic parity validation against Seurat and UCell reference implementations.

Initial public beta snapshot.
