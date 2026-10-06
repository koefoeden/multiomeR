# multiomeR 1.2.0 (unreleased)

Adds the public configuration, figures and results of the multiomeR manuscript,
lowers the memory use of CellBender input, and lets targets stores move between
disks. Each configuration directory now carries its own crew controllers.

## Configuration

- Read `crew_controllers.R` from the selected configuration directory, so each
  directory carries controllers that suit its aggregations; the default file
  moves to `configuration/crew_controllers.R`.
- Use CellBender counts whenever a GEM well gives a CellBender H5 file. The
  `GEM_well_add_cellbender` column is no longer needed; a leftover column is
  ignored.
- Move `mixed_human_31x` and the nested `comparison_1x` to `comparison_20x`
  aggregations into the new `configuration_dev/`, which builds the manuscript's
  results in its own targets store (`TAR_PROJECT=dev`, store `outputs_dev`),
  with Slurm workers on one node. `configuration/` keeps the demo and the small
  examples, and validation now covers whichever validation aggregations the
  selected directory defines.

## Fixes and performance

- Read CellBender matrices in blocks of 2,000 barcodes and keep only the gene
  features, so a 34,498-nucleus GEM well peaks at 5.2 GB instead of 17.1 GB.
  The written matrices are byte-identical.
- Load Roadmap chromHMM annotations per aggregation, so activating another
  aggregation or selecting another configuration directory no longer rebuilds
  them.
- Record unresolved, project-relative paths for opened BPCells matrices and
  fragments and for downloaded reference files, so a targets store keeps
  working after it moves, for example from node-local disk to shared storage.

## Manuscript

- Add `pipeline_manuscript_figures/`, which builds Figure 1, Supplementary
  Figure S2 and Supplementary Table S4 of the manuscript from public data in
  `configuration_dev/`, with the results and their data under `results/`. It
  replaces `manuscript_figures/`.
- In Figure 1, spell cell types alike in panels C to E and show protein-coding
  genes in panel D with labels that no longer overlap.

## Documentation

- Remove the dropped HC3 scan from the workflow overview figure.
- Complete and case-protect bibliography entries.
- Add `CITATION.cff`, which GitHub shows as "Cite this repository" and Zenodo
  uses for the authors of archived releases.
- Ask users, in the README and on the manual's home page, to cite multiomeR,
  targets and BPCells, and the methods behind the results they report;
  `CITATION.cff` lists targets and BPCells as references.

## Migration

- A custom configuration directory needs its own `crew_controllers.R`: copy
  `configuration/crew_controllers.R` and adapt it. Move any changes you made to
  the former root `crew_controllers.R` into it.
- Select `configuration_dev/` to rebuild `mixed_human_31x` or the comparison
  aggregations.
- Target invalidation: targets that open BPCells matrices or fragments
  (`GEX_counts_BPCells_matrix`, `fragments_w_prefix_bpcells`,
  `aggregated_counts_BPCells_matrix.GEX`, `consensus_peak_BPCells_matrix.ATAC`,
  `motif_family_accessibility_BPCells_matrix.ATAC`,
  `gene_score_archr_BPCells_matrix.ATAC` and the pseudobulk count matrices)
  rerun with unchanged data, and their consumers rerun with them. Smaller
  reruns that reproduce their results: `GEM_well_config_tsv`,
  `GEM_well_metadata_tibble`, `GEM_well_annotation_metadata_tibble`,
  `GEX_counts_BPcells_matrix_dir`, `chromHMMs_list_proj.ATAC` and the five
  downloaded reference-file targets. The Figure 1 targets affect only
  `configuration_dev/`.

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
[RELEASES.md](RELEASES.md).

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
