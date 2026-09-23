# Algorithmic implementations, deviations and validation

multiomeR reimplements a small number of reference algorithms so they can operate on the workflow's native matrices and graph state. The algorithms themselves are described on the primary-module and optional-module pages. This page records why each was reimplemented, where it deliberately differs from its reference, and what the executable validation establishes.

The evidence labels are intentionally narrow:

- **Reference-parity tested** means the repository and named reference implementation run on the same deterministic fixture and their returned values are compared directly.
- **Reference-similarity tested** means exact equality is not an appropriate contract, so predefined similarity thresholds are checked against the named reference implementation.
- **Algorithmically derived** means the implementation is checked against an independent mathematical result, not against another software implementation.

Passing these fixtures does not validate every dataset, parameter regime, approximate-neighbour realization, biological interpretation, or downstream target. The test suite contains only such reference comparisons. Each test asserts the reference version it was written against; versions are locked by `pixi.lock`. The [CI workflow](https://github.com/koefoeden/multiomeR/blob/main/.github/workflows/algorithm-validation.yaml) runs the complete suite when tests, relevant helpers, or the Pixi environment change.

Run the complete suite with `pixi run --use-environment-activation-cache test`. The narrower `pixi run --use-environment-activation-cache test-algorithm-validation` task runs only the slow UCell, AMULET, WNN, and SCAVENGE tests.

| Implementation | Evidence status | Reference | Validation contract |
|---|---|---|---|
| BPCells-native UCell | Reference-parity tested | UCell | Identical values, dimensions, and dimnames |
| BPCells-native AMULET | Reference-parity tested | scDblFinder | Identical metrics and loci, including order |
| Native WNN | Reference-similarity tested | Seurat | Modality-weight Spearman and neighbour-overlap thresholds per fixture |
| Sparse SCAVENGE propagation | Algorithmically derived and reference-parity tested | SCAVENGE source at `8ee8b173d965` | Closed-form propagation; identical seed samples, exceedance counts and significant-cell calls |
| Peak–gene donor-slope REML and Kenward–Roger kernels | Reference-parity tested | lme4 and pbkrtest | Coefficients, df and P-values within 1e-6; identical fit statuses |
| Peak–gene HC3 statistics and compact BH breakpoints | Reference-parity tested | sandwich; `stats::p.adjust()` | HC3 coefficients, errors and P-values within 1e-10; identical FDR per chromosome slice |

The peak–gene rows belong to the analyses in [Peak–gene correlation](methods_peak_gene_correlation.md); their tests are `test-peak-gene-hierarchical-parity.R` and `test-peak-gene-correlation-parity.R`.

## BPCells-native UCell scoring

**Reference algorithm.** [`UCell::ScoreSignatures_UCell()`](https://bioconductor.org/packages/release/bioc/html/UCell.html) calculates per-cell signature scores from descending feature ranks, caps ranks at `maxRank`, combines positive and negative signatures, and clips negative combined scores to zero.

**Reason for reimplementation.** The workflow keeps gene-by-cell counts in BPCells-backed matrices. Materializing the complete matrix or building a Seurat object solely for marker scoring would discard that storage contract, so multiomeR ranks bounded cell chunks and returns metadata-ready scores directly.

**Deliberate deviations and consequences.** Only one cell chunk is materialized at a time, and optional fork workers operate across chunks; this changes memory and execution behaviour but not the tested values. The helper returns a data frame instead of mutating a Seurat object. The target-level marker validator rejects configured genes missing from the reference, whereas the lower-level helper still exposes UCell's impute and skip modes. The production annotation reuses the chunked ranking helper and scores signed signatures per cell before aggregation, which costs more computation than a positive-only rank-summary shortcut.

**Implementation.** `calculate_BPCells_UCell_scores_from_matrix()` and `rank_UCell_count_chunk()` in `R/processing_GEX_helpers.R`; the cluster annotation built on them is in `R/cluster_annotation_helpers.R` and described in [Cell-type annotation and motif accessibility](methods_annotation_and_motifs.md#cell-type-annotation).

**Validation.** `tests/testthat/test-scoring-parity.R` compares signed signatures, with imputed and skipped missing genes, on a deterministic BPCells fixture and requires `identical()` values, dimensions, and dimnames. It also runs the production annotation path on unsigned, signed and negative-only signatures and compares per-cell scores, cluster means and matched-control summaries with reference scores averaged within clusters, within 1e-12, and checks that chunking and fork workers leave them unchanged. The production cell-cycle scorer is compared with `Seurat::CellCycleScoring()`: phases are identical and scores agree within 1e-6, because BPCells normalizes counts at lower floating-point precision.

```bash
pixi run --use-environment-activation-cache test-scoring-parity
```

## BPCells-native AMULET

**Reference algorithm.** [`scDblFinder::amulet()`](https://bioconductor.org/packages/release/bioc/html/scDblFinder.html) detects likely scATAC-seq doublets from the number of genomic loci covered by more than two fragments, after filtering fragment sizes and excluded regions and removing loci recurrently covered across many cells, and derives Poisson P-values with Benjamini–Hochberg correction.

**Reason for reimplementation.** The per-GEM-well workflow already stores Cell Ranger ATAC fragments as compressed BPCells directories. Passing the fragment file to scDblFinder materializes chromosome-scale `GRanges` objects; the local implementation streams the BPCells fragments and retains only one chromosome's selected fragments while calculating coverage runs.

**Deliberate deviations and consequences.** Only unique-fragment operation is supported: BPCells fragments do not retain Cell Ranger's PCR-duplicate counts, so requesting non-unique expansion fails explicitly. Cell Ranger's inclusive end coordinate is shifted back by one base to match scDblFinder's BED import. BPCells does not export its fragment iterator header, so the native helper mirrors that private interface, verifies the pinned BPCells commit before use, and compiles a small shared library in each worker. A BPCells upgrade must revalidate this interface and the parity fixture before updating the pin.

**Implementation.** `src/amulet_bpcells.cpp` iterates fragments and calculates coverage runs; `calculate_amulet_metrics_BPCells()` in `R/amulet_BPCells_helpers.R` holds the interface check, high-overlap filtering and statistics; the `amulet_metrics_tibble` target consumes the prefixed BPCells fragments.

**Validation.** `tests/testthat/test-amulet-parity.R` requires `identical()` results against scDblFinder for the bundled fragment-file metrics, the production call with prefixed Cell Ranger barcodes, a deterministic multi-chromosome loci fixture, and the corresponding full AMULET metrics, and requires an explicit error for PCR-duplicate expansion.

```bash
pixi run --use-environment-activation-cache test-amulet-parity
```

## BPCells-backed ATAC scDblFinder feature aggregation

**Reference algorithm.** With `aggregateFeatures = TRUE`, [`scDblFinder::scDblFinder()`](https://bioconductor.org/packages/release/bioc/html/scDblFinder.html) performs its own TF-IDF-based feature clustering and sums peaks into feature groups before artificial-doublet classification.

**Reason for adaptation.** Materializing and transforming every peak within each GEM well can exhaust worker memory before classification. The pipeline instead derives feature groups once from the aggregation's LSI loadings, sums the disk-backed peak matrix with BPCells, and passes the compact matrix to scDblFinder with `aggregateFeatures = FALSE`.

**Deliberate deviations and consequences.** This is not a reimplementation of the classifier, which runs unchanged. Global LSI-derived groups replace scDblFinder's per-GEM-well feature groups, so the aggregated matrix, and therefore scores and calls, can differ; exact parity is not expected. The GEX path calls `scDblFinder::scDblFinder()` directly through the same memory-bounding per-GEM-well wrapper.

**Implementation.** `get_feature_groups_from_LSI_loadings()` and `aggregate_BPCells_rows_by_group()` in `R/processing_GEX_helpers.R`, called from `extra_targets/ATAC_targets.R`.

**Validation scope.** The reference-parity fixtures do not establish equivalence for this path. Its settings are described in [Preprocessing and nucleus QC](methods_preprocessing_and_QC.md).

## Native weighted nearest neighbors

**Reference algorithm.** [`Seurat::FindMultiModalNeighbors()`](https://satijalab.org/seurat/reference/findmultimodalneighbors) constructs cell-specific modality weights from within- and cross-modality neighbourhood prediction, collects candidate neighbours across modalities, and selects a weighted multimodal neighbour set.

**Reason for reimplementation.** The pipeline already holds aligned RNA and ATAC embeddings and needs reusable neighbour indices, distances, and modality weights without creating a Seurat object. Native graph state also feeds UMAP, Leiden clustering, SCAVENGE, and the Seurat/Signac export.

**Deliberate deviations and consequences.** BPCells HNSW replaces Seurat's Annoy search, so candidate sets need not be identical. The helper does not expose Seurat's optional smoothing or cross-constant list, and BPCells builds the downstream SNN graph rather than Seurat `Neighbor` and `Graph` objects. The helper's `seed` argument is not consulted by the HNSW calls and does not control neighbour-search randomness. These choices can change weights, selected neighbours, SNN edges, clusters, and UMAP coordinates, so correlation and overlap, not exact equality, are the validation contract.

**Implementation.** `weighted_nearest_neighbors_BPCells()` in `R/processing_multimodal_helpers.R` with the small-SNN bandwidth kernel in `src/wnn_snn_bandwidth.cpp`; `extra_targets/WNN_targets.R` aligns the embeddings and wires the result into clustering, UMAP and metadata targets.

**Validation.** `tests/testthat/test-wnn-parity.R` compares deterministic RNA/ATAC fixtures with Seurat at production search settings: a 400-cell fixture with candidate range 200 and `k` of 20 and 50, and a 160-cell stress fixture with `k = 15` and candidate range 50. Each case must exceed its thresholds for modality-weight Spearman correlation and mean and lower-quartile neighbour-set overlap, and one- and two-thread results must be identical. This does not assert equality of selected neighbours, SNN weights, clustering, or UMAP.

```bash
pixi run --use-environment-activation-cache test-algorithm-validation
```

## Sparse SCAVENGE propagation and significance

**Reference algorithm.** [SCAVENGE at commit `8ee8b173d965`](https://github.com/sankaranlab/SCAVENGE/tree/8ee8b173d965009a696b2a590d5b17b28b7cf851) selects high chromVAR z-score seed cells, constructs a binary mutual-nearest-neighbour graph, performs a column-normalized random walk with restart, caps and rescales the propagation score into a trait relevance score (TRS), and uses degree-matched seed permutations to identify significant cells.

**Reason for reimplementation.** The reference package's dependency stack predates the pipeline's R/Bioconductor environment. multiomeR needs sparse propagation over its native SNN graphs and must avoid materializing a cell-by-permutation score matrix for large cell sets.

**Deliberate deviations and consequences.** The reference builds a mutual-kNN graph, whereas multiomeR uses the binary support of its BPCells-derived SNN graph; edge weights are discarded, but topology can differ. Only per-cell exceedance counts and the cluster medians needed downstream are retained from the permutations, and random walks rather than random-number generation are parallelized, so the sampled null does not depend on the core count. Seed and scale-factor helpers guarantee at least one selected cell for small inputs, the degree sampler handles one-cell strata explicitly, and the random walk has a maximum-iteration guard. Cluster-level permutation medians, add-one P-values, and Benjamini–Hochberg adjustment within each grouping column are pipeline extensions.

**Implementation.** `R/SCAVENGE_helpers.R` and `src/scavenge_random_walk.cpp`; `module_genetic_enrichment/SCAVENGE_graph_targets.R` builds the graph and result records, and `SCAVENGE_group_targets.R` combines summaries and plots.

**Validation.** `tests/testthat/test-scavenge-parity.R` uses a deterministic fixture with heterogeneous-degree graph blocks and nonuniform edge weights, so it also tests conversion to binary adjacency. The iterative random walk must match the closed-form solution

\[ s = r\left(I - (1-r)P\right)^{-1}p_0 \]

within 1e-10. Compact local reference functions reproduce the relevant SCAVENGE code at the pinned commit without installing its dependency stack; propagation and transformed scores must agree within 1e-12, and the degree-matched seed samples, streamed exceedance counts, empirical P-values, significant-cell calls and one- versus two-core results must be identical. This does not establish parity of graph construction, chromVAR inputs, or biological interpretation.

```bash
pixi run --use-environment-activation-cache test-algorithm-validation
```
