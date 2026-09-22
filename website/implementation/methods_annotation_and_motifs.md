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