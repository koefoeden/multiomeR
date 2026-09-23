# Cell-type annotation and motif accessibility

This chapter covers marker-signature cluster annotation and motif-family accessibility. The BPCells-native UCell scorer is compared with UCell in [Algorithmic implementations](algorithm_validation.md#bpcells-native-ucell-scoring). The tables follow the layout in [Methods and parameter tables](implementation_conventions.md#methods-and-parameter-tables).

## Cell-type annotation

{{< include _shared_methods/cell_type_annotation.md >}}

### Signature scoring

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Signatures | Marker genes per label | Configurable |  | [`aggregation_GEX_marker_genes`](../parameters.html#aggregation_GEX_marker_genes) |
| Ranks | Rank cap | Fixed | `min(1500, number of genes)` | `prepare_cluster_UCell_controls()` |
| Ranks | Direction and ties | Fixed | descending raw counts per nucleus; ties averaged | `rank_UCell_count_chunk()` |
| Scores | Signed signatures | Fixed | positive-gene score minus negative-gene score, clipped with `pmax(0, score)` per nucleus before averaging | `score_signed_UCell_cells()` |
| Scores | Partitions annotated | Fixed | GEX `PCA_harmony_SNN_cluster` before doublet removal; ATAC `LSI_harmony_SNN_cluster` before ATAC doublet removal; WNN `WNN_harmony_SNN_cluster` | `cluster_UCell_evidence.GEX`, `cluster_UCell_evidence.ATAC`, `cluster_UCell_evidence.WNN` |

### Matched-control cluster annotation

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Reference | Nuclei and counts | Fixed | up to 50 nuclei per GEM well, sampled from the pre-doublet-filter GEX nuclei; full aggregated GEX counts | `cluster_UCell_controls.GEX`, `prepare_cluster_UCell_controls()` |
| Controls | Random mappings | Fixed | 999 | `build_UCell_controls()` |
| Controls | Matching coordinates | Fixed | `log1p(abundance × 1e4)` and `asin(sqrt(detection))`, each standardized; Euclidean distance | `build_UCell_controls()` |
| Controls | Candidate pool and draw | Fixed | `neighbours = 50`; pool of the max(4 × neighbours, 2 × marker genes) nearest eligible genes, which is 200 for panels of up to 100 marker genes; each marker drawn uniformly from its `neighbours` nearest candidates not yet used in the replicate | `build_UCell_controls()` |
| Controls | Random seed | Fixed | 20260910 | `prepare_cluster_UCell_controls()`, `build_UCell_controls()` |
| Adjusted score | Background | Fixed | 0.95 quantile of the label's matched-control scores | `score_UCell_group_evidence()` |
| Assignment | Advantage rule | Fixed | `min(best, best − second)`; requires best \> 0 and no exact tie | `assign_UCell_cluster_evidence()` |
| Assignment | Minimum advantage | Configurable |  | [`aggregation_cluster_annotation_min_advantage`](../parameters.html#aggregation_cluster_annotation_min_advantage) |
| Diagnostics | Marker deletion | Fixed | each marker of the leading label, with its matched control, omitted in turn (labels with at least 2 markers); stability = fraction of deletions keeping a positive advantage at least the minimum | `score_cluster_UCell_summaries()`, `evaluate_cluster_UCell_evidence()` |
| Diagnostics | Cell-deletion blocks | Fixed | 10 blocks stratified by cluster and GEM well, each left out once; stability reported when at least 2 blocks are assessable | `make_annotation_blocks()`, `summarize_cluster_UCell_counts()`, `evaluate_cluster_UCell_evidence()` |
| Diagnostics | Marker detection | Fixed | positive markers detected in at least 10% of a cluster's nuclei | `score_UCell_group_evidence()` |
| Diagnostics | GEM-well agreement | Fixed | GEM wells with at least 25 nuclei in the cluster; agreement reported when at least 2 wells are assessed | `score_cluster_UCell_summaries()`, `evaluate_cluster_UCell_evidence()` |
| Plots | Marker-set order | Fixed | `hclust(method = "ward.D2")` on Euclidean distances between adjusted-score profiles | `plot_UCell_annotation_dot()` |

## Motif families and motif accessibility

{{< include _shared_methods/motif_accessibility.md >}}

| Step | Setting | Status | Value | Source |
|---|---|---|---|---|
| Families | JASPAR 2026 CORE vertebrate files | Fixed | 233 familial root motifs; family membership of 1,019 motifs | `resources/JASPAR2026_vertebrate_familial_root_motifs.tf`, `resources/JASPAR2026_vertebrate_motif_families.tsv` |
| Families | Configured transcription factors | Configurable |  | [`aggregation_ATAC_marker_TFs`](../parameters.html#aggregation_ATAC_marker_TFs) |
| Families | Symbol resolution | Fixed | case-insensitive TF name, or `TF__motifID` for symbols in several families | `resolve_marker_motif_families()` |
| Scanning | Motif matrices | Fixed | PFMs with a uniform 0.25 background on the `+` strand; all-zero columns dropped | `read_JASPAR_familial_root_PFMatrixList()` |
| Scanning | Match thresholds | Fixed | package defaults: `p.cutoff = 5e-5`, `bg = "subject"`, `w = 7` | `get_motif_matrix_from_peak_ranges()`, `motifmatchr::matchMotifs()` |
| Scanning | Genome | Fixed | BSgenome UCSC hg38, mm10 or mm39, matching the reference | `get_chromVAR_genome_obj()` |
| chromVAR | Nuclei | Fixed | ATAC nuclei after doublet removal | `chromVAR_obj.ATAC` |
| chromVAR | Peaks, bias and expectation | Fixed | zero-count peaks removed; `betterChromVAR::addGCBias()` with missing values set to 0; expectation = mean count per peak over nuclei | `get_chromVAR_obj_from_peak_matrix()`, `get_chromVAR_peak_expectation()` |
| chromVAR | Background | Fixed | package defaults: default background bins, no shrinkage (`shrinkage = "none"`) | `betterChromVAR::getBackgroundBins()`, `betterChromVAR::computeBackgrounds()` |
| chromVAR | Deviations | Fixed | analytic deviations and z-scores; package default `denominator = "global"` | `compute_chromVAR_annotation_chunk_result()`, `betterChromVAR::computeDeviationsAnalytic()` |
| Summaries | Per-cell-type test | Fixed | `BPCells::marker_features(method = "wilcoxon")` on z-scores by ATAC cell type; BH across all family-by-group tests; effect = difference in mean z-score | `get_marker_motif_family_accessibility_from_chromVAR_BPCells_z_scores()` |
