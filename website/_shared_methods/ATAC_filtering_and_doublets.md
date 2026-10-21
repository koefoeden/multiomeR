ATAC processing starts from the nuclei kept after GEX processing.
The peak-calling groups and the consensus peak matrix are built on this cell set.
Per-nucleus quality-control metrics are checked against configured aggregation-level exclusions.
The fragments in peaks and their fraction of all fragments are counted from the stored fragments, each fragment once, even when its two cut sites fall in different peaks.
After ATAC clustering, a separate ATAC-based `scDblFinder` analysis follows per GEM well, with `dbr.sd = 1` [@germain2022_scdblfinder].
Instead of scDblFinder's internal per-well feature aggregation, peaks are grouped once into 50 feature groups.
The groups come from *k*-means clustering of the peaks' aggregation-wide LSI coordinates on dimensions 2–20 (loadings scaled by the singular values, as in scDblFinder).
The peaks of each group are summed with BPCells before scDblFinder's unchanged classifier runs.
This change limits memory use, but its scores are not expected to equal those from the internal aggregation.
That aggregation transforms every peak within each well, which can exhaust memory before classification.
ATAC doublets are removed at the nucleus and cluster level as for GEX.
