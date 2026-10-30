ATAC processing starts from the nuclei kept after GEX processing.
The peak-calling groups and the consensus peak matrix are built on this cell set.
Per-nucleus quality-control metrics are checked against configured aggregation-level exclusions.
The fragments in peaks and their fraction of all fragments are counted from the stored fragments, each fragment once, even when its two cut sites fall in different peaks.
A separate ATAC-based `scDblFinder` analysis then follows per GEM well, with `dbr.sd = 1`, for all nuclei that pass these exclusions [@germain2022_scdblfinder].
It groups the nuclei by the WNN cell types of a joint graph that still holds the ATAC doublets, described with the WNN cell set.
In our tests, these groups gave as many calls as ATAC cell types, with slightly more agreement between seeds, and ATAC annotation no longer delays the joint graph.
Instead of scDblFinder's internal per-well feature aggregation, peaks are grouped once into 50 feature groups.
The groups come from *k*-means clustering of the peaks' aggregation-wide LSI coordinates on the configured dimensions within 2–20 (loadings scaled by the singular values, as in scDblFinder).
The peaks of each group are summed with BPCells before scDblFinder's unchanged classifier runs.
This change limits memory use, but its scores are not expected to equal those from the internal aggregation.
That aggregation transforms every peak within each well, which can exhaust memory before classification.
Nuclei called as ATAC doublets are removed.
All nuclei of an ATAC Leiden cluster are removed when more than a configured share of them are called doublets, but only if RNA depth or AMULET also points to doublets.
The cluster's median RNA depth, relative to its GEM wells' nuclei not called doublets, must reach a configured ratio, or AMULET must call more than that share of its nuclei doublets.
A rare population that the ATAC calls flag without doublet-like RNA, such as mast cells, is thus kept.
