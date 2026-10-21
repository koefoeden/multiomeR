Fragments on the standard chromosomes are used for peak calling.
Peak-calling groups are defined by a configurable metadata column of the nuclei kept after GEX processing, by default their GEX clusters.
GEX clusters serve because ATAC clustering needs the peak matrix.
Calling peaks per group keeps those of rare populations, which a pooled call would dilute [@stuart2021_signac; @granja2021_archr].
Groups larger than 50,000 nuclei are randomly downsampled for peak discovery, which limits the runtime and memory of peak calling.
Each group's discovery nuclei are split into two random halves.
Both draws are seeded by the group's name, so a group's halves do not change when other groups are added or removed.
Peaks are called in each half with either MACS3 or the BPCells tile-based caller [@zhang2008_macs; @parks2025_bpcells].
MACS3 runs on each half's fragment left ends, written in one pass per GEM well, with a −75-bp shift, a 150-bp extension, a 150-bp maximum gap and summit calling.
Only one of each fragment's two insertions thus enters the pileup, as when MACS3 reads a fragment file.
The gap is fixed because MACS3's default, the mean length of the first ten input lines, would be 1 bp for these single-base cut sites.
The tile-based caller runs on the stored fragments, using 500-bp tiles and an FDR of 0.01.
Consensus peaks follow the ArchR strategy of iterative fixed-width peak construction [@granja2021_archr].
Summits are extended by 250 bp to 500-bp peaks, and peaks overlapping any blacklisted region of the reference genome are removed [@amemiya2019_encode_blacklist; @ogata2023_excluderanges].
Overlapping peaks are collapsed by summit significance, first within each half and then across the two halves.
Only peaks called in both halves are kept.
Like ArchR's pseudo-replicate rule [@granja2021_archr], this removes weak peaks called in only one half.
Each group keeps its strongest peaks, up to 500 per discovery nucleus and 150,000 in all, as ArchR does by default [@granja2021_archr].
The group-level sets are then combined into one non-overlapping consensus set.
Here, overlapping peaks are ranked by the quantile of their summit significance within their group.
The caps and within-group quantiles keep large or deeply sequenced groups from dominating the consensus [@corces2018_cancer_atac].
For all nuclei kept after GEX processing, BPCells counts each fragment once in every peak that holds at least one of its two Tn5 cut sites.
This paired-insertion counting, also Signac's default, counts a fragment with both cut sites in a peak once rather than twice, and ignores fragments that only span a peak [@miao2024_pic; @stuart2021_signac].
Signals that describe positions, such as peak-calling pileups, TSS enrichment and coverage tracks, still use cut sites.

ATAC dimensional reduction starts from the consensus peak matrix of the nuclei that pass the aggregation-level ATAC exclusions.
Signac's TF-IDF method 1, with a scale factor of 10,000, normalizes this matrix [@stuart2021_signac].
A BPCells singular-value decomposition of the result gives the latent semantic indexing (LSI) embeddings [@parks2025_bpcells].
It reads the matrix in a fixed number of blocks, so its result does not depend on the number of threads.
When memory allows, the transformed matrix or its counts are held in memory during the decomposition.
At least 50 components are computed, of which the configured ones are used downstream.
The pipeline does not test them for correlation with sequencing depth.
Instead, the default dimension settings leave out the first component, which often reflects depth [@stuart2021_signac].
