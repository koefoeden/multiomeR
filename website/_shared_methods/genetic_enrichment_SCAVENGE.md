SCAVENGE-style trait-relevance scores are calculated from the nucleus-level z-scores on the WNN graph [@yu2022_scavenge].
They use multiomeR's own sparse-matrix implementation of the SCAVENGE propagation strategy rather than the reference package.
The reference package's dependencies predate the pipeline's R environment, and the local version propagates directly over the sparse WNN graph.
Nuclei whose one-sided normal-tail probability is at most 0.05 are seeds.
When they exceed the maximum seed fraction, 5% of nuclei by default, those with the highest z-scores are kept.
When no nucleus qualifies, all scores are zero.
A shared-nearest-neighbor graph is built from the WNN neighbors and made unweighted, with every remaining edge set to one, and nuclei without edges are excluded.
Seed signal is propagated by a random walk with restart probability 0.05 by default until convergence.
Nuclei that no seed reaches remain unscored.
Scores are capped at their 95th percentile, min–max scaled and multiplied by the mean z-score of the top 1% of nuclei.
Unlike the reference, the pipeline computes no permutation P-values.
The reference's degree-matched null scatters seeds across the graph, whereas observed seeds cluster within cell types.
Nearly every seeded group would therefore appear significant.
Scores are summarized per metadata group.
