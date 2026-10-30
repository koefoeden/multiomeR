Optional gene activity scores summarize the accessibility around each gene in each nucleus.
They follow ArchR's gene-score model as BPCells implements it [@granja2021_archr; @parks2025_bpcells].
The genome is cut into 500 bp tiles, and blacklisted tiles are left out.
Each gene is extended 5 kb upstream of its start, and it takes the tiles of this extended body and those up to 100 kb on either side, stopping at a neighbouring gene.
A tile's weight is exp(−*d*/5 kb) + exp(−1), where *d* is its distance from the extended body.
Shorter genes get up to five times the weight, as in ArchR, so that short and long genes score on a similar scale.
A nucleus's fragment counts are capped at 4 per tile, and its weighted gene sums are scaled to total 10,000.
The weights come from the Ensembl genes of the aggregation's reference, and the scores cover the nuclei that pass ATAC filtering.
They are used for one check: whether the GEX-derived cell types have accessibility at their marker genes.
Accessibility around a gene does not show that it is transcribed, and its window can hold the regulatory regions of neighbouring genes.
