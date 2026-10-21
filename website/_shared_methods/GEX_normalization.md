For each aggregation, the quality-controlled gene-expression matrices are combined with their GEM-well-prefixed barcodes.
Genes with 50 or fewer total counts are excluded before dimensional reduction.
Normalization and dimensional reduction use either the disk-backed BPCells Pearson-residual workflow or Seurat `SCTransform` [@parks2025_bpcells; @hafemeister2019_sctransform].
With the BPCells backend, analytic Pearson residuals are computed with one overdispersion parameter (θ = 100) for all genes and clipped to ±√*N* for *N* nuclei [@lause2021_pearson_residuals].
Per-gene estimates of θ would absorb differences between cell types and flatten marker-gene residuals, and the √*N* clip keeps markers of rare populations [@lause2021_pearson_residuals].
Configured cell-level covariates can be regressed out.
The 3,000 genes with the highest residual variance are kept, and a truncated singular-value decomposition gives at least 50 principal components, of which the configured ones are used downstream.
It reads these genes' counts from memory in a fixed number of blocks, so its result does not depend on the number of threads.
The Seurat backend applies SCTransform with its v2 regularization to the same cells and genes [@choudhary2022_sctransform_v2].
It keeps the same number of variable features and computes the principal components from the residual matrix.
A cell-cycle difference score can be regressed out but is not by default.
Its S and G2M scores are computed from log-normalized counts with the 2019 Seurat cell-cycle gene sets [@tirosh2016_melanoma].
