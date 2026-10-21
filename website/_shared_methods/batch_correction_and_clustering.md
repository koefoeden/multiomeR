Harmony correction can be applied separately to the selected GEX principal components and ATAC LSI dimensions, with at most 25 iterations and `lambda = 1` [@korsunsky2019_harmony; @patikas2026_harmony2].
It is off by default, since correcting a variable can also remove biology confounded with it [@luecken2022_integration_benchmark].
The configured metadata columns are combined into one interaction batch factor, and the ATAC correction can add its own columns.
Nuclei with missing covariate values are removed with a warning.

For each modality, the nearest neighbors of every nucleus are found once on the selected dimensions, by hierarchical navigable small world (HNSW) search with cosine distance [@malkov2020_hnsw; @melville2024_rcpphnsw].
The search index is built on one thread, so a rebuild finds the same neighbors, while the search itself runs in parallel.
Clustering and UMAP share these neighbors.
They are converted to a shared-nearest-neighbor (SNN) graph with Jaccard weights, removing edges with weights below 1/15.
The graph is clustered with the modularity objective of the Leiden algorithm in ten iterations [@parks2025_bpcells; @traag2019_leiden].
Clusters are renumbered by size, and clusters of fewer than 100 nuclei by default are removed from the GEX, ATAC and WNN clusterings.
UMAP embeddings are computed with uwot from the same neighbors [@mcinnes2018_umap; @melville2026_uwot].
GEX and ATAC layouts start from the first two principal components of the clustering dimensions, and WNN layouts from a spectral embedding of the neighbor graph.
uwot's batch optimization gives the same layout for any number of threads.

Most clusters below the size threshold hold a single nucleus or a few nuclei that share too few neighbors with any larger cluster, so the filter mainly removes such outlier nuclei.
Rare populations smaller than the threshold are removed as well.
