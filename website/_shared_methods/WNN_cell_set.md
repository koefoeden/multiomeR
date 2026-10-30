The WNN graph is first built on all nuclei that pass the ATAC exclusions and are present in both embeddings, ATAC doublets included.
Its Leiden clusters of at least the minimum size are labelled from the WNN marker signatures, and these labels group the nuclei for ATAC doublet detection, with the nuclei of smaller clusters as one more group.
The graph is then restricted to the nuclei kept after ATAC processing, and clustered again.
When small joint clusters are removed in turn, the graph is restricted to the remaining nuclei, while the cluster labels are kept.
At each restriction, every kept nucleus keeps its highest-scoring kept neighbors and its original modality weights.
A nucleus left with too few neighbors has them selected again among the kept nuclei, with the same weights.
The UMAP is computed once, from the final neighbors.
The GEX, ATAC and WNN metadata therefore describe successively smaller sets of nuclei, not interchangeable versions of the Cell Ranger cell calls.
