WNN analysis uses the nuclei kept after ATAC processing that are present in both embeddings.
When small joint clusters are removed, each kept nucleus keeps its highest-scoring kept neighbors and its original modality weights, while the cluster labels are kept.
A nucleus left with too few neighbors has them selected again among the kept nuclei, with the same weights.
The UMAP is computed once, from the final neighbors.
The GEX, ATAC and WNN metadata therefore describe successively smaller sets of nuclei, not interchangeable versions of the Cell Ranger cell calls.
