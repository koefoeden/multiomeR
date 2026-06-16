Approach to generating graphs that are both human readable and still somewhat reflect the actual graph structure:

1) Select a "topic" of interest, i.e. parallel processing, GEX-aggregation, ATAC, WNN, differential analyses, psbulked genetic enrichment and single-cell genetic-enrichment.
2) Find downstream target(s) that capture this topic meaningfully
3) Use tar_mermaid() to generate the draft .mmd file
4) use 'python mermaid_gen_patterns.py' to draft a list of all current nodes
5) Iteratively edit the pattern-nodes-list while running 'python mermaid_bypass_patterns.py' to discard uninformative notes
6) Perform final clean-up, i.e. remove aggregations-specific suffixes etc.

## Source-level graph labels

Curated graph membership is recorded in target descriptions with repeated
single-value tags:

```r
description = "Harmony-corrected SCTransform GEX PCA embeddings. [part_of_graph:GEX] [part_of_graph:WNN] [part_of_graph:seurat_export]"
```

`[part_of_graph:<graph_id>]` means that the target should remain visible as an
explicit node in that named graph after graph-pruning code bypasses
uninformative dependencies. It is not a target-run selector, a checkpoint, or a
global importance label.

Current graph IDs are:

- `parallel`
- `GEX`
- `ATAC`
- `WNN`
- `seurat_export`
- `differential_analyses`
- `genetic_enrichment_single_nucleus`
- `genetic_enrichment_pseudobulk`
