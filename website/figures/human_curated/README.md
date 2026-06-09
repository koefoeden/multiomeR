Approach to generating graphs that are both human readable and still somewhat reflect the actual graph structure:

1) Select a "topic" of interest, i.e. parallel processing, GEX-aggregation, ATAC, WNN, differential analyses, psbulked genetic enrichment and single-cell genetic-enrichment.
2) Find downstream target(s) that capture this topic meaningfully
3) Use tar_mermaid() to generate the draft .mmd file
4) use 'python mermaid_gen_patterns.py' to draft a list of all current nodes
5) Iteratively edit the pattern-nodes-list while running 'python mermaid_bypass_patterns.py' to discard uninformative notes
6) Perform final clean-up, i.e. remove aggregations-specific suffixes etc.
