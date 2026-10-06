# Pipeline manuscript figures

`targets.R` extends the pipeline graph with the figures and tables of the
multiomeR manuscript that come from public data. It uses the aggregations in
[`configuration_dev/`](../configuration_dev/README.md) and writes its results
below `results/`:

- `1.png`: Figure 1, selected outputs of `mixed_human_31x`.
- `S2.png` and `data/seurat_signac_comparison_resources.tsv`: Supplementary
  Figure S2, resource use of multiomeR and a conventional Seurat/Signac workflow
  on `comparison_1x` to `comparison_20x`.
- `S4.md` and `data/algorithm_parity/`: Supplementary Table S4, the
  reimplemented algorithms against their references on `comparison_5x` objects.

## Run

Select `configuration_dev/` and its store, then build the results from the
repository root inside a Slurm allocation (see the configuration's README):

```r
targets::tar_make(
  script = "pipeline_manuscript_figures/targets.R",
  names = tidyselect::any_of(c("figure_1.mixed_human_31x", "supplementary_figure_S2",
    "supplementary_table_S4.comparison_5x"))
)
```

The conventional Seurat/Signac chain needs fragtk; install it into the
environment with `pixi run install-fragtk`. The graph builds only the
dependencies of the selected results. S2 reads the recorded metadata of both
workflows rather than depending on their targets, so build it after they have
finished; before that it stops with missing runtimes. Its measurements come
from the run that builds the workflows, so build the results afterwards from a
copy of that store if needed. S2 also needs the per-job resource history of
[slurm-monitor](https://github.com/koefoeden/slurm-monitor), read from
`SLURM_MONITOR_HISTORY_DIR` (default `~/slurm_monitor_history`): schema 2
files (`<job>.v2.tsv`) with the columns `name`, `timestamp`, `target`,
`ram_cur_gib`, `cpu_cur` and `scope`, one row per sample of a crew worker job.

## Figure 1

Five panel targets read pipeline data, never saved pipeline plots; a sixth
combines them and the file target exports an 8 × 11 inch page at 300 dpi.
Selections are deterministic but follow the data.

- **A, B:** the eligible, estimable, non-promoter peak–gene link with the lowest
  hierarchical P-value (ties: cell type, gene, peak). A shows the peak, gene
  bodies, focal-cell-type coverage and a schematic link with 25 kb flanks; B the
  aggregate residuals coloured by donor with a descriptive regression line. The
  ranking is not an FDR claim, and aggregates are not independent donors.
- **C:** up to six traits with a raw chromVAR Z ≥ qnorm(0.95), ordered by their
  strongest cell type, across all cell types; stars mark unadjusted upper-tail
  P ≤ 0.05 and P ≤ 0.01.
- **D:** for the two strongest cell-type–trait combinations, the locus with the
  largest absolute contribution to the deviation, as aligned facets of variant
  contributions (point area: PIP), protein-coding genes, focal coverage and
  consensus peaks.
- **E:** SCAVENGE trait-relevance scores of the same two traits on the WNN UMAP,
  with a shared colour scale, unclipped scores and cell-type labels.

## Supplementary Figure S2

The conventional chain is a compact `{targets}` implementation of a standard
Seurat v5/Signac v2 analysis on BPCells-backed assays, in the same graph, store
and controllers as multiomeR: it writes each GEM well's gene-expression counts
to a BPCells directory and runs SCTransform, PCA and RNA clustering on their
joined matrix, repeats the same MACS3, fixed-width, blacklist and consensus-peak
algorithm in one target, quantifies the peaks per GEM well with
`Signac::FeatureMatrix()` and its default fragtk backend into BPCells
directories, and runs Signac's TF-IDF/LSI on them before building the WNN graph,
UMAP and clusters. Seurat's kernels use the cores of their Slurm worker, and
peak calling and quantification run one `future` worker per core.

`comparison_1x` to `comparison_20x` take the first 1, 2, 5, 10 and 20 GEM wells
of `mixed_human_31x` without per-well QC filters, CellBender or doublet removal,
so both workflows keep every called nucleus and yield identical consensus peaks.
Shared cluster-fragment preparation is excluded from both workflows, and
cell-type annotation, which the conventional chain does not perform, from
multiomeR. Critical paths come from recorded target runtimes, with dynamic
branches treated as concurrent; CPU time and RAM use multiply each target's
runtime by its job's mean sampled CPU use and unreclaimable memory, which leaves
out page cache, and jobs too short to be sampled get the median of the sampled
jobs; disk space counts retained objects and files inside the store.

## Supplementary Table S4

Each native and reference implementation runs as its own six-core target in a
fresh R process, which records the wall time and peak process-tree RSS, on
fixtures built from `comparison_5x` targets: UCell, WNN, ATAC scDblFinder with
and without clusters, AMULET and SCAVENGE propagation. The table target stops
when a parity threshold fails.
