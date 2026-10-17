# Pipeline manuscript figures

`targets.R` extends the pipeline graph with the figures and tables of the
multiomeR manuscript that come from public data. It uses the aggregations in
[`configuration_dev/`](../configuration_dev/README.md) and writes its results
below `results/`:

- `1.png`: Figure 1, selected outputs of `mixed_human_31x`, and
  `data/figure_1_aggregation.tsv`, `data/figure_1_link.tsv` and
  `data/figure_1_traits.tsv`, the values the manuscript quotes from it.
- `S3.png` and `data/cell_type_trait_enrichment.tsv`: Supplementary Figure S3,
  the cell-type deviations of every GWAS study of `mixed_human_31x`, of which
  Figure 1D shows a selection.
- `data/table_S3.tsv`: Supplementary Table S3, the GWAS studies of
  `mixed_human_31x` with their Open Targets records and credible sets, and the
  credible-set loci and variants that overlap consensus peaks and the peaks they
  weight.
- `S2.png` and its data, Supplementary Table S6, on what multiomeR costs and
  saves: `data/seurat_signac_comparison_resources.tsv` compares multiomeR with a
  conventional Seurat/Signac workflow on `comparison_1x` to `comparison_31x`;
  `data/benchmark_sizing.tsv` and `data/benchmark_sizing_stages.tsv` measure
  multiomeR with its defaults on `sizing_1x` to `sizing_91x`;
  `data/benchmark_revisions.tsv` measures what five revisions of `revision_31x`
  rerun; `data/benchmark_simulation_check.tsv` compares the simulated wall
  time with measured runs; and `data/benchmark_worker_attribution.tsv` gives the
  worker time that no target claims while each measurement ran. Its inputs are `data/benchmark_runs/`, the targets
  each measured run built, and `data/benchmark_allocations.tsv`, the worker
  cores and RAM that each counted target requested in the run that measured it.
- `S4.md` and `data/algorithm_parity/`: Supplementary Table S4, the
  reimplemented algorithms against their references on `comparison_5x` objects.

## Run

Select `configuration_dev/` and its store, then build the results from the
repository root inside a Slurm allocation (see the configuration's README):

```r
targets::tar_make(
  script = "pipeline_manuscript_figures/targets.R",
  names = tidyselect::any_of(c("figure_1.mixed_human_31x", "supplementary_figure_S3.mixed_human_31x",
    "supplementary_table_S3.mixed_human_31x", "supplementary_figure_S2", "supplementary_table_S4.comparison_5x"))
)
```

The conventional Seurat/Signac chain needs fragtk; install it into the
environment with `pixi run install-fragtk`. The graph builds only the
dependencies of the selected results.

S2 and Table S6 read recorded runtimes, not the targets they measure, so build
them after the benchmark runs have finished. Before that they stop with a
missing run or runtime. `targets` allows a target to read its own store's
metadata only as a local `format = "file"` target, and the measured
aggregations must be active, because inactive ones are absent from
`tar_network()`. S2 also needs the per-job resource history of
[slurm-monitor](https://github.com/koefoeden/slurm-monitor) in
`SLURM_MONITOR_HISTORY_DIR` (default `~/slurm_monitor_history`), which it reads
with `slurm_monitor history samples`: the schema 2 `timestamp`, `raw_job_id`,
`name`, `elapsed`, `target`, `ram_cur_gib`, `cpu_cur`, `ram_lim_gib`,
`cpus_alloc` and `scope` of each sample of a crew worker job.

Each benchmark run builds its targets alone on one node, from a clean checkout
with the store on node-local disk. `record_run.py` then records the targets the
run built, with their runtimes, completion times, bytes and errors, in
`data/benchmark_runs/<run>.tsv`. A later run's
record of a target replaces an earlier one. Later runs overwrite the store's
metadata, so the records, not the store, hold the measurements.

## Figure 1

Six panel targets read pipeline data, never saved pipeline plots; a seventh
combines them and the file target exports an 8 × 9.5 inch page at 300 dpi.
Selections are deterministic but follow the data.

- **A:** for each WNN cell type, the configured marker gene whose expression
  and promoter accessibility, in cell-type pseudobulks summed over donors and
  scaled across cell types, are most specific to it.
- **B, C:** the eligible, estimable, non-promoter peak–gene link with the lowest
  hierarchical P-value (ties: cell type, gene, peak). B shows the peak, gene
  bodies, focal-cell-type coverage and a schematic link with 25 kb flanks; C the
  aggregate residuals coloured by donor with a descriptive regression line. The
  ranking is not an FDR claim, and aggregates are not independent donors.
- **D:** broadly supported cell-type–study enrichments across all assigned
  cell types: chromVAR q < 0.05 and a top locus below half of the summed
  contribution, so that no pair rests on one or a few loci; the configured
  studies already have at least 20 effective peaks, (Σw)²/Σw². In decreasing z, a pair is shown unless its
  study is shown already or its cell type shows a study of the same category.
  Stars mark q < 0.05.
- **E:** for one of these pairs, colorectal cancer in colonocytes, the locus
  with the largest absolute contribution to the deviation, as aligned facets of variant contributions (point area: PIP),
  protein-coding genes, focal coverage and consensus peaks.
- **F:** SCAVENGE trait-relevance scores of the same trait on the WNN UMAP,
  with unclipped scores and the WNN cell-type labels of D.

Supplementary Figure S3 shows every study across the assigned cell types with
the module's unscaled deviation heatmap and GWAS metadata tracks, on Figure 1D's
colour scale and with its stars; larger deviations take the scale's end colours.

## Supplementary Figure S2 and Table S6

S2 answers three questions: what the pipeline machinery costs or saves on the
same algorithms (A-E), whether multiomeR with its defaults runs on a given
data size and what it costs (F-H), and what revising an analysis costs (I).
Sizes and resources use logarithmic axes, so equal gaps are equal ratios.

### Like for like

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

`comparison_1x` to `comparison_31x` take the first 1, 2, 5, 10, 20 and all 31
GEM wells of `mixed_human_31x` without per-well QC filters, CellBender or
doublet removal, so both workflows keep every called nucleus and yield identical
consensus peaks. multiomeR runs its Seurat SCTransform backend there, so both
workflows use the same algorithms. It is measured twice: up to its WNN cell
metadata, and with the multimodal Seurat/Signac object it exports from it,
which the figure shows because it matches the conventional chain's output.
Shared cluster-fragment preparation is excluded from all workflows, and analyses
that the conventional chain does not perform from both multiomeR measurements:
cell-type annotation, doublet detection, QC metrics, gene and peak annotation,
chromatin-state projection, motif analysis, and the ATAC-only and uncorrected
RNA embeddings. These measurements come from the store's executions before the
benchmark runs (`store_before_runs.tsv.gz`), so later rebuilds of shared
targets leave them unchanged.

A third line measures the same multiomeR result with its default backend,
BPCells_native, which computes Pearson residuals on disk: `sizing_1x` to
`sizing_31x` repeat the comparison aggregations with only that setting changed,
and `sizing_91x` adds 60 ENCODE multiome series, for 1.16 million nuclei. Its
per-well steps are those of the comparison aggregations. SCTransform cannot run
at that size in either workflow: Seurat converts the joined count matrix to a
`dgCMatrix`, which holds at most 2^31 - 1 non-zero entries, and `sizing_91x` has
2.17 billion. F shows each size's RNA non-zero entries against that limit, the
genes detected in each called nucleus summed from Cell Ranger's per-barcode
metrics. They equal the count in Seurat's error.

### Sizing

The sizing measures cover multiomeR's full analysis of each `sizing_`
aggregation: its own targets, every target of its GEM wells, such as their QC
plots, and every upstream target, including shared setup such as reference
downloads, which a first run builds once. Rows of
`benchmark_sizing.tsv` give the like-for-like subsets of both multiomeR
workflows, the primary module and, for `sizing_31x` and `sizing_91x`, the
primary module with the genetic-enrichment, differential-analysis and peak-gene
modules. `benchmark_sizing_stages.tsv` splits the full analysis into setup,
per-well steps, the primary module's aggregation steps and each optional module,
which takes the upstream targets that only it needs, such as the Open Targets
datasets. CellBender runs upstream of the pipeline, in cellranger-wrapper, and
is left out.

### Revisions

`revision_31x` copies `mixed_human_31x`, the configured study with QC,
CellBender, doublet removal, tuned markers and all three optional modules. Its
full build is the `revision_base` run plus the setup and per-well targets it
reused from the store. Five runs then change it one step further each, so each
run measures one change on top of the previous ones: the WNN clustering
resolution from 0.4 to 0.6, an added mast-cell marker set, an added dilated
cardiomyopathy GWAS, age added to the cardiomyocyte sex model, and an added
liver GEM well with its donor. The analysed configuration lacks that well, so
its recorded targets take the dependencies of the same targets of another well.
Each run's targets, critical path, CPU time and simulated wall time are
compared with the full build's. The revision runs and `sizing_91x` ran on
node-local copies of the store that were not copied back, so their records hold
their only measurements. At 1.16 million nuclei, 14 QC plots exceeded ragg's
50,000-pixel limit, so the last `sizing_91x` driver run rebuilt the 17 failed
targets with plot sizes clamped in `save_plots_structured()`; that clamp is a
benchmark-only change, not on the main branch. A script without saved intermediates reruns
everything, and careful manual checkpointing narrows that gap.

### Measures

- **Critical path:** the longest chain of recorded runtimes through the
  dependency graph, as with unlimited workers. Dynamic branches run
  concurrently, and excluded or unchanged targets keep their place with no
  runtime. So do input files outside the store, such as configuration tables:
  `targets` records a file's modification time as its completion time, so a run
  that finds the file unchanged records no build of it. The conventional chain's critical path is the sum of its runtimes, as
  when a script runs its steps one after another.
- **Simulated wall time:** the same runtimes scheduled on the workers of
  `configuration_dev/crew_controllers.R`: 16 light, 6 heavy, one large and one
  huge worker, each target routed to the tier its RAM request selects. A target
  is ready once its upstream targets have finished, and each tier's workers take
  ready targets in the order they became ready. Each target also takes 2 s
  for dispatching, loading and saving (`BENCHMARK_TASK_OVERHEAD_SECONDS`), fitted
  to the single-run records, which it brings within 10% of their measured
  elapsed times, against 9-16% short without it.
  `benchmark_simulation_check.tsv` compares the simulation, with and without
  that overhead, with each recorded run's elapsed time: the periods in which its
  targets ran, from the first recorded start to the last completion, joined
  across gaps of up to ten minutes. It also gives each run's median start gap,
  the time from the last completion of a target's upstream targets in the run to
  its start, which in large runs includes waiting for a busy worker.
  `sizing_small`, `sizing_91x` and `revision_5` are left out of the check,
  because each spans several driver runs: after a failed download and a rebuild
  of a truncated annotation, after fixes for empty marker levels, nuclei without
  counts and nuclei without chromVAR scores, and after the fix for empty marker
  levels.
- **CPU time and memory:** each target's runtime times its job's mean sampled
  CPU use and unreclaimable memory, which leaves out page cache. Jobs too short
  to be sampled get the median of the sampled jobs measured with them. Memory
  counts what the processes held, including memory a worker kept from its
  earlier targets, as in any run of either workflow.
- **Peak memory of a worker:** the largest sampled memory of a worker while
  it ran a counted step, which sets the machine size, and that step. It
  includes memory the worker kept from its earlier steps, so the step that
  reached it need not be the one that needed it.
- **Reserved CPU and RAM:** each target's runtime plus its share of its
  worker's idle time, times the cores and memory that Slurm booked for that
  worker: the tier's `cpus_per_task` and `memory_gigabytes_required`, which
  crew.cluster books in MiB, so a 48 GB limit is 48 GiB or 51.5 GB. A worker's
  idle time is the lifetime of its Slurm job, from its start to its last
  sample, less the runtimes of the targets it ran: start-up, loading and saving
  data, waiting for work, and the two idle minutes before it exits. It is shared
  among those targets in proportion to their runtimes. A target belongs to the
  job whose samples name it within its runtime. A target too short to be
  sampled is assigned to a worker of its tier that was alive and free
  throughout its runtime, preferring workers that no sample names, then the
  latest started, as a new worker takes the target that was waiting for it.
  `requested_*` columns multiply the same times by the cores and RAM that the
  targets requested in the measured runs (`benchmark_allocations.tsv`), and
  `requested_runtime_*` the runtimes alone, the definition before idle time
  was counted. Memory limits leave headroom above the requests, so reserved
  RAM exceeds requested RAM. The figure shows the measured CPU time and memory,
  not reservations.
- **Disk space:** retained objects and files inside the store plus the Cell
  Ranger inputs that every workflow reads: filtered matrices, fragment files
  and their indexes.

### Accounting choices

- Unlimited workers favour the workflow with more parallel steps, multiomeR.
  The simulated wall time shows the configured budget instead. The shorter
  critical path comes from running steps in parallel, while CPU time stays
  similar.
- Crediting the conventional chain with a separate allocation per step favours
  it: a single script would hold its largest allocation, up to 500 GB for peak
  calling, for its whole run.
- Charging worker idle time weighs on multiomeR, whose many short steps leave
  workers idle between them, while the conventional chain's long steps keep its
  workers busy. Worker jobs that ran no assignable target are charged to no
  workflow, which favours multiomeR, whose short steps are the ones that cannot
  be sampled. `benchmark_worker_attribution.tsv` gives that unattributed share
  while each measurement ran.
- Leaving out the analyses that only multiomeR performs keeps its
  like-for-like measures to the shared result; the sizing measures include them.

`data/benchmark_runs/store_before_runs.tsv.gz` holds the stems and current
branches of every target upstream of the comparison, sizing and revision
aggregations' targets, as `tar_meta()` recorded them on 16 October 2026 before
the benchmark runs. The sizing and revision measures read the store's stems
first, for targets that no run rebuilt and the snapshot predates, such as the
comparison GEM wells' QC plots, then the snapshot and the records. Targets built
only on node-local store copies are absent from the store's metadata and take
the location, inside or outside the store, of their family's targets in other
aggregations or GEM wells. `benchmark_allocations.tsv` records the requests of
multiomeR 1.2.0 for the comparison runs' targets, and those at commit 2e8cbe9f
for the other targets the benchmark counts, which the sizing and revision runs
requested.

## Supplementary Table S4

Each native and reference implementation runs as its own six-core target in a
fresh R process, which records the wall time and peak process-tree RSS, on
fixtures built from `comparison_5x` targets: UCell, WNN, ATAC scDblFinder with
and without clusters, AMULET and SCAVENGE propagation. The table target stops
when a parity threshold fails.
