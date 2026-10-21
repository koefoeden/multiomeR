# Pipeline manuscript figures

`targets.R` extends the pipeline graph with the figures and tables of the
multiomeR manuscript that come from public data. It uses the aggregations in
[`configuration_dev/`](../configuration_dev/README.md) and writes its results
below `results/`:

- `1.png`: Figure 1, selected outputs of `mixed_human_31x`, and
  `data/figure_1_aggregation.tsv`, `data/figure_1_link.tsv` and
  `data/figure_1_traits.tsv`, the values the manuscript quotes from it.
- `S5.png` and `data/cell_type_trait_enrichment.tsv`: Supplementary Figure S5,
  the cell-type deviations of every GWAS study of `mixed_human_31x`, of which
  Figure 1D shows a selection.
- `S6.png` and `data/cell_type_heritability_enrichment.tsv`: Supplementary
  Figure S6, the S-LDSC heritability enrichment of the cell types' peaks for
  every GWAS study of `mixed_human_31x` with genome-wide summary statistics.
- `data/table_S3.tsv`: Supplementary Table S3, the GWAS studies of
  `mixed_human_31x` with their Open Targets records and credible sets, and the
  credible-set loci and variants that overlap consensus peaks and the peaks they
  weight.
- `S2.png`, `S3.png` and `S4.png` and their data, Supplementary Table S6, on what multiomeR costs and
  saves: `data/seurat_signac_comparison_resources.tsv` compares multiomeR with a
  conventional Seurat/Signac workflow, normalizing with SCTransform and with
  Seurat's default log-normalization, on `comparison_1x` to `comparison_31x`;
  `data/benchmark_sizing.tsv` and `data/benchmark_sizing_stages.tsv` measure
  multiomeR with its defaults on `sizing_1x` to `sizing_91x`, whose nuclei, RNA
  non-zero entries, ATAC fragments and input files `data/benchmark_sizes.tsv`
  gives;
  `data/benchmark_revisions.tsv` measures what five revisions of `revision_31x`
  rerun; `data/benchmark_simulation_check.tsv` compares the wall time
  simulated on the configured workers with measured runs; and `data/benchmark_worker_attribution.tsv` gives the
  worker time that no target claims while each run ran. Its inputs are `data/benchmark_runs/`, the targets
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
  names = tidyselect::any_of(c("figure_1.mixed_human_31x", "supplementary_figure_S5.mixed_human_31x",
    "supplementary_figure_S6.mixed_human_31x", "supplementary_table_S3.mixed_human_31x", "supplementary_figures_S2_S4", "supplementary_table_S4.comparison_5x"))
)
```

The conventional Seurat/Signac chain needs fragtk; install it into the
environment with `pixi run install-fragtk`. The graph builds only the
dependencies of the selected results.

Supplementary Figures S2-S4 and Table S6 read recorded runtimes, not the targets they measure or any
store, so build them after the benchmark runs have finished. Before that they
stop with a missing run or runtime. The measured aggregations must be active,
because inactive ones are absent from `tar_network()`. They also need the per-job
resource history of
[slurm-monitor](https://github.com/koefoeden/slurm-monitor) in
`SLURM_MONITOR_HISTORY_DIR` (default `~/slurm_monitor_history`), which it reads
with `slurm_monitor history samples`: the schema 2 `timestamp`, `raw_job_id`,
`node`, `name`, `elapsed`, `target`, `ram_cur_gib`, `cpu_cur`, `ram_lim_gib`,
`cpus_alloc` and `scope` of each sample of a crew worker job.

The benchmark runs form four series, each on its own node, from a clean
checkout with its own store on node-local disk; a series' first run starts from
an empty store. After each run, `record_run.py` records the targets it built,
with their runtimes, completion times, bytes, errors and node, and the input
files outside the store that it tracked, in `data/benchmark_runs/<run>.tsv`.
The records, not the stores, hold the measurements. Monitor samples and worker
jobs are matched to a recorded target by its name and node, as the series ran
at the same time and built targets of the same names. The committed records
name the nodes `node_1` to `node_3`. A table outside the repository, named by
`BENCHMARK_NODE_LABELS`, gives each node's hostname (`node`) and label
(`label`), and the build maps the hostnames in the monitor's samples and in new
records to those labels.

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

Supplementary Figure S5 shows every study across the assigned cell types with
the module's unscaled deviation heatmap and GWAS metadata tracks, on Figure 1D's
colour scale and with its stars; larger deviations take the scale's end colours.

## Supplementary Figures S2-S4 and Table S6

Three figures answer three questions: what the pipeline machinery costs or
saves on the same algorithms (S2), what multiomeR's full analysis with its
defaults needs at a given data size (S3), and what revising an analysis costs
(S4). All axes are linear and start at zero, so values and differences read
directly, and the smallest aggregations bunch near the origin.
Titles and panel tags use DejaVu Sans Bold, which the environment ships inside
matplotlib; without it they fall back to a regular face.

### Like for like

The conventional chain is a compact `{targets}` implementation of a standard
Seurat v5/Signac v2 analysis on BPCells-backed assays, in the same graph, store
and controllers as multiomeR: it writes each GEM well's gene-expression counts
to a BPCells directory and runs SCTransform, PCA and RNA clustering on their
joined matrix, repeats the same MACS3, fixed-width, blacklist and consensus-peak
algorithm in one target on the cut sites that multiomeR writes for its
peak-calling clusters, quantifies the peaks per GEM well with
`Signac::FeatureMatrix()` and its default fragtk backend into BPCells
directories, and runs Signac's TF-IDF/LSI on them before building the WNN graph,
UMAP and clusters. Each step runs alone in a 32-core, 256-GB Slurm job, as on
one workstation: Seurat's kernels use its cores, and peak calling and
quantification run one `future` worker per core, in the fixed chunks of
future.apply's default scheduling, as Signac's `CallPeaks()` does. Both
clusterings leave singleton communities unmerged (`group.singletons = FALSE`),
as multiomeR does: Seurat's default merges each into a cluster in a
single-threaded loop over singletons and clusters, which took 2.3 of the
log-normalization chain's 6.9 hours at 31 samples in a first run.

A second conventional chain follows Seurat's default normalization instead of
SCTransform: `NormalizeData()` log-normalizes the joined counts on disk,
`FindVariableFeatures()` selects 2,000 variable genes, and `ScaleData()` scales
only those into an in-memory matrix before the same PCA. Its other steps and
settings are the SCTransform chain's, and it is measured on the same store, so
both chains share their peak calling and count matrices. At 91 samples, where
SCTransform cannot run, the
`comparison_lognorm_91x` series measures it on `sizing_91x`'s store after the
sizing run.

`comparison_1x` to `comparison_31x` take the first 1, 2, 5, 10, 20 and all 31
GEM wells of `mixed_human_31x` without per-well QC filters, CellBender or
doublet removal, so both workflows keep every called nucleus and yield identical
consensus peaks. multiomeR runs its Seurat SCTransform backend there and, like
Signac, counts each fragment once per peak, so both workflows use the same
algorithms. It is measured twice: up to its WNN cell metadata, and with the
multimodal Seurat/Signac object it exports from it, which S2 shows
because it matches the conventional chain's output. The shared preparation of
peak-calling clusters and their cut sites is excluded from all workflows, and
analyses that the conventional chain does not perform from both multiomeR
measurements: cell-type annotation, doublet detection, QC metrics, gene and peak
annotation, chromatin-state projection, motif analysis, and the ATAC-only and
uncorrected RNA embeddings. multiomeR and the SCTransform chain are measured in
one run, the `comparison` series, which builds all six aggregations from an
empty store.

Another line measures the same multiomeR result with its default normalization,
the BPCells_native backend, which computes Pearson residuals on disk: `sizing_1x` to
`sizing_31x` repeat the comparison aggregations with only that setting changed,
and `sizing_91x` adds 60 ENCODE multiome series, for 1.16 million nuclei. Its
per-well steps are those of the comparison aggregations. SCTransform cannot run
at that size in either workflow: Seurat converts the joined count matrix to a
`dgCMatrix`, which holds at most 2^31 - 1 non-zero entries, and `sizing_91x` has
2.17 billion. Log-normalization keeps the counts on disk and runs at 91
samples, so both tools' default normalizations reach that size.
`data/benchmark_sizes.tsv` gives each size's RNA non-zero entries, the genes
detected in each called nucleus summed from Cell Ranger's per-barcode metrics,
which equal the count in Seurat's error. In S2, colour gives the tool, blue for multiomeR and
orange for Seurat/Signac, and the darker shade each tool's default
normalization. Each line is named in its colour by tool and normalization
method: to the right of its last point for lines that end at 31 samples, and
above and left of it for the two default normalizations, which reach 91. On
linear axes, the differences read directly from the lines.

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
is left out. S3 shows the full analysis's simulated wall time, CPU time, peak
memory of a worker and disk space, with the primary module alone and, as
unconnected points at the two sizes that ran them, with all three optional
modules. Both series take multiomeR's two blues from S2, the primary module the
darker one of its default normalization, and are named in their colours above
and left of their last points. The optional modules' name gives the value of the
largest analysis. Two plot targets of `sizing_91x`, with nothing
downstream of them, are left out of its totals rather than counted as zero, as
its `excluded_targets` column states. The ATAC motif-family marker volcano plot
facets by ATAC cell type, and at about a hundred cell types its automatic
canvas, 42,300 pixels square, is more than the PNG device can allocate: its
worker crashed twice, each time about a minute after it started, against 66 s
for the whole plot at 31 wells. The genetic-enrichment variant detail plots were
not run: on earlier code they ran over 1.5 h and up to 110 GB per GWAS without
finishing, against at most 2.3 min per GWAS at 31 wells now. The plot records
they draw from finished and are counted. The sizing aggregations keep clusters
of any size, so `sizing_91x` had 1,699 cell-type groups, most of them unassigned
clusters of one or a few nuclei. Genetic enrichment now profiles only cell
types with at least 100 nuclei. S-LDSC's cell-type fits and matched-size
annotation records crashed at 48 GB on `sizing_91x` before that change and now
request 60 GB workers. No target needed a benchmark-only change.

### Revisions

`revision_31x` copies `mixed_human_31x`, the configured study with QC,
CellBender, doublet removal, tuned markers and all three optional modules. Its
full build is the `revision_base` run, from an empty store. Five runs on the
same store then change it one step further each, so each
run measures one change on top of the previous ones: the WNN clustering
resolution from 0.4 to 0.6, an added mast-cell marker set, an added dilated
cardiomyopathy GWAS, age added to the cardiomyocyte sex model, and an added
liver GEM well with its donor. The analysed configuration lacks that well, so
its recorded targets take the dependencies of the same targets of another well.
Each run's targets, critical path, CPU time and simulated wall time are
compared with the full build's. S4 shows each run's simulated wall time as a
share of the full build's, with the static targets it reran and its measured
elapsed time. A script without saved intermediates reruns everything, the
dashed line at 100%, and careful manual checkpointing narrows that gap.

### Measures

- **Critical path:** the longest chain of recorded runtimes through the
  dependency graph, as with unlimited workers. Dynamic branches run
  concurrently, and excluded or unchanged targets keep their place with no
  runtime. So do input files outside the store, such as configuration tables:
  `targets` records a file's modification time as its completion time, so a run
  that finds the file unchanged records no build of it. The conventional
  chain's critical path is the sum of its runtimes, as when a script runs its
  steps one after another. Table S6 gives it; the figures show the simulated
  wall time.
- **Simulated wall time,** which S2A, S3A and S4 show: the same runtimes
  scheduled on one workstation with 32 cores and 256 GB (`BENCHMARK_BUDGET`),
  with as many steps at once as fit. Each execution holds the cores and memory
  that its target requests, as a local worker would, such as one core and 16 GB
  for a light target. Memory, not cores, limits it: at most 16 light or 4 heavy
  steps fit at once. A target is ready once its upstream targets have finished. Ready executions
  start in the order they became ready, each as soon as it fits in the free
  cores and memory, so a smaller step can start while an earlier, larger one
  waits. Each multiomeR target also takes 2 s for dispatching, loading and
  saving (`BENCHMARK_TASK_OVERHEAD_SECONDS`). The estimate is optimistic: the
  recorded runtimes were measured with fewer steps running at once, on at most
  16 light, 6 heavy, one large and one huge worker, and some steps use more
  memory than they request. The conventional chain's simulated wall time is the
  sum of its runtimes, measured with each step alone in the workstation's job,
  as is its critical path, so it takes no per-target overhead.
  `benchmark_simulation_check.tsv` checks the recorded runtimes and that
  overhead against measured runs. It schedules each run on the workers that
  ran it, the 16 light, 6 heavy, one large and one huge worker of
  `configuration_dev/crew_controllers.R`, each execution holding one worker of
  its tier. With the overhead, the simulation came to 0.77-1.07 times the
  measured elapsed time of each single-run record, against 0.67-1.00 times
  without it, which falls short most for runs of many short targets, such as `revision_1`. The full revision build
  is the exception, at 0.53: other users' jobs on its node held up its workers,
  so 191 of its targets waited over a minute between dispatch and start, against
  none in an earlier build of the same aggregation. The check compares the
  simulation, with and without that overhead, with each recorded run's elapsed
  time: the periods in which its targets ran, from the first recorded start to
  the last completion, joined across gaps of up to ten minutes. It also gives each run's median start gap,
  the time from the last completion of a target's upstream targets in the run to
  its start, which in large runs includes waiting for a busy worker.
  Records that span several driver runs are left out: `comparison`,
  `comparison_lognorm_91x`, `sizing_31x` and `sizing_91x`.
- **CPU time and memory:** each target's runtime times its job's mean sampled
  CPU use and unreclaimable memory, which leaves out page cache. Jobs too short
  to be sampled get the median of the sampled jobs measured with them. Memory
  counts what the processes held, including memory a worker kept from its
  earlier targets, as in any run of either workflow.
- **Peak memory of a worker:** the largest sampled memory of a worker while
  it ran a counted step, which sets the machine size, and that step. It
  includes memory the worker kept from its earlier steps, so the step that
  reached it need not be the one that needed it. The Seurat/Signac steps
  therefore run garbage collection first: their workstation worker runs the
  steps of every aggregation one after another, as no script would. Without
  it, the five-sample peak calling, which ran right after the 31-sample joint
  object, started at 59 GB and peaked at 63 GB, against 28 GB with it.
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
  RAM exceeds requested RAM. The figures show the measured CPU time and memory,
  not reservations.
- **Disk space:** retained objects and files inside the store plus the Cell
  Ranger inputs that every workflow reads: filtered matrices, fragment files
  and their indexes.

### Accounting choices

- Unlimited workers favour the workflow with more parallel steps, multiomeR,
  so the figures show the simulated wall time within one server's cores and
  memory instead, and Table S6 also gives the critical path. As many workers
  as fit, with runtimes measured on fewer, still favour multiomeR. The shorter
  run time comes from running steps in parallel, while CPU time stays similar.
- Crediting the conventional chain with a separate allocation per step favours
  it: a single script would hold its largest allocation, up to 500 GB for peak
  calling, for its whole run.
- Charging worker idle time weighs on multiomeR, whose many short steps leave
  workers idle between them, while the conventional chain's long steps keep its
  workers busy. Worker jobs that ran no assignable target are charged to no
  workflow, which favours multiomeR, whose short steps are the ones that cannot
  be sampled. `benchmark_worker_attribution.tsv` gives that unattributed share
  while each run ran.
- Leaving out the analyses that only multiomeR performs keeps its
  like-for-like measures to the shared result; the sizing measures include them.

The series ran from 18 to 21 October 2026 on three nodes of one model (AMD EPYC
7543, 128 cores, 2 TB RAM), each with its workers on its own node and its own
store on that node's disk:

- `comparison` at commit f9a109d5, with both Seurat/Signac chains measured at
  594af40e, each step alone on the 32-core, 256-GB workstation worker;
- `comparison_lognorm_91x` on `sizing_91x`'s store after its sizing run,
  measured again in the same way;
- `sizing_small`, which builds `sizing_1x` to `sizing_20x`, at 950472fd, on the
  comparison's node after it had finished;
- `sizing_31x` and `sizing_91x`, each from an empty store at 113d5633 with all
  three optional modules, topped up after failed Open Targets downloads and
  after S-LDSC fixes up to 9d6c026e, from which genetic enrichment profiles
  only cell types with at least 100 nuclei;
- `revision_base` at 113d5633, then `revision_1` to `revision_5` on its store,
  at the five configuration commits on top of it.

Each target counts with its latest build. Each series' first run builds the
setup, such as reference downloads, and the per-well steps itself, and compiles
the native kernels into its own cache, so every measured aggregation is a first
build. Before 950472fd, `save_plots_structured()` drew only the first part of a
patchwork plot, which made those plots cheaper, so the sizing and revision
series were measured at 950472fd or later; the comparison builds no plot
targets. The nodes also ran other users' jobs. Their median load averages,
sampled from 12:00 to 02:00 on 19-20 October and counting the series' own
workers, were about 50 cores on node_4 and 28 on node_2. node_1
was not sampled, and a spot check showed about 29. Times may therefore be
somewhat longer than on idle nodes. The Seurat/Signac chains were measured
again on 20 and 21 October, finally with garbage collection at 594af40e. A first log-normalization series on its own node
read the largest Cell Ranger outputs from network storage 4 to 10 times as
slowly as the comparison's node, twice, with most cores idle, so that chain was
moved to the comparison's store. On `comparison_lognorm_91x`'s node, reading
the 91 Cell Ranger gene-expression matrices took 43 and then 39 minutes while
its access to network storage was slow, against 8 minutes in the final
measurement. Each record keeps its latest build. S-LDSC's cell-type fits and matched-size
annotation records moved from 16 GB to 60 GB workers during the sizing top-ups,
so `benchmark_allocations.tsv` gives each recorded target's request in its
latest build.

## Supplementary Table S4

Each native and reference implementation runs as its own six-core target in a
fresh R process, which records the wall time and peak process-tree RSS, on
fixtures built from `comparison_5x` targets: UCell, WNN, ATAC scDblFinder with
and without clusters, AMULET and SCAVENGE propagation. The table target stops
when a parity threshold fails.
