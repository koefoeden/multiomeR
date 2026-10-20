# Supplementary Figures S2-S4 and Table S6: what multiomeR costs and what it saves, one figure per
# question (README.md). Like for like: multiomeR and a conventional Seurat/Signac workflow
# on the same algorithms. Sizing: multiomeR with its defaults up to 1.16 million nuclei.
# Revisions: what changing a configured analysis reruns. Runtimes and bytes come from recorded
# executions, dependencies from the targets graph, and CPU and memory from slurm_monitor.

BENCHMARK_SCRIPT <- "pipeline_manuscript_figures/targets.R"
BENCHMARK_DATA_DIR <- "pipeline_manuscript_figures/results/data"
# One TSV per run, written by record_run.py.
BENCHMARK_RUNS_DIR <- file.path(BENCHMARK_DATA_DIR, "benchmark_runs")
# The records of each aggregation's full build, in build order. Each series of runs built its
# aggregations from an empty store on its own node, so a measurement reads only its own series.
BENCHMARK_AGGREGATION_RUNS <- c(
  stats::setNames(rep(list("comparison"), 6), paste0("comparison_", c(1, 2, 5, 10, 20, 31), "x")),
  stats::setNames(rep(list("sizing_small"), 5), paste0("sizing_", c(1, 2, 5, 10, 20), "x")),
  list(sizing_31x = "sizing_31x", sizing_91x = "sizing_91x", revision_31x = "revision_base"))
# Records that span several driver runs, after failures or top-ups, have no single elapsed time.
BENCHMARK_MULTIPLE_DRIVER_RUNS <- c("comparison", "comparison_lognorm", "comparison_lognorm_91x", "sizing_31x", "sizing_91x")
# The worker cores and RAM that each counted target requested in the run that measured it.
BENCHMARK_ALLOCATIONS <- file.path(BENCHMARK_DATA_DIR, "benchmark_allocations.tsv")
BENCHMARK_MODULES <- c("genetic_enrichment", "differential_analyses", "peak_gene_correlation")
BENCHMARK_STAGES <- c("setup", "per-well", "aggregation", BENCHMARK_MODULES)
BENCHMARK_SIZING_AGGREGATIONS <- paste0("sizing_", c(1, 2, 5, 10, 20, 31, 91), "x")
# Targets left out of the sizing measures, with the reason that benchmark_sizing.tsv records.
BENCHMARK_EXCLUDED_TARGETS <- c(
  "motif_family_accessibility_marker_volcano_plots.7_ATAC_QC.sizing_91x" = paste(
    "QC plot faceted by ATAC cell type, whose canvas at about a hundred cell types, 42,300 pixels square, exceeds",
    "what the PNG device can allocate: its worker crashed twice, about a minute after each start (66 s at 31 wells)"),
  "chromVAR_variant_contribution_detail_plots.cell_type_pseudobulk.genetic_enrichment.sizing_91x" = paste(
    "not run: on earlier code they took over 1.5 h and up to 110 GB per GWAS at 91 samples without finishing",
    "(at most 2.3 min per GWAS at 31 wells)"))
# Seconds the simulated wall time adds to each target for dispatching, loading and saving. With it,
# the simulation on the configured workers comes to 0.77-1.07 times the measured elapsed time of the
# single-run records, apart from the full revision build (0.53), whose workers waited for other users'
# jobs on its node.
BENCHMARK_TASK_OVERHEAD_SECONDS <- 2
# The cores and GB of the server for which the figures estimate run time, with as many workers of
# each tier as fit. That is optimistic: the recorded runtimes were measured with fewer steps at once.
BENCHMARK_BUDGET <- c(cores = 128, RAM_GB = 2000)
# Each revision run changes revision_31x one step further, in this order.
BENCHMARK_REVISIONS <- c(
  revision_1 = "WNN clustering resolution 0.4 to 0.6",
  revision_2 = "Mast-cell marker set added",
  revision_3 = "Dilated cardiomyopathy GWAS added",
  revision_4 = "Age added to the cardiomyocyte sex model",
  revision_5 = "Liver sample added")

# Like for like: excluded targets keep their graph position with no runtime and are left out
# of every resource sum.
SEURAT_SIGNAC_COMPARISON_WORKFLOWS <- local({
  multiomeR_exclude <- c(
    # Peak-calling cluster and cut-site preparation is shared by both workflows.
    "^blacklist_GRanges([.]ATAC|_by_genome)(\\.|$)", "^BCs_per_peak_cluster_list[.]ATAC(\\.|$)",
    "^peak_calling_cluster_(names|discovery_tibble)[.]ATAC(\\.|$)",
    "^peak_calling_cut_sites_per_GEM_well[.]fragments[.]ATAC(\\.|$)",
    # The conventional chain performs no cell-type annotation, doublet detection,
    # QC metrics, gene or peak annotation, chromatin-state projection or motif
    # analysis, and no ATAC-only or uncorrected RNA embedding.
    "^cluster_UCell_", "^(UCell_(GEX|WNN)_marker_genes_list|GEX_marker_genes_vec)(\\.|$)",
    "^metadata_w_cell_types(_unfiltered|_annotation)?_tibble[.]",
    "^(amulet|scDblFinder)_", "^(ATAC_qc_metrics|ATAC_blacklist_metrics|per_barcode_metrics|peak_fragment_counts)_tibble[.]",
    "^metadata_w_QC", "Ensembl_(gene_)?annotations?_(GRanges_list|by_hub_ID)(\\.|$)",
    "^(consensus_peak_annotated|signac_annotation)_GRanges(\\.|$)", "^(roadmap_)?chromHMMs_",
    "^(chromVAR|motif_family|peak_TF_motif|JASPAR)",
    "^UMAP_embeddings_tibble[.](ATAC|GEX_non_harmony)[.]", "^LSI_clusters[.]")
  list(
    # multiomeR's own result is its WNN cell metadata; the Seurat/Signac object
    # is an optional export built from it.
    multiomeR = list(endpoint = "metadata_w_clusters_tibble.WNN.", exclude = multiomeR_exclude),
    "multiomeR + Seurat exports" = list(endpoint = "multimodal_Seurat_object.8_multimodal_QC.",
      exclude = multiomeR_exclude),
    # Everything outside the conventional chain is imported from multiomeR. Its
    # peak calling uses multiomeR's clusters, so in the graph it can overlap the
    # RNA steps; a script runs the steps one after another.
    "Seurat/Signac" = list(endpoint = "comparison_seurat_signac_object.",
      exclude = "^(?!comparison_seurat_signac_)", sequential = TRUE),
    # The same chain with Seurat's default log-normalization, measured in its own run from an empty
    # store, which builds the steps it shares with the SCTransform chain again, and on sizing_91x's
    # store after its sizing run.
    "Seurat/Signac log-normalization" = list(endpoint = "comparison_seurat_signac_lognorm_object.",
      exclude = "^(?!comparison_seurat_signac_)", sequential = TRUE,
      run = c(stats::setNames(rep("comparison_lognorm", 6), paste0("comparison_", c(1, 2, 5, 10, 20, 31), "x")),
        sizing_91x = "comparison_lognorm_91x")))
})

benchmark_aggregation_tibble <- function(aggregations) {
  GEM_well_tibble <- build_GEM_well_tibble()
  aggregation_tibble_all_from_yaml <- read_aggregation_config_tibble(config_file = configuration_path("cfg_aggregations.yaml"))
  aggregation_tibble <- build_aggregation_tibble(
    aggregation_tibble_all_from_yaml = aggregation_tibble_all_from_yaml,
    GEM_well_tibble = GEM_well_tibble
  )

  missing_aggregations <- setdiff(aggregations, aggregation_tibble$aggregation)
  if (length(missing_aggregations)) {
    stop(
      "Aggregation(s) not found in cfg_aggregations.yaml: ",
      paste(missing_aggregations, collapse = ", "),
      call. = FALSE
    )
  }

  aggregation_tibble |>
    dplyr::mutate(
      benchmark_cellranger_count_dirs = purrr::map(
        .data$aggregation_GEM_well_IDs,
        \(GEM_well_IDs) GEM_well_tibble$GEM_well_cellranger_arc_count_dir[match(GEM_well_IDs, GEM_well_tibble$GEM_well_ID)]
      )
    ) |>
    dplyr::slice(match(aggregations, .data$aggregation))
}

# Size measures from Cell Ranger's per-barcode metrics of the called nuclei: their RNA non-zero
# entries are their detected genes. Every workflow also reads the filtered matrices, fragment
# files and their indexes from outside the store, which count towards its disk space.
benchmark_cellranger_inputs <- function(aggregations) {
  aggregation_tibble <- benchmark_aggregation_tibble(aggregations)
  GEM_wells <- purrr::map_dfr(unique(unlist(aggregation_tibble$benchmark_cellranger_count_dirs)), \(dir) {
    metrics <- data.table::fread(file.path(dir, "outs", "per_barcode_metrics.csv"),
      select = c("is_cell", "gex_genes_count", "atac_fragments"))
    metrics <- metrics[metrics$is_cell == 1, ]
    tibble::tibble(dir = dir, nuclei = nrow(metrics),
      RNA_nonzero_entries = sum(as.numeric(metrics$gex_genes_count)),
      ATAC_fragments = sum(as.numeric(metrics$atac_fragments)),
      input_disk_space_GB = sum(file.size(file.path(dir, "outs",
        c("filtered_feature_bc_matrix.h5", "atac_fragments.tsv.gz", "atac_fragments.tsv.gz.tbi")))) / 1e9)
  })
  aggregation_tibble |>
    dplyr::select("aggregation", dir = "benchmark_cellranger_count_dirs") |>
    tidyr::unnest("dir") |>
    dplyr::left_join(GEM_wells, by = "dir") |>
    dplyr::summarise(GEM_well_count = dplyr::n(), cellranger_input_nuclei = sum(.data$nuclei),
      RNA_nonzero_entries = sum(.data$RNA_nonzero_entries), ATAC_fragments = sum(.data$ATAC_fragments),
      input_disk_space_GB = sum(.data$input_disk_space_GB), .by = "aggregation")
}

# The committed records name each node by a neutral label. BENCHMARK_NODE_LABELS names a table kept
# outside the repository, with each node's hostname (node) and label (label), which maps the hostnames
# in the monitor's samples and in new records to those labels.
benchmark_node_labels <- function(nodes) {
  path <- Sys.getenv("BENCHMARK_NODE_LABELS")
  if (!nzchar(path)) stop("Set BENCHMARK_NODE_LABELS to the table of benchmark node labels.", call. = FALSE)
  labels <- readr::read_tsv(path, col_types = "cc")
  dplyr::coalesce(labels$label[match(nodes, labels$node)], nodes)
}

# The recorded rows of the named runs, one per stem or branch: the targets each run built and,
# with in_store false, the input files outside the store that it tracked.
read_benchmark_runs <- function(labels) {
  purrr::map_dfr(labels, \(label) {
    path <- file.path(BENCHMARK_RUNS_DIR, paste0(label, ".tsv"))
    if (!file.exists(path)) stop("No recorded benchmark run ", label, " in ", BENCHMARK_RUNS_DIR, call. = FALSE)
    readr::read_tsv(path, col_types = "ccdTdccl", na = "") |> dplyr::mutate(run = label, node = benchmark_node_labels(.data$node))
  }) |>
    dplyr::mutate(static_target = dplyr::coalesce(.data$parent, .data$name))
}

# The latest recorded build of each target in the named runs.
read_benchmark_executions <- function(labels) {
  read_benchmark_runs(labels) |>
    dplyr::filter(.data$in_store) |>
    dplyr::slice_tail(n = 1, by = "name")
}

# Measured use comes from the per-job history of slurm_monitor, which samples each
# crew worker's cgroup every 10-15 s together with its running target. RAM is the
# memory the kernel cannot reclaim, so page cache from reading and writing files is
# left out; shared pages count once and child processes are included.
read_slurm_monitor_samples <- function(since) {
  # The command reads the compacted history and the files of jobs not yet folded into
  # it. Only schema 2 holds unreclaimable RAM in ram_cur_gib; rows from sstat hold RSS
  # and are left out.
  cmd <- paste("slurm_monitor history samples --schema 2 --name '^crew-worker-' --since",
    shQuote(format(since, "%Y-%m-%d %H:%M:%S", tz = "Europe/Copenhagen")),
    "--columns timestamp,raw_job_id,node,name,elapsed,target,ram_cur_gib,cpu_cur,ram_lim_gib,cpus_alloc,scope")
  samples <- data.table::fread(cmd = cmd, sep = "\t", na.strings = "NA",
    colClasses = list(character = c("timestamp", "raw_job_id", "elapsed"))) |>
    dplyr::filter(!is.na(.data$ram_cur_gib), !.data$scope %in% "sstat")
  tibble::tibble(job = samples$raw_job_id, node = benchmark_node_labels(samples$node), worker = samples$name, target = samples$target,
    time = .POSIXct(as.POSIXct(samples$timestamp, tz = "Europe/Copenhagen"), tz = "UTC"), elapsed = samples$elapsed,
    RAM_GB = samples$ram_cur_gib * 2^30 / 1e9, cpu_cores = samples$cpu_cur,
    booked_RAM_GB = samples$ram_lim_gib * 2^30 / 1e9, booked_cores = samples$cpus_alloc)
}

# Seconds of a Slurm elapsed time, [D-][HH:]MM:SS.
slurm_elapsed_seconds <- function(elapsed) {
  days <- as.numeric(ifelse(grepl("-", elapsed), sub("-.*", "", elapsed), "0"))
  clock <- sub(".*-", "", elapsed)
  clock <- ifelse(lengths(regmatches(clock, gregexpr(":", clock))) == 1L, paste0("0:", clock), clock)
  days * 86400 + as.vector(matrix(as.numeric(unlist(strsplit(clock, ":"))), ncol = 3, byrow = TRUE) %*% c(3600, 60, 1))
}

# Mean measured RAM and CPU of each execution within its recorded runtime on its node. Executions
# too short to be sampled get the median use of the sampled ones.
add_measured_resources <- function(executions, samples) {
  matched <- executions |>
    dplyr::select("name", "node", "seconds", "completed_utc") |>
    dplyr::inner_join(samples, by = c(name = "target", "node"), relationship = "many-to-many") |>
    dplyr::filter(.data$time >= .data$completed_utc - .data$seconds - 20, .data$time <= .data$completed_utc + 20)
  measured <- dplyr::summarise(matched, used_RAM_GB = mean(.data$RAM_GB), peak_RAM_GB = max(.data$RAM_GB),
    used_cores = mean(.data$cpu_cores, na.rm = TRUE), .by = "name")
  executions |>
    dplyr::left_join(measured, by = "name") |>
    dplyr::mutate(sampled = !is.na(.data$used_RAM_GB),
      used_RAM_GB = dplyr::coalesce(.data$used_RAM_GB, stats::median(.data$used_RAM_GB, na.rm = TRUE)),
      used_cores = dplyr::coalesce(.data$used_cores, stats::median(.data$used_cores, na.rm = TRUE)))
}

# The crew worker jobs of the configured tiers: each job's tier from its name, its lifetime from
# the start that its first sample's elapsed time gives to its last sample, and the cores and
# memory that Slurm booked for it. crew.cluster books memory in MiB, 1024 per configured GB.
benchmark_worker_jobs <- function(samples, tiers) {
  samples |>
    dplyr::arrange(.data$time) |>
    dplyr::summarise(node = dplyr::first(.data$node), worker = dplyr::first(.data$worker), first_sample = dplyr::first(.data$time),
      first_elapsed = dplyr::first(.data$elapsed), ended = dplyr::last(.data$time),
      booked_cores = max(.data$booked_cores), booked_RAM_GB = max(.data$booked_RAM_GB), .by = "job") |>
    dplyr::mutate(tier = sub("^crew-worker-(.*)-[0-9a-f]+$", "\\1", .data$worker),
      started = .data$first_sample - slurm_elapsed_seconds(.data$first_elapsed)) |>
    dplyr::filter(.data$tier %in% tiers) |>
    dplyr::select("job", "node", "tier", "started", "ended", "booked_cores", "booked_RAM_GB")
}

# Each execution's worker job on its node: the job whose samples name it within its runtime, or
# else, for executions too short to be sampled, a job of its tier that was alive and free
# throughout its runtime. Unnamed jobs are preferred, then the latest started, as a newly started
# worker takes the target that was waiting for it. Idle workers name no target, "default" or "-".
benchmark_assign_workers <- function(executions, samples, jobs) {
  labelled <- dplyr::filter(samples, !.data$target %in% c(NA, "default", "-", "__AUTOMETRIC__"), .data$job %in% jobs$job)
  named <- executions |>
    dplyr::select("name", "node", "seconds", "completed_utc") |>
    dplyr::inner_join(labelled[c("job", "node", "target", "time")], by = c(name = "target", "node"), relationship = "many-to-many") |>
    dplyr::filter(.data$time >= .data$completed_utc - .data$seconds - 20, .data$time <= .data$completed_utc + 20) |>
    dplyr::summarise(samples = dplyr::n(), .by = c("name", "node", "completed_utc", "job")) |>
    dplyr::slice_max(.data$samples, n = 1, with_ties = FALSE, by = c("name", "node", "completed_utc"))
  executions <- executions |>
    dplyr::left_join(named[c("name", "node", "completed_utc", "job")], by = c("name", "node", "completed_utc")) |>
    dplyr::mutate(start = .data$completed_utc - .data$seconds)
  busy <- split(executions[!is.na(executions$job), c("start", "completed_utc")], executions$job[!is.na(executions$job)])
  named_jobs <- unique(named$job)
  waiting <- which(is.na(executions$job) & !is.na(executions$tier))
  for (i in waiting[order(executions$start[waiting])]) {
    candidates <- jobs$job[jobs$tier == executions$tier[[i]] & jobs$node == executions$node[[i]] &
      jobs$started <= executions$start[[i]] & jobs$ended + 30 >= executions$completed_utc[[i]]]
    free <- candidates[vapply(candidates, \(job) {
      intervals <- busy[[job]]
      is.null(intervals) || !any(intervals$start < executions$completed_utc[[i]] & intervals$completed_utc > executions$start[[i]])
    }, logical(1))]
    if (!length(free)) next
    job <- free[order(free %in% named_jobs, -as.numeric(jobs$started[match(free, jobs$job)]))][[1]]
    executions$job[[i]] <- job
    busy[[job]] <- rbind(busy[[job]], executions[i, c("start", "completed_utc")])
  }
  # A job's samples can name a target that no execution matches, such as a build that a later
  # one replaced; that target counts with the span of its samples.
  unmatched <- labelled |>
    dplyr::anti_join(named, by = c("job", target = "name")) |>
    dplyr::summarise(seconds = as.numeric(max(.data$time) - min(.data$time), units = "secs"), .by = c("job", "target"))
  list(executions = dplyr::select(executions, "name", "node", "completed_second", "job"),
    jobs = dplyr::bind_rows(dplyr::filter(executions, !is.na(.data$job)), unmatched) |>
      dplyr::summarise(busy_seconds = sum(.data$seconds), targets_run = dplyr::n(), .by = "job") |>
      dplyr::right_join(jobs, by = "job") |>
      dplyr::mutate(idle_seconds = as.numeric(.data$ended - .data$started, units = "secs") - .data$busy_seconds))
}

upstream_targets <- function(graph, names) {
  as.character(unique(names(unlist(igraph::ego(graph, order = igraph::vcount(graph), nodes = names, mode = "in")))))
}

# Static targets of an aggregation's full analysis by stage: its own primary-module targets and
# every target of its GEM wells, with their upstream setup targets, then each optional module's
# own targets and the upstream targets that only it needs. The conventional Seurat/Signac chain,
# which sizing_91x also maps, is not part of it.
benchmark_stages <- function(graph, aggregation, GEM_well_IDs) {
  vertices <- igraph::V(graph)$name
  own <- vertices[endsWith(vertices, paste0(".", aggregation)) & !startsWith(vertices, "comparison_seurat_signac_")]
  module <- purrr::map_chr(own, \(name)
    c(BENCHMARK_MODULES[endsWith(name, paste0(".", BENCHMARK_MODULES, ".", aggregation))], NA)[[1]])
  GEM_well_targets <- vertices[sub("^.*[.]", "", vertices) %in% GEM_well_IDs]
  stages <- tibble::tibble(name = upstream_targets(graph, c(own[is.na(module)], GEM_well_targets))) |>
    dplyr::mutate(stage = dplyr::case_when(endsWith(.data$name, paste0(".", aggregation)) ~ "aggregation",
      sub("^.*[.]", "", .data$name) %in% GEM_well_IDs ~ "per-well", .default = "setup"))
  for (module_name in BENCHMARK_MODULES) {
    stages <- dplyr::bind_rows(stages, tibble::tibble(stage = module_name,
      name = setdiff(upstream_targets(graph, own[module %in% module_name]), stages$name)))
  }
  stages
}

# Revision 5 adds a GEM well that the analysed configuration lacks: its targets take the
# dependencies of the same targets of a reference GEM well.
add_GEM_well_targets <- function(graph, GEM_well_IDs, reference_GEM_well) {
  edges <- igraph::as_data_frame(graph)
  pattern <- paste0("[.]", reference_GEM_well, "$")
  reference <- edges[grepl(pattern, edges$from) | grepl(pattern, edges$to), ]
  cloned <- purrr::map_dfr(GEM_well_IDs, \(GEM_well)
    dplyr::mutate(reference, dplyr::across(c("from", "to"), \(x) sub(pattern, paste0(".", GEM_well), x))))
  edges <- dplyr::bind_rows(edges, cloned)
  igraph::graph_from_data_frame(edges, vertices = data.frame(name = union(igraph::V(graph)$name, c(edges$from, edges$to))))
}

# Longest path through the graph, weighting each target with its slowest execution, so dynamic
# branches run concurrently; targets without runtime keep their place in the graph.
critical_path_seconds <- function(graph, weights) {
  order <- igraph::as_ids(igraph::topo_sort(graph))
  parents <- lapply(igraph::adjacent_vertices(graph, order, mode = "in"), igraph::as_ids)
  finish <- stats::setNames(numeric(length(order)), order)
  for (i in seq_along(order)) finish[[i]] <- weights[[order[[i]]]] + max(0, finish[parents[[i]]])
  max(finish)
}

# Wall time by list scheduling: a target is ready once its upstream targets have finished, and
# whenever capacity is free, the earliest ready execution that fits in it starts, so a smaller one
# can start in a gap that an earlier one does not fit. Each execution holds its row of `demand`,
# one column per entry of `capacity`, for its recorded runtime plus `overhead_seconds`, which can
# stand for dispatching, loading and saving.
simulate_wall_time_seconds <- function(graph, executions, demand, capacity, overhead_seconds = 0) {
  order <- igraph::as_ids(igraph::topo_sort(graph))
  children <- stats::setNames(lapply(igraph::adjacent_vertices(graph, order, mode = "out"), igraph::as_ids), order)
  waiting <- igraph::degree(graph, order, mode = "in")
  rows <- split(seq_len(nrow(executions)), factor(executions$static_target, levels = order))
  unfinished <- lengths(rows)
  seconds <- executions$seconds + overhead_seconds
  demand <- t(demand)
  queued <- running <- integer()
  ends <- numeric()
  now <- 0
  ready <- order[waiting == 0]
  repeat {
    while (length(ready)) {
      target <- ready[[1]]
      ready <- ready[-1]
      if (unfinished[[target]] == 0L) {
        released <- children[[target]]
        waiting[released] <- waiting[released] - 1
        ready <- c(ready, released[waiting[released] == 0])
      } else {
        queued <- c(queued, rows[[target]])
      }
    }
    repeat {
      fits <- queued[colSums(demand[, queued, drop = FALSE] <= capacity) == length(capacity)]
      if (!length(fits)) break
      capacity <- capacity - demand[, fits[[1]]]
      running <- c(running, fits[[1]])
      ends <- c(ends, now + seconds[[fits[[1]]]])
      queued <- queued[queued != fits[[1]]]
    }
    if (!length(ends)) break
    finished <- which.min(ends)
    now <- ends[[finished]]
    capacity <- capacity + demand[, running[[finished]]]
    target <- executions$static_target[[running[[finished]]]]
    ends <- ends[-finished]
    running <- running[-finished]
    unfinished[[target]] <- unfinished[[target]] - 1L
    # A finished target releases its downstream targets on the next pass.
    if (unfinished[[target]] == 0L) ready <- target
  }
  stopifnot(all(unfinished == 0L), all(waiting == 0))
  now
}

# Everything the measures share: the graph; worker tiers with their configured worker counts,
# requests and Slurm bookings; each counted target's request and tier; the recorded runs and the
# input files they tracked; the monitor's samples; and every recorded build's worker job.
benchmark_context <- function() {
  network <- targets::tar_network(script = BENCHMARK_SCRIPT, targets_only = TRUE, outdated = FALSE, reporter = "silent")
  controllers <- source(configuration_path("crew_controllers.R"))$value
  tiers <- tibble::tibble(tier = vapply(controllers$controller_list, \(x) x$launcher$name, character(1)),
    workers = vapply(controllers$controller_list, \(x) x$launcher$workers, numeric(1)),
    tier_booked_cores = vapply(controllers$controller_list, \(x) x$launcher$options_cluster$cpus_per_task, numeric(1)),
    tier_booked_RAM_GB = vapply(controllers$controller_list,
      \(x) x$launcher$options_cluster$memory_gigabytes_required, numeric(1)) * 2^30 / 1e9) |>
    dplyr::inner_join(controllers$controller_resources_tibble, by = c(tier = "controller_name"))
  allocations <- readr::read_tsv(BENCHMARK_ALLOCATIONS, show_col_types = FALSE) |>
    dplyr::left_join(tiers[c("tier", "cores", "RAM_GB")], by = c("cores", "RAM_GB"))
  runs <- tools::file_path_sans_ext(basename(fs::dir_ls(BENCHMARK_RUNS_DIR, glob = "*.tsv")))
  records <- read_benchmark_runs(runs)
  # Every recorded build, so that each worker's busy time counts the targets of every run.
  executions <- records |>
    dplyr::filter(.data$in_store) |>
    dplyr::mutate(completed_second = floor(as.numeric(.data$completed_utc))) |>
    dplyr::left_join(allocations[c("static_target", "tier")], by = "static_target")
  timed <- executions$seconds >= 1
  samples <- read_slurm_monitor_samples(since = min(executions$completed_utc[timed] - executions$seconds[timed]) - 3600)
  jobs <- benchmark_worker_jobs(samples, tiers$tier)
  # Other pipelines' workers on the same node can run targets of the same names.
  samples <- dplyr::filter(samples, .data$job %in% jobs$job)
  workers <- benchmark_assign_workers(executions, samples, jobs)
  list(graph = igraph::graph_from_data_frame(network$edges, vertices = network$vertices["name"]),
    workers = stats::setNames(tiers$workers, tiers$tier), tiers = tiers, allocations = allocations, runs = runs,
    input_files = unique(records$static_target[!records$in_store]), samples = samples,
    execution_jobs = workers$executions, worker_jobs = workers$jobs)
}

# Measured CPU and memory; reserved time, the runtime plus a share of the worker job's idle time
# in proportion to runtime, or equal if no target of the job took measurable time; the cores and
# memory that Slurm booked for the job, or for the tier where no job is known; and the request and
# tier. Executions too short to be sampled take the median use of the others passed with them.
enrich_executions <- function(context, executions) {
  executions |>
    add_measured_resources(context$samples) |>
    dplyr::mutate(completed_second = floor(as.numeric(.data$completed_utc))) |>
    dplyr::left_join(context$execution_jobs, by = c("name", "node", "completed_second")) |>
    dplyr::left_join(context$worker_jobs[c("job", "idle_seconds", "busy_seconds", "targets_run", "booked_cores", "booked_RAM_GB")],
      by = "job") |>
    dplyr::left_join(context$allocations, by = "static_target") |>
    dplyr::left_join(context$tiers[c("tier", "tier_booked_cores", "tier_booked_RAM_GB")], by = "tier") |>
    dplyr::mutate(reserved_seconds = .data$seconds + dplyr::coalesce(.data$idle_seconds *
        dplyr::if_else(.data$busy_seconds > 0, .data$seconds / .data$busy_seconds, 1 / .data$targets_run), 0),
      booked_cores = dplyr::coalesce(.data$booked_cores, .data$tier_booked_cores),
      booked_RAM_GB = dplyr::coalesce(.data$booked_RAM_GB, .data$tier_booked_RAM_GB))
}

# Sums over enriched executions. Reserved resources are what Slurm booked; requested ones what
# the targets asked for, with and without idle time. Disk space counts the recorded builds'
# objects and files, which are inside the store.
summarise_executions <- function(executions) {
  peak <- dplyr::slice_max(executions, .data$peak_RAM_GB, n = 1, with_ties = FALSE, na_rm = FALSE)
  tibble::tibble(
    runtime_hours = sum(executions$seconds) / 3600,
    reserved_CPU_hours = sum(executions$reserved_seconds * executions$booked_cores) / 3600,
    reserved_RAM_GB_hours = sum(executions$reserved_seconds * executions$booked_RAM_GB) / 3600,
    requested_CPU_hours = sum(executions$reserved_seconds * executions$cores) / 3600,
    requested_RAM_GB_hours = sum(executions$reserved_seconds * executions$RAM_GB) / 3600,
    requested_runtime_CPU_hours = sum(executions$seconds * executions$cores) / 3600,
    requested_runtime_RAM_GB_hours = sum(executions$seconds * executions$RAM_GB) / 3600,
    used_CPU_hours = sum(executions$seconds * executions$used_cores) / 3600,
    used_RAM_GB_hours = sum(executions$seconds * executions$used_RAM_GB) / 3600,
    peak_RAM_GB = peak$peak_RAM_GB, peak_RAM_target = peak$static_target,
    sampled_runtime_fraction = sum(executions$seconds[executions$sampled]) / sum(executions$seconds),
    worker_assigned_runtime_fraction = sum(executions$seconds[!is.na(executions$job)]) / sum(executions$seconds),
    disk_space_GB = sum(executions$bytes) / 1e9,
    targets = dplyr::n_distinct(executions$static_target), executions = nrow(executions))
}

# The share of worker time that no target claims, for the jobs of the configured tiers that ran
# on the executions' node while they ran. Shares are of job lifetimes and of the CPU-hours that
# Slurm booked.
benchmark_worker_attribution <- function(context, executions) {
  periods <- benchmark_active_periods(executions)
  jobs <- context$worker_jobs |>
    dplyr::filter(.data$node %in% executions$node, purrr::map2_lgl(as.numeric(.data$started), as.numeric(.data$ended),
      \(started, ended) any(started < periods$to & ended > periods$from))) |>
    dplyr::mutate(lifetime_hours = as.numeric(.data$ended - .data$started, units = "hours"))
  unattributed <- is.na(jobs$targets_run)
  tibble::tibble(worker_jobs = nrow(jobs), worker_hours = sum(jobs$lifetime_hours),
    unattributed_worker_hours = sum(jobs$lifetime_hours[unattributed]),
    unattributed_share = sum(jobs$lifetime_hours[unattributed]) / sum(jobs$lifetime_hours),
    unattributed_booked_CPU_share = sum((jobs$lifetime_hours * jobs$booked_cores)[unattributed]) /
      sum(jobs$lifetime_hours * jobs$booked_cores))
}

# Worker attribution while each recorded run ran.
benchmark_worker_attribution_tibble <- function(context) {
  stats::setNames(lapply(context$runs, read_benchmark_executions), context$runs) |>
    purrr::map_dfr(\(executions) benchmark_worker_attribution(context, executions), .id = "run")
}

# Resource measures of the counted targets among `scope`, the static targets that a result
# needs, from their enriched executions; uncounted targets keep their graph position with no
# runtime. So do input files outside the store, such as configuration tables: targets records a
# file's modification time as its completion time, so a run records no build of it. The simulated
# wall time gives each execution the cores and memory that Slurm books for a worker of its tier,
# within the budget. A sequential workflow's critical path and simulated wall time are the sum of its
# runtimes, as when a script runs its steps one after another.
benchmark_measures <- function(context, executions, scope, counted = rep(TRUE, length(scope)), sequential = FALSE,
                               graph = context$graph) {
  graph <- igraph::induced_subgraph(graph, scope)
  executions <- dplyr::filter(executions, .data$static_target %in% scope[counted])
  missing <- c(setdiff(scope[counted], c(executions$static_target, context$input_files)),
    unique(executions$static_target[is.na(executions$tier)]))
  if (length(missing)) {
    stop("No recorded runtime or allocation for target(s): ", paste(utils::head(missing, 20), collapse = ", "), call. = FALSE)
  }
  if (any(!is.na(executions$error))) {
    stop("Recorded runs hold errored target(s): ", paste(utils::head(executions$name[!is.na(executions$error)], 20),
      collapse = ", "), call. = FALSE)
  }
  weights <- stats::setNames(numeric(length(scope)), scope)
  slowest <- tapply(executions$seconds, executions$static_target, max)
  weights[names(slowest)] <- slowest
  tibble::tibble(
    critical_path_hours = if (sequential) sum(executions$seconds) / 3600 else critical_path_seconds(graph, weights) / 3600,
    simulated_wall_time_hours = if (sequential) sum(executions$seconds) / 3600 else
      simulate_wall_time_seconds(graph, executions, cbind(executions$tier_booked_cores, executions$tier_booked_RAM_GB),
        BENCHMARK_BUDGET, BENCHMARK_TASK_OVERHEAD_SECONDS) / 3600,
    summarise_executions(executions))
}

# The static targets upstream of a workflow's endpoint in an aggregation, and which count.
benchmark_workflow_scope <- function(graph, workflow, aggregation) {
  spec <- SEURAT_SIGNAC_COMPARISON_WORKFLOWS[[workflow]]
  scope <- upstream_targets(graph, paste0(spec$endpoint, aggregation))
  list(scope = scope, counted = !Reduce(`|`, lapply(spec$exclude, grepl, x = scope, perl = TRUE), FALSE),
    sequential = isTRUE(spec$sequential))
}

# Measures of a workflow's counted executions, which alone set the median use of unsampled ones.
benchmark_workflow_measures <- function(context, executions, workflow, aggregation) {
  workflow_scope <- benchmark_workflow_scope(context$graph, workflow, aggregation)
  counted <- workflow_scope$scope[workflow_scope$counted]
  benchmark_measures(context, enrich_executions(context, dplyr::filter(executions, .data$static_target %in% counted)),
    workflow_scope$scope, workflow_scope$counted, workflow_scope$sequential)
}

# Like for like: multiomeR, with and without its Seurat exports, and the Seurat/Signac chain
# on comparison_1x to comparison_31x, all built in one run, and the Seurat/Signac chain with
# log-normalization, built in its own run and also on sizing_91x, where SCTransform cannot run.
seurat_signac_comparison_resource_tibble <- function(context) {
  tidyr::expand_grid(aggregation = comparison_aggregation_names(), workflow = names(SEURAT_SIGNAC_COMPARISON_WORKFLOWS)) |>
    tibble::add_row(aggregation = "sizing_91x", workflow = "Seurat/Signac log-normalization") |>
    purrr::pmap_dfr(\(aggregation, workflow) tibble::tibble(aggregation = aggregation, workflow = workflow,
      benchmark_workflow_measures(context, read_benchmark_executions(
        SEURAT_SIGNAC_COMPARISON_WORKFLOWS[[workflow]]$run[[aggregation]] %||% BENCHMARK_AGGREGATION_RUNS[[aggregation]]),
        workflow, aggregation))) |>
    dplyr::left_join(benchmark_cellranger_inputs(c(comparison_aggregation_names(), "sizing_91x")), by = "aggregation") |>
    dplyr::mutate(disk_space_GB = .data$disk_space_GB + .data$input_disk_space_GB) |>
    dplyr::relocate("GEM_well_count", "cellranger_input_nuclei", .after = "aggregation")
}

# Sizing: multiomeR with its default backend on sizing_1x to sizing_91x, each built in full from
# an empty store by its series. `sizing` has one row per aggregation and scope: the like-for-like
# subsets of the two multiomeR workflows, the primary module and, where the aggregation runs
# them, the optional modules; `stages` splits the full analysis by stage.
benchmark_sizing_tibbles <- function(context) {
  aggregation_tibble <- benchmark_aggregation_tibble(BENCHMARK_SIZING_AGGREGATIONS)
  inputs <- benchmark_cellranger_inputs(BENCHMARK_SIZING_AGGREGATIONS)
  results <- purrr::map2(aggregation_tibble$aggregation, aggregation_tibble$aggregation_GEM_well_IDs, \(aggregation, GEM_well_IDs) {
    executions <- read_benchmark_executions(BENCHMARK_AGGREGATION_RUNS[[aggregation]])
    stages <- benchmark_stages(context$graph, aggregation, GEM_well_IDs)
    excluded <- intersect(names(BENCHMARK_EXCLUDED_TARGETS), stages$name)
    stages <- dplyr::filter(stages, !.data$name %in% excluded)
    analysis <- enrich_executions(context, dplyr::filter(executions, .data$static_target %in% stages$name)) |>
      dplyr::left_join(stages, by = c(static_target = "name"))
    primary <- stages$name[stages$stage %in% BENCHMARK_STAGES[1:3]]
    sizing <- dplyr::bind_rows(
      purrr::map_dfr(c("multiomeR", "multiomeR + Seurat exports"), \(workflow) tibble::tibble(
        scope = paste("like for like:", workflow), benchmark_workflow_measures(context, executions, workflow, aggregation))),
      tibble::tibble(scope = "primary module", benchmark_measures(context, analysis, primary),
        excluded_targets = stringr::str_c(excluded, ": ", BENCHMARK_EXCLUDED_TARGETS[excluded], collapse = "; ")),
      if (length(primary) < nrow(stages)) {
        tibble::tibble(scope = "primary and optional modules", benchmark_measures(context, analysis, stages$name),
          excluded_targets = stringr::str_c(excluded, ": ", BENCHMARK_EXCLUDED_TARGETS[excluded], collapse = "; "))
      })
    list(sizing = dplyr::mutate(sizing, aggregation = aggregation, .before = 1),
      stages = purrr::map_dfr(split(analysis, factor(analysis$stage, BENCHMARK_STAGES), drop = TRUE),
        summarise_executions, .id = "stage") |> dplyr::mutate(aggregation = aggregation, .before = 1))
  })
  list(
    sizing = purrr::map_dfr(results, "sizing") |>
      dplyr::left_join(inputs, by = "aggregation") |>
      dplyr::mutate(disk_space_GB = .data$disk_space_GB + .data$input_disk_space_GB),
    stages = purrr::map_dfr(results, "stages") |>
      dplyr::left_join(inputs[c("aggregation", "GEM_well_count", "cellranger_input_nuclei")], by = "aggregation"))
}

# A recorded run's graph: the targets of GEM wells that the analysed configuration lacks take the
# dependencies of a reference GEM well's targets. Any other target missing from the graph stops.
benchmark_run_graph <- function(graph, static_targets, reference_GEM_well) {
  graph <- add_GEM_well_targets(graph, unique(sub("^.*[.]", "", setdiff(static_targets, igraph::V(graph)$name))),
    reference_GEM_well)
  unknown <- setdiff(static_targets, igraph::V(graph)$name)
  if (length(unknown)) stop("Recorded target(s) missing from the graph: ", paste(utils::head(unknown, 20), collapse = ", "), call. = FALSE)
  graph
}

# The periods in which recorded targets ran, joined across gaps of up to ten minutes, which
# leaves out the time between the driver runs of one record.
benchmark_active_periods <- function(executions) {
  executions <- dplyr::arrange(executions, .data$completed_utc - .data$seconds)
  start <- as.numeric(executions$completed_utc - executions$seconds)
  end <- cummax(as.numeric(executions$completed_utc))
  period <- cumsum(c(TRUE, start[-1] > end[-length(end)] + 600))
  tibble::tibble(from = as.numeric(tapply(start, period, min)), to = as.numeric(tapply(end, period, max)))
}

# A run's elapsed time: its active periods, from the first recorded start to the last completion.
benchmark_elapsed_hours <- function(executions) {
  periods <- benchmark_active_periods(executions)
  sum(periods$to - periods$from) / 3600
}

# Revisions: revision_31x built in full from an empty store by revision_base, then changed one step
# further by each revision run, which records what it rebuilt. Shares compare each rebuild with
# the full build.
benchmark_revision_tibble <- function(context) {
  aggregation <- "revision_31x"
  GEM_well_IDs <- benchmark_aggregation_tibble(aggregation)$aggregation_GEM_well_IDs[[1]]
  stage_counts <- function(stages) {
    counts <- table(factor(stages, BENCHMARK_STAGES))
    tibble::as_tibble_row(stats::setNames(as.list(as.integer(counts)), paste0("targets_", gsub("-", "_", BENCHMARK_STAGES))))
  }
  stages <- benchmark_stages(context$graph, aggregation, GEM_well_IDs)
  full <- benchmark_measures(context, enrich_executions(context, dplyr::filter(
    read_benchmark_executions(BENCHMARK_AGGREGATION_RUNS[[aggregation]]), .data$static_target %in% stages$name)), stages$name)
  revisions <- purrr::map_dfr(names(BENCHMARK_REVISIONS), \(revision) {
    executions <- read_benchmark_executions(revision)
    graph <- benchmark_run_graph(context$graph, unique(executions$static_target), GEM_well_IDs[[1]])
    revision_stages <- benchmark_stages(graph, aggregation,
      c(GEM_well_IDs, setdiff(sub("^.*[.]", "", igraph::V(graph)$name), sub("^.*[.]", "", igraph::V(context$graph)$name))))
    outside <- setdiff(executions$static_target, revision_stages$name)
    if (length(outside)) stop(revision, " rebuilt targets outside ", aggregation, ": ", paste(utils::head(outside, 20), collapse = ", "), call. = FALSE)
    rerun <- revision_stages$name %in% executions$static_target
    tibble::tibble(revision = revision, change = BENCHMARK_REVISIONS[[revision]], targets_total = nrow(revision_stages),
      measured_elapsed_hours = benchmark_elapsed_hours(executions),
      benchmark_measures(context, enrich_executions(context, executions), revision_stages$name, rerun, graph = graph),
      stage_counts(revision_stages$stage[rerun]))
  })
  dplyr::bind_rows(
    tibble::tibble(revision = "full build", change = "revision_31x from scratch", targets_total = nrow(stages), full,
      stage_counts(stages$stage)),
    revisions) |>
    dplyr::mutate(critical_path_share = .data$critical_path_hours / full$critical_path_hours,
      CPU_share = .data$used_CPU_hours / full$used_CPU_hours,
      wall_time_share = .data$simulated_wall_time_hours / full$simulated_wall_time_hours)
}

# The gap between each recorded target's start and the last completion of its upstream targets
# built in the same run: the main process's dispatching, loading data and any wait for a worker.
benchmark_start_gaps <- function(graph, executions) {
  edges <- igraph::as_data_frame(graph)
  edges <- edges[edges$from %in% executions$static_target & edges$to %in% executions$static_target, ]
  ready <- dplyr::summarise(executions, done = max(.data$completed_utc), .by = "static_target") |>
    dplyr::inner_join(edges, by = c(static_target = "from")) |>
    dplyr::summarise(ready = max(.data$done), .by = "to")
  executions |>
    dplyr::inner_join(ready, by = c(static_target = "to")) |>
    dplyr::transmute(gap = as.numeric(.data$completed_utc - .data$seconds - .data$ready, units = "secs")) |>
    dplyr::pull("gap")
}

# Each recorded run's measured elapsed time against the wall time simulated for its executions on
# the configured workers that ran it, each execution holding one worker of its tier, with the
# targets it did not build counted as finished, with and without the per-target overhead, and the
# run's median start gap. Records of several driver runs are left out.
benchmark_simulation_check <- function(context) {
  reference_GEM_well <- benchmark_aggregation_tibble("revision_31x")$aggregation_GEM_well_IDs[[1]][[1]]
  purrr::map_dfr(setdiff(context$runs, BENCHMARK_MULTIPLE_DRIVER_RUNS), \(label) {
    executions <- read_benchmark_executions(label)
    graph <- benchmark_run_graph(context$graph, unique(executions$static_target), reference_GEM_well)
    scope <- upstream_targets(graph, unique(executions$static_target))
    enriched <- enrich_executions(context, executions)
    measures <- benchmark_measures(context, enriched, scope, scope %in% executions$static_target, graph = graph)
    on_workers <- \(overhead_seconds) simulate_wall_time_seconds(igraph::induced_subgraph(graph, scope), enriched,
      outer(enriched$tier, names(context$workers), "=="), context$workers, overhead_seconds) / 3600
    start_gap <- stats::median(benchmark_start_gaps(graph, executions))
    tibble::tibble(run = label, measured_elapsed_hours = benchmark_elapsed_hours(executions),
      simulated_wall_time_hours = on_workers(BENCHMARK_TASK_OVERHEAD_SECONDS),
      measures[c("critical_path_hours", "runtime_hours", "targets", "executions")],
      simulated_without_overhead_hours = on_workers(0), median_start_gap_seconds = start_gap)
  }) |>
    dplyr::mutate(simulated_to_measured = .data$simulated_wall_time_hours / .data$measured_elapsed_hours,
      simulated_without_overhead_to_measured = .data$simulated_without_overhead_hours / .data$measured_elapsed_hours)
}

# Supplementary Figures S2-S4, one per question. S2, same methods: multiomeR on SCTransform with
# its Seurat exports against the Seurat/Signac chain on SCTransform, each tool with its default
# normalization, and the RNA non-zero entries against the dgCMatrix limit that stops SCTransform. S3,
# sizing: multiomeR's full analysis with its defaults. S4, revisions. All axes are linear from zero, so
# values and differences read directly, and the smallest aggregations bunch near the origin.
# Hue gives the tool, orange for Seurat/Signac and blue for multiomeR, and the darker shade each tool's
# default normalization. Every pair of the four stays apart under normal and colour-deficient vision.
BENCHMARK_WORKFLOW_COLOURS <- c("Seurat/Signac (SCTransform)" = "#eb6834", "Seurat/Signac (log-normalization)" = "#983e19",
  "multiomeR (SCTransform)" = "#3987e5", "multiomeR (default normalization)" = "#184f95")
BENCHMARK_DEFAULT_COLOUR <- BENCHMARK_WORKFLOW_COLOURS[["multiomeR (default normalization)"]]
# S3's series, both multiomeR with its default normalization, in multiomeR's two shades.
BENCHMARK_SIZING_COLOURS <- c("Primary module" = BENCHMARK_DEFAULT_COLOUR,
  "With optional modules" = BENCHMARK_WORKFLOW_COLOURS[["multiomeR (SCTransform)"]])
# The names written beside each workflow's line, by normalization method.
BENCHMARK_LINE_LABELS <- c("Seurat/Signac (SCTransform)" = "Seurat/Signac (SCT)", "multiomeR (SCTransform)" = "multiomeR (SCT)",
  "multiomeR (default normalization)" = "multiomeR (default)", "Seurat/Signac (log-normalization)" = "Seurat/Signac (default)")
# Width and height in inches.
BENCHMARK_FIGURE_SIZES <- list(S2 = c(8, 6.4), S3 = c(8, 5.6), S4 = c(8, 3.3))

# DejaVu Sans, the font that the other figures fall back to, with its bold face, which the
# environment ships only inside matplotlib. Without it, titles and tags fall back to a regular face.
BENCHMARK_FONT <- "BenchmarkSans"
register_benchmark_font <- function() {
  if (BENCHMARK_FONT %in% systemfonts::registry_fonts()$family) return(invisible())
  dir <- Sys.glob(file.path(R.home(), "..", "python3*", "site-packages", "matplotlib", "mpl-data", "fonts", "ttf"))
  files <- file.path(dir[1], c("DejaVuSans.ttf", "DejaVuSans-Bold.ttf"))
  if (length(dir) && all(file.exists(files))) systemfonts::register_font(BENCHMARK_FONT, plain = files[[1]], bold = files[[2]])
}

benchmark_panel <- function(plot, title, y) {
  plot + figure_1_theme() + ggplot2::labs(title = title, y = y) +
    ggplot2::theme(text = ggplot2::element_text(family = BENCHMARK_FONT),
      plot.title = ggplot2::element_text(size = 8.5, face = "bold", hjust = 0, margin = ggplot2::margin(b = 3)))
}

benchmark_nuclei_scale <- function() {
  ggplot2::scale_x_continuous("Input nuclei", breaks = seq(0, 1e6, by = 250e3),
    labels = scales::label_number(scale_cut = scales::cut_short_scale()), limits = c(0, NA),
    expand = ggplot2::expansion(mult = c(.02, .05)))
}

# `top` leaves room above the data for labels.
benchmark_linear_scale <- function(..., top = .1) {
  ggplot2::scale_y_continuous(..., limits = c(0, NA), expand = ggplot2::expansion(mult = c(0, top)))
}

benchmark_text <- function(...) ggplot2::geom_text(..., family = BENCHMARK_FONT, size = 2.1, lineheight = .9, colour = "grey25")

# Names each workflow's line in its colour, on a white box that masks lines passing behind it. Lines
# that end before the largest size are named to the right of their last point, and names closer than
# 1.3 lines move apart around their mean. Lines that reach the largest size, such as multiomeR's default,
# are named above and left of their last point, in the order of those points, and above any other such
# point or name to the right that they would overlap.
benchmark_line_labels <- function(values, labels = BENCHMARK_LINE_LABELS) {
  # One text line's height in data units, 4.6% of a panel whose y axis runs 22% beyond the data, and
  # one character's width, 2.4% of an x axis that runs 7% beyond the data.
  line <- .046 * 1.22 * max(values$value)
  char <- .024 * 1.07 * max(values$cellranger_input_nuclei)
  ends <- dplyr::slice_max(values, .data$cellranger_input_nuclei, by = "workflow") |>
    dplyr::mutate(label = labels[as.character(.data$workflow)], width = nchar(.data$label) * char)
  largest <- max(ends$cellranger_input_nuclei)
  right <- dplyr::filter(ends, .data$cellranger_input_nuclei < largest) |> dplyr::arrange(.data$value)
  spread <- right$value
  for (i in seq_along(spread)[-1]) spread[[i]] <- max(spread[[i]], spread[[i - 1]] + 1.3 * line)
  right$value <- spread - mean(spread) + mean(right$value)
  left <- dplyr::filter(ends, .data$cellranger_input_nuclei == largest) |> dplyr::arrange(.data$value)
  points <- left$value
  left$value <- points + .3 * line
  for (i in seq_len(nrow(left))) {
    if (i > 1) left$value[[i]] <- max(left$value[[i]], left$value[[i - 1]] + 1.3 * line)
    crossed <- points[-i][abs(left$value[[i]] + .5 * line - points[-i]) < .8 * line]
    if (length(crossed)) left$value[[i]] <- max(crossed) + .6 * line
    repeat {
      overlap <- right$cellranger_input_nuclei + 2.5e4 + right$width > largest - left$width[[i]] &
        abs(left$value[[i]] + .5 * line - right$value) < line
      if (!any(overlap)) break
      left$value[[i]] <- max(right$value[overlap]) + .8 * line
    }
  }
  list(
    ggplot2::geom_label(data = right, ggplot2::aes(label = .data$label), hjust = 0, nudge_x = 2.5e4, family = BENCHMARK_FONT,
      size = 2.3, fill = "white", linewidth = 0, label.padding = grid::unit(1, "pt"), show.legend = FALSE),
    ggplot2::geom_label(data = left, ggplot2::aes(label = .data$label), hjust = 1, vjust = 0, family = BENCHMARK_FONT,
      size = 2.3, fill = "white", linewidth = 0, label.padding = grid::unit(1, "pt"), show.legend = FALSE))
}

# A figure's panels, tagged from A.
benchmark_figure <- function(plots, ncol) {
  patchwork::wrap_plots(plots, ncol = ncol, axis_titles = "collect") + patchwork::plot_annotation(tag_levels = "A") &
    ggplot2::theme(plot.tag = ggplot2::element_text(family = BENCHMARK_FONT, face = "bold", size = 11))
}

# Supplementary Figures S2, S3 and S4 from the Table S6 data, named by figure.
plot_benchmark_figures <- function(tables) {
  like_for_like <- dplyr::bind_rows(
    tables$seurat_signac_comparison_resources |>
      dplyr::filter(.data$workflow != "multiomeR") |>
      dplyr::mutate(workflow = dplyr::recode(.data$workflow, "multiomeR + Seurat exports" = "multiomeR (SCTransform)",
        "Seurat/Signac" = "Seurat/Signac (SCTransform)", "Seurat/Signac log-normalization" = "Seurat/Signac (log-normalization)")),
    tables$benchmark_sizing |>
      dplyr::filter(.data$scope == "like for like: multiomeR + Seurat exports") |>
      dplyr::mutate(workflow = "multiomeR (default normalization)")) |>
    dplyr::mutate(workflow = factor(.data$workflow, names(BENCHMARK_WORKFLOW_COLOURS)))
  # Title and axis title of each measure.
  measures <- list(simulated_wall_time_hours = c("Run time on 128 cores, 2 TB", "Estimated hours"),
    used_CPU_hours = c("CPU time", "CPU-hours"),
    used_RAM_GB_hours = c("Memory use summed over the run", "GB-hours"),
    peak_RAM_GB = c("Memory of the largest worker", "Peak GB"),
    disk_space_GB = c("Disk space", "GB, with Cell Ranger inputs"))
  like_panels <- purrr::imap(measures, \(words, measure) {
    values <- dplyr::transmute(like_for_like, .data$cellranger_input_nuclei, .data$workflow, value = .data[[measure]])
    plot <- ggplot2::ggplot(values, ggplot2::aes(.data$cellranger_input_nuclei, .data$value, colour = .data$workflow)) +
      ggplot2::geom_line(linewidth = .5) + ggplot2::geom_point(size = 1.2) + benchmark_line_labels(values) +
      benchmark_nuclei_scale() + benchmark_linear_scale(labels = scales::label_number(drop0trailing = TRUE), top = .22) +
      ggplot2::scale_colour_manual(values = BENCHMARK_WORKFLOW_COLOURS, guide = "none")
    benchmark_panel(plot, words[[1]], words[[2]])
  })
  # The joined RNA count matrix that SCTransform converts holds at most 2^31 - 1 non-zero entries. The
  # label of the largest size, just above that limit, sits a line higher to clear it.
  sizes <- dplyr::mutate(tables$benchmark_sizes, label = sprintf("%d samples:\n%.2f billion", .data$GEM_well_count, .data$RNA_nonzero_entries / 1e9),
    vjust = dplyr::if_else(.data$GEM_well_count == max(.data$GEM_well_count), -1.1, -.2))
  ceiling <- ggplot2::ggplot(sizes, ggplot2::aes(.data$cellranger_input_nuclei, .data$RNA_nonzero_entries)) +
    ggplot2::geom_hline(yintercept = 2^31 - 1, linetype = "dashed", linewidth = .3, colour = "grey30") +
    ggplot2::annotate("text", x = min(sizes$cellranger_input_nuclei), y = 2^31 - 1, hjust = 0, vjust = -.4, family = BENCHMARK_FONT,
      size = 2.1, lineheight = .9, colour = "grey25",
      label = sprintf("In-memory matrix limit:\n%.2f billion", (2^31 - 1) / 1e9)) +
    ggplot2::geom_line(linewidth = .5, colour = "grey55") + ggplot2::geom_point(size = 1.2, colour = "grey30") +
    benchmark_text(data = dplyr::filter(sizes, .data$GEM_well_count %in% c(31, max(.data$GEM_well_count))),
      ggplot2::aes(label = .data$label, vjust = .data$vjust), hjust = 1.12) +
    benchmark_nuclei_scale() +
    benchmark_linear_scale(labels = scales::label_number(scale_cut = scales::cut_short_scale()), top = .3)
  ceiling <- benchmark_panel(ceiling, "SCTransform's size limit", "Non-zero RNA counts")

  # Sizing: multiomeR's full analysis with its defaults. The primary module ran at every size, all three
  # optional modules only at 31 and 91 samples, so their two points stay unconnected; their name gives the
  # value of the largest analysis.
  sizing <- tables$benchmark_sizing |>
    dplyr::filter(!startsWith(.data$scope, "like for like")) |>
    dplyr::mutate(workflow = factor(dplyr::if_else(.data$scope == "primary module", "Primary module", "With optional modules"),
      names(BENCHMARK_SIZING_COLOURS)))
  sizing_measures <- list(simulated_wall_time_hours = c("Run time on 128 cores, 2 TB", "Estimated hours", "%.1f h"),
    used_CPU_hours = c("CPU time", "CPU-hours", "%.0f CPU-hours"),
    peak_RAM_GB = c("Memory of the largest worker", "Peak GB", "%.0f GB"),
    disk_space_GB = c("Disk space", "GB, with Cell Ranger inputs", "%.0f GB"))
  sizing_panels <- purrr::imap(sizing_measures, \(words, measure) {
    values <- dplyr::transmute(sizing, .data$cellranger_input_nuclei, .data$workflow, value = .data[[measure]])
    largest <- dplyr::filter(values, .data$workflow == "With optional modules") |> dplyr::slice_max(.data$cellranger_input_nuclei)
    labels <- c("Primary module" = "Primary module", "With optional modules" = sprintf(paste("With optional modules:", words[[3]]), largest$value))
    plot <- ggplot2::ggplot(values, ggplot2::aes(.data$cellranger_input_nuclei, .data$value, colour = .data$workflow)) +
      ggplot2::geom_line(data = dplyr::filter(values, .data$workflow == "Primary module"), linewidth = .5) +
      ggplot2::geom_point(ggplot2::aes(shape = .data$workflow), size = 1.6) + benchmark_line_labels(values, labels) +
      ggplot2::scale_colour_manual(values = BENCHMARK_SIZING_COLOURS, guide = "none") +
      ggplot2::scale_shape_manual(values = c(16, 17), guide = "none") +
      benchmark_nuclei_scale() + benchmark_linear_scale(labels = scales::label_number(drop0trailing = TRUE), top = .22)
    benchmark_panel(plot, words[[1]], words[[2]])
  })

  # Revisions, in the order applied: each change's estimated wall time as a share of the full build's,
  # with the steps it reran and its measured elapsed time. A script without saved intermediates reruns
  # everything.
  full <- dplyr::filter(tables$benchmark_revisions, .data$revision == "full build")
  revisions <- tables$benchmark_revisions |>
    dplyr::filter(.data$revision != "full build") |>
    dplyr::mutate(label = sprintf("%d. %s\n%s of %s steps rerun", dplyr::row_number(), .data$change, scales::comma(.data$targets),
        scales::comma(.data$targets_total)),
      label = factor(.data$label, rev(.data$label)),
      value = sprintf("%s, %.0f min", dplyr::if_else(.data$wall_time_share < .01, "<1%",
        scales::label_percent(accuracy = 1)(.data$wall_time_share)), 60 * .data$measured_elapsed_hours))
  revision_plot <- ggplot2::ggplot(revisions, ggplot2::aes(.data$wall_time_share, .data$label)) +
    ggplot2::geom_col(width = .62, fill = BENCHMARK_DEFAULT_COLOUR) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", linewidth = .3, colour = "grey30") +
    ggplot2::annotate("text", x = 1, y = nrow(revisions) + .45, hjust = -.04, vjust = 0, family = BENCHMARK_FONT, size = 2.1,
      colour = "grey25", label = sprintf("Full build, %.1f h", full$simulated_wall_time_hours)) +
    # Long bars hold their value inside, clear of the full build's line.
    benchmark_text(data = dplyr::filter(revisions, .data$wall_time_share < .5), ggplot2::aes(label = .data$value), hjust = -.08) +
    ggplot2::geom_text(data = dplyr::filter(revisions, .data$wall_time_share >= .5), ggplot2::aes(label = .data$value),
      hjust = 1.08, family = BENCHMARK_FONT, size = 2.1, colour = "white") +
    ggplot2::scale_x_continuous("Share of a full build's estimated run time", labels = scales::label_percent(),
      limits = c(0, 1.32), breaks = c(0, .25, .5, .75, 1), expand = ggplot2::expansion(mult = c(0, .02))) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = c(.6, 1)))
  revision_plot <- benchmark_panel(revision_plot, "Run time of each change after a full build, in the order applied", NULL) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(), axis.text.y = ggplot2::element_text(size = 7, hjust = 1))

  list(S2 = benchmark_figure(c(like_panels, list(ceiling)), ncol = 3), S3 = benchmark_figure(sizing_panels, ncol = 2),
    S4 = revision_plot)
}

# Saves each figure as <name>.png at its size, with the bold font where the environment has it.
save_benchmark_figures <- function(plots, output_dir) {
  register_benchmark_font()
  paths <- file.path(output_dir, paste0(names(plots), ".png"))
  purrr::walk2(names(plots), paths, \(name, path) ggplot2::ggsave(path, plots[[name]], device = ragg::agg_png,
    width = BENCHMARK_FIGURE_SIZES[[name]][[1]], height = BENCHMARK_FIGURE_SIZES[[name]][[2]], dpi = 300, create.dir = TRUE))
  paths
}
