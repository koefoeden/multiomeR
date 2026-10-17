# Supplementary Figure S2 and Table S6: what multiomeR costs and what it saves, for three
# questions (README.md). Like for like: multiomeR and a conventional Seurat/Signac workflow
# on the same algorithms. Sizing: multiomeR with its defaults up to 1.16 million nuclei.
# Revisions: what changing a configured analysis reruns. Runtimes and bytes come from recorded
# executions, dependencies from the targets graph, and CPU and memory from slurm_monitor.

BENCHMARK_SCRIPT <- "pipeline_manuscript_figures/targets.R"
BENCHMARK_DATA_DIR <- "pipeline_manuscript_figures/results/data"
# One TSV per run, written by record_run.py; store_before_runs holds the store's executions
# of every target upstream of the benchmark aggregations before the benchmark runs.
BENCHMARK_RUNS_DIR <- file.path(BENCHMARK_DATA_DIR, "benchmark_runs")
BENCHMARK_SIZING_RUNS <- c("store_before_runs", "sizing_small", "sizing_31x", "sizing_91x")
BENCHMARK_RUNS <- c(BENCHMARK_SIZING_RUNS, "revision_base")
# Records that span several driver runs, after failures, have no single elapsed time.
BENCHMARK_MULTIPLE_DRIVER_RUNS <- c("sizing_small", "sizing_91x", "revision_5")
# The worker cores and RAM that each counted target requested in the run that measured it.
BENCHMARK_ALLOCATIONS <- file.path(BENCHMARK_DATA_DIR, "benchmark_allocations.tsv")
BENCHMARK_MODULES <- c("genetic_enrichment", "differential_analyses", "peak_gene_correlation")
BENCHMARK_STAGES <- c("setup", "per-well", "aggregation", BENCHMARK_MODULES)
BENCHMARK_SIZING_AGGREGATIONS <- paste0("sizing_", c(1, 2, 5, 10, 20, 31, 91), "x")
# Targets left out of the sizing measures, with the reason that benchmark_sizing.tsv records.
BENCHMARK_EXCLUDED_TARGETS <- c(
  "motif_family_accessibility_marker_volcano_plots.7_ATAC_QC.sizing_91x" = paste(
    "QC plot that does not scale to about a hundred ATAC cell types, faceted with up to 20 repelled labels",
    "per direction and facet: its worker crashed three times and a retry ran 3.9 h unfinished (62 s at 31 wells)"))
# Seconds the simulated wall time adds to each target for dispatching, loading and saving, fitted
# to the single-run records: it brings each within 10% of its measured elapsed time.
BENCHMARK_TASK_OVERHEAD_SECONDS <- 2
# Each revision run changes revision_31x one step further, in this order.
BENCHMARK_REVISIONS <- c(
  revision_1 = "WNN clustering resolution 0.4 to 0.6",
  revision_2 = "Mast-cell marker set added",
  revision_3 = "Dilated cardiomyopathy GWAS added",
  revision_4 = "Age added to the cardiomyocyte sex model",
  revision_5 = "Liver GEM well added")

# Like for like: excluded targets keep their graph position with no runtime and are left out
# of every resource sum.
SEURAT_SIGNAC_COMPARISON_WORKFLOWS <- local({
  multiomeR_exclude <- c(
    # Cluster-fragment preparation is shared by both workflows.
    "^blacklist_GRanges[.]ATAC(\\.|$)", "^BCs_per_peak_cluster_list[.]ATAC(\\.|$)",
    "^peak_calling_cluster_(names|discovery_tibble)[.]ATAC(\\.|$)",
    "^fragments_per_(cluster|peak_calling_cluster_discovery)[.]fragments[.]ATAC(\\.|$)",
    # The conventional chain performs no cell-type annotation, doublet detection,
    # QC metrics, gene or peak annotation, chromatin-state projection or motif
    # analysis, and no ATAC-only or uncorrected RNA embedding.
    "^cluster_UCell_", "^(UCell_GEX_marker_genes_list|GEX_marker_genes_vec)(\\.|$)",
    "^metadata_w_cell_types(_unfiltered|_annotation)?_tibble[.]",
    "^(amulet|scDblFinder)_", "^(ATAC_qc_metrics|per_barcode_metrics|blacklist_counts)_tibble[.]",
    "^metadata_w_QC", "Ensembl_(gene_)?annotations?_GRanges_list(\\.|$)",
    "^(consensus_peak_annotated|signac_annotation)_GRanges(\\.|$)", "^chromHMMs_",
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
      exclude = "^(?!comparison_seurat_signac_)", sequential = TRUE))
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

# The recorded executions of the named runs, one row per stem or branch.
read_benchmark_runs <- function(labels) {
  purrr::map_dfr(labels, \(label) {
    path <- Sys.glob(file.path(BENCHMARK_RUNS_DIR, paste0(label, ".tsv*")))
    if (length(path) != 1L) stop("No recorded benchmark run ", label, " in ", BENCHMARK_RUNS_DIR, call. = FALSE)
    readr::read_tsv(path, col_types = "ccdTdc", na = "") |> dplyr::mutate(run = label)
  })
}

# The latest execution of each target in the store's stems and the named runs, in that order.
# The store adds targets that no run rebuilt and the snapshot predates, such as the per-well QC
# plots of the comparison GEM wells; branches come from the records alone, which leaves out
# branches that a later build replaced.
benchmark_executions <- function(context, labels) {
  dplyr::bind_rows(context$store_stems, read_benchmark_runs(labels)) |>
    dplyr::slice_tail(n = 1, by = "name") |>
    dplyr::mutate(static_target = dplyr::coalesce(.data$parent, .data$name))
}

# The latest recorded execution of each target in the named runs.
read_benchmark_executions <- function(labels) {
  read_benchmark_runs(labels) |>
    dplyr::slice_tail(n = 1, by = "name") |>
    dplyr::mutate(static_target = dplyr::coalesce(.data$parent, .data$name))
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
    "--columns timestamp,raw_job_id,name,elapsed,target,ram_cur_gib,cpu_cur,ram_lim_gib,cpus_alloc,scope")
  samples <- data.table::fread(cmd = cmd, sep = "\t", na.strings = "NA",
    colClasses = list(character = c("timestamp", "raw_job_id", "elapsed"))) |>
    dplyr::filter(!is.na(.data$ram_cur_gib), !.data$scope %in% "sstat")
  tibble::tibble(job = samples$raw_job_id, worker = samples$name, target = samples$target,
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

# Mean measured RAM and CPU of each execution within its recorded runtime. Executions too
# short to be sampled get the median use of the sampled ones.
add_measured_resources <- function(executions, samples) {
  matched <- executions |>
    dplyr::select("name", "seconds", "completed_utc") |>
    dplyr::inner_join(samples, by = c(name = "target"), relationship = "many-to-many") |>
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
    dplyr::summarise(worker = dplyr::first(.data$worker), first_sample = dplyr::first(.data$time),
      first_elapsed = dplyr::first(.data$elapsed), ended = dplyr::last(.data$time),
      booked_cores = max(.data$booked_cores), booked_RAM_GB = max(.data$booked_RAM_GB), .by = "job") |>
    dplyr::mutate(tier = sub("^crew-worker-(.*)-[0-9a-f]+$", "\\1", .data$worker),
      started = .data$first_sample - slurm_elapsed_seconds(.data$first_elapsed)) |>
    dplyr::filter(.data$tier %in% tiers) |>
    dplyr::select("job", "tier", "started", "ended", "booked_cores", "booked_RAM_GB")
}

# Each execution's worker job: the job whose samples name it within its runtime, or else, for
# executions too short to be sampled, a job of its tier that was alive and free throughout its
# runtime. Unnamed jobs are preferred, then the latest started, as a newly started worker takes
# the target that was waiting for it. Idle workers name no target, "default" or "-".
benchmark_assign_workers <- function(executions, samples, jobs) {
  labelled <- dplyr::filter(samples, !.data$target %in% c(NA, "default", "-", "__AUTOMETRIC__"), .data$job %in% jobs$job)
  named <- executions |>
    dplyr::select("name", "seconds", "completed_utc") |>
    dplyr::inner_join(labelled[c("job", "target", "time")], by = c(name = "target"), relationship = "many-to-many") |>
    dplyr::filter(.data$time >= .data$completed_utc - .data$seconds - 20, .data$time <= .data$completed_utc + 20) |>
    dplyr::summarise(samples = dplyr::n(), .by = c("name", "completed_utc", "job")) |>
    dplyr::slice_max(.data$samples, n = 1, with_ties = FALSE, by = c("name", "completed_utc"))
  executions <- executions |>
    dplyr::left_join(named[c("name", "completed_utc", "job")], by = c("name", "completed_utc")) |>
    dplyr::mutate(start = .data$completed_utc - .data$seconds)
  busy <- split(executions[!is.na(executions$job), c("start", "completed_utc")], executions$job[!is.na(executions$job)])
  named_jobs <- unique(named$job)
  waiting <- which(is.na(executions$job) & !is.na(executions$tier))
  for (i in waiting[order(executions$start[waiting])]) {
    candidates <- jobs$job[jobs$tier == executions$tier[[i]] & jobs$started <= executions$start[[i]] &
      jobs$ended + 30 >= executions$completed_utc[[i]]]
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
  list(executions = dplyr::select(executions, "name", "completed_second", "job"),
    jobs = dplyr::bind_rows(dplyr::filter(executions, !is.na(.data$job)), unmatched) |>
      dplyr::summarise(busy_seconds = sum(.data$seconds), targets_run = dplyr::n(), .by = "job") |>
      dplyr::right_join(jobs, by = "job") |>
      dplyr::mutate(idle_seconds = as.numeric(.data$ended - .data$started, units = "secs") - .data$busy_seconds))
}

# A static target's family: its name without the aggregation or GEM well it was mapped to.
benchmark_target_family <- function(static_target) sub("[.][^.]*$", "", static_target)

upstream_targets <- function(graph, names) {
  as.character(unique(names(unlist(igraph::ego(graph, order = igraph::vcount(graph), nodes = names, mode = "in")))))
}

# Static targets of an aggregation's full analysis by stage: its own primary-module targets and
# every target of its GEM wells, with their upstream setup targets, then each optional module's
# own targets and the upstream targets that only it needs.
benchmark_stages <- function(graph, aggregation, GEM_well_IDs) {
  vertices <- igraph::V(graph)$name
  own <- vertices[endsWith(vertices, paste0(".", aggregation))]
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

# Wall time on a worker budget by list scheduling: a target is ready once its upstream targets
# have finished, and the workers of its tier take ready executions in the order they became
# ready, each for its recorded runtime plus `overhead_seconds`, which can stand for dispatching,
# loading and saving. benchmark_simulation_check() compares it with measured runs.
simulate_wall_time_seconds <- function(graph, executions, workers, overhead_seconds = 0) {
  order <- igraph::as_ids(igraph::topo_sort(graph))
  children <- stats::setNames(lapply(igraph::adjacent_vertices(graph, order, mode = "out"), igraph::as_ids), order)
  waiting <- igraph::degree(graph, order, mode = "in")
  seconds <- split(executions$seconds + overhead_seconds, factor(executions$static_target, levels = order))
  tier <- stats::setNames(executions$tier[match(order, executions$static_target)], order)
  unfinished <- lengths(seconds)
  queued_targets <- stats::setNames(rep(list(character()), length(workers)), names(workers))
  queued_seconds <- stats::setNames(rep(list(numeric()), length(workers)), names(workers))
  next_task <- stats::setNames(rep(1L, length(workers)), names(workers))
  ends <- numeric()
  running_targets <- running_tiers <- character()
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
        queued_targets[[tier[[target]]]] <- c(queued_targets[[tier[[target]]]], rep(target, unfinished[[target]]))
        queued_seconds[[tier[[target]]]] <- c(queued_seconds[[tier[[target]]]], seconds[[target]])
      }
    }
    for (worker_tier in names(workers)) {
      while (workers[[worker_tier]] > 0 && next_task[[worker_tier]] <= length(queued_targets[[worker_tier]])) {
        ends <- c(ends, now + queued_seconds[[worker_tier]][[next_task[[worker_tier]]]])
        running_targets <- c(running_targets, queued_targets[[worker_tier]][[next_task[[worker_tier]]]])
        running_tiers <- c(running_tiers, worker_tier)
        next_task[[worker_tier]] <- next_task[[worker_tier]] + 1L
        workers[[worker_tier]] <- workers[[worker_tier]] - 1
      }
    }
    if (!length(ends)) break
    finished <- which.min(ends)
    now <- ends[[finished]]
    target <- running_targets[[finished]]
    workers[[running_tiers[[finished]]]] <- workers[[running_tiers[[finished]]]] + 1
    ends <- ends[-finished]
    running_targets <- running_targets[-finished]
    running_tiers <- running_tiers[-finished]
    unfinished[[target]] <- unfinished[[target]] - 1L
    # A finished target releases its downstream targets on the next pass.
    if (unfinished[[target]] == 0L) ready <- target
  }
  stopifnot(all(unfinished == 0L), all(waiting == 0))
  now
}

# Everything the measures share: the graph; worker tiers with their configured worker counts,
# requests and Slurm bookings; each counted target's request and tier; whether its outputs are
# inside the store; the monitor's samples; and every known execution's worker job.
benchmark_context <- function(store = targets::tar_config_get("store")) {
  network <- targets::tar_network(script = BENCHMARK_SCRIPT, store = store,
    targets_only = TRUE, outdated = FALSE, reporter = "silent")
  meta <- targets::tar_meta(store = store, fields = c("type", "parent", "time", "seconds", "bytes", "format", "path", "error")) |>
    dplyr::filter(.data$type %in% c("stem", "branch"))
  controllers <- source(configuration_path("crew_controllers.R"))$value
  tiers <- tibble::tibble(tier = vapply(controllers$controller_list, \(x) x$launcher$name, character(1)),
    workers = vapply(controllers$controller_list, \(x) x$launcher$workers, numeric(1)),
    tier_booked_cores = vapply(controllers$controller_list, \(x) x$launcher$options_cluster$cpus_per_task, numeric(1)),
    tier_booked_RAM_GB = vapply(controllers$controller_list,
      \(x) x$launcher$options_cluster$memory_gigabytes_required, numeric(1)) * 2^30 / 1e9) |>
    dplyr::inner_join(controllers$controller_resources_tibble, by = c(tier = "controller_name"))
  allocations <- readr::read_tsv(BENCHMARK_ALLOCATIONS, show_col_types = FALSE) |>
    dplyr::left_join(tiers[c("tier", "cores", "RAM_GB")], by = c("cores", "RAM_GB"))
  store_prefix <- paste0(fs::path_abs(store), "/")
  runs <- tools::file_path_sans_ext(basename(fs::dir_ls(BENCHMARK_RUNS_DIR)), compression = TRUE)
  # Every known execution, so that each worker's busy time counts the targets of other runs.
  all_executions <- dplyr::bind_rows(read_benchmark_runs(runs),
    dplyr::transmute(meta, .data$name, parent = dplyr::na_if(.data$parent, ""), .data$seconds, completed_utc = .data$time)) |>
    dplyr::mutate(static_target = dplyr::coalesce(.data$parent, .data$name),
      completed_second = floor(as.numeric(.data$completed_utc))) |>
    dplyr::distinct(dplyr::pick("name", "completed_second"), .keep_all = TRUE) |>
    dplyr::select("name", "static_target", "seconds", "completed_utc", "completed_second") |>
    dplyr::left_join(allocations[c("static_target", "tier")], by = "static_target")
  timed <- all_executions$seconds >= 1
  samples <- read_slurm_monitor_samples(since = min(all_executions$completed_utc[timed] - all_executions$seconds[timed]) - 3600)
  workers <- benchmark_assign_workers(all_executions, samples, benchmark_worker_jobs(samples, tiers$tier))
  list(graph = igraph::graph_from_data_frame(network$edges, vertices = network$vertices["name"]),
    workers = stats::setNames(tiers$workers, tiers$tier), tiers = tiers, allocations = allocations,
    in_store = dplyr::transmute(meta, .data$name, family = benchmark_target_family(dplyr::coalesce(dplyr::na_if(.data$parent, ""), .data$name)),
      in_store = .data$format != "file" | purrr::map_lgl(.data$path, \(path) isTRUE(all(startsWith(fs::path_abs(path), store_prefix))))),
    samples = samples,
    store_stems = meta |> dplyr::filter(.data$type == "stem") |>
      dplyr::transmute(.data$name, parent = NA_character_, .data$seconds, completed_utc = .POSIXct(.data$time, tz = "UTC"), .data$bytes,
        error = dplyr::na_if(.data$error, ""), run = "store"),
    input_file_families = dplyr::transmute(meta, family = benchmark_target_family(dplyr::coalesce(dplyr::na_if(.data$parent, ""), .data$name)),
      input_file = .data$format == "file" & !purrr::map_lgl(.data$path, \(path) isTRUE(all(startsWith(fs::path_abs(path), store_prefix))))) |>
      dplyr::summarise(input_file = all(.data$input_file), .by = "family") |> dplyr::filter(.data$input_file) |> dplyr::pull("family"),
    execution_jobs = workers$executions, worker_jobs = workers$jobs)
}

# Measured CPU and memory; reserved time, the runtime plus a share of the worker job's idle time
# in proportion to runtime, or equal if no target of the job took measurable time; the cores and
# memory that Slurm booked for the job, or for the tier where no job is known; the request and
# tier; and whether the outputs are inside the store. Executions too short to be sampled take the
# median use of the others passed with them.
enrich_executions <- function(context, executions) {
  executions |>
    add_measured_resources(context$samples) |>
    dplyr::mutate(completed_second = floor(as.numeric(.data$completed_utc))) |>
    dplyr::left_join(context$execution_jobs, by = c("name", "completed_second")) |>
    dplyr::left_join(context$worker_jobs[c("job", "idle_seconds", "busy_seconds", "targets_run", "booked_cores", "booked_RAM_GB")],
      by = "job") |>
    dplyr::left_join(context$allocations, by = "static_target") |>
    dplyr::left_join(context$tiers[c("tier", "tier_booked_cores", "tier_booked_RAM_GB")], by = "tier") |>
    dplyr::mutate(reserved_seconds = .data$seconds + dplyr::coalesce(.data$idle_seconds *
        dplyr::if_else(.data$busy_seconds > 0, .data$seconds / .data$busy_seconds, 1 / .data$targets_run), 0),
      booked_cores = dplyr::coalesce(.data$booked_cores, .data$tier_booked_cores),
      booked_RAM_GB = dplyr::coalesce(.data$booked_RAM_GB, .data$tier_booked_RAM_GB)) |>
    # Targets built only on node-local store copies are absent from the store's metadata and take
    # the location of their family's targets in other aggregations or GEM wells.
    dplyr::left_join(context$in_store[c("name", "in_store")], by = "name") |>
    dplyr::mutate(family = benchmark_target_family(.data$static_target)) |>
    dplyr::left_join(dplyr::summarise(context$in_store, family_in_store = all(.data$in_store), .by = "family"), by = "family") |>
    dplyr::mutate(in_store = dplyr::coalesce(.data$in_store, .data$family_in_store))
}

# Sums over enriched executions. Reserved resources are what Slurm booked; requested ones what
# the targets asked for, with and without idle time. Disk space counts objects and files inside
# the store.
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
    disk_space_GB = sum(executions$bytes[executions$in_store]) / 1e9,
    targets = dplyr::n_distinct(executions$static_target), executions = nrow(executions))
}

# The share of worker time that no target claims, for the jobs of the configured tiers that ran
# while a set of executions ran. Shares are of job lifetimes and of the CPU-hours that Slurm booked.
benchmark_worker_attribution <- function(context, executions) {
  periods <- benchmark_active_periods(executions)
  jobs <- context$worker_jobs |>
    dplyr::filter(purrr::map2_lgl(as.numeric(.data$started), as.numeric(.data$ended),
      \(started, ended) any(started < periods$to & ended > periods$from))) |>
    dplyr::mutate(lifetime_hours = as.numeric(.data$ended - .data$started, units = "hours"))
  unattributed <- is.na(jobs$targets_run)
  tibble::tibble(worker_jobs = nrow(jobs), worker_hours = sum(jobs$lifetime_hours),
    unattributed_worker_hours = sum(jobs$lifetime_hours[unattributed]),
    unattributed_share = sum(jobs$lifetime_hours[unattributed]) / sum(jobs$lifetime_hours),
    unattributed_booked_CPU_share = sum((jobs$lifetime_hours * jobs$booked_cores)[unattributed]) /
      sum(jobs$lifetime_hours * jobs$booked_cores))
}

# Worker attribution while the like-for-like executions and each recorded run ran.
benchmark_worker_attribution_tibble <- function(context) {
  executions <- read_benchmark_executions("store_before_runs")
  counted <- unique(unlist(purrr::pmap(tidyr::expand_grid(workflow = names(SEURAT_SIGNAC_COMPARISON_WORKFLOWS),
    aggregation = comparison_aggregation_names()), \(workflow, aggregation) {
      workflow_scope <- benchmark_workflow_scope(context$graph, workflow, aggregation)
      workflow_scope$scope[workflow_scope$counted]
    })))
  runs <- c(setdiff(BENCHMARK_RUNS, "store_before_runs"), names(BENCHMARK_REVISIONS))
  c(list("like for like" = dplyr::filter(executions, .data$static_target %in% counted)),
    stats::setNames(lapply(runs, read_benchmark_executions), runs)) |>
    purrr::map_dfr(\(executions) benchmark_worker_attribution(context, executions), .id = "measurement")
}

# Resource measures of the counted targets among `scope`, the static targets that a result
# needs, from their enriched executions; uncounted targets keep their graph position with no
# runtime. So do input files outside the store, such as configuration tables: targets records a
# file's modification time as its completion time, so a run that finds the file unchanged records
# no build of it. A sequential workflow's critical path is the sum of its runtimes, as when a
# script runs its steps one after another.
benchmark_measures <- function(context, executions, scope, counted = rep(TRUE, length(scope)), sequential = FALSE,
                               graph = context$graph) {
  graph <- igraph::induced_subgraph(graph, scope)
  executions <- dplyr::filter(executions, .data$static_target %in% scope[counted])
  inputs <- scope[benchmark_target_family(scope) %in% context$input_file_families]
  missing <- c(setdiff(scope[counted], c(executions$static_target, inputs)), unique(executions$static_target[is.na(executions$tier)]))
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
    simulated_wall_time_hours = if (sequential) NA_real_ else
      simulate_wall_time_seconds(graph, executions, context$workers, BENCHMARK_TASK_OVERHEAD_SECONDS) / 3600,
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
# on comparison_1x to comparison_31x, measured before the benchmark runs.
seurat_signac_comparison_resource_tibble <- function(context) {
  executions <- read_benchmark_executions("store_before_runs")
  tidyr::expand_grid(aggregation = comparison_aggregation_names(), workflow = names(SEURAT_SIGNAC_COMPARISON_WORKFLOWS)) |>
    purrr::pmap_dfr(\(aggregation, workflow) tibble::tibble(aggregation = aggregation, workflow = workflow,
      benchmark_workflow_measures(context, executions, workflow, aggregation))) |>
    dplyr::left_join(benchmark_cellranger_inputs(comparison_aggregation_names()), by = "aggregation") |>
    dplyr::mutate(disk_space_GB = .data$disk_space_GB + .data$input_disk_space_GB) |>
    dplyr::relocate("GEM_well_count", "cellranger_input_nuclei", .after = "aggregation")
}

# Sizing: multiomeR with its default backend on sizing_1x to sizing_91x, whose runs built the
# aggregations' own targets on the setup and per-well targets already in the store. `sizing`
# has one row per aggregation and scope: the like-for-like subsets of the two multiomeR
# workflows, the primary module and, where the aggregation runs them, the optional modules;
# `stages` splits the full analysis by stage.
benchmark_sizing_tibbles <- function(context) {
  executions <- benchmark_executions(context, BENCHMARK_SIZING_RUNS)
  aggregation_tibble <- benchmark_aggregation_tibble(BENCHMARK_SIZING_AGGREGATIONS)
  inputs <- benchmark_cellranger_inputs(BENCHMARK_SIZING_AGGREGATIONS)
  results <- purrr::map2(aggregation_tibble$aggregation, aggregation_tibble$aggregation_GEM_well_IDs, \(aggregation, GEM_well_IDs) {
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

# Revisions: revision_31x built in full, revision_base on the setup and per-well targets already
# in the store, then changed one step further by each revision run, which records what it
# rebuilt. Shares compare each rebuild with the full build.
benchmark_revision_tibble <- function(context) {
  aggregation <- "revision_31x"
  GEM_well_IDs <- benchmark_aggregation_tibble(aggregation)$aggregation_GEM_well_IDs[[1]]
  stage_counts <- function(stages) {
    counts <- table(factor(stages, BENCHMARK_STAGES))
    tibble::as_tibble_row(stats::setNames(as.list(as.integer(counts)), paste0("targets_", gsub("-", "_", BENCHMARK_STAGES))))
  }
  stages <- benchmark_stages(context$graph, aggregation, GEM_well_IDs)
  full <- benchmark_measures(context, enrich_executions(context,
    dplyr::filter(benchmark_executions(context, BENCHMARK_RUNS), .data$static_target %in% stages$name)), stages$name)
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
# the configured workers, with the targets it did not build counted as finished, with and without
# the per-target overhead, and the run's median start gap. Records of several driver runs are
# left out.
benchmark_simulation_check <- function(context) {
  reference_GEM_well <- benchmark_aggregation_tibble("revision_31x")$aggregation_GEM_well_IDs[[1]][[1]]
  runs <- setdiff(c(BENCHMARK_RUNS, names(BENCHMARK_REVISIONS)), c("store_before_runs", BENCHMARK_MULTIPLE_DRIVER_RUNS))
  purrr::map_dfr(runs, \(label) {
    executions <- read_benchmark_executions(label)
    graph <- benchmark_run_graph(context$graph, unique(executions$static_target), reference_GEM_well)
    scope <- upstream_targets(graph, unique(executions$static_target))
    enriched <- enrich_executions(context, executions)
    measures <- benchmark_measures(context, enriched, scope, scope %in% executions$static_target, graph = graph)
    without_overhead <- simulate_wall_time_seconds(igraph::induced_subgraph(graph, scope),
      dplyr::filter(enriched, .data$static_target %in% scope), context$workers) / 3600
    start_gap <- stats::median(benchmark_start_gaps(graph, executions))
    tibble::tibble(run = label, measured_elapsed_hours = benchmark_elapsed_hours(executions),
      measures[c("simulated_wall_time_hours", "critical_path_hours", "runtime_hours", "targets", "executions")],
      simulated_without_overhead_hours = without_overhead, median_start_gap_seconds = start_gap)
  }) |>
    dplyr::mutate(simulated_to_measured = .data$simulated_wall_time_hours / .data$measured_elapsed_hours,
      simulated_without_overhead_to_measured = .data$simulated_without_overhead_hours / .data$measured_elapsed_hours)
}

# Supplementary Figure S2. Rows A-F: like for like, multiomeR with its Seurat exports on
# SCTransform and on its default backend against the Seurat/Signac chain, and the RNA non-zero
# entries against the dgCMatrix limit that stops SCTransform. G-H: sizing with defaults. I:
# revisions. Axes of sizes and resources are logarithmic, so equal gaps are equal ratios.
BENCHMARK_WORKFLOW_COLOURS <- c("multiomeR, SCTransform" = "#2a78d6", "Seurat/Signac" = "#eb6834",
  "multiomeR, default backend" = "#1baf7a")
BENCHMARK_STAGE_COLOURS <- c(setup = "#6250d6", "per-well" = "#eda100", aggregation = "#1baf7a",
  genetic_enrichment = "#e34948", differential_analyses = "#008300", peak_gene_correlation = "#e87ba4")

benchmark_nuclei_scale <- function() {
  ggplot2::scale_x_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale()),
    expand = ggplot2::expansion(mult = .06))
}

plot_benchmark_figure <- function(tables) {
  like_for_like <- dplyr::bind_rows(
    tables$seurat_signac_comparison_resources |>
      dplyr::filter(.data$workflow != "multiomeR") |>
      dplyr::mutate(workflow = dplyr::recode(.data$workflow, "multiomeR + Seurat exports" = "multiomeR, SCTransform")),
    tables$benchmark_sizing |>
      dplyr::filter(.data$scope == "like for like: multiomeR + Seurat exports") |>
      dplyr::mutate(workflow = "multiomeR, default backend")) |>
    dplyr::mutate(workflow = factor(.data$workflow, names(BENCHMARK_WORKFLOW_COLOURS)))
  measures <- c(critical_path_hours = "Critical path (h)", used_CPU_hours = "CPU time (CPU-hours)",
    used_RAM_GB_hours = "Memory (GB-hours)", peak_RAM_GB = "Peak memory of a worker (GB)", disk_space_GB = "Disk space (GB)")
  # Fold change of the Seurat/Signac chain over multiomeR on SCTransform at the largest shared size.
  largest_shared <- max(tables$seurat_signac_comparison_resources$cellranger_input_nuclei)
  like_panels <- lapply(names(measures), \(measure) {
    values <- dplyr::transmute(like_for_like, .data$cellranger_input_nuclei, .data$workflow, value = .data[[measure]])
    shared <- dplyr::filter(values, .data$cellranger_input_nuclei == largest_shared,
      .data$workflow %in% c("multiomeR, SCTransform", "Seurat/Signac"))
    ratio <- shared$value[shared$workflow == "Seurat/Signac"] / shared$value[shared$workflow == "multiomeR, SCTransform"]
    # Differences under 10% get no label.
    fold <- tibble::tibble(x = largest_shared, y = sqrt(prod(shared$value)),
      label = sprintf("%.1f\u00d7 %s", pmax(ratio, 1 / ratio), if (ratio > 1) "lower" else "higher"))[abs(log(ratio)) >= log(1.1), ]
    ggplot2::ggplot(values, ggplot2::aes(.data$cellranger_input_nuclei, .data$value, colour = .data$workflow)) +
      ggplot2::geom_line(linewidth = .5) + ggplot2::geom_point(size = 1.2) +
      ggplot2::geom_segment(data = tibble::tibble(x = largest_shared, y = min(shared$value), yend = max(shared$value))[seq_len(nrow(fold)), ],
        ggplot2::aes(.data$x, .data$y, xend = .data$x, yend = .data$yend), inherit.aes = FALSE,
        linetype = "dashed", linewidth = .3, colour = "grey45") +
      ggplot2::geom_label(data = fold, ggplot2::aes(.data$x, .data$y, label = .data$label), inherit.aes = FALSE,
        family = "NimbusSan", size = 2, colour = "grey25", linewidth = 0, label.padding = grid::unit(.6, "pt"), hjust = 1.08) +
      benchmark_nuclei_scale() +
      ggplot2::scale_y_log10(labels = scales::label_number(drop0trailing = TRUE)) +
      ggplot2::scale_colour_manual(values = BENCHMARK_WORKFLOW_COLOURS, name = NULL, drop = FALSE) +
      figure_1_theme() + ggplot2::labs(x = "Input nuclei", y = measures[[measure]])
  })
  # The joined RNA count matrix that SCTransform converts holds at most 2^31 - 1 non-zero entries.
  sizes <- dplyr::distinct(tables$benchmark_sizing, .data$cellranger_input_nuclei, .data$RNA_nonzero_entries)
  ceiling <- ggplot2::ggplot(sizes, ggplot2::aes(.data$cellranger_input_nuclei, .data$RNA_nonzero_entries)) +
    ggplot2::geom_hline(yintercept = 2^31 - 1, linetype = "dashed", linewidth = .3, colour = "#e34948") +
    ggplot2::annotate("text", x = min(sizes$cellranger_input_nuclei), y = 2^31 - 1, vjust = 1.5, hjust = 0,
      label = "dgCMatrix limit:\nSCTransform stops", family = "NimbusSan", size = 2.2, lineheight = .9, colour = "grey25") +
    ggplot2::geom_line(linewidth = .5, colour = "grey30") + ggplot2::geom_point(size = 1.2, colour = "grey30") +
    ggplot2::geom_text(data = dplyr::slice_max(sizes, .data$RNA_nonzero_entries, n = 1),
      ggplot2::aes(label = sprintf("%.2f billion", .data$RNA_nonzero_entries / 1e9)),
      hjust = 1.15, vjust = 1.6, family = "NimbusSan", size = 2.2, colour = "grey25") +
    benchmark_nuclei_scale() +
    ggplot2::scale_y_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    figure_1_theme() + ggplot2::labs(x = "Input nuclei", y = "RNA non-zero entries")

  # Sizing: wall time of the full analysis on unlimited workers and on the configured workers.
  full <- tables$benchmark_sizing |>
    dplyr::filter(!startsWith(.data$scope, "like for like")) |>
    dplyr::slice_tail(n = 1, by = "aggregation") |>
    tidyr::pivot_longer(c("critical_path_hours", "simulated_wall_time_hours"), names_to = "measure") |>
    dplyr::mutate(measure = dplyr::recode(.data$measure, critical_path_hours = "Unlimited workers",
      simulated_wall_time_hours = "Configured workers (simulated)"))
  wall_time <- ggplot2::ggplot(full, ggplot2::aes(.data$cellranger_input_nuclei, .data$value, linetype = .data$measure)) +
    ggplot2::geom_line(linewidth = .5, colour = BENCHMARK_WORKFLOW_COLOURS[["multiomeR, default backend"]]) +
    ggplot2::geom_point(ggplot2::aes(shape = .data$scope), size = 1.4, colour = BENCHMARK_WORKFLOW_COLOURS[["multiomeR, default backend"]]) +
    ggplot2::scale_shape_manual(values = c("primary module" = 16, "primary and optional modules" = 17), name = NULL) +
    ggplot2::scale_linetype_manual(values = c("Unlimited workers" = "dashed", "Configured workers (simulated)" = "solid"), name = NULL) +
    benchmark_nuclei_scale() + ggplot2::scale_y_log10(labels = scales::label_number(drop0trailing = TRUE)) +
    figure_1_theme() + ggplot2::theme(legend.position = "inside", legend.position.inside = c(.02, .98),
      legend.justification = c(0, 1), legend.box.just = "left", legend.spacing.y = grid::unit(0, "pt"),
      legend.key.height = grid::unit(8, "pt"), legend.background = ggplot2::element_rect(fill = "white", colour = NA)) +
    ggplot2::labs(x = "Input nuclei", y = "Wall time (h)")
  stages <- tables$benchmark_sizing_stages |>
    dplyr::mutate(stage = factor(.data$stage, rev(names(BENCHMARK_STAGE_COLOURS))),
      GEM_wells = factor(.data$GEM_well_count, sort(unique(.data$GEM_well_count))))
  totals <- dplyr::summarise(stages, CPU_hours = sum(.data$used_CPU_hours), .by = "GEM_wells")
  stage_shares <- ggplot2::ggplot(stages, ggplot2::aes(.data$GEM_wells, .data$used_CPU_hours, fill = .data$stage)) +
    ggplot2::geom_col(position = "fill", width = .7, colour = "white", linewidth = .3) +
    ggplot2::geom_text(data = totals, ggplot2::aes(.data$GEM_wells, 1, label = scales::label_number(accuracy = .1)(.data$CPU_hours)),
      inherit.aes = FALSE, vjust = -.4, family = "NimbusSan", size = 2) +
    ggplot2::scale_fill_manual(values = BENCHMARK_STAGE_COLOURS, name = NULL, breaks = names(BENCHMARK_STAGE_COLOURS),
      labels = \(x) gsub("_", " ", x)) +
    ggplot2::scale_y_continuous(labels = scales::label_percent(), expand = ggplot2::expansion(mult = c(0, .1))) +
    figure_1_theme() + ggplot2::theme(panel.grid.major.x = ggplot2::element_blank()) +
    ggplot2::labs(x = "GEM wells (CPU-hours above)", y = "Share of CPU time")

  # Revisions: each rebuild as shares of the full build, with the targets it reran and its measured elapsed time.
  revisions <- tables$benchmark_revisions |>
    dplyr::filter(.data$revision != "full build") |>
    dplyr::mutate(label = sprintf("%s (%s of %s targets, %.0f min)", .data$change, .data$targets, .data$targets_total,
      60 * .data$measured_elapsed_hours),
      label = factor(.data$label, rev(.data$label))) |>
    tidyr::pivot_longer(c("critical_path_share", "CPU_share", "wall_time_share"), names_to = "measure", values_to = "share") |>
    dplyr::mutate(measure = factor(dplyr::recode(.data$measure, critical_path_share = "Critical path",
      CPU_share = "CPU time", wall_time_share = "Wall time"), c("Critical path", "CPU time", "Wall time")))
  revision_plot <- ggplot2::ggplot(revisions, ggplot2::aes(.data$share, .data$label)) +
    ggplot2::geom_col(width = .6, fill = BENCHMARK_WORKFLOW_COLOURS[["multiomeR, default backend"]]) +
    ggplot2::geom_text(ggplot2::aes(label = dplyr::if_else(.data$share < .01, "<1%", scales::label_percent(accuracy = 1)(.data$share))),
      hjust = -.15, family = "NimbusSan", size = 2) +
    ggplot2::facet_wrap(ggplot2::vars(.data$measure), nrow = 1) +
    ggplot2::scale_x_continuous(labels = scales::label_percent(), limits = c(0, 1.3), breaks = c(0, .5, 1)) +
    figure_1_theme() + ggplot2::theme(panel.grid.major.y = ggplot2::element_blank()) +
    ggplot2::labs(x = "Share of the full build", y = NULL)

  like_for_like_panels <- patchwork::wrap_plots(c(like_panels, list(ceiling)), ncol = 3, guides = "collect") &
    ggplot2::theme(legend.position = "top")
  # The revision labels are long; freeing that panel keeps them from widening every row's margin.
  patchwork::wrap_plots(like_for_like_panels, patchwork::wrap_plots(wall_time, stage_shares, widths = c(1, 1.6)),
    patchwork::free(revision_plot), ncol = 1, heights = c(2, 1.1, .9)) +
    patchwork::plot_annotation(tag_levels = "A") &
    ggplot2::theme(plot.tag = ggplot2::element_text(family = "NimbusSan", face = "bold", size = 11))
}
