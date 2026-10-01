# One-node SLURM controllers for the aggregations in this directory.

# Every worker runs on the compute node of the Slurm allocation that runs the
# pipeline and connects to it through the loopback interface, as local workers
# would. Runtimes are therefore measured on one CPU model, and each worker's
# Slurm log carries crew's resource metrics with the target it runs. The worker
# counts fit a 128-CPU, 2-TB node; adjust them and the wall time for your cluster.

node <- Sys.getenv("SLURMD_NODENAME")
if (!nzchar(node)) stop("Run configuration_dev inside a Slurm allocation.", call. = FALSE)
# Worker logs stay outside the store, which may be on node-local disk, so that a
# monitor on another node can read the target each worker runs.
log_file <- file.path(fs::dir_create(file.path(getwd(), "logs", "scheduler")), "%j.txt")

tiers <- tibble::tribble(
  ~controller_name, ~cores, ~RAM_GB, ~gpus, ~workers,
  "slurm-light",         1,      16,     0,       28,
  "slurm-heavy",         6,      60,     0,       10,
  "slurm-large",         6,     200,     0,        2,
  "slurm-huge",          6,     500,     0,        1
)

controller_list <- purrr::pmap(tiers, function(controller_name, cores, RAM_GB, gpus, workers) {
  crew.cluster::crew_controller_slurm(
    name = controller_name,
    workers = workers,
    host = "127.0.0.1",
    seconds_idle = 120,
    crashes_max = 1,
    options_metrics = crew::crew_options_metrics(path = "/dev/stdout", seconds_interval = 30),
    options_cluster = crew.cluster::crew_options_slurm(
      script_lines = paste0("#SBATCH --nodelist=", node),
      log_output = log_file,
      log_error = log_file,
      cpus_per_task = cores,
      memory_gigabytes_required = RAM_GB,
      time_minutes = 24 * 60
    )
  )
})

list(
  controller_resources_tibble = tiers[c("controller_name", "cores", "RAM_GB", "gpus")],
  controller_list = controller_list
)
