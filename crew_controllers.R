# Default local Crew controller setup.

# Simple so-called local crew-controller setup, that works with the get_tar_resources() helper, and
# fits within a 16 CPU, 256 GB RAM machine.

controller_list <- list(
  crew::crew_controller_local(
    name = "local-light",
    workers = 4
  ),
  crew::crew_controller_local(
    name = "local-heavy",
    workers = 2,
    crashes_max = 1
  )
)


controller_resources_tibble <- tibble::tribble(
  ~controller_name , ~cores , ~RAM_GB , ~gpus ,
  "local-light"    ,      1 ,      16 ,     0 ,
  "local-heavy"    ,      6 ,      60 ,     0
)


# Example HPC-scheduler (here SLURM) setup with tiered controllers.
# The resources table controls routing; the matching controller options below
# are where CPU, memory, and GPU requests are submitted to SLURM.
#
# controller_resources_tibble <- tibble::tribble(
#   ~controller_name, ~cores, ~RAM_GB, ~gpus,
#   "light",              1,      16,     0,
#   "heavy",             15,     240,     0,
#   "gpu",               15,     120,     1
# )
#
# controller_list <- list(
#   crew.cluster::crew_controller_slurm(
#     name = "light",
#     workers = 30,
#     options_cluster = crew.cluster::crew_options_slurm(
#       cpus_per_task = 1,
#       memory_gigabytes_required = 16
#     )
#   ),
#   crew.cluster::crew_controller_slurm(
#     name = "heavy",
#     workers = 10,
#     options_cluster = crew.cluster::crew_options_slurm(
#       cpus_per_task = 15,
#       memory_gigabytes_required = 240
#     )
#   ),
#   crew.cluster::crew_controller_slurm(
#     name = "gpu",
#     workers = 2,
#     options_cluster = crew.cluster::crew_options_slurm(
#       cpus_per_task = 15,
#       memory_gigabytes_required = 120,
#       script_lines = "#SBATCH --gres=gpu:1"
#     )
#   )
# )

list(
  controller_resources_tibble = controller_resources_tibble,
  controller_list = controller_list
)
