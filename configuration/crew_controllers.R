# Default local Crew controller setup.

# Simple so-called local crew-controller setup, that works with the get_tar_resources() helper, and
# fits within a 16 CPU, 256 GB RAM machine. configuration_dev/crew_controllers.R shows SLURM
# controllers.

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

list(
  controller_resources_tibble = controller_resources_tibble,
  controller_list = controller_list
)
