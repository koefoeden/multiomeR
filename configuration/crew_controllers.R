# Default local Crew controller setup.

# Simple so-called local crew-controller setup, that works with the get_tar_resources() helper, and
# runs the demo within a 16 CPU, 32 GB RAM machine. Only the optional Seurat_SCT GEX backend uses
# local-large. configuration_dev/crew_controllers.R shows SLURM controllers.

controller_list <- list(
  crew::crew_controller_local(
    name = "local-light",
    workers = 4
  ),
  crew::crew_controller_local(
    name = "local-heavy",
    workers = 2,
    crashes_max = 1
  ),
  crew::crew_controller_local(
    name = "local-large",
    workers = 1,
    crashes_max = 1
  )
)


controller_resources_tibble <- tibble::tribble(
  ~controller_name , ~cores , ~RAM_GB , ~gpus ,
  "local-light"    ,      1 ,      16 ,     0 ,
  "local-heavy"    ,      6 ,      32 ,     0 ,
  "local-large"    ,      6 ,      60 ,     0
)

list(
  controller_resources_tibble = controller_resources_tibble,
  controller_list = controller_list
)
