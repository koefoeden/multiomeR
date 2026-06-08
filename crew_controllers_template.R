# Crew controller setup template.
#
# This file enables a simple, so-called local crew-controller, that works with the get_tar_resources() helper.
# o crew_controllers.R for a local, non-SLURM setup, then
# adjust the resource table and controller list for your machine. The generous
# resource row lets the public target graph construct before site-specific
# scheduler tiers have been configured.

controller_resources_tibble <- tibble::tribble(
  ~controller_name , ~cores , ~RAM_GB , ~gpus ,
  "local"          , Inf    , Inf     , Inf
)

controller_list <- list(
  crew::crew_controller_local(
    name = "local",
    workers = 3,
    crashes_max = 1
  )
)

list(
  controller_resources_tibble = controller_resources_tibble,
  controller_list = controller_list
)
