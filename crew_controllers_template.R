# Local controller setup template.
#
# Copy this file to crew_controllers.R for a local, non-SLURM setup, then
# adjust the resource table and controller list for your machine. The generous
# resource row lets the public target graph construct before site-specific
# scheduler tiers have been configured.

controller_resources_tibble <- tibble::tribble(
  ~controller_name , ~cores , ~RAM_GB , ~gpus ,
  "local"          ,     15 ,     200 ,     0
)

controller_list <- list(
  crew::crew_controller_local(
    name = "local",
    workers = 15,
    crashes_max = 1
  )
)

list(
  controller_resources_tibble = controller_resources_tibble,
  controller_list = controller_list
)
