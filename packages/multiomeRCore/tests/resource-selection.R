library(multiomeRCore)

# Match the label-sorted Esrum table that previously selected 15 CPUs/200 GB
# for a 1 CPU/32 GB request.
resources <- expand.grid(
  cores = c(1, 6, 15, 63),
  RAM_GB = c(16, 32, 60, 200)
)
resources$controller_name <- paste0(resources$cores, "core-", resources$RAM_GB, "GB")
resources$gpus <- 0
resources <- resources[order(resources$controller_name),
  c("controller_name", "cores", "RAM_GB", "gpus")]
resources <- rbind(resources, data.frame(
  controller_name = "GPU", cores = 15, RAM_GB = 200, gpus = 1
))

for (rows in list(seq_len(nrow(resources)), rev(seq_len(nrow(resources))))) {
  available <- resources[rows, ]
  for (index in seq_len(nrow(resources))) {
    requested <- resources[index, ]
    selected <- get_tar_resources(
      cores_req = requested$cores,
      RAM_GB_req = requested$RAM_GB,
      gpus_req = requested$gpus,
      controller_resources = available
    )
    stopifnot(identical(selected$crew$controller, requested$controller_name))
  }
  stopifnot(
    identical(get_tar_resources(controller_resources = available)$crew$controller, "1core-16GB"),
    identical(get_tar_resources(4, 16, controller_resources = available)$crew$controller, "6core-16GB"),
    identical(get_tar_resources(1, 8, controller_resources = available)$crew$controller, "1core-16GB"),
    identical(get_tar_resources(6, 16, RAM_GB_per_extra_core = 4,
      controller_resources = available)$crew$controller, "6core-60GB"),
    inherits(try(get_tar_resources(64, 16, controller_resources = available), silent = TRUE), "try-error")
  )

  apply_crew_controller_options(list(
    controller_resources_tibble = available,
    controller_list = lapply(available$controller_name, function(name) {
      crew::crew_controller_local(name = name, workers = 1)
    })
  ))
  stopifnot(
    identical(targets::tar_option_get("resources")$crew$controller, "1core-16GB"),
    identical(get_tar_resources(1, 32)$crew$controller, "1core-32GB")
  )
}
