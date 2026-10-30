# R CMD tools need no project runtime. Pixi's post-link scripts run R CMD INSTALL in
# this folder to install the BSgenome data packages, before the GitHub packages exist.
if (!identical(Sys.getenv("MULTIOMER_SKIP_PROJECT_BOOTSTRAP"), "1") && !nzchar(Sys.getenv("R_CMD"))) {
  source("R/bootstrap_helpers.R")
  load_project_runtime()
}
