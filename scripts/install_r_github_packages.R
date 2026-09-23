source("scripts/github_packages.R")

is_installed_ref <- function(package, ref) {
  identical(tryCatch(utils::packageDescription(package)[["RemoteSha"]], error = function(error) NULL), ref)
}

install_github_package <- function(package) {
  if (is_installed_ref(package$package, package$ref)) {
    message(package$package, " is already installed at ", package$ref)
    return(invisible(NULL))
  }

  message("Installing ", package$package, " from ", package$repo, "@", package$ref)
  remotes::install_github(
    package$repo,
    subdir = package$subdir,
    ref = package$ref,
    upgrade = "never",
    dependencies = FALSE,
    build_vignettes = FALSE
  )
}

invisible(lapply(github_packages, install_github_package))
