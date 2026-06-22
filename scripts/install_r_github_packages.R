github_packages <- list(
  list(
    package = "BPCells",
    repo = "bnprks/BPCells",
    subdir = "r",
    ref = "28759cdd512578b6cbe549e226e1cd52a2d2308c"
  ),
  list(
    package = "Signac",
    repo = "stuart-lab/signac",
    ref = "5d66cf6c34322309d5b8aece3fa294f5531b8eee"
  ),
  list(
    package = "betterChromVAR",
    repo = "plger/betterChromVAR",
    ref = "82ae1e4ada8a43c721aa2ae1c3f2c6dfc3d5a637"
  )
)

is_installed_ref <- function(package, ref) {
  description <- tryCatch(
    utils::packageDescription(package),
    error = function(error) NULL
  )

  if (!is.list(description)) {
    return(FALSE)
  }

  remote_sha <- description[["RemoteSha"]]
  !is.null(remote_sha) && identical(remote_sha, ref)
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
