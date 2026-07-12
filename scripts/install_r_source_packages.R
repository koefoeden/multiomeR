is_macos <- identical(Sys.info()[["sysname"]], "Darwin")

source_packages <- if (is_macos) list(
  list(
    package = "InteractionSet",
    version = "1.38.0",
    url = "https://bioconductor.org/packages/3.22/bioc/src/contrib/InteractionSet_1.38.0.tar.gz"
  ),
  list(
    package = "GenomicInteractions",
    version = "1.44.0",
    url = "https://bioconductor.org/packages/3.22/bioc/src/contrib/GenomicInteractions_1.44.0.tar.gz"
  ),
  list(
    package = "CAGEfightR",
    version = "1.30.0",
    url = "https://bioconductor.org/packages/3.22/bioc/src/contrib/CAGEfightR_1.30.0.tar.gz"
  )
) else list()

macos_github_packages <- if (is_macos) list(
  list(
    package = "base64url",
    repo = "mllg/base64url",
    ref = "65e251ea703c9c9b4fd3f78b00a9ca46d4b044cb"
  ),
  list(
    package = "secretbase",
    repo = "shikokuchuo/secretbase",
    ref = "f010c1fb4c0107a2ba4bbe5601c4b3bcb76b9c50"
  ),
  list(
    package = "graphql",
    repo = "ropensci/graphql",
    ref = "16035c85cd400d11d99011cfd06d66024b9025ba"
  ),
  list(
    package = "ghql",
    repo = "ropensci/ghql",
    ref = "ebbf30d87ae8eacd79ba69977572c28998c519f8"
  ),
  list(
    package = "targets",
    repo = "ropensci/targets",
    ref = "e65b01c27bc38f45760e7bb52ef021f1656b9846"
  ),
  list(
    package = "tarchetypes",
    repo = "ropensci/tarchetypes",
    ref = "020909b152e81f5891242672a3c2ef42312cfc2b"
  )
) else list()

github_packages <- c(macos_github_packages, list(
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
))

is_installed_version <- function(package, version) {
  tryCatch(
    identical(as.character(utils::packageVersion(package)), version),
    error = function(error) FALSE
  )
}

install_source_package <- function(package) {
  if (is_installed_version(package$package, package$version)) {
    message(package$package, " ", package$version, " is already installed")
    return(invisible(NULL))
  }

  message("Installing ", package$package, " ", package$version)
  remotes::install_url(
    package$url,
    upgrade = "never",
    dependencies = FALSE,
    build_vignettes = FALSE
  )
  stopifnot(is_installed_version(package$package, package$version))
}

is_installed_ref <- function(package, ref) {
  description <- tryCatch(
    utils::packageDescription(package),
    error = function(error) NULL
  )

  is.list(description) && identical(description[["RemoteSha"]], ref)
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
  stopifnot(is_installed_ref(package$package, package$ref))
}

invisible(lapply(source_packages, install_source_package))
invisible(lapply(github_packages, install_github_package))
