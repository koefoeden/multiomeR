# Bounded live GET requests; no resource is written into the targets cache.
fetch_validation_resource <- function(url, destination, max_bytes, range = NULL) {
  args <- c('--fail', '--location', '--silent', '--show-error',
            '--connect-timeout', '15', '--max-time', '60', '--retry', '1',
            '--max-filesize', sprintf('%.0f', max_bytes), '--output', destination,
            '--write-out', '%{http_code}')
  if (!is.null(range)) args <- c(args, '--range', range)
  response <- suppressWarnings(system2('curl', shQuote(c(args, url)), stdout = TRUE))
  status <- attr(response, 'status')
  if ((!is.null(status) && status != 0L) ||
      !identical(response, if (is.null(range)) '200' else '206')) {
    stop('Live resource check failed: ', url, ' (HTTP ', paste(response, collapse = ''), ')')
  }
  if (!file.exists(destination) || file.info(destination)$size == 0 ||
      file.info(destination)$size > max_bytes) stop('Invalid response size: ', url)
  invisible(destination)
}

run_validation_preflight <- function() {
  directory <- tempfile('multiomer-preflight-')
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  results <- list()
  check <- function(name, expr) {
    message('Preflight: ', name)
    result <- tryCatch({ force(expr); 'PASS' }, error = function(error) conditionMessage(error))
    results[[name]] <<- data.frame(check = name, result = result)
  }
  source('scripts/github_packages.R', local = TRUE)
  for (package in github_packages) {
    check(paste0('package/', package$package), {
      description <- utils::packageDescription(package$package)
      stopifnot(identical(description$RemoteSha, package$ref))
      loadNamespace(package$package)
      path <- if (is.null(package$subdir)) 'DESCRIPTION' else paste0(package$subdir, '/DESCRIPTION')
      url <- paste0('https://raw.githubusercontent.com/', package$repo, '/', package$ref, '/', path)
      file <- file.path(directory, package$package)
      fetch_validation_resource(url, file, 100000)
      stopifnot(read.dcf(file)[1, 'Package'] == package$package)
    })
  }
  urls <- open_targets_dataset_urls()
  for (dataset in names(urls)) {
    check(paste0('OpenTargets/', dataset), {
      index <- file.path(directory, 'index')
      fetch_validation_resource(urls[[dataset]], index, 2e6)
      shard <- parse_open_targets_listing(readLines(index, warn = FALSE))[[1]]
      file <- file.path(directory, 'parquet-header')
      fetch_validation_resource(paste0(urls[[dataset]], shard), file, 4, '0-3')
      stopifnot(identical(readBin(file, 'raw', n = 4), charToRaw('PAR1')))
    })
  }
  check('CollecTRI/checksum-and-schema', {
    file <- file.path(directory, 'CollecTRI.csv')
    resource <- collectri_source()
    fetch_validation_resource(resource$url, file, 20e6)
    stopifnot(identical(digest::digest(file = file, algo = 'sha256'), resource$sha256))
    stopifnot(nrow(read_CollecTRI_human_network(file)) > 0L)
  })
  check('AnnotationHub/metadata-and-Ensembl', {
    hub_url <- AnnotationHub::getAnnotationHubOption('URL')
    file <- file.path(directory, 'sqlite-header')
    fetch_validation_resource(paste0(hub_url, '/metadata/annotationhub.sqlite3'), file, 16, '0-15')
    stopifnot(identical(readBin(file, 'raw', n = 16), c(charToRaw('SQLite format 3'), as.raw(0))))
    # Only Hub metadata can be refreshed here, never the large annotation files.
    hub <- with_annotation_hub_cache_lock({
      withCallingHandlers(
        AnnotationHub::AnnotationHub(
          cache = file.path(targets::tar_config_get('store'), 'files', 'AnnotationHub'),
          localHub = FALSE, ask = FALSE
        ), warning = function(warning) stop(conditionMessage(warning))
      )
    })
    ids <- annotation_hub_ensembl_ids()
    info <- AnnotationHub::getInfoOnIds(hub, ids)
    stopifnot(setequal(info$ah_id, ids), nrow(info) == length(ids))
    for (fetch_id in info$fetch_id) {
      fetch_validation_resource(paste0(hub_url, '/fetch/', fetch_id), file, 16, '0-15')
      stopifnot(identical(readBin(file, 'raw', n = 16), c(charToRaw('SQLite format 3'), as.raw(0))))
    }
  })
  dplyr::bind_rows(results)
}
