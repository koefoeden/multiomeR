# Shared pinned sources for pipeline downloads and live validation probes.
open_targets_dataset_urls <- function() {
  datasets <- c('credible_set', 'study', 'evidence_gwas_credible_sets', 'target')
  stats::setNames(paste0('https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.03/output/', datasets, '/'), datasets)
}

collectri_source <- function() {
  list(url = 'https://rescued.omnipathdb.org/CollecTRI.csv',
       sha256 = '86c90b30f2cc75c189da1f0a8c353d1547287cd656a9fac1c678634285bcb4e0')
}

amulet_repeat_list_source <- function() {
  list(url = 'https://github.com/UcarLab/AMULET/releases/download/v1.0/RestrictionRepeatLists.zip',
       sha256 = '97b9df0f375f2ec6a10d72fdd67fa9d7ec93ff0acd1bb9dd76b47cf1e014681c')
}

parse_open_targets_listing <- function(dataset_index) {
  parquet_files <- stringr::str_match_all(paste(dataset_index, collapse = '\n'), 'href="([^"/]+\\.parquet)"')[[1]][, 2]
  if (!length(parquet_files)) stop('No Parquet files found in Open Targets listing.')
  unique(parquet_files)
}

annotation_hub_ensembl_ids <- function() {
  names(ENSEMBL_ANNOTATION_HUB_MD5)
}
