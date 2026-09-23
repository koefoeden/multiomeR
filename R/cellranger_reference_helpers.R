CELLRANGER_REFERENCE_ANNOTATION_HUB_IDS <- c(
  "gencode.v44" = "AH113665",
  "gencode.vM33" = "AH113713",
  "gencode.v32" = "AH75011",
  "gencode.vM23" = "AH75036"
)

resolve_cellranger_reference_json <- function(fragment_file, reference_json_files) {
  connection <- gzfile(fragment_file, open = "rt")
  on.exit(close(connection))
  header <- character()
  repeat {
    line <- readLines(connection, n = 1L, warn = FALSE)
    if (length(line) == 0L || !startsWith(line, "#")) break
    header <- c(header, line)
  }
  hash_keys <- c("reference_fasta_hash", "reference_gtf_hash")
  hashes <- vapply(hash_keys, function(key) {
    prefix <- paste0("# ", key, "=")
    value <- substring(header[startsWith(header, prefix)], nchar(prefix) + 1L)
    if (length(value) != 1L || !nzchar(value)) {
      stop("Missing or duplicated ", key, " in fragment header: ", fragment_file, call. = FALSE)
    }
    value
  }, character(1))
  matches <- vapply(reference_json_files, function(path) {
    reference <- jsonlite::read_json(path)
    identical(reference$fasta_hash, unname(hashes[[1]])) &&
      identical(reference[["gtf_hash.gz"]], unname(hashes[[2]]))
  }, logical(1))
  if (sum(matches) != 1L) {
    stop(
      "Expected exactly one supplied reference JSON matching the FASTA/GTF hashes in ",
      fragment_file, "; found ", sum(matches),
      ". Supply the matching reference.json under reference_metadata/ and remove duplicate matches.",
      call. = FALSE
    )
  }
  reference_json_files[matches]
}

resolve_aggregation_gene_features <- function(gene_features, GEM_well_IDs, aggregation) {
  shared_features <- gene_features[[1]]
  matching_features <- vapply(gene_features, identical, logical(1), y = shared_features)
  if (!all(matching_features)) {
    stop(
      "Aggregation '", aggregation,
      "' has incompatible ordered gene definitions in the actual Cell Ranger inputs. ",
      "Compared with GEM well '", GEM_well_IDs[[1]], "', mismatching wells: ",
      paste(GEM_well_IDs[!matching_features], collapse = ", "),
      ". Use inputs processed with the same reference; changing the configured ",
      "reference path does not change the count matrices.",
      call. = FALSE
    )
  }
  shared_features
}

read_cellranger_reference_json <- function(reference_json_file) {
  required_fields <- c("genomes", "input_gtf_files", "organism")
  cellranger_reference <- jsonlite::read_json(reference_json_file)
  missing_fields <- setdiff(required_fields, names(cellranger_reference))
  if (length(missing_fields) > 0L) {
    stop(
      "Cell Ranger reference JSON is missing required field(s): ",
      paste(missing_fields, collapse = ", "),
      ". File: ",
      reference_json_file,
      call. = FALSE
    )
  }

  if (length(cellranger_reference$genomes) != 1L || length(cellranger_reference$input_gtf_files) != 1L) {
    stop(
      "Cell Ranger reference JSON must define exactly one genome and input GTF. File: ",
      reference_json_file,
      call. = FALSE
    )
  }

  gencode_version <- stringr::str_remove(
    cellranger_reference$input_gtf_files[[1]],
    "[.]primary_assembly.*$"
  )
  annotation_hub_id <- CELLRANGER_REFERENCE_ANNOTATION_HUB_IDS[[gencode_version]]
  if (is.null(annotation_hub_id)) {
    stop(
      "Unsupported Gencode version '",
      gencode_version,
      "' in Cell Ranger reference JSON: ",
      reference_json_file,
      call. = FALSE
    )
  }

  cellranger_reference$gencode_version <- gencode_version
  cellranger_reference$annot_hub_code <- annotation_hub_id
  cellranger_reference$reference_json_file <- reference_json_file
  cellranger_reference
}

assert_cellranger_reference_matches_features <- function(
  cellranger_reference,
  feature_genomes,
  GEM_well_ID
) {
  feature_genomes <- unique(feature_genomes[!is.na(feature_genomes) & nzchar(feature_genomes)])
  reference_genome <- cellranger_reference$genomes[[1]]

  if (!identical(feature_genomes, reference_genome)) {
    stop(
      "Cell Ranger feature genome and reference JSON genome differ for GEM well '",
      GEM_well_ID,
      "': feature HDF5 = ",
      paste(feature_genomes, collapse = ", "),
      "; reference JSON = ",
      reference_genome,
      ".",
      call. = FALSE
    )
  }

  cellranger_reference
}

# Each GEM well's reference JSON is the unique match of its FASTA/GTF hashes,
# so wells share a reference exactly when they resolve to the same JSON file.
resolve_aggregation_cellranger_reference <- function(
  cellranger_references,
  GEM_well_IDs,
  aggregation
) {
  reference_json_files <- purrr::map_chr(cellranger_references, "reference_json_file")
  if (dplyr::n_distinct(reference_json_files) > 1L) {
    stop(
      "Aggregation '",
      aggregation,
      "' combines GEM wells assigned to different Cell Ranger references:\n",
      paste0("- ", GEM_well_IDs, ": ", reference_json_files, collapse = "\n"),
      call. = FALSE
    )
  }
  cellranger_references[[1]]
}
