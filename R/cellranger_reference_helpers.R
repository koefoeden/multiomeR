CELLRANGER_REFERENCE_ANNOTATION_HUB_IDS <- c(
  "gencode.v44" = "AH113665",
  "gencode.vM33" = "AH113713",
  "gencode.v32" = "AH75011",
  "gencode.vM23" = "AH75036"
)

read_cellranger_reference_json <- function(reference_json_file) {
  required_fields <- c(
    "fasta_hash",
    "genomes",
    "gtf_hash.gz",
    "input_fasta_files",
    "input_gtf_files",
    "mkref_version",
    "non_nuclear_contigs",
    "organism",
    "primary_contigs",
    "version"
  )

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

  if (
    length(cellranger_reference$genomes) != 1L ||
      length(cellranger_reference$input_fasta_files) != 1L ||
      length(cellranger_reference$input_gtf_files) != 1L
  ) {
    stop(
      "Cell Ranger reference JSON must define exactly one genome, input FASTA, and input GTF. File: ",
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

cellranger_reference_identity <- function(cellranger_reference) {
  identity_fields <- c(
    "fasta_hash",
    "genomes",
    "gtf_hash.gz",
    "input_fasta_files",
    "input_gtf_files",
    "mkref_version",
    "non_nuclear_contigs",
    "organism",
    "primary_contigs",
    "version"
  )
  cellranger_reference[identity_fields]
}

cellranger_reference_label <- function(cellranger_reference) {
  paste0(
    cellranger_reference$genomes[[1]],
    " ",
    cellranger_reference$version,
    " / ",
    cellranger_reference$gencode_version,
    " (",
    cellranger_reference$reference_json_file,
    ")"
  )
}

resolve_aggregation_cellranger_reference <- function(
  cellranger_references,
  GEM_well_IDs,
  aggregation
) {
  reference_identities <- purrr::map(
    cellranger_references,
    cellranger_reference_identity
  )
  shared_reference <- reference_identities[[1]]
  matching_reference <- purrr::map_lgl(
    reference_identities,
    identical,
    y = shared_reference
  )

  if (!all(matching_reference)) {
    reference_assignments <- paste0(
      "- ",
      GEM_well_IDs,
      ": ",
      purrr::map_chr(cellranger_references, cellranger_reference_label),
      collapse = "\n"
    )
    stop(
      "Aggregation '",
      aggregation,
      "' combines GEM wells assigned to different Cell Ranger references:\n",
      reference_assignments,
      call. = FALSE
    )
  }

  cellranger_references[[1]]
}
