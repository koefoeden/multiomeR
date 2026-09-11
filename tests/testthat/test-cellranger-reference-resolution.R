source_project_file("R/cellranger_reference_helpers.R")
source_project_file("R/cellranger_arc_subset_helpers.R")

testthat::test_that("reference resolution uses both hashes and requires a unique match", {
  directory <- tempfile()
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  fragment_file <- file.path(directory, "fragments.tsv.gz")
  write_header <- function(lines) {
    connection <- gzfile(fragment_file, "wt")
    writeLines(lines, connection)
    close(connection)
  }
  reference_file <- file.path(directory, "reference.json")
  jsonlite::write_json(list(fasta_hash = "fasta", "gtf_hash.gz" = "gtf"), reference_file, auto_unbox = TRUE)
  header <- c(
    "# reference_path=/irrelevant/renamed-reference",
    "# reference_fasta_hash=fasta", "# reference_gtf_hash=gtf"
  )
  write_header(c(header, "chr1\t1\t2\tbarcode\t1"))
  testthat::expect_identical(
    resolve_cellranger_reference_json(fragment_file, reference_file), reference_file
  )
  testthat::expect_error(resolve_cellranger_reference_json(fragment_file, character()), "found 0")
  duplicate <- file.path(directory, "duplicate.json")
  file.copy(reference_file, duplicate)
  testthat::expect_error(
    resolve_cellranger_reference_json(fragment_file, c(reference_file, duplicate)), "found 2"
  )
  jsonlite::write_json(list(fasta_hash = "fasta", "gtf_hash.gz" = "other"), reference_file, auto_unbox = TRUE)
  testthat::expect_error(resolve_cellranger_reference_json(fragment_file, reference_file), "found 0")
  write_header(header[-3])
  testthat::expect_error(resolve_cellranger_reference_json(fragment_file, reference_file), "reference_gtf_hash")
  write_header(c(header, header[2]))
  testthat::expect_error(resolve_cellranger_reference_json(fragment_file, reference_file), "duplicated")
})

testthat::test_that("fragment subsets preserve reference headers while filtering barcodes", {
  directory <- tempfile()
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  input <- file.path(directory, "input.tsv.gz")
  output <- file.path(directory, "output.tsv.gz")
  barcodes <- file.path(directory, "barcodes.txt")
  header <- c("# reference_fasta_hash=fasta", "# reference_gtf_hash=gtf")
  connection <- gzfile(input, "wt")
  writeLines(c(header, "chr1\t1\t2\tkeep\t1", "chr1\t3\t4\tdrop\t1"), connection)
  close(connection)
  writeLines("keep", barcodes)
  cellranger_arc_subset_write_fragments(input, output, barcodes)
  connection <- gzfile(output, "rt")
  observed <- readLines(connection)
  close(connection)
  testthat::expect_identical(observed, c(header, "chr1\t1\t2\tkeep\t1"))
  testthat::expect_true(file.exists(paste0(output, ".tbi")))
})
