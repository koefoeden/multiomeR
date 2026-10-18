# Demultiplexing by genotype ----
#' Get cellsnp dir
#'
#' Run cellsnp-lite on one GEM well BAM using CellRanger-called barcodes and donor VCF variants.
#'
#' @param bam_file Input BAM path for one 10x Genomics GEM well.
#' @param cellranger_barcodes_tsv Path to a TSV containing CellRanger-called barcodes, used as the cellsnp-lite sample whitelist.
#' @param GEM_well_donors_VCF_file VCF path containing donor genotypes for the GEM well.
#' @param cores Number of CPU cores requested for external tools or parallel work.
#' @param minMAF Minimum minor allele frequency threshold passed to cellsnp-lite (`GEM_well_cellsnp_min_MAF`).
#' @param minCOUNT Minimum SNP read-count threshold passed to cellsnp-lite (`GEM_well_cellsnp_min_count`).
#' @return Path to the cellsnp-lite output directory containing `cellSNP.base.vcf.gz`.
#' @keywords internal

get_cellsnp_dir <- function(
  bam_file,
  cellranger_barcodes_tsv,
  GEM_well_donors_VCF_file,
  cores,
  minMAF,
  minCOUNT
) {
  cellsnp_out_dir <- get_structured_file_path(filetype = NULL)

  # -T streams each chromosome once instead of seeking to every SNP (-R), with identical counts, but at a
  # position the VCF lists more than once, such as a split multi-allelic site, it genotypes only the first
  # record. So -T genotypes the positions listed once and -R the others, merged back into VCF order.
  # The awk program writes the two site lists and prints the contigs in VCF order.
  split_program <- '
    /^#/ { print > once; print > repeated; next }
    $1 != contig { contig = $1; print contig }
    $1 FS $2 != position { flush(); position = $1 FS $2 }
    { records[++n] = $0 }
    END { flush() }
    function flush(  i) { for (i = 1; i <= n; i++) print records[i] > (n > 1 ? repeated : once); n = 0 }'
  VCF_files <- c(once = tempfile(fileext = ".vcf"), repeated = tempfile(fileext = ".vcf"))
  repeated_out_dir <- tempfile()
  on.exit(unlink(c(VCF_files, repeated_out_dir), recursive = TRUE))
  VCF_contigs <- processx::run(
    "bash",
    c("-o", "pipefail", "-c", "gzip -cdf \"$1\" | cut -f 1-8 | awk -F '\t' -v once=\"$2\" -v repeated=\"$3\" \"$4\"",
      "bash", GEM_well_donors_VCF_file, VCF_files[["once"]], VCF_files[["repeated"]], split_program)
  )$stdout |>
    strsplit("\n") |>
    unlist()

  run_cellsnp <- function(VCF_args, out_dir, threads) {
    run_w_error_check(
      command_string = "cellsnp-lite",
      arguments_chr = c(
        c("-s", bam_file),
        c("-b", cellranger_barcodes_tsv),
        c("-O", out_dir),
        VCF_args,
        c("-p", threads),
        c("--minMAF", minMAF),
        c("--minCOUNT", minCOUNT),
        c("--UMItag", "None"),
        c("--gzip")
      )
    )
  }
  run_cellsnp(c("-T", VCF_files[["once"]], "--chrom", paste(VCF_contigs, collapse = ",")), cellsnp_out_dir, cores)
  if (!all(startsWith(readLines(VCF_files[["repeated"]]), "#"))) {
    run_cellsnp(c("-R", VCF_files[["repeated"]]), repeated_out_dir, 1)
    merge_cellsnp_dirs(cellsnp_out_dir, repeated_out_dir, VCF_contigs)
  }

  if (all(startsWith(readLines(file.path(cellsnp_out_dir, "cellSNP.base.vcf.gz")), "#"))) {
    stop("cellsnp-lite found no donor VCF SNP passing its filters in ", bam_file,
      "; check that the VCF and BAM share contig names and genome build.", call. = FALSE)
  }
  return(cellsnp_out_dir)
}


#' Merge cellsnp dirs
#'
#' Merge the SNPs of a second cellsnp-lite output for the same barcodes into a first, ordering them by VCF
#' contig and position as one cellsnp-lite -R run over both SNP lists would. The two outputs must cover
#' different positions.
#'
#' @param cellsnp_dir cellsnp-lite output directory, overwritten with the merged output.
#' @param other_cellsnp_dir cellsnp-lite output directory for other positions and the same barcodes.
#' @param VCF_contigs Contigs in input VCF order.
#' @return `cellsnp_dir`, invisibly.
#' @keywords internal

merge_cellsnp_dirs <- function(cellsnp_dir, other_cellsnp_dir, VCF_contigs) {
  base_VCF_lines <- readLines(file.path(cellsnp_dir, "cellSNP.base.vcf.gz"))
  is_header <- startsWith(base_VCF_lines, "#")
  other_records <- grep("^#", readLines(file.path(other_cellsnp_dir, "cellSNP.base.vcf.gz")), value = TRUE, invert = TRUE)
  records <- c(base_VCF_lines[!is_header], other_records)
  merged_order <- order(match(sub("\t.*", "", records), VCF_contigs), as.integer(sub("^[^\t]*\t([^\t]*).*", "\\1", records)))
  new_SNP_index <- order(merged_order)
  base_VCF_connection <- gzfile(file.path(cellsnp_dir, "cellSNP.base.vcf.gz"), "w")
  writeLines(c(base_VCF_lines[is_header], records[merged_order]), base_VCF_connection)
  close(base_VCF_connection)

  n_barcodes <- length(readLines(file.path(cellsnp_dir, "cellSNP.samples.tsv")))
  read_entries <- function(dir, tag) {
    matrix(scan(file.path(dir, paste0("cellSNP.tag.", tag, ".mtx")), what = integer(), skip = 3, quiet = TRUE), ncol = 3, byrow = TRUE)
  }
  for (tag in c("AD", "DP", "OTH")) {
    other_entries <- read_entries(other_cellsnp_dir, tag)
    other_entries[, 1] <- other_entries[, 1] + sum(!is_header)
    entries <- rbind(read_entries(cellsnp_dir, tag), other_entries)
    entries[, 1] <- new_SNP_index[entries[, 1]]
    entries <- entries[order(entries[, 1], entries[, 2]), , drop = FALSE]
    mtx_file <- file.path(cellsnp_dir, paste0("cellSNP.tag.", tag, ".mtx"))
    writeLines(c("%%MatrixMarket matrix coordinate integer general", "%", paste(length(records), n_barcodes, nrow(entries), sep = "\t")), mtx_file)
    data.table::fwrite(as.data.frame(entries), mtx_file, append = TRUE, sep = "\t", col.names = FALSE)
  }
  invisible(cellsnp_dir)
}


#' Get vireo donor ids tibble
#'
#' Run vireo on cellsnp-lite output and return donor-probability metadata keyed by barcode.
#'
#' @param cellsnp_dir Directory produced by cellsnp-lite, containing genotype likelihood files consumed by vireo.
#' @param GEM_well_donors_VCF_file VCF path containing donor genotypes for the GEM well.
#' @param GEM_well_n_donors Expected donor count passed to vireo.
#' @param GEM_well_donor_id Configured donor IDs used to name vireo probability columns.
#' @param cellranger_barcodes_tsv Path to a TSV containing CellRanger-called barcodes, used as the cellsnp-lite sample whitelist.
#' @param cores Number of CPU cores requested for external tools or parallel work.
#' @param do_learn Logical; when TRUE, let vireo learn donor genotypes instead of only using supplied priors.
#' @param learn_iterations Number of vireo genotype-learning iterations when `do_learn` is TRUE.
#' @return Tibble of barcodes, donor calls, donor probabilities, and configured donor labels.
#' @keywords internal

get_vireo_donor_ids_tibble <- function(
  cellsnp_dir,
  GEM_well_donors_VCF_file,
  GEM_well_n_donors,
  GEM_well_donor_id,
  cellranger_barcodes_tsv,
  cores,
  do_learn = FALSE,
  learn_iterations = NULL
) {
  vireo_out_dir <- get_structured_file_path(filetype = NULL)

  if (length(cellsnp_dir) == 0L) {
    # skip demultiplexing and generate dummy donor_id tibble
    vireo_donor_ids_tibble <- tibble::tibble(
      cell = cellranger_barcodes_tsv %>% readr::read_tsv(col_names = FALSE) %>% dplyr::pull(),
      donor_id = GEM_well_donor_id,
      prob_max = NA,
      prob_doublet = NA,
      n_vars = NA,
      best_doublet = NA,
      doublet_logLikRatio = NA,
      vireo_type = "not_demultiplexed"
    )
  } else {
    # setup flags and run vireo
    learn_args <- if (do_learn) {
      c("--forceLearnGT", "-M", as.character(learn_iterations))
    } else {
      NULL
    }

    genotype_args <- if (!is.na(GEM_well_donors_VCF_file)) {
      c("-d", GEM_well_donors_VCF_file, "-t", "GT")
    } else {
      NULL
    }

    run_w_error_check(
      command_string = "vireo",
      arguments_chr = c("-p", cores, "-c", cellsnp_dir, "-N", GEM_well_n_donors, "-o", vireo_out_dir, genotype_args, learn_args)
    )

    # read vireo output
    vireo_donor_ids_tibble <- file.path(vireo_out_dir, "donor_ids.tsv") %>%
      readr::read_tsv(col_types = "ccdddccd") %>%
      dplyr::mutate(
        vireo_type = dplyr::case_when(!donor_id %in% c("doublet", "unassigned") ~ "singlet", .default = donor_id),
        donor_id = best_singlet # effectively discard original donor_id, and create a new one based on the best singlet
      )
  }

  # format
  vireo_donor_ids_tibble %>%
    dplyr::rename(
      barcode = cell,
      vireo_max_prob_singlet = prob_max,
      vireo_max_prob_doublet = prob_doublet,
      vireo_n_vars = n_vars,
      vireo_best_doublet = best_doublet,
      vireo_doublet_logLikRatio = doublet_logLikRatio
    )
}


#' Map donors to samples
#'
#' Replace each nucleus's donor with that donor's sample in this GEM well, so
#' that a donor who gives several samples, such as one per condition, enters the
#' aggregation as separate biological replicates.
#'
#' @param vireo_donor_ids_tibble Per-barcode donor assignments from `get_vireo_donor_ids_tibble()`.
#' @param sample_ids Named character vector giving each donor's sample in this GEM well.
#' @return `vireo_donor_ids_tibble` with each `donor_id` replaced by its sample.
#' @keywords internal

map_donors_to_samples <- function(vireo_donor_ids_tibble, sample_ids) {
  unmapped_donor_ids <- setdiff(stats::na.omit(vireo_donor_ids_tibble$donor_id), names(sample_ids))
  if (length(unmapped_donor_ids) > 0L) {
    stop("cfg_GEM_well_samples.tsv gives no sample for donor(s): ", paste(unmapped_donor_ids, collapse = ", "), call. = FALSE)
  }
  dplyr::mutate(vireo_donor_ids_tibble, donor_id = unname(sample_ids[.data$donor_id]))
}
