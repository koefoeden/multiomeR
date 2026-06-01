# Demultiplexing by genotype ----
#' Get cellsnp dir
#'
#' Run cellsnp-lite on one reaction BAM using CellRanger-called barcodes and donor VCF variants.
#'
#' @param bam_file Input BAM path for one 10x reaction.
#' @param cellranger_barcodes_tsv Path to a TSV containing CellRanger-called barcodes, used as the cellsnp-lite sample whitelist.
#' @param reaction_donors_VCF_file VCF path containing donor genotypes for the reaction.
#' @param cores Number of CPU cores requested for external tools or parallel work.
#' @param minMAF Minimum minor allele frequency threshold passed to cellsnp-lite.
#' @param minCOUNT Minimum SNP read-count threshold passed to cellsnp-lite.
#' @return Path to the cellsnp-lite output directory containing `cellSNP.base.vcf.gz`.
#' @keywords internal

get_cellsnp_dir <- function(
  bam_file,
  cellranger_barcodes_tsv,
  reaction_donors_VCF_file,
  cores,
  minMAF = 0.1,
  minCOUNT = 20
) {
  cellsnp_out_dir <- get_structured_file_path(filetype = NULL)

  run_w_error_check(
    command_string = "cellsnp-lite",
    arguments_chr = c(
      c("-s", bam_file),
      c("-b", cellranger_barcodes_tsv),
      c("-O", cellsnp_out_dir),
      c("-R", reaction_donors_VCF_file),
      c("-p", cores),
      c("--minMAF", minMAF),
      c("--minCOUNT", minCOUNT),
      c("--UMItag", "None"),
      c("--gzip")
    )
  )

  vcf_out_lines <- cellsnp_out_dir %>% file.path("cellSNP.base.vcf.gz") %>% readLines() %>% length()

  if (vcf_out_lines <= 1) {
    stop("cellsnp-lite failed to run. Check the output directory for errors. VCF file contains ", vcf_out_lines, " line().")
  }
  return(cellsnp_out_dir)
}


#' Get vireo donor ids tibble
#'
#' Run vireo on cellsnp-lite output and return donor-probability metadata keyed by barcode.
#'
#' @param cellsnp_dir Directory produced by cellsnp-lite, containing genotype likelihood files consumed by vireo.
#' @param reaction_donors_VCF_file VCF path containing donor genotypes for the reaction.
#' @param reaction_n_donors Expected donor count passed to vireo.
#' @param reaction_donor_id Configured donor IDs used to name vireo probability columns.
#' @param cellranger_barcodes_tsv Path to a TSV containing CellRanger-called barcodes, used as the cellsnp-lite sample whitelist.
#' @param cores Number of CPU cores requested for external tools or parallel work.
#' @param do_learn Logical; when TRUE, let vireo learn donor genotypes instead of only using supplied priors.
#' @param learn_iterations Number of vireo genotype-learning iterations when `do_learn` is TRUE.
#' @return Tibble of barcodes, donor calls, donor probabilities, and configured donor labels.
#' @keywords internal

get_vireo_donor_ids_tibble <- function(
  cellsnp_dir,
  reaction_donors_VCF_file,
  reaction_n_donors,
  reaction_donor_id,
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
      donor_id = reaction_donor_id,
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

    genotype_args <- if (!is.na(reaction_donors_VCF_file)) {
      c("-d", reaction_donors_VCF_file, "-t", "GT")
    } else {
      NULL
    }

    run_w_error_check(
      command_string = "vireo",
      arguments_chr = c("-p", cores, "-c", cellsnp_dir, "-N", reaction_n_donors, "-o", vireo_out_dir, genotype_args, learn_args)
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
