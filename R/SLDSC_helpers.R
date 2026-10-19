#' Download the S-LDSC GRCh38 reference
#'
#' Extract the 1000 Genomes phase 3 EUR PLINK files, baseline-LD v2.2 model and
#' HapMap3 regression weights that ldsc's authors lifted to GRCh38.
#'
#' @param url Zenodo URL of the `GRCh38.tgz` archive.
#' @param md5 Expected MD5 checksum of the archive.
#' @return The reference directory.
#' @keywords internal

download_SLDSC_reference <- function(url, md5) {
  reference_dir <- file.path(targets::tar_config_get("store"), "files", "S-LDSC", "GRCh38")
  archive_dir <- withr::local_tempdir()
  archive <- file.path(archive_dir, "GRCh38.tgz")
  curl::curl_download(url, archive)
  if (!identical(unname(tools::md5sum(archive)), md5)) {
    stop("S-LDSC reference archive checksum mismatch: ", url)
  }
  utils::untar(archive, exdir = archive_dir)
  unlink(reference_dir, recursive = TRUE)
  for (part in c("baselineLD_v2.2", "plink_files", "weights")) {
    utils::untar(file.path(archive_dir, "GRCh38", paste0(part, ".tgz")), exdir = reference_dir)
  }
  reference_dir
}

#' Read the baseline-LD model at the S-LDSC regression SNPs
#'
#' Regression SNPs are the HapMap3 SNPs outside the MHC that carry both
#' baseline-LD scores and regression weights, kept in reference order so that
#' jackknife blocks match ldsc's.
#'
#' @param reference_dir Directory returned by `download_SLDSC_reference()`.
#' @return A list with the regression SNP tibble (including reference alleles
#'   and weight LD scores), the baseline LD-score matrix and the per-annotation
#'   counts of common (MAF > 5%) reference SNPs.
#' @keywords internal

read_SLDSC_baseline_record <- function(reference_dir) {
  LD_score_tibble <- purrr::map(1:22, \(chromosome) {
    LD_score_tibble <- data.table::fread(file.path(reference_dir, "baselineLD_v2.2", paste0("baselineLD.", chromosome, ".l2.ldscore.gz")))
    weight_tibble <- data.table::fread(file.path(reference_dir, "weights", paste0("weights.hm3_noMHC.", chromosome, ".l2.ldscore.gz")))
    bim <- read_SLDSC_bim(reference_dir, chromosome)
    bim_idx <- match(LD_score_tibble$SNP, bim$SNP)
    LD_score_tibble[, `:=`(
      weight_L2 = weight_tibble$L2[match(SNP, weight_tibble$SNP)],
      A1 = bim$A1[bim_idx],
      A2 = bim$A2[bim_idx]
    )][!is.na(weight_L2)]
  }) |>
    data.table::rbindlist()
  annotation_cols <- setdiff(names(LD_score_tibble), c("CHR", "SNP", "BP", "weight_L2", "A1", "A2"))
  M_5_50 <- purrr::map(1:22, \(chromosome) scan(
    file.path(reference_dir, "baselineLD_v2.2", paste0("baselineLD.", chromosome, ".l2.M_5_50")),
    quiet = TRUE
  )) |>
    purrr::reduce(`+`)
  list(
    SNP_tibble = tibble::as_tibble(LD_score_tibble[, .(CHR, SNP, BP, A1, A2, weight_L2)]),
    L2_matrix = as.matrix(LD_score_tibble[, ..annotation_cols]),
    M_5_50 = rlang::set_names(M_5_50, annotation_cols)
  )
}

read_SLDSC_bim <- function(reference_dir, chromosome) {
  data.table::fread(
    file.path(reference_dir, "plink_files", paste0("1000G.EUR.hg38.", chromosome, ".bim")),
    header = FALSE,
    col.names = c("CHR", "SNP", "CM", "BP", "A1", "A2")
  )
}

#' Find the GWAS Catalog summary statistics of a study
#'
#' Prefer the GWAS Catalog's GRCh38 harmonised file, in GWAS-SSF format when
#' the study also has an older-format one, and fall back to a submitted GWAS-SSF
#' file.
#'
#' @param study_id GWAS Catalog study accession.
#' @return The file URL, or `NA` when the study has no genome-wide summary
#'   statistics.
#' @keywords internal

get_GWAS_catalog_sumstats_url <- function(study_id) {
  study_number <- as.integer(stringr::str_remove(study_id, "^GCST"))
  digits <- nchar(study_id) - 4L
  first_study <- (study_number - 1L) %/% 1000L * 1000L + 1L
  study_url <- sprintf(
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST%0*d-GCST%0*d/%s/",
    digits, first_study, digits, first_study + 999L, study_id
  )
  # A study without summary statistics has no directory.
  list_files <- \(url) stringr::str_match_all(read_URL_text(url) %||% "", 'href="([^"?/][^"]*)"')[[1]][, 2]
  harmonised_files <- stringr::str_subset(list_files(paste0(study_url, "harmonised/")), "\\.h\\.tsv\\.gz$")
  # Older-format files, named after the PubMed ID, lack per-variant sample sizes.
  harmonised_files <- c(intersect(paste0(study_id, ".h.tsv.gz"), harmonised_files), harmonised_files)
  if (length(harmonised_files) > 0L) {
    return(paste0(study_url, "harmonised/", harmonised_files[[1]]))
  }
  submitted_files <- stringr::str_subset(list_files(study_url), "\\.tsv(\\.gz)?$")
  if (length(submitted_files) > 0L) paste0(study_url, submitted_files[[1]]) else NA_character_
}

#' Read a text file from a URL
#'
#' @param url File URL.
#' @return The file's text, or `NULL` when the server reports it missing.
#' @keywords internal

read_URL_text <- function(url) {
  # Server errors are retried; they must not pass as missing files.
  response <- retry_EBI_request(\() {
    response <- curl::curl_fetch_memory(url)
    if (!response$status_code %in% c(200L, 404L)) {
      stop("HTTP status ", response$status_code, ": ", url)
    }
    response
  })
  if (response$status_code == 404L) NULL else rawToChar(response$content)
}

#' Retry a request to EBI's FTP server
#'
#' The server refuses connections while busy, so failed requests are retried
#' after randomized pauses that double from about half a minute up to about ten
#' minutes.
#'
#' @param request Function performing the request.
#' @param attempts Maximum number of attempts.
#' @return The request's value.
#' @keywords internal

retry_EBI_request <- function(request, attempts = 8L) {
  for (attempt in seq_len(attempts)) {
    result <- tryCatch(request(), error = identity)
    if (!inherits(result, "error") || attempt == attempts) {
      break
    }
    Sys.sleep(min(30 * 2^(attempt - 1), 600) * stats::runif(1, 0.5, 1.5))
  }
  if (inherits(result, "error")) stop(result) else result
}

#' Download a file over parallel byte ranges
#'
#' EBI's FTP server serves each connection at a few MB/s, so a server that
#' accepts byte ranges sends the file over several connections at once.
#'
#' @param url File URL.
#' @param path Destination file.
#' @param n_ranges Number of connections, at most the six per host of curl's
#'   connection pool.
#' @keywords internal

download_in_ranges <- function(url, path, n_ranges = 6L) {
  headers <- curl::parse_headers_list(curl::curl_fetch_memory(url, curl::new_handle(nobody = TRUE))$headers)
  if (!identical(headers[["accept-ranges"]], "bytes")) {
    return(curl::curl_download(url, path))
  }
  size <- as.numeric(headers[["content-length"]])
  bounds <- floor(seq(0, size, length.out = n_ranges + 1L))
  parts <- paste0(path, ".part", seq_len(n_ranges))
  on.exit(unlink(parts))
  handles <- purrr::map(seq_len(n_ranges), \(i) {
    curl::new_handle(url = url, range = sprintf("%.0f-%.0f", bounds[[i]], bounds[[i + 1L]] - 1))
  })
  # Multiplexing would queue the ranges on one HTTP/1.1 connection.
  result <- curl::multi_download(handles, parts, progress = FALSE, multiplex = FALSE)
  if (!all(result$status_code == 206L) || sum(file.size(parts)) != size) {
    stop("Ranged download failed: ", url)
  }
  unlink(path)
  file.append(path, parts)
}

#' Summarize the samples of GWAS-SSF metadata
#'
#' @param metadata_text Text of a GWAS-SSF `-meta.yaml` file, or `NULL`.
#' @return One row with the total sample size and the shares of the ancestry
#'   groups that Open Targets study metadata also report, or `NULL` when the
#'   metadata lists no samples.
#' @keywords internal

get_GWAS_SSF_metadata_tibble <- function(metadata_text) {
  # The GWAS Catalog lists some sample groups several times.
  samples <- unique(yaml::yaml.load(metadata_text %||% "")$samples)
  if (length(samples) == 0L) {
    return(NULL)
  }
  sample_sizes <- purrr::map_dbl(samples, \(sample) as.numeric(sample$sample_size %||% NA))
  categories <- purrr::map(samples, \(sample) as.character(unlist(sample$sample_ancestry_category)))
  # A sample group with several ancestry categories counts equally towards each.
  category_sizes <- rep(sample_sizes / lengths(categories), lengths(categories))
  ancestry_groups <- dplyr::recode_values(
    unlist(categories),
    c("European", "Finnish") ~ "EUR",
    "East Asian" ~ "EAS",
    c("African American or Afro-Caribbean", "African unspecified", "Sub-Saharan African") ~ "AFR",
    "Hispanic or Latin American" ~ "AMR",
    "South Asian" ~ "SAS",
    default = "OTH"
  )
  ancestry_shares <- tapply(
    category_sizes, factor(ancestry_groups, c("EUR", "EAS", "AFR", "AMR", "SAS", "OTH")), sum, default = 0
  ) / sum(category_sizes)
  tibble::tibble(
    sample_size = sum(sample_sizes),
    !!!rlang::set_names(as.list(ancestry_shares), paste0("ancestry_", names(ancestry_shares)))
  )
}

#' Fetch one GWAS source's summary statistics at the S-LDSC regression SNPs
#'
#' Use the configured `summaryStatistics` file, URL or GWAS Catalog accession
#' and otherwise, for a GWAS Catalog `sourceId`, the Catalog's summary
#' statistics. Only their S-LDSC columns at the regression SNPs are kept, as
#' text, so that a munging change reuses them instead of downloading again.
#'
#' @param sumstats_source One row with `sourceId` and `summaryStatistics`.
#' @param sumstats_files Tracked local summary-statistics file and its
#'   GWAS-SSF metadata, or `character()`.
#' @param regression_SNPs Regression SNP identifiers.
#' @return A record with the source, the summary statistics location and, as
#'   available, the text of their GWAS-SSF metadata and the kept rows, which
#'   are missing for a source without genome-wide summary statistics or a file
#'   without rsIDs.
#' @keywords internal

fetch_SLDSC_sumstats_rows <- function(sumstats_source, sumstats_files, regression_SNPs) {
  location <- dplyr::coalesce(sumstats_source$summaryStatistics, sumstats_source$sourceId)
  sumstats_url <- if (stringr::str_detect(location, "^GCST[0-9]+$")) {
    get_GWAS_catalog_sumstats_url(location)
  } else {
    sumstats_source$summaryStatistics
  }
  record <- list(sourceId = sumstats_source$sourceId, sumstats_url = sumstats_url)
  if (is.na(sumstats_url)) {
    return(record)
  }
  if (length(sumstats_files) > 0L) {
    sumstats_path <- sumstats_files[[1]]
    metadata_text <- if (length(sumstats_files) == 2L) readr::read_file(sumstats_files[[2]])
  } else {
    sumstats_path <- withr::local_tempfile(fileext = paste0(".", basename(sumstats_url)))
    retry_EBI_request(\() download_in_ranges(sumstats_url, sumstats_path))
    metadata_text <- read_URL_text(paste0(sumstats_url, "-meta.yaml"))
  }
  record$metadata_text <- metadata_text
  header <- strsplit(readLines(sumstats_path, n = 1L), "\t")[[1]]
  pick <- \(candidates) candidates[candidates %in% header][1]
  columns <- c(
    SNP = pick(c("hm_rsid", "rsid", "rs_id")),
    effect_allele = pick(c("hm_effect_allele", "effect_allele")),
    other_allele = pick(c("hm_other_allele", "other_allele")),
    P = "p_value",
    N = pick(c("n", "sample_size")),
    FRQ = pick(c("hm_effect_allele_frequency", "effect_allele_frequency"))
  )
  # GWAS-SSF makes rsIDs optional in submitted files.
  if (is.na(columns[["SNP"]])) {
    return(record)
  }
  columns <- columns[!is.na(columns)]
  # Streaming keeps memory low; values stay text so that munging decides how to
  # parse them. A single thread keeps the file's row order, which decides the
  # duplicate SNPs that munging keeps.
  withr::local_options(arrow.use_threads = FALSE)
  record$sumstats_tibble <- arrow::open_delim_dataset(
    sumstats_path,
    delim = "\t",
    convert_options = arrow::csv_convert_options(
      col_types = arrow::schema(!!!rlang::set_names(rep(list(arrow::utf8()), length(columns)), columns)),
      null_values = c("", "NA", "#NA"),
      strings_can_be_null = TRUE
    )
  ) |>
    dplyr::select(dplyr::all_of(columns)) |>
    dplyr::filter(SNP %in% regression_SNPs) |>
    dplyr::collect()
  record
}

#' Munge GWAS summary statistics for S-LDSC
#'
#' Follow ldsc's `munge_sumstats.py --merge-alleles` filtering order on the
#' text columns of a GWAS Catalog harmonised or GWAS-SSF file at the
#' regression SNPs: keep complete, unambiguous, reference-matching alleles, MAF
#' > 1% when reported and, for per-SNP sample sizes, N of at least two thirds of
#' the 90th percentile.
#'
#' @param sumstats_tibble Text columns `SNP`, `effect_allele`, `other_allele`,
#'   `P` and optionally `N` and `FRQ` from `fetch_SLDSC_sumstats_rows()`.
#' @param regression_SNP_tibble Regression SNPs with reference alleles.
#' @return One row per retained SNP with its chi-square statistic and sample
#'   size, which is missing when the file has no `n` or `sample_size` column.
#' @keywords internal

munge_GWAS_sumstats <- function(sumstats_tibble, regression_SNP_tibble) {
  sumstats_tibble <- sumstats_tibble |>
    # Some harmonised files carry optional columns without any values.
    dplyr::select(dplyr::where(\(x) !all(is.na(x)))) |>
    stats::na.omit() |>
    dplyr::mutate(
      # P-values below the double range, such as 3.35E-1003, parse to zero and
      # are removed with other P-values outside (0, 1], as ldsc removes them.
      dplyr::across(dplyr::any_of(c("P", "N", "FRQ")), \(x) suppressWarnings(as.numeric(x))),
      effect_allele = toupper(effect_allele),
      other_allele = toupper(other_allele)
    )
  complement <- c(A = "T", C = "G", G = "C", T = "A")
  sumstats_tibble <- sumstats_tibble |>
    dplyr::filter(
      P > 0, P <= 1,
      effect_allele %in% names(complement), other_allele %in% names(complement),
      effect_allele != other_allele, complement[effect_allele] != other_allele,
      if ("FRQ" %in% names(sumstats_tibble)) pmin(FRQ, 1 - FRQ) > 0.01 else TRUE
    ) |>
    dplyr::distinct(SNP, .keep_all = TRUE)
  sumstats_tibble <- if ("N" %in% names(sumstats_tibble)) {
    dplyr::filter(sumstats_tibble, N >= stats::quantile(N, 0.9) / 1.5)
  } else {
    dplyr::mutate(sumstats_tibble, N = NA_real_)
  }
  reference <- regression_SNP_tibble[match(sumstats_tibble$SNP, regression_SNP_tibble$SNP), ]
  allele_pair <- \(a, b) paste(pmin(a, b), pmax(a, b))
  reference_pair <- allele_pair(reference$A1, reference$A2)
  matched <- allele_pair(sumstats_tibble$effect_allele, sumstats_tibble$other_allele) == reference_pair |
    allele_pair(complement[sumstats_tibble$effect_allele], complement[sumstats_tibble$other_allele]) == reference_pair
  tibble::tibble(
    SNP = sumstats_tibble$SNP,
    chi2 = stats::qchisq(sumstats_tibble$P, df = 1, lower.tail = FALSE),
    N = sumstats_tibble$N
  )[matched, ]
}

#' Resolve a configured local summary-statistics file
#'
#' @param summaryStatistics Configured `summaryStatistics` path, URL or GWAS
#'   Catalog accession, or `NA`.
#' @return The normalized local filename, followed by its GWAS-SSF metadata
#'   file when present, or `character()` for a URL, an accession or no setting.
#' @keywords internal

resolve_SLDSC_sumstats_files <- function(summaryStatistics) {
  if (is.na(summaryStatistics) || stringr::str_detect(summaryStatistics, "^([a-z]+://|GCST[0-9]+$)")) {
    return(character())
  }
  files <- paste0(summaryStatistics, c("", "-meta.yaml"))
  normalizePath(files[c(TRUE, file.exists(files[[2]]))], mustWork = TRUE)
}

#' Munge one GWAS source's fetched S-LDSC summary statistics
#'
#' SNPs without a per-SNP sample size take the total sample size of the
#' GWAS-SSF metadata, so that N always describes the summary statistics, never
#' a credible-set study.
#'
#' @param sumstats_rows_record Record from `fetch_SLDSC_sumstats_rows()`.
#' @param regression_SNP_tibble Regression SNPs with reference alleles.
#' @return A record with the source, the summary statistics location, the
#'   sample size and ancestry shares of their metadata, and the munged
#'   statistics, which are `NULL` when the source has no usable genome-wide
#'   summary statistics.
#' @keywords internal

get_SLDSC_sumstats_record <- function(sumstats_rows_record, regression_SNP_tibble) {
  record <- c(
    sumstats_rows_record[c("sourceId", "sumstats_url")],
    list(metadata_tibble = get_GWAS_SSF_metadata_tibble(sumstats_rows_record$metadata_text))
  )
  if (is.null(sumstats_rows_record$sumstats_tibble)) {
    return(record)
  }
  sumstats_tibble <- munge_GWAS_sumstats(sumstats_rows_record$sumstats_tibble, regression_SNP_tibble)
  if (anyNA(sumstats_tibble$N)) {
    sample_size <- record$metadata_tibble$sample_size %||% NA_real_
    if (is.na(sample_size)) {
      stop("Summary statistics without per-variant n need GWAS-SSF metadata with sample sizes: ", record$sumstats_url)
    }
    sumstats_tibble$N <- sample_size
  }
  record$sumstats_tibble <- sumstats_tibble
  record
}

#' Annotate peaks with cell-type accessibility specificity
#'
#' Give each peak an all-peaks indicator and, per cell type, its share of the
#' depth-normalized pseudobulk accessibility summed over cell types. This share
#' is chromVAR's observed-over-expected accessibility ratio under equal cell-type
#' weights, divided by the number of cell types.
#'
#' @param pseudobulk_counts_matrix Peak-by-cell-type ATAC counts.
#' @return Dense peak-by-annotation matrix for peaks with any counts.
#' @keywords internal

get_SLDSC_peak_annotation_matrix <- function(pseudobulk_counts_matrix) {
  counts_matrix <- methods::as(pseudobulk_counts_matrix, "dgCMatrix")
  counts_matrix <- as.matrix(counts_matrix[Matrix::rowSums(counts_matrix) > 0, , drop = FALSE])
  normalized_matrix <- sweep(counts_matrix, 2, colSums(counts_matrix), "/")
  cbind(all_peaks = 1, normalized_matrix / rowSums(normalized_matrix))
}

#' Score gene windows by cell-type expression specificity
#'
#' Sum each cell type's GEX pseudobulks, divide them by the cell type's total
#' counts and give each gene its share of these depth-normalized profiles
#' summed over cell types, as for the peaks. A gene counts as expressed in a
#' cell type with at least 10 counts and 1 count per million there; its
#' specificity is zero elsewhere. Only autosomal protein-coding genes expressed
#' in any cell type are kept.
#'
#' @param pseudobulk_counts_matrix Gene-by-pseudobulk GEX counts with columns
#'   named `<cell type>_<donor>`.
#' @param Ensembl_annotations_GRanges_list List with the Ensembl `genes`.
#' @param window Bases added to each side of the gene body.
#' @return A list with the gene-by-cell-type specificity matrix and the gene
#'   windows, both named by GEX feature.
#' @keywords internal

get_SLDSC_gene_specificity_record <- function(pseudobulk_counts_matrix, Ensembl_annotations_GRanges_list, window = 100000L) {
  counts_matrix <- methods::as(pseudobulk_counts_matrix, "dgCMatrix")
  clusters <- get_pseudobulk_sample_tibble(counts_matrix)$cluster
  counts_matrix <- as.matrix(counts_matrix %*% Matrix::t(Matrix::fac2sparse(factor(clusters))))
  normalized_matrix <- sweep(counts_matrix, 2, colSums(counts_matrix), "/")
  gene_tibble <- make_peak_gene_correlation_gene_TSS_tibble(Ensembl_annotations_GRanges_list, counts_matrix) |>
    dplyr::filter(gene_biotype == "protein_coding", chr %in% paste0("chr", 1:22))
  expressed <- counts_matrix[gene_tibble$gene_matrix_feature, ] >= 10 &
    normalized_matrix[gene_tibble$gene_matrix_feature, ] >= 1e-6
  gene_tibble <- gene_tibble[rowSums(expressed) > 0, ]
  expressed <- expressed[gene_tibble$gene_matrix_feature, ]
  normalized_matrix <- normalized_matrix[gene_tibble$gene_matrix_feature, ]
  list(
    specificity_matrix = normalized_matrix / rowSums(normalized_matrix) * expressed,
    ranges = GenomicRanges::GRanges(gene_tibble$chr, IRanges::IRanges(
      pmax(gene_tibble$gene_start - window, 1L), gene_tibble$gene_end + window, names = gene_tibble$gene_matrix_feature
    ))
  )
}

#' Annotate gene windows of each cell type's most specific genes
#'
#' Each cell type's annotation marks the windows of its most specific expressed
#' genes, 10% of the genes scored (LDSC-SEG); an all-genes annotation marks the
#' windows of all of them.
#'
#' @param gene_specificity_record Record from
#'   `get_SLDSC_gene_specificity_record()`.
#' @return A list with the binary gene-by-annotation matrix and the gene
#'   windows, both named by GEX feature.
#' @keywords internal

get_SLDSC_gene_annotation_record <- function(gene_specificity_record) {
  specificity_matrix <- gene_specificity_record$specificity_matrix
  n_top <- round(0.1 * nrow(specificity_matrix))
  top_matrix <- apply(specificity_matrix, 2, \(share) rank(-share, ties.method = "first") <= n_top & share > 0)
  list(annotation_matrix = cbind(all_genes = 1, top_matrix * 1), ranges = gene_specificity_record$ranges)
}

#' Read PLINK genotypes as allele dosages
#'
#' @param bed_path SNP-major PLINK `.bed` file.
#' @param n_samples Number of samples in the `.fam` file.
#' @param n_SNPs Number of SNPs in the `.bim` file.
#' @param SNP_idx Indices of the SNPs to read.
#' @return A sample-by-SNP dosage matrix with missing genotypes as `NA`.
#' @keywords internal

read_plink_bed_dosages <- function(bed_path, n_samples, n_SNPs, SNP_idx) {
  bytes_per_SNP <- ceiling(n_samples / 4)
  bed_raw <- readBin(bed_path, "raw", n = 3L + bytes_per_SNP * n_SNPs)
  byte_idx <- 3L + outer(seq_len(bytes_per_SNP), (SNP_idx - 1L) * bytes_per_SNP, "+")
  # Each byte packs four genotypes, low bits first: 00 is two A1 alleles, 10 one,
  # 11 none and 01 missing.
  dosage_lookup <- vapply(
    0:255,
    \(byte) c(2, NA, 1, 0)[bitwAnd(bitwShiftR(byte, c(0L, 2L, 4L, 6L)), 3L) + 1L],
    numeric(4)
  )
  dosages <- dosage_lookup[, as.integer(bed_raw[byte_idx]) + 1L]
  dim(dosages) <- c(4L * bytes_per_SNP, length(SNP_idx))
  dosages[seq_len(n_samples), , drop = FALSE]
}

#' Compute one chromosome's LD scores of range annotations
#'
#' Reproduce `ldsc.py --l2 --ld-wind-cm 1` for annotations that are non-zero
#' only at reference SNPs inside the ranges, computing only the regression
#' SNPs' scores: ldsc's bias-corrected r^2 summed over its window of whole
#' 50-SNP chunks that cover 1 cM, after dropping monomorphic SNPs.
#'
#' @param chromosome Autosome number.
#' @param annotation_matrix Range-by-annotation values.
#' @param annotation_ranges Non-overlapping ranges in the order of the
#'   annotation rows.
#' @param reference_dir Directory returned by `download_SLDSC_reference()`.
#' @param regression_SNPs Regression SNP identifiers.
#' @param block_size Regression SNPs per cross-product block.
#' @return A list with the regression SNPs' LD-score matrix and each
#'   annotation's sum and sum of squares over common (MAF > 5%) reference SNPs.
#' @keywords internal

get_chromosome_SLDSC_LD_score_record <- function(chromosome, annotation_matrix, annotation_ranges,
  reference_dir, regression_SNPs, block_size = 500L) {
  prefix <- file.path(reference_dir, "plink_files", paste0("1000G.EUR.hg38.", chromosome))
  bim <- read_SLDSC_bim(reference_dir, chromosome)
  n_bim_SNPs <- nrow(bim)
  bim <- bim[, bim_idx := .I][data.table::fread(paste0(prefix, ".frq"))$MAF > 0]
  n_samples <- length(readLines(paste0(prefix, ".fam")))
  # ldsc pairs each 50-SNP chunk with whole chunks back to the 1 cM window of
  # the chunk's first SNP; the SNPs within 1 cM of the chromosome start form one
  # block. Each SNP's window is therefore the contiguous range of SNPs from its
  # left limit to the last SNP whose left limit it reaches.
  chunk_size <- 50L
  block_left <- findInterval(bim$CM - 1, bim$CM, left.open = TRUE)
  first_block_end <- min(ceiling(match(TRUE, block_left > 0, nomatch = nrow(bim)) / chunk_size) * chunk_size, nrow(bim))
  chunk_start <- (seq_len(nrow(bim)) - 1L) %/% chunk_size * chunk_size
  left_limit <- 1L + ifelse(
    seq_len(nrow(bim)) <= first_block_end,
    0L,
    chunk_start - ceiling((chunk_start - block_left[chunk_start + 1L]) / chunk_size) * chunk_size
  )
  right_limit <- findInterval(seq_len(nrow(bim)), left_limit)

  SNP_range_idx <- GenomicRanges::findOverlaps(
    GenomicRanges::GRanges(paste0("chr", chromosome), IRanges::IRanges(bim$BP, width = 1L)),
    annotation_ranges,
    select = "first"
  )
  annotated_idx <- which(!is.na(SNP_range_idx))
  regression_idx <- which(bim$SNP %in% regression_SNPs)
  genotype_idx <- sort(union(annotated_idx, regression_idx))
  dosages <- read_plink_bed_dosages(paste0(prefix, ".bed"), n_samples, n_bim_SNPs, bim$bim_idx[genotype_idx])
  allele_mean <- colMeans(dosages, na.rm = TRUE)
  missing_idx <- which(is.na(dosages), arr.ind = TRUE)
  dosages[missing_idx] <- allele_mean[missing_idx[, 2]]
  genotype_z <- scale(dosages, center = allele_mean, scale = sqrt(colMeans(dosages^2) - allele_mean^2))
  rm(dosages)
  SNP_annotation_matrix <- annotation_matrix[SNP_range_idx[annotated_idx], , drop = FALSE]
  annotated_cols <- match(annotated_idx, genotype_idx)
  regression_cols <- match(regression_idx, genotype_idx)

  L2_matrix <- matrix(0, length(regression_idx), ncol(SNP_annotation_matrix),
    dimnames = list(bim$SNP[regression_idx], colnames(SNP_annotation_matrix)))
  for (block_start in seq(1L, length(regression_idx), by = block_size)) {
    block <- block_start:min(block_start + block_size - 1L, length(regression_idx))
    block_idx <- regression_idx[block]
    window <- which(annotated_idx >= min(left_limit[block_idx]) & annotated_idx <= max(right_limit[block_idx]))
    in_window <- outer(left_limit[block_idx], annotated_idx[window], "<=") &
      outer(right_limit[block_idx], annotated_idx[window], ">=")
    r2 <- (crossprod(genotype_z[, regression_cols[block], drop = FALSE], genotype_z[, annotated_cols[window], drop = FALSE]) / n_samples)^2
    L2_matrix[block, ] <- ((r2 - (1 - r2) / (n_samples - 2)) * in_window) %*% SNP_annotation_matrix[window, , drop = FALSE]
  }
  MAF <- pmin(allele_mean, 2 - allele_mean)[annotated_cols] / 2
  common_annotation_matrix <- SNP_annotation_matrix[MAF > 0.05, , drop = FALSE]
  list(
    L2_matrix = L2_matrix,
    M_5_50 = colSums(common_annotation_matrix),
    sum_sq_5_50 = colSums(common_annotation_matrix^2)
  )
}

#' Compute genome-wide LD scores of range annotations
#'
#' A reference SNP takes each annotation's largest value over the ranges that
#' contain it, which for binary annotations of overlapping ranges, such as gene
#' windows, marks their union.
#'
#' @param annotation_matrix Range-by-annotation values named by range.
#' @param annotation_ranges Ranges named like the annotation rows.
#' @inheritParams get_chromosome_SLDSC_LD_score_record
#' @param regression_SNP_tibble Regression SNPs in reference order.
#' @param threads BLAS threads for the genotype cross-products.
#' @return A list with the LD-score matrix aligned to the regression SNPs and
#'   each annotation's sum and sum of squares over common reference SNPs.
#' @keywords internal

get_SLDSC_LD_score_record <- function(annotation_matrix, annotation_ranges, reference_dir, regression_SNP_tibble, threads = 1L) {
  previous_threads <- RhpcBLASctl::blas_get_num_procs()
  RhpcBLASctl::blas_set_num_threads(threads)
  on.exit(RhpcBLASctl::blas_set_num_threads(previous_threads), add = TRUE)
  pieces <- GenomicRanges::disjoin(annotation_ranges[rownames(annotation_matrix)], with.revmap = TRUE)
  chromosome_records <- purrr::map(1:22, get_chromosome_SLDSC_LD_score_record,
    annotation_matrix = apply(annotation_matrix, 2, \(values) max(IRanges::extractList(values, pieces$revmap))),
    annotation_ranges = pieces,
    reference_dir = reference_dir,
    regression_SNPs = regression_SNP_tibble$SNP
  )
  L2_matrix <- do.call(rbind, purrr::map(chromosome_records, "L2_matrix"))
  list(
    L2_matrix = L2_matrix[regression_SNP_tibble$SNP, , drop = FALSE],
    M_5_50 = purrr::reduce(purrr::map(chromosome_records, "M_5_50"), `+`),
    sum_sq_5_50 = purrr::reduce(purrr::map(chromosome_records, "sum_sq_5_50"), `+`)
  )
}

#' Fit S-LDSC models that add candidate annotations to shared annotations
#'
#' Reproduce ldsc's partitioned-heritability regression (`ldsc.py --h2` with
#' several annotations): weights from its aggregate heritability guess without
#' reweighting iterations, a free intercept, chi-square < max(80, 0.001 max N),
#' and a 200-block jackknife over SNPs in reference order. Unlike ldsc, the
#' weights use the total LD score of fixed weight annotations, so that all
#' models of a GWAS share them and one set of block cross-products.
#'
#' Each candidate's robustness is measured by the blocks whose deletion lowers
#' its coefficient most: the share of the coefficient that deleting the 1, 5 or
#' 10 most influential blocks together removes, and the coefficient's z-score
#' without them, from a jackknife over the remaining blocks.
#'
#' @param chi2,N Chi-square statistic and sample size per regression SNP.
#' @param weight_L2 Regression-weight LD score per SNP.
#' @param L2_matrix,M LD scores and common-SNP sums of the annotations.
#' @param weight_annotations Annotations whose total LD score sets the
#'   regression weights.
#' @param shared_annotations Annotations in every model.
#' @param candidate_groups List of the annotations each model adds. A model
#'   with an annotation without common SNPs, such as the genes of a cell type
#'   that expresses none, cannot be fitted and has no rows.
#' @param n_blocks Jackknife blocks.
#' @param n_influential_blocks Numbers of most influential blocks to delete.
#' @return One row per candidate with its coefficient, its delete-one-block
#'   jackknife estimates and its robustness, and the model's total heritability
#'   and intercept, each with a jackknife standard error.
#' @keywords internal

fit_SLDSC_candidate_regressions <- function(chi2, N, weight_L2, L2_matrix, M, weight_annotations,
  shared_annotations, candidate_groups, n_blocks = 200L, n_influential_blocks = c(1L, 5L, 10L)) {
  kept <- chi2 < max(0.001 * max(N), 80)
  chi2 <- chi2[kept]
  N <- N[kept]
  weight_L2_total <- rowSums(L2_matrix[kept, weight_annotations, drop = FALSE])
  weight_M <- sum(M[weight_annotations])
  h2_guess <- min(max(weight_M * (mean(chi2) - 1) / mean(weight_L2_total * N), 0), 1)
  sqrt_weights <- sqrt(1 / (2 * (1 + h2_guess * N / weight_M * pmax(weight_L2_total, 1))^2) / pmax(weight_L2[kept], 1))
  N_bar <- mean(N)
  # Each model family gets its own cross-products: BLAS sums the same products
  # in another order when the design matrix has other columns.
  annotations <- c(shared_annotations, unique(unlist(candidate_groups)))
  design_matrix <- cbind(N * L2_matrix[kept, annotations, drop = FALSE] / N_bar, 1, chi2) * sqrt_weights
  dimnames(design_matrix) <- NULL
  separators <- floor(seq(0, nrow(design_matrix), length.out = n_blocks + 1L))
  block_crossproducts <- vapply(
    seq_len(n_blocks),
    \(b) crossprod(design_matrix[(separators[[b]] + 1L):separators[[b + 1L]], , drop = FALSE]),
    matrix(0, ncol(design_matrix), ncol(design_matrix))
  )
  intercept_col <- ncol(design_matrix) - 1L
  purrr::map_dfr(candidate_groups, \(candidates) {
    if (any(M[candidates] == 0)) return(NULL)
    model_annotations <- c(shared_annotations, candidates)
    model_cols <- c(match(model_annotations, annotations), intercept_col)
    XtX_blocks <- block_crossproducts[model_cols, model_cols, , drop = FALSE]
    Xty_blocks <- block_crossproducts[model_cols, ncol(design_matrix), ]
    fit <- jackknife_LDSC_blocks(XtX_blocks, Xty_blocks)
    annotation_idx <- seq_along(model_annotations)
    candidate_idx <- length(shared_annotations) + seq_along(candidates)
    model_M <- M[model_annotations]
    # Robustness is needed only where q < 0.05 is possible.
    robustness_matrix <- vapply(candidate_idx, \(i) {
      if (fit$estimate[[i]] / sqrt(fit$cov[i, i]) < stats::qnorm(0.95)) {
        return(rep(NA_real_, 2L * length(n_influential_blocks)))
      }
      influential_blocks <- order(fit$estimate[[i]] - fit$delete_values[i, ], decreasing = TRUE)
      reduced_estimates <- vapply(n_influential_blocks, \(n) {
        deleted <- influential_blocks[seq_len(n)]
        reduced_fit <- jackknife_LDSC_blocks(XtX_blocks[, , -deleted, drop = FALSE], Xty_blocks[, -deleted, drop = FALSE])
        c(reduced_fit$estimate[[i]], sqrt(reduced_fit$cov[i, i]))
      }, numeric(2))
      c(1 - reduced_estimates[1, ] / fit$estimate[[i]], reduced_estimates[1, ] / reduced_estimates[2, ])
    }, rlang::set_names(numeric(2L * length(n_influential_blocks)), c(
      paste0("share_in_top_", n_influential_blocks, "_blocks"),
      paste0("z_without_top_", n_influential_blocks, "_blocks")
    )))
    tibble::tibble(
      annotation = candidates,
      tau = fit$estimate[candidate_idx] / N_bar,
      tau_se = sqrt(diag(fit$cov)[candidate_idx]) / N_bar,
      tau_delete_values = purrr::map(candidate_idx, \(i) fit$delete_values[i, ] / N_bar),
      h2 = sum(model_M * fit$estimate[annotation_idx]) / N_bar,
      h2_se = sqrt(drop(model_M %*% fit$cov[annotation_idx, annotation_idx] %*% model_M)) / N_bar,
      intercept = fit$estimate[[length(model_cols)]],
      intercept_se = sqrt(fit$cov[length(model_cols), length(model_cols)]),
      n_SNPs = nrow(design_matrix),
      tibble::as_tibble(t(robustness_matrix))
    )
  })
}

#' Sum S-LDSC regression cross-products over jackknife blocks
#'
#' Weight the regression as ldsc's `--h2` does: weights from its aggregate
#' heritability guess without reweighting iterations, after removing SNPs with
#' chi-square >= max(80, 0.001 max N). Unlike ldsc, the weights use the total
#' LD score of fixed weight annotations, so that all models of a GWAS can share
#' them. Blocks are contiguous SNPs in reference order.
#'
#' @inheritParams fit_SLDSC_candidate_regressions
#' @param annotations Annotations of the design matrix, followed there by the
#'   intercept and the chi-square statistic.
#' @return A list with the per-block cross-product array of the design matrix,
#'   the mean sample size and the number of SNPs kept.
#' @keywords internal

get_SLDSC_block_crossproducts <- function(chi2, N, weight_L2, L2_matrix, M, weight_annotations, annotations, n_blocks) {
  kept <- chi2 < max(0.001 * max(N), 80)
  chi2 <- chi2[kept]
  N <- N[kept]
  weight_L2_total <- rowSums(L2_matrix[kept, weight_annotations, drop = FALSE])
  weight_M <- sum(M[weight_annotations])
  h2_guess <- min(max(weight_M * (mean(chi2) - 1) / mean(weight_L2_total * N), 0), 1)
  sqrt_weights <- sqrt(1 / (2 * (1 + h2_guess * N / weight_M * pmax(weight_L2_total, 1))^2) / pmax(weight_L2[kept], 1))
  N_bar <- mean(N)
  design_matrix <- cbind(N * L2_matrix[kept, annotations, drop = FALSE] / N_bar, 1, chi2) * sqrt_weights
  dimnames(design_matrix) <- NULL
  separators <- floor(seq(0, nrow(design_matrix), length.out = n_blocks + 1L))
  list(
    block_crossproducts = vapply(
      seq_len(n_blocks),
      \(b) crossprod(design_matrix[(separators[[b]] + 1L):separators[[b + 1L]], , drop = FALSE]),
      matrix(0, ncol(design_matrix), ncol(design_matrix))
    ),
    N_bar = N_bar,
    n_SNPs = nrow(design_matrix)
  )
}

#' Solve block cross-products with ldsc's delete-one-block jackknife
#'
#' @param XtX_blocks Array of per-block X'X matrices.
#' @param Xty_blocks Matrix of per-block X'y columns.
#' @return The full-data estimate, the jackknife covariance and the estimates
#'   without each block.
#' @keywords internal

jackknife_LDSC_blocks <- function(XtX_blocks, Xty_blocks) {
  n_blocks <- dim(XtX_blocks)[[3]]
  XtX <- rowSums(XtX_blocks, dims = 2)
  Xty <- rowSums(Xty_blocks)
  estimate <- solve(XtX, Xty)
  delete_values <- vapply(
    seq_len(n_blocks),
    \(b) solve(XtX - XtX_blocks[, , b], Xty - Xty_blocks[, b]),
    numeric(length(Xty))
  )
  pseudovalues <- n_blocks * estimate - (n_blocks - 1) * delete_values
  list(estimate = estimate, cov = stats::cov(t(pseudovalues)) / n_blocks, delete_values = delete_values)
}

#' Test one GWAS for heritability enrichment in each cell type's peaks and genes
#'
#' Fit three models per cell type, each with the baseline-LD model, LDSC-SEG
#' style: its peak annotation with all peaks, its gene annotation with all
#' genes, and both jointly with both controls. All models share the regression
#' weights of the baseline-LD model and all peaks.
#'
#' @param GWAS_input_record One normalized GWAS input record, whose study
#'   sample size and ancestry describe summary statistics without metadata.
#' @param sumstats_record The record of the same source from
#'   `get_SLDSC_sumstats_record()`.
#' @param baseline_record Record from `read_SLDSC_baseline_record()`.
#' @param peak_LD_score_record,gene_LD_score_record Records from
#'   `get_SLDSC_LD_score_record()` of the peak and gene annotations.
#' @return One row per model, tested annotation and cell type with the
#'   coefficient, its z-score and one-sided P-value, the standardized effect
#'   size tau* (missing when total heritability is not positive), the
#'   coefficient's delete-one-block jackknife estimates and its robustness at
#'   z >= 1.645, the model's total heritability and intercept, and the summary
#'   statistics' sample size and ancestry shares; estimates are missing without
#'   summary statistics.
#' @keywords internal

get_SLDSC_cell_type_tibble <- function(GWAS_input_record, sumstats_record, baseline_record, peak_LD_score_record, gene_LD_score_record) {
  if (!identical(sumstats_record$sourceId, GWAS_input_record$sourceId)) {
    stop("S-LDSC summary statistics of ", sumstats_record$sourceId, " were paired with ", GWAS_input_record$GWAS_ID, ".")
  }
  clusters <- setdiff(colnames(peak_LD_score_record$L2_matrix), "all_peaks")
  # Each model adds its controls and one cell type's tested annotations to the baseline-LD model.
  models <- list(
    peaks = list(controls = "all_peaks", tested = "peaks"),
    genes = list(controls = "all_genes", tested = "genes"),
    joint = list(controls = c("all_peaks", "all_genes"), tested = c("peaks", "genes"))
  )
  model_tibble <- purrr::imap_dfr(models, \(definition, name) tidyr::expand_grid(model = name, annotation = definition$tested, cluster = clusters))
  metadata_tibble <- sumstats_record$metadata_tibble %||%
    dplyr::select(GWAS_input_record$metadata_tibble, sample_size, dplyr::starts_with("ancestry_"))
  sumstats_tibble <- sumstats_record$sumstats_tibble
  if (is.null(sumstats_tibble)) {
    return(tibble::tibble(
      GWAS_ID = GWAS_input_record$GWAS_ID,
      model_tibble,
      tau = NA_real_, tau_se = NA_real_, tau_delete_values = list(NULL), z = NA_real_, p = NA_real_, tau_star = NA_real_,
      h2 = NA_real_, h2_se = NA_real_, intercept = NA_real_, intercept_se = NA_real_, n_SNPs = NA_integer_,
      sumstats_url = sumstats_record$sumstats_url, metadata_tibble
    ))
  }
  kept <- baseline_record$SNP_tibble$SNP %in% sumstats_tibble$SNP
  sumstats_idx <- match(baseline_record$SNP_tibble$SNP[kept], sumstats_tibble$SNP)
  # Cell-type annotations are named <modality>:<cluster>.
  LD_score_records <- purrr::imap(list(peaks = peak_LD_score_record, genes = gene_LD_score_record), \(record, modality) {
    columns <- c(paste0("all_", modality), clusters)
    annotations <- c(paste0("all_", modality), paste0(modality, ":", clusters))
    L2_matrix <- record$L2_matrix[kept, columns, drop = FALSE]
    colnames(L2_matrix) <- annotations
    list(
      L2_matrix = L2_matrix,
      M_5_50 = rlang::set_names(record$M_5_50[columns], annotations),
      sum_sq_5_50 = rlang::set_names(record$sum_sq_5_50[columns], annotations)
    )
  })
  L2_matrix <- do.call(cbind, c(list(baseline_record$L2_matrix[kept, ]), purrr::map(LD_score_records, "L2_matrix")))
  M <- c(baseline_record$M_5_50, unlist(unname(purrr::map(LD_score_records, "M_5_50"))))
  sum_sq <- unlist(unname(purrr::map(LD_score_records, "sum_sq_5_50")))
  M_common <- baseline_record$M_5_50[["baseL2"]]
  annotation_sd <- sqrt(sum_sq / M_common - (M[names(sum_sq)] / M_common)^2)
  baseline_annotations <- colnames(baseline_record$L2_matrix)
  purrr::imap_dfr(models, \(model, name) {
    fit_SLDSC_candidate_regressions(
      chi2 = sumstats_tibble$chi2[sumstats_idx],
      N = sumstats_tibble$N[sumstats_idx],
      weight_L2 = baseline_record$SNP_tibble$weight_L2[kept],
      L2_matrix = L2_matrix,
      M = M,
      weight_annotations = c(baseline_annotations, "all_peaks"),
      shared_annotations = c(baseline_annotations, model$controls),
      candidate_groups = purrr::map(clusters, \(cluster) paste0(model$tested, ":", cluster))
    ) |>
      dplyr::mutate(model = name)
  }) |>
    tidyr::separate_wider_delim(annotation, ":", names = c("annotation", "cluster")) |>
    dplyr::mutate(
      GWAS_ID = GWAS_input_record$GWAS_ID,
      z = tau / tau_se,
      p = stats::pnorm(z, lower.tail = FALSE),
      tau_star = dplyr::if_else(h2 > 0, tau * unname(annotation_sd[paste0(annotation, ":", cluster)]) * M_common / h2, NA_real_),
      sumstats_url = sumstats_record$sumstats_url,
      metadata_tibble
    ) |>
    dplyr::relocate(GWAS_ID, model, annotation, cluster, tau, tau_se, z, p, tau_star)
}

#' Plot S-LDSC cell-type enrichment for all GWAS
#'
#' @param SLDSC_tibble Combined cell-type S-LDSC results of one model and
#'   annotation, with q-values and top sets.
#' @param GWAS_inputs_tibble Metadata of the same GWAS, whose sample size and
#'   ancestry are replaced by those of the S-LDSC results.
#' @param compartments_patterns Optional cell-type compartment patterns.
#' @return Annotated heatmap.
#' @keywords internal

plot_SLDSC_cell_type_heatmap <- function(SLDSC_tibble, GWAS_inputs_tibble, compartments_patterns = NULL) {
  list_GWAS <- \(rows) paste(sort(unique(SLDSC_tibble$GWAS_ID[rows])), collapse = ", ")
  without_sumstats <- is.na(SLDSC_tibble$z)
  low_h2 <- !without_sumstats & SLDSC_tibble$h2 / SLDSC_tibble$h2_se < 7
  # A few extreme tau* would otherwise wash out the colours of all others.
  colour_limit <- signif(stats::quantile(abs(SLDSC_tibble$tau_star), 0.99, na.rm = TRUE, names = FALSE), 2)
  SLDSC_tibble <- SLDSC_tibble |>
    dplyr::mutate(
      tau_star_colour = pmax(pmin(tau_star, colour_limit), -colour_limit),
      significance_label = get_significance_stars(p, q)
    )
  GWAS_metadata_tracks_plot <- GWAS_inputs_tibble |>
    dplyr::rows_update(
      dplyr::distinct(dplyr::select(SLDSC_tibble, GWAS_ID, sample_size, dplyr::starts_with("ancestry_"))),
      by = "GWAS_ID"
    ) |>
    dplyr::arrange(Category, GWAS_ID) |>
    dplyr::mutate(GWAS_ID = factor(GWAS_ID, levels = rev(GWAS_ID))) |>
    plot_GWAS_metadata_tracks()
  plot_GWAS_by_cluster_heatmap(
    SLDSC_tibble,
    GWAS_metadata_tracks_plot = GWAS_metadata_tracks_plot,
    compartments_patterns = compartments_patterns,
    fill_col = "tau_star_colour",
    fill_label = "Standardized effect size tau*",
    support_label_col = "significance_label",
    top_set_col = "in_top_set"
  ) + patchwork::plot_annotation(
    title = "Heritability enrichment in cell-type-specific accessibility (S-LDSC)",
    subtitle = paste(
      "Red marks cell types whose specifically accessible peaks carry more per-SNP heritability than peaks in general.",
      "Stars mark nominal (*) and FDR-significant (**) enrichment; a black outline marks the top significant cell type and grey outlines the cell types not significantly below it; together they contain the most relevant cell type.",
      sep = "\n"
    ),
    caption = stringr::str_wrap(paste(
      "Each SNP in a consensus peak carries the peak's share of depth-normalized pseudobulk accessibility in the cell type.",
      "Each cell type is fitted jointly with an all-peaks annotation and the 97 baseline-LD v2.2 annotations to HapMap3 SNP summary statistics from the GWAS Catalog or a configured release, using GRCh38 1000 Genomes European LD.",
      sprintf("Colour: tau*, the change in per-SNP heritability per standard deviation of the annotation, relative to the mean per-SNP heritability; values beyond ±%s take the end colours.", colour_limit),
      "Stars: * one-sided coefficient P < 0.05 from its 200-block jackknife; ** also Benjamini-Hochberg q < 0.05 over all GWAS and cell types.",
      "Outlines: in GWAS with a cell type at q < 0.05, black marks the one of these with the largest tau* and grey every cell type whose shortfall from it is not significant at one-sided P < 0.05 divided by the number of other cell types",
      "(Bonferroni, so the set holds the cell type with the largest true tau* with at least 95% probability); shortfall standard errors come from paired jackknife differences over the same blocks.",
      if (any(without_sumstats)) paste0("Grey rows have no usable genome-wide summary statistics: ", list_GWAS(without_sumstats), "."),
      if (any(low_h2)) paste0("Heritability z < 7, so partitioning is less reliable: ", list_GWAS(low_h2), "."),
      "Sample and ancestry tracks describe the summary statistics, from their GWAS-SSF metadata or else the credible-set study; method and loci tracks describe the credible sets.",
      "GWAS with non-European samples are approximated with European LD."
    ), width = 150),
    theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0))
  )
}

#' Pair S-LDSC and pseudobulk chromVAR cell-type z-scores and P-values
#'
#' @param SLDSC_tibble Combined cell-type S-LDSC results.
#' @param chromVAR_deviation_tibble Cell-type pseudobulk chromVAR results.
#' @return One row per GWAS with summary statistics and cell type.
#' @keywords internal

pair_SLDSC_chromVAR_z <- function(SLDSC_tibble, chromVAR_deviation_tibble) {
  SLDSC_tibble |>
    dplyr::filter(!is.na(z)) |>
    dplyr::select(GWAS_ID, Category, cluster, SLDSC_z = z, SLDSC_p = p, h2, h2_se, sample_size, ancestry_EUR, n_cells) |>
    dplyr::inner_join(
      dplyr::select(chromVAR_deviation_tibble, GWAS_ID, cluster, chromVAR_z = z, chromVAR_p = z_p),
      by = c("GWAS_ID", "cluster")
    )
}

#' Compare S-LDSC with pseudobulk chromVAR cell-type enrichment
#'
#' @inheritParams pair_SLDSC_chromVAR_z
#' @return Faceted scatter plot with one panel per GWAS with summary statistics.
#' @keywords internal

plot_SLDSC_chromVAR_comparison <- function(SLDSC_tibble, chromVAR_deviation_tibble) {
  comparison_tibble <- pair_SLDSC_chromVAR_z(SLDSC_tibble, chromVAR_deviation_tibble) |>
    dplyr::mutate(
      panel = sprintf(
        "%s\nρ = %.2f; h² z = %.0f",
        GWAS_ID,
        stats::cor(SLDSC_z, chromVAR_z, method = "spearman"),
        dplyr::first(h2 / h2_se)
      ),
      .by = GWAS_ID
    ) |>
    dplyr::arrange(Category, GWAS_ID) |>
    dplyr::mutate(panel = forcats::fct_inorder(panel))
  if (nrow(comparison_tibble) == 0L) {
    return(ggplot2::ggplot() + ggplot2::theme_void() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No GWAS with genome-wide summary statistics"))
  }

  comparison_tibble |>
    ggplot2::ggplot(ggplot2::aes(chromVAR_z, SLDSC_z, colour = cluster)) +
    ggplot2::geom_hline(yintercept = stats::qnorm(0.95), linetype = "dashed", colour = "grey60") +
    ggplot2::geom_point(size = 1.2) +
    ggrepel::geom_text_repel(
      data = \(data) dplyr::filter(data, pmin(SLDSC_p, chromVAR_p) < 0.05),
      ggplot2::aes(label = cluster),
      size = 2.2,
      max.overlaps = 20,
      show.legend = FALSE
    ) +
    ggplot2::facet_wrap(ggplot2::vars(panel), ncol = 5) +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 4, override.aes = list(size = 2))) +
    ggplot2::labs(
      title = "Cell-type enrichment: S-LDSC versus pseudobulk chromVAR",
      subtitle = paste(
        "Cell types supported by both genome-wide heritability and fine-mapped variant accessibility fall in the upper right;",
        "labelled cell types elsewhere are disagreements to investigate. Shared axes compare signal strength between GWAS.",
        sep = "\n"
      ),
      caption = stringr::str_wrap(paste(
        "x: pseudobulk chromVAR z-score: the background z-score of fine-mapped posterior-probability-weighted peaks, standardized by those of 10,000 copies of the GWAS with its credible sets moved 0.5-2 Mb along their chromosomes.",
        "y: S-LDSC coefficient z-score of the cell type's peak-specificity annotation, conditional on all peaks and the baseline-LD model.",
        "Cell types at unadjusted one-sided P < 0.05 in either method are labelled; the dashed line marks this level for S-LDSC (z = 1.645).",
        "chromVAR P-values are the share of moved copies that score at least as high; their null is skewed, so no single chromVAR z-score marks P = 0.05.",
        "Panel titles give the Spearman correlation across cell types and the S-LDSC total-heritability z-score.",
        "Only GWAS with genome-wide summary statistics are shown."
      ), width = 150),
      x = "Pseudobulk chromVAR z-score",
      y = "S-LDSC coefficient z-score",
      colour = "Cell type"
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      legend.position = "bottom",
      strip.text = ggplot2::element_text(size = 7),
      plot.caption = ggplot2::element_text(hjust = 0),
      plot.caption.position = "plot"
    )
}

#' Plot S-LDSC and pseudobulk chromVAR agreement across GWAS and cell types
#'
#' @inheritParams pair_SLDSC_chromVAR_z
#' @param GWAS_inputs_tibble GWAS metadata with fine-mapping method and
#'   credible-set loci.
#' @return Patchwork of per-GWAS and per-cell-type agreement against predictors.
#' @keywords internal

plot_SLDSC_chromVAR_agreement <- function(SLDSC_tibble, chromVAR_deviation_tibble, GWAS_inputs_tibble) {
  comparison_tibble <- pair_SLDSC_chromVAR_z(SLDSC_tibble, chromVAR_deviation_tibble)
  if (nrow(comparison_tibble) == 0L) {
    return(ggplot2::ggplot() + ggplot2::theme_void() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No GWAS with genome-wide summary statistics"))
  }
  n_cell_types <- dplyr::n_distinct(comparison_tibble$cluster)
  # Fieller's approximate standard error of Fisher-transformed Spearman correlations.
  rho_margin <- stats::qnorm(0.975) * sqrt(1.06 / (n_cell_types - 3))
  GWAS_tibble <- comparison_tibble |>
    dplyr::summarise(
      rho = stats::cor(SLDSC_z, chromVAR_z, method = "spearman"),
      chromVAR_spread = stats::sd(chromVAR_z),
      SLDSC_spread = stats::sd(SLDSC_z),
      h2_z = stats::median(h2 / h2_se),
      sample_size = dplyr::first(sample_size),
      ancestry_EUR = dplyr::first(ancestry_EUR),
      supported_in = dplyr::case_when(
        any(SLDSC_p < 0.05) & any(chromVAR_p < 0.05) ~ "Both methods",
        any(SLDSC_p < 0.05) ~ "S-LDSC only",
        any(chromVAR_p < 0.05) ~ "chromVAR only",
        .default = "Neither"
      ),
      .by = c(GWAS_ID, Category)
    ) |>
    dplyr::left_join(
      dplyr::select(GWAS_inputs_tibble, GWAS_ID, finemappingMethod, n_credible_set_loci),
      by = "GWAS_ID", relationship = "one-to-one"
    ) |>
    dplyr::arrange(Category, GWAS_ID)
  # Standardizing within GWAS removes power differences before comparing GWAS per cell type.
  cell_type_tibble <- comparison_tibble |>
    dplyr::mutate(dplyr::across(c(SLDSC_z, chromVAR_z), \(z) (z - mean(z)) / stats::sd(z)), .by = GWAS_ID) |>
    dplyr::summarise(
      rho = stats::cor(SLDSC_z, chromVAR_z, method = "spearman"),
      chromVAR_spread = stats::sd(chromVAR_z),
      SLDSC_spread = stats::sd(SLDSC_z),
      log10_nuclei = log10(dplyr::first(n_cells)),
      .by = cluster
    )
  # One facet per predictor, titled with its Spearman correlation with agreement.
  pivot_predictors <- \(tibble, predictors) {
    tibble |>
      tidyr::pivot_longer(dplyr::all_of(unname(predictors)), names_to = "predictor", values_to = "value") |>
      dplyr::filter(!is.na(value)) |>
      dplyr::mutate(predictor = factor(predictor, predictors, names(predictors))) |>
      dplyr::mutate(
        strip = sprintf(
          "%s\nρ = %.2f, P = %.2g (n = %d)",
          predictor,
          stats::cor(value, rho, method = "spearman"),
          stats::cor.test(value, rho, method = "spearman", exact = FALSE)$p.value,
          dplyr::n()
        ),
        .by = predictor
      ) |>
      dplyr::mutate(strip = forcats::fct_reorder(strip, as.integer(predictor)))
  }
  method_colors <- c("SuSie" = "#238B45", "SuSiE-inf" = "#41B6C4", "PICS" = "#F16913")
  method_colors <- c(method_colors, make_named_heatmap_palette(setdiff(GWAS_tibble$finemappingMethod, names(method_colors)), palette = "Dark2"))
  GWAS_scales <- list(
    ggplot2::scale_colour_manual(values = method_colors, name = "Fine-mapping"),
    ggplot2::scale_shape_manual(values = c("Both methods" = 16, "S-LDSC only" = 2, "chromVAR only" = 6, "Neither" = 4),
      name = "Cell type at P < 0.05 in")
  )
  zero_line <- ggplot2::geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3)
  base_theme <- ggplot2::theme_bw(base_size = 9)

  GWAS_plot <- GWAS_tibble |>
    dplyr::mutate(GWAS_ID = forcats::fct_rev(forcats::fct_inorder(GWAS_ID)), Category = forcats::fct_inorder(Category)) |>
    ggplot2::ggplot(ggplot2::aes(rho, GWAS_ID, colour = finemappingMethod, shape = supported_in)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
    ggplot2::geom_linerange(ggplot2::aes(xmin = tanh(atanh(rho) - rho_margin), xmax = tanh(atanh(rho) + rho_margin)), show.legend = FALSE) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_grid(ggplot2::vars(Category), scales = "free_y", space = "free_y", switch = "y",
      labeller = ggplot2::as_labeller(\(x) stringr::str_replace_all(x, "_", " "))) +
    ggplot2::scale_y_discrete(position = "right") +
    GWAS_scales +
    ggplot2::labs(subtitle = "Agreement per GWAS", x = "Spearman ρ across cell types", y = NULL) +
    base_theme +
    ggplot2::theme(strip.text.y.left = ggplot2::element_text(angle = 0), strip.placement = "outside",
      panel.spacing.y = ggplot2::unit(1, "pt"))

  support_levels <- c("P >= 0.05", "P < 0.05")
  support_plot <- comparison_tibble |>
    dplyr::count(
      chromVAR = factor(support_levels[(chromVAR_p < 0.05) + 1L], support_levels),
      SLDSC = factor(support_levels[(SLDSC_p < 0.05) + 1L], support_levels)
    ) |>
    ggplot2::ggplot(ggplot2::aes(chromVAR, SLDSC)) +
    ggplot2::geom_tile(ggplot2::aes(fill = n), colour = "white", show.legend = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = n), size = 3.5) +
    ggplot2::scale_fill_gradient(low = "grey95", high = "grey60", transform = "log10") +
    ggplot2::labs(subtitle = "GWAS-cell-type pairs", x = "chromVAR", y = "S-LDSC") +
    base_theme +
    ggplot2::theme(panel.grid = ggplot2::element_blank())

  GWAS_predictor_plot <- GWAS_tibble |>
    dplyr::mutate(log10_sample_size = log10(sample_size), log10_loci = log10(n_credible_set_loci)) |>
    pivot_predictors(c(
      "chromVAR z spread across cell types (SD)" = "chromVAR_spread",
      "S-LDSC z spread across cell types (SD)" = "SLDSC_spread",
      "S-LDSC total heritability z-score" = "h2_z",
      "Credible-set loci (log10)" = "log10_loci",
      "Summary-statistics sample size (log10)" = "log10_sample_size",
      "Summary-statistics European share" = "ancestry_EUR"
    )) |>
    ggplot2::ggplot(ggplot2::aes(value, rho, colour = finemappingMethod, shape = supported_in)) +
    zero_line +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(ggplot2::vars(strip), scales = "free_x", ncol = 3) +
    GWAS_scales +
    ggplot2::labs(subtitle = "Predictors of per-GWAS agreement", x = NULL, y = "Spearman ρ across cell types") +
    base_theme

  cell_type_plot <- cell_type_tibble |>
    pivot_predictors(c(
      "chromVAR standardized z spread across GWAS (SD)" = "chromVAR_spread",
      "S-LDSC standardized z spread across GWAS (SD)" = "SLDSC_spread",
      "Nuclei (log10)" = "log10_nuclei"
    )) |>
    ggplot2::ggplot(ggplot2::aes(value, rho, label = cluster)) +
    zero_line +
    ggplot2::geom_point(size = 1.2, colour = "grey30") +
    ggrepel::geom_text_repel(size = 2, colour = "grey30", max.overlaps = Inf, segment.size = 0.2, seed = 1) +
    ggplot2::facet_wrap(ggplot2::vars(strip), scales = "free_x", ncol = 3) +
    ggplot2::labs(subtitle = "Agreement per cell type across GWAS", x = NULL, y = "Spearman ρ across GWAS") +
    base_theme

  patchwork::wrap_plots(
    A = GWAS_plot, B = support_plot, C = GWAS_predictor_plot, D = cell_type_plot,
    design = "AC\nAC\nAD\nBD", widths = c(1, 2), heights = c(2.6, 2.6, 1.6, 2.4), guides = "collect"
  ) +
    patchwork::plot_annotation(
      title = "Agreement between S-LDSC and pseudobulk chromVAR cell-type enrichment",
      subtitle = paste(
        "Check which GWAS and cell types the two methods rank alike, and whether agreement follows each method's cell-type signal,",
        "GWAS power, fine-mapping, ancestry, trait category or cell-type size. Agreement is limited where one method barely separates cell types.",
        sep = "\n"
      ),
      caption = stringr::str_wrap(paste(
        sprintf("A: Spearman correlation of S-LDSC coefficient z-scores and pseudobulk chromVAR z-scores across %d cell types;", n_cell_types),
        "bars are approximate 95% intervals (Fisher transformation) that treat cell types as independent.",
        "B: GWAS-cell-type pairs by whether each method reaches unadjusted one-sided P < 0.05; chromVAR P-values come from 10,000 copies of the GWAS with its credible sets moved 0.5-2 Mb along their chromosomes.",
        "C: Spread is the standard deviation of a method's z-scores across cell types; heritability z-score is S-LDSC total SNP heritability over its jackknife standard error.",
        "Sample size and European share describe the summary statistics, from their GWAS-SSF metadata or else the credible-set study.",
        "D: Spearman correlation across GWAS after standardizing each method's z-scores across cell types within each GWAS, which removes differences in GWAS power;",
        "spread is the standard deviation of these standardized z-scores across GWAS.",
        "Strips in C and D give each predictor's Spearman correlation with agreement and its unadjusted P; GWAS share cohorts and traits, so P-values are descriptive.",
        sprintf("Only the %d GWAS with genome-wide summary statistics are shown.", nrow(GWAS_tibble))
      ), width = 200),
      tag_levels = "A",
      theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0))
    ) &
    ggplot2::theme(legend.position = "bottom")
}

#' Compare S-LDSC cell-type enrichment of accessibility, expression and both
#'
#' @param SLDSC_tibble Combined cell-type S-LDSC results of all models.
#' @return Patchwork of the enriched GWAS-cell-type pairs by the modality that
#'   explains them, and of the single-modality and joint z-scores.
#' @keywords internal

plot_SLDSC_modality_comparison <- function(SLDSC_tibble) {
  if (!any(SLDSC_tibble$q < 0.05 & SLDSC_tibble$model != "joint", na.rm = TRUE)) {
    return(ggplot2::ggplot() + ggplot2::theme_void() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No GWAS-cell-type pair reaches q < 0.05"))
  }
  pair_tibble <- SLDSC_tibble |>
    dplyr::filter(!is.na(z)) |>
    tidyr::pivot_wider(
      id_cols = c(GWAS_ID, Category, cluster),
      names_from = c(model, annotation),
      values_from = c(z, q)
    ) |>
    dplyr::mutate(
      single = dplyr::case_when(
        q_peaks_peaks < 0.05 & q_genes_genes < 0.05 ~ "Both",
        q_peaks_peaks < 0.05 ~ "Accessibility only",
        q_genes_genes < 0.05 ~ "Expression only"
      ),
      joint = dplyr::case_when(
        is.na(single) ~ NA_character_,
        q_joint_peaks < 0.05 & q_joint_genes < 0.05 ~ "Both add",
        q_joint_peaks < 0.05 ~ "Accessibility adds",
        q_joint_genes < 0.05 ~ "Expression adds",
        .default = "Neither adds"
      ) |>
        factor(c("Accessibility adds", "Expression adds", "Both add", "Neither adds")),
      single = factor(single, c("Both", "Accessibility only", "Expression only"))
    )
  enriched_tibble <- dplyr::filter(pair_tibble, !is.na(single))
  joint_colors <- c("Accessibility adds" = "#D95F02", "Expression adds" = "#1B9E77", "Both add" = "#7570B3", "Neither adds" = "grey60")
  joint_scale <- ggplot2::scale_colour_manual(values = joint_colors, name = "Joint model: q < 0.05 for", drop = FALSE,
    guide = ggplot2::guide_legend(override.aes = list(size = 3, shape = 16)))
  base_theme <- ggplot2::theme_bw(base_size = 9)
  z_threshold <- \(model_annotation) min(pair_tibble[[paste0("z_", model_annotation)]][pair_tibble[[paste0("q_", model_annotation)]] < 0.05], Inf)

  matrix_plot <- enriched_tibble |>
    dplyr::arrange(Category, GWAS_ID) |>
    dplyr::mutate(
      GWAS_ID = forcats::fct_rev(forcats::fct_inorder(GWAS_ID)),
      cluster = forcats::fct_reorder(cluster, as.integer(GWAS_ID), .fun = max, .desc = TRUE)
    ) |>
    ggplot2::ggplot(ggplot2::aes(cluster, GWAS_ID, colour = joint, shape = single)) +
    ggplot2::geom_point(ggplot2::aes(size = pmax(z_peaks_peaks, z_genes_genes)), stroke = 0.8) +
    ggplot2::facet_grid(ggplot2::vars(Category), scales = "free_y", space = "free_y", switch = "y",
      labeller = ggplot2::as_labeller(\(x) stringr::str_replace_all(x, "_", " "))) +
    joint_scale +
    ggplot2::scale_shape_manual(values = c("Both" = 16, "Accessibility only" = 15, "Expression only" = 17),
      name = "Single-modality model: q < 0.05 for", drop = FALSE, guide = ggplot2::guide_legend(override.aes = list(size = 3))) +
    ggplot2::scale_size_area(max_size = 4, breaks = c(3, 5, 7), name = "Larger single-modality z") +
    ggplot2::scale_y_discrete(position = "right") +
    ggplot2::labs(subtitle = "Enriched GWAS-cell-type pairs", x = NULL, y = NULL) +
    base_theme +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.text.y.left = ggplot2::element_text(angle = 0),
      strip.placement = "outside",
      panel.spacing.y = ggplot2::unit(1, "pt")
    )

  scatter_plot <- \(x, y, x_label, y_label, subtitle, diagonal) {
    ggplot2::ggplot(pair_tibble, ggplot2::aes(.data[[paste0("z_", x)]], .data[[paste0("z_", y)]])) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
      (if (diagonal) ggplot2::geom_abline(colour = "grey60", linewidth = 0.3)) +
      ggplot2::geom_hline(yintercept = z_threshold(y), linetype = "dashed", colour = "grey60") +
      ggplot2::geom_vline(xintercept = z_threshold(x), linetype = "dashed", colour = "grey60") +
      ggplot2::geom_point(data = \(data) dplyr::filter(data, is.na(joint)), colour = "grey85", size = 0.6) +
      ggplot2::geom_point(data = \(data) dplyr::filter(data, !is.na(joint)), ggplot2::aes(colour = joint, shape = single), size = 1.4) +
      ggplot2::scale_colour_manual(values = joint_colors, guide = "none") +
      ggplot2::scale_shape_manual(values = c("Both" = 16, "Accessibility only" = 15, "Expression only" = 17), guide = "none") +
      ggplot2::labs(subtitle = subtitle, x = x_label, y = y_label) +
      base_theme
  }
  single_plot <- scatter_plot("peaks_peaks", "genes_genes", "Accessibility z (peak model)", "Expression z (gene model)",
    "Single-modality models", diagonal = TRUE)
  joint_plot <- scatter_plot("joint_peaks", "joint_genes", "Accessibility z beyond expression", "Expression z beyond accessibility",
    "Joint model", diagonal = FALSE)

  count_plot <- enriched_tibble |>
    dplyr::count(Category, joint) |>
    dplyr::mutate(Category = stringr::str_replace_all(Category, "_", " ")) |>
    ggplot2::ggplot(ggplot2::aes(n, forcats::fct_rev(Category), fill = joint)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::scale_fill_manual(values = joint_colors, guide = "none") +
    ggplot2::labs(subtitle = "Enriched pairs per trait category", x = "GWAS-cell-type pairs", y = NULL) +
    base_theme

  patchwork::wrap_plots(
    A = matrix_plot, B = single_plot, C = joint_plot, D = count_plot,
    design = "AB\nAC\nAD", widths = c(2.2, 1), guides = "collect"
  ) +
    patchwork::plot_annotation(
      title = "Heritability enrichment in cell-type-specific accessibility and expression (S-LDSC)",
      subtitle = paste(
        "Check which GWAS-cell-type enrichments specific peaks and specific genes each detect, and which modality still explains",
        "heritability when both are fitted together. Pairs where neither adds share their signal between the two annotations.",
        sep = "\n"
      ),
      caption = stringr::str_wrap(paste(
        "Peak annotation: each SNP in a consensus peak carries the peak's share of depth-normalized pseudobulk accessibility in the cell type.",
        "Gene annotation: SNPs within 100 kb of the bodies of the cell type's most specifically expressed genes, 10% of the autosomal protein-coding genes expressed in any cell type",
        "(at least 10 counts and 1 count per million in a cell type; specificity is a gene's share of depth-normalized pseudobulk expression across cell types).",
        "Single-modality models fit one cell-type annotation with all peaks or all genes; the joint model fits both with both controls; all include the 97 baseline-LD v2.2 annotations and share regression weights.",
        "z: coefficient over its 200-block jackknife standard error. q: Benjamini-Hochberg-adjusted one-sided P over GWAS and cell types within each model and annotation; dashed lines mark the smallest z with q < 0.05.",
        "Pairs are shown when a single-modality model reaches q < 0.05 (A, coloured points in B and C); others are grey.",
        "Point size in A: the larger single-modality z. GWAS with non-European samples are approximated with European LD."
      ), width = 180),
      tag_levels = "A",
      theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0))
    ) &
    ggplot2::theme(legend.position = "bottom", legend.box = "horizontal")
}

#' Plot how concentrated significant S-LDSC enrichments are in jackknife blocks
#'
#' @param SLDSC_tibble Combined cell-type S-LDSC results of all models, with
#'   robustness columns for coefficients at q < 0.05.
#' @param n_blocks Jackknife blocks of the fits.
#' @return Faceted scatter plot of the coefficient share carried by the most
#'   influential blocks against the coefficient z-score.
#' @keywords internal

plot_SLDSC_block_influence <- function(SLDSC_tibble, n_blocks = 200L) {
  if (!any(SLDSC_tibble$q < 0.05, na.rm = TRUE)) {
    return(ggplot2::ggplot() + ggplot2::theme_void() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No coefficient reaches q < 0.05"))
  }
  n_influential_blocks <- as.integer(stringr::str_extract(grep("^share_in_top_", names(SLDSC_tibble), value = TRUE), "[0-9]+"))
  influence_tibble <- SLDSC_tibble |>
    dplyr::filter(!is.na(share_in_top_1_blocks)) |>
    dplyr::mutate(term = factor(
      paste0(c(peaks = "Accessibility", genes = "Expression")[annotation], dplyr::if_else(model == "joint", " (joint)", "")),
      c("Accessibility", "Expression", "Accessibility (joint)", "Expression (joint)")
    )) |>
    tidyr::pivot_longer(
      dplyr::matches("^(share_in|z_without)_top_[0-9]+_blocks$"),
      names_to = c(".value", "n_blocks_deleted"),
      names_pattern = "^(share_in|z_without)_top_([0-9]+)_blocks$"
    ) |>
    dplyr::mutate(
      n_blocks_deleted = as.integer(n_blocks_deleted),
      robust = factor(dplyr::if_else(z_without >= stats::qnorm(0.95), "z >= 1.645", "z < 1.645"), c("z >= 1.645", "z < 1.645")),
      facet = sprintf("Top %d block%s deleted: %d of %d keep z >= 1.645", n_blocks_deleted,
        dplyr::if_else(n_blocks_deleted == 1L, "", "s"), sum(robust == "z >= 1.645"), dplyr::n()),
      .by = n_blocks_deleted
    ) |>
    dplyr::mutate(facet = forcats::fct_reorder(facet, n_blocks_deleted))
  # Blocks of equal expected contribution with normal noise: deleting the K
  # largest removes K / n of the coefficient plus the sum of the K largest of n
  # standard normal order statistics (Blom), scaled by the coefficient's z.
  z_range <- range(influence_tibble$z)
  even_tibble <- tidyr::expand_grid(n_blocks_deleted = n_influential_blocks, z = seq(z_range[[1]], z_range[[2]], length.out = 100)) |>
    dplyr::mutate(
      share_in = n_blocks_deleted / n_blocks + purrr::map_dbl(n_blocks_deleted, \(k) sum(stats::qnorm((n_blocks - seq_len(k) + 0.625) / (n_blocks + 0.25)))) /
        (sqrt(n_blocks) * z)
    ) |>
    dplyr::left_join(dplyr::distinct(influence_tibble, n_blocks_deleted, facet), by = "n_blocks_deleted")
  influence_tibble |>
    ggplot2::ggplot(ggplot2::aes(z, share_in)) +
    ggplot2::geom_hline(yintercept = 1, colour = "grey80", linewidth = 0.3) +
    ggplot2::geom_line(data = even_tibble, colour = "grey40", linetype = "dashed") +
    ggplot2::geom_point(ggplot2::aes(colour = term, shape = robust), size = 1.4, alpha = 0.8) +
    ggplot2::facet_wrap(ggplot2::vars(facet), nrow = 1) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_colour_manual(values = c("Accessibility" = "#D95F02", "Expression" = "#1B9E77",
      "Accessibility (joint)" = "#FDB863", "Expression (joint)" = "#80CDC1"), name = "Coefficient", drop = FALSE) +
    ggplot2::scale_shape_manual(values = c("z >= 1.645" = 16, "z < 1.645" = 4), name = "Without these blocks") +
    ggplot2::labs(
      title = "Concentration of significant S-LDSC enrichments in the genome",
      subtitle = paste(
        "Points far above the dashed curve rest on a few genome regions; points near it spread their enrichment widely.",
        "Crosses lose nominal support once those regions are left out.",
        sep = "\n"
      ),
      caption = stringr::str_wrap(paste(
        sprintf("Each point is a cell-type coefficient at q < 0.05 in one model. Blocks are the %d contiguous jackknife blocks of regression SNPs (about 15 Mb each).", n_blocks),
        "The most influential blocks are those whose single deletion lowers the coefficient most; they are deleted together and the model refitted.",
        "y: share of the coefficient they carry; x: the coefficient's full z-score; shape: z-score without them, from a jackknife over the remaining blocks.",
        "Dashed curve: expected share if every block contributed equally with normally distributed noise; heavier-tailed block noise also lifts points above it.",
        "SNPs with chi-square above 80, or above 0.001 times the largest sample size, are excluded from every fit, so the largest-effect loci never enter."
      ), width = 160),
      x = "Coefficient z-score (log scale)",
      y = "Share of the coefficient in the deleted blocks"
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(legend.position = "bottom", plot.caption = ggplot2::element_text(hjust = 0))
}
