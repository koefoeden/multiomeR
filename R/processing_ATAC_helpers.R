get_gene_annotation_GRanges_from_EnsDb <- function(ensdb) {
  genes <- ensdb |>
    ensembldb::genes() |>
    GenomeInfoDb::`seqlevelsStyle<-`("UCSC") |>
    GenomeInfoDb::keepStandardChromosomes(pruning.mode = "coarse")
  names(genes) <- genes$gene_name %||% genes$gene_id
  genes
}

#' Get blacklist GRanges
#'
#' Load genome-specific ENCODE/excluderanges blacklist intervals as GRanges.
#'
#' @param genome Genome build key used to choose chromosome sizes, blacklist resources, and external-tool parameters.
#' @param annot_hub_interface Optional pre-created `AnnotationHub::AnnotationHub()`
#'   object, useful for tests or shared cache handling. `GRCm39` is loaded from
#'   the repo BED resource instead.
#' @return GRanges of blacklist intervals in UCSC chromosome style.
#' @keywords internal

get_blacklist_GRanges <- function(genome, annot_hub_interface = NULL) {
  if (identical(genome, "GRCm39")) {
    return(
      readr::read_tsv(
        "resources/mm39.excluderanges.bed",
        col_names = c("seqnames", "start", "end", "width", "strand", "reason"),
        show_col_types = FALSE
      ) |>
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    )
  }

  title <- switch(
    genome,
    "GRCh38" = "hg38.Kundaje.GRCh38_unified_Excludable",
    "mm10" = "mm10.Boyle.mm10-Excludable.v2",
    stop("Invalid genome for blacklist selection.")
  )
  annot_hub_interface <- annot_hub_interface %||% AnnotationHub::AnnotationHub(ask = FALSE)
  hits <- AnnotationHub::query(annot_hub_interface, c("excluderanges", genome))
  hit_idx <- which(S4Vectors::mcols(hits)$title == title)
  if (length(hit_idx) != 1L) {
    stop("Could not resolve one AnnotationHub blacklist resource for ", genome, ".")
  }
  hits[[hit_idx]]
}

# Peaks ----

get_standard_chroms <- function(genome) {
  switch(
    genome,
    "GRCh38" = paste0("chr", c(1:22, "X", "Y")),
    "mm10" = paste0("chr", c(1:19, "X", "Y")),
    "GRCm39" = paste0("chr", c(1:19, "X", "Y")),
    stop("Invalid genome. Must be 'GRCh38', 'mm10', or 'GRCm39'.")
  )
}

get_effective_genome_size <- function(genome) {
  switch(
    genome,
    "GRCh38" = 2.913e9,
    "GRCm39" = 2.65e9,
    "mm10" = 2.65e9,
    stop("Invalid genome. Must be 'mm10','GRCh38' or 'GRCm39'.")
  )
}

get_chrom_sizes_for_BPCells_tile_calling <- function(genome) {
  chrom_sizes_file <- switch(
    genome,
    "mm10" = "resources/mm10.chrom.sizes",
    "GRCh38" = "resources/GRCh38.chrom.sizes",
    "GRCm39" = "resources/GRCm39.chrom.sizes",
    stop("Invalid genome. Must be 'mm10','GRCh38' or 'GRCm39'.")
  )

  readr::read_tsv(
    chrom_sizes_file,
    col_names = c("chr", "end"),
    show_col_types = FALSE
  ) |>
    dplyr::mutate(start = 0L) |>
    dplyr::select(chr, start, end) |>
    dplyr::filter(chr %in% get_standard_chroms(genome))
}

#' Order genes for BPCells gene score archr
#'
#' Reproduce BPCells/ArchR gene-score tile ordering for downstream matrix construction.
#'
#' @param genes Gene annotation GRanges used by ArchR-style gene-score ordering.
#' @param chromosome_sizes Data frame/GRanges-like chromosome ranges accepted by
#'   BPCells range normalization.
#' @param blacklist Optional blacklist ranges excluded from gene-score tiles.
#' @param tile_width Width in bases of the tiles used by the ArchR gene-score
#'   algorithm.
#' @param addArchRBug Logical passed through to BPCells' ArchR-compatible tile
#'   generation to reproduce ArchR behavior when required.
#' @return A character vector of gene names ordered by first appearance in the
#'   generated gene-score tile table.
#' @keywords internal

order_genes_for_BPCells_gene_score_archr <- function(
  genes,
  chromosome_sizes,
  blacklist = NULL,
  tile_width = 500,
  addArchRBug = FALSE
) {
  chromosome_sizes <- BPCells:::normalize_ranges(chromosome_sizes)
  if (!is.null(blacklist)) {
    blacklist <- BPCells:::normalize_ranges(blacklist)
  }

  tiles <- BPCells:::gene_score_tiles_archr(
    genes = genes,
    chromosome_sizes = chromosome_sizes,
    tile_width = tile_width,
    addArchRBug = addArchRBug
  )

  if (!is.null(blacklist) && length(blacklist) > 0) {
    blacklist_tiles <- unique(BPCells:::range_overlaps(tiles, blacklist)$from)
    if (length(blacklist_tiles) > 0) {
      tiles <- tiles[-blacklist_tiles, , drop = FALSE]
    }
  }
  if (nrow(tiles) == 0) {
    stop("No ArchR gene-score tiles retained after blacklist filtering.", call. = FALSE)
  }

  # BPCells sizes the gene-weight matrix from max(tiles$gene_idx) after
  # blacklist filtering. Genes without retained tiles therefore need to be
  # before scored genes, otherwise trailing zero-tile genes make BPCells assign
  # too many row names to too few matrix rows.
  scored_gene_idx <- unique(tiles$gene_idx)
  genes[c(setdiff(seq_along(genes), scored_gene_idx), scored_gene_idx)]
}

#' Get roadmap chromHMMs from annotation hub
#'
#' Fetch Roadmap ChromHMM annotation resources from AnnotationHub by EDACC name.
#'
#' @param annot_hub_interface `AnnotationHub::AnnotationHub()` object to query.
#' @param roadmap_EDACC_names Character vector/list of Roadmap EDACC tissue IDs
#'   to retrieve.
#' @return A named list of AnnotationHub ChromHMM resources keyed by requested
#'   EDACC name.
#' @keywords internal

get_roadmap_chromHMMs_from_annotation_hub <- function(annot_hub_interface, roadmap_EDACC_names) {
  roadmap_EDACC_names <- roadmap_EDACC_names |>
    unlist(use.names = FALSE) |>
    as.character() |>
    purrr::discard(is.na) |>
    purrr::discard(\(x) x == "") |>
    unique()
  if (length(roadmap_EDACC_names) == 0) {
    return(list())
  }

  roadmap_df <- annot_hub_interface[["AH41830"]]
  missing_cols <- setdiff(c("EID", "EDACC_NAME"), colnames(roadmap_df))
  if (length(missing_cols) > 0) {
    stop("Roadmap metadata is missing required column(s): ", paste(missing_cols, collapse = ", "))
  }

  roadmap_lookup <- roadmap_df |>
    dplyr::as_tibble() |>
    dplyr::select(EID, EDACC_NAME) |>
    dplyr::filter(.data$EDACC_NAME %in% roadmap_EDACC_names) |>
    dplyr::distinct(.data$EDACC_NAME, .keep_all = TRUE)

  missing_EDACC_names <- setdiff(roadmap_EDACC_names, roadmap_lookup$EDACC_NAME)
  if (length(missing_EDACC_names) > 0) {
    stop("Roadmap EDACC name(s) not found: ", paste(missing_EDACC_names, collapse = ", "))
  }

  core_marks_metadata <- annot_hub_interface |>
    AnnotationHub::query("coreMarks_mnemonics.bed.gz") |>
    S4Vectors::mcols() |>
    as.data.frame() |>
    tibble::rownames_to_column("ah_id") |>
    dplyr::transmute(
      ah_id = .data$ah_id,
      EID = stringr::str_match(.data$title, "^(E[0-9]+)_.*coreMarks_mnemonics\\.bed\\.gz$")[, 2]
    ) |>
    dplyr::filter(!is.na(.data$EID)) |>
    dplyr::distinct(.data$EID, .keep_all = TRUE)

  chromHMM_lookup <- roadmap_lookup |>
    dplyr::left_join(core_marks_metadata, by = "EID") |>
    dplyr::arrange(match(.data$EDACC_NAME, roadmap_EDACC_names))
  missing_chromHMMs <- chromHMM_lookup$EDACC_NAME[is.na(chromHMM_lookup$ah_id)]
  if (length(missing_chromHMMs) > 0) {
    stop("Roadmap chromHMM resource(s) not found for EDACC name(s): ", paste(missing_chromHMMs, collapse = ", "))
  }

  hg19_to_hg38_chain <- annot_hub_interface[["AH14150"]]
  chromHMM_lookup$ah_id |>
    purrr::set_names(chromHMM_lookup$EDACC_NAME) |>
    purrr::map(\(annotation_hub_id) {
      marks_hg38_granges <- annot_hub_interface[[annotation_hub_id]] |>
        rtracklayer::liftOver(chain = hg19_to_hg38_chain) |>
        unlist(use.names = FALSE)
      marks_hg38_granges$abbr <- factor(
        marks_hg38_granges$abbr,
        levels = stringr::str_sort(unique(marks_hg38_granges$abbr), numeric = TRUE)
      )
      split(marks_hg38_granges, marks_hg38_granges$abbr)
    })
}

filter_GRanges_to_chrom_sizes <- function(ranges, chromosome_sizes, allow_empty = FALSE) {
  valid_chroms <- chromosome_sizes$chr
  filtered_ranges <- GenomeInfoDb::keepSeqlevels(
    ranges,
    intersect(GenomeInfoDb::seqlevels(ranges), valid_chroms),
    pruning.mode = "coarse"
  )

  filtered_ranges <- filtered_ranges[
    as.character(GenomeInfoDb::seqnames(filtered_ranges)) %in% valid_chroms
  ]
  if (!allow_empty && length(filtered_ranges) == 0) {
    stop("No ranges remain after filtering to chromosome sizes.")
  }

  filtered_ranges
}

empty_peak_GRanges <- function() {
  gr <- GenomicRanges::GRanges()
  S4Vectors::mcols(gr)$name <- character()
  S4Vectors::mcols(gr)$score <- integer()
  S4Vectors::mcols(gr)$fold_change <- numeric()
  S4Vectors::mcols(gr)$neg_log10pvalue_summit <- numeric()
  S4Vectors::mcols(gr)$neg_log10qvalue_summit <- numeric()
  S4Vectors::mcols(gr)$relative_summit_position <- integer()
  gr
}

#' Build peak calling cluster discovery tibble
#'
#' Build peak-calling branch metadata with optional barcode downsampling.
#'
#' @param BCs_per_peak_cluster_list Named list of barcode vectors, one element
#'   per peak-calling cluster.
#' @param peak_calling_cluster_names Cluster names to include; defaults to the
#'   names of `BCs_per_peak_cluster_list`.
#' @param max_cells_per_cluster Maximum number of barcodes retained per cluster;
#'   larger clusters are downsampled for peak calling.
#' @param seed Random seed passed to stochastic clustering, sampling, or embedding code for reproducibility.
#' @return A tibble with identifiers and derived columns consumed by downstream targets.
#' @keywords internal

build_peak_calling_cluster_discovery_tibble <- function(
  BCs_per_peak_cluster_list,
  peak_calling_cluster_names = names(BCs_per_peak_cluster_list),
  max_cells_per_cluster = 50000,
  seed = 1
) {
  max_cells_per_cluster <- max_cells_per_cluster %||% Inf
  cluster_barcodes <- BCs_per_peak_cluster_list[peak_calling_cluster_names]

  purrr::imap_dfr(cluster_barcodes, \(barcodes, cluster_name) {
    barcodes <- sort(unique(as.character(barcodes)))
    n_cells <- length(barcodes)
    set.seed(seed + match(cluster_name, peak_calling_cluster_names))

    discovery_barcodes <- if (n_cells <= max_cells_per_cluster) {
      barcodes
    } else {
      sample(barcodes, max_cells_per_cluster)
    }

    tibble::tibble(
      peak_calling_cluster_name = cluster_name,
      BCs_per_peak_cluster = list(barcodes),
      BCs_for_peak_discovery = list(sort(discovery_barcodes)),
      n_cells_in_cluster = n_cells,
      n_cells_for_peak_discovery = length(discovery_barcodes),
      was_downsampled = length(discovery_barcodes) < n_cells,
      downsample_fraction = length(discovery_barcodes) / n_cells
    )
  })
}

#' Write ATAC fragments for peak calling cluster
#'
#' Write a per-cluster fragment BED from BPCells fragments and barcode groups.
#'
#' @param ATAC_combined_BPCells_fragment_obj Combined BPCells fragment object for one aggregation, with prefixed cell names.
#' @param BCs_per_peak_cluster Character vector of barcodes assigned to one peak-calling cluster.
#' @param peak_calling_cluster_name Cluster label used in output filenames and peak names.
#' @param output_suffix Suffix appended to command outputs so dynamic branches do not overwrite each other.
#' @param allow_empty Logical; when TRUE, write an empty fragments output instead of failing on zero selected cells/fragments.
#' @return The created file path, submitted job ID, or command result needed by the caller.
#' @keywords internal

write_ATAC_fragments_for_peak_calling_cluster <- function(
  ATAC_combined_BPCells_fragment_obj,
  BCs_per_peak_cluster,
  peak_calling_cluster_name,
  output_suffix = peak_calling_cluster_name,
  allow_empty = FALSE
) {
  fragment_obj <- BPCells::select_cells(ATAC_combined_BPCells_fragment_obj, BCs_per_peak_cluster)

  Sys.setenv(LC_ALL = "C")
  zipped_frag_file_out_path <- get_structured_file_path(
    filetype = "tsv.gz",
    override_suffix = output_suffix
  )
  tmp_frag_path <- tempfile(fileext = ".tsv")
  on.exit(fs::file_delete(tmp_frag_path), add = TRUE)

  BPCells::write_fragments_10x(fragments = fragment_obj, path = tmp_frag_path, append_5th_column = TRUE)

  run_w_error_check(
    command_string = "bgzip",
    arguments_chr = c("-f", "-o", zipped_frag_file_out_path, "--threads", "1", tmp_frag_path)
  )

  if (length(readr::read_lines(zipped_frag_file_out_path, n_max = 1)) == 0) {
    if (!allow_empty) {
      stop("No fragments produced for cluster ", peak_calling_cluster_name)
    }
    return(zipped_frag_file_out_path)
  }

  run_w_error_check(
    command_string = "tabix",
    arguments_chr = c("-f", "-p", "bed", zipped_frag_file_out_path)
  )

  zipped_frag_file_out_path
}

#' Call peaks w MACS3
#'
#' Run MACS3 peak calling for one cluster-specific fragment file.
#'
#' @param ATAC_fragments_per_cluster Path to a BED-like fragment file containing
#'   only the cells for one peak-calling cluster.
#' @param ATAC_peak_calling_cluster_names Cluster label used in errors and, by
#'   default, the output filename suffix.
#' @param genome Genome build key used to choose chromosome sizes, blacklist resources, and external-tool parameters.
#' @param output_suffix Suffix used for the structured output directory and
#'   MACS3 `--name`; defaults to the cluster label.
#' @param allow_no_peaks Logical; when `TRUE`, empty input or empty MACS3 output
#'   is represented as an empty narrowPeak file instead of an error.
#' @return Path to the cluster narrowPeak file. MACS3 summit and XLS sidecar
#'   files are deleted because downstream code only consumes narrowPeak output.
#' @keywords internal

call_peaks_w_MACS3 <- function(ATAC_fragments_per_cluster,
                               ATAC_peak_calling_cluster_names,
                               genome,
                               output_suffix = ATAC_peak_calling_cluster_names,
                               allow_no_peaks = FALSE) {
  out_dir <- get_structured_file_path(override_suffix = output_suffix)
  prefix_path <- file.path(out_dir, output_suffix)
  out_file <- stringr::str_glue("{prefix_path}_peaks.narrowPeak")

  if (allow_no_peaks && length(readr::read_lines(ATAC_fragments_per_cluster, n_max = 1)) == 0) {
    if (fs::file_exists(out_file)) {
      fs::file_delete(out_file)
    }
    fs::file_create(out_file)
    return(out_file)
  }

  genome_size <- as.character(get_effective_genome_size(genome))
  run_w_error_check(
    command_string = "macs3",
    arguments_chr = c(
      "callpeak",
      c("-t", ATAC_fragments_per_cluster),
      c("-g", genome_size),
      c("-f", "BED"),
      c("--nomodel"),
      c("--shift", "-75"),
      c("--extsize", "150"),
      c("--call-summits"),
      c("--name", output_suffix),
      c("--outdir", out_dir),
      c("--keep-dup", "all")
    )
  )

  fs::file_delete(stringr::str_glue("{prefix_path}_summits.bed")) # remove summit file, not needed
  fs::file_delete(stringr::str_glue("{prefix_path}_peaks.xls")) # remove xls file, not needed
  if (!fs::file_exists(out_file)) {
    fs::file_create(out_file)
  }
  if (readr::read_lines(out_file) %>% length() == 0 && !allow_no_peaks) {
    stop("No peaks produced for cluster ", ATAC_peak_calling_cluster_names)
  }
  return(out_file)
}

#' Call peaks w BPCells tile
#'
#' Call peaks with BPCells tile enrichment for one cluster and optional chromosome.
#'
#' @param ATAC_combined_BPCells_fragment_obj Combined BPCells fragment object for one aggregation, with prefixed cell names.
#' @param ATAC_BCs_per_peak_cluster Character vector of barcode names selecting
#'   the cells in the cluster to peak-call.
#' @param ATAC_peak_calling_cluster_names Cluster label written into peak names
#'   and used in error messages.
#' @param genome Genome build key used to choose chromosome sizes, blacklist resources, and external-tool parameters.
#' @param chromosome Chromosome name processed by this branch, matching the naming style used by the input fragments/ranges.
#' @param output_suffix Suffix used for the structured narrowPeak output path.
#' @param allow_no_peaks Logical; when `TRUE`, no called peaks yields an empty
#'   narrowPeak file instead of an error.
#' @param peak_width Width in bases for candidate BPCells tile peaks.
#' @param peak_tiling Number of shifted tilings used by BPCells to scan candidate peaks.
#' @param fdr_cutoff FDR threshold passed to `BPCells::call_peaks_tile()`.
#' @return Path to a narrowPeak-format file written from BPCells peak calls.
#' @keywords internal

call_peaks_w_BPCells_tile <- function(
  ATAC_combined_BPCells_fragment_obj,
  ATAC_BCs_per_peak_cluster,
  ATAC_peak_calling_cluster_names,
  genome,
  chromosome = NULL,
  output_suffix = ATAC_peak_calling_cluster_names,
  allow_no_peaks = FALSE,
  peak_width = 500,
  peak_tiling = 3,
  fdr_cutoff = 0.01
) {
  fragment_obj <- BPCells::select_cells(ATAC_combined_BPCells_fragment_obj, ATAC_BCs_per_peak_cluster)
  if (!is.null(chromosome)) {
    fragment_obj <- BPCells::select_chromosomes(fragment_obj, chromosome)
  }
  cell_groups <- rep.int(ATAC_peak_calling_cluster_names, length(BPCells::cellNames(fragment_obj)))
  chromosome_sizes <- get_chrom_sizes_for_BPCells_tile_calling(genome)
  if (!is.null(chromosome)) {
    chromosome_sizes <- dplyr::filter(chromosome_sizes, chr == chromosome)
  }

  BPCells_tile_peaks <- BPCells::call_peaks_tile(
    fragments = fragment_obj,
    chromosome_sizes = chromosome_sizes,
    cell_groups = cell_groups,
    effective_genome_size = get_effective_genome_size(genome),
    peak_width = peak_width,
    peak_tiling = peak_tiling,
    fdr_cutoff = fdr_cutoff,
    merge_peaks = "none"
  )

  if (nrow(BPCells_tile_peaks) == 0) {
    if (!allow_no_peaks) {
      stop("No BPCells tile peaks produced for cluster ", ATAC_peak_calling_cluster_names)
    }
    out_file <- get_structured_file_path(
      filetype = "narrowPeak",
      override_suffix = output_suffix
    )
    if (fs::file_exists(out_file)) {
      fs::file_delete(out_file)
    }
    fs::file_create(out_file)
    return(out_file)
  }

  out_file <- get_structured_file_path(
    filetype = "narrowPeak",
    override_suffix = output_suffix
  )

  BPCells_tile_peaks |>
    dplyr::mutate(
      name = paste(ATAC_peak_calling_cluster_names, "peak", dplyr::row_number(), sep = "_"),
      score = pmin(1000L, round(-10 * log10(pmax(q_val, .Machine$double.xmin)))),
      strand = ".",
      fold_change = enrichment,
      neg_log10pvalue_summit = -log10(pmax(p_val, .Machine$double.xmin)),
      neg_log10qvalue_summit = -log10(pmax(q_val, .Machine$double.xmin)),
      relative_summit_position = as.integer(floor((end - start) / 2))
    ) |>
    dplyr::select(
      chr,
      start,
      end,
      name,
      score,
      strand,
      fold_change,
      neg_log10pvalue_summit,
      neg_log10qvalue_summit,
      relative_summit_position
    ) |>
    readr::write_tsv(out_file, col_names = FALSE)

  out_file
}

#' Get peak GRanges w fixed width
#'
#' Construct or transform genomic ranges with coordinates compatible with downstream ATAC helpers.
#'
#' @param narrow_peak_path Path to a MACS3/BPCells narrowPeak file. Missing or
#'   empty files return an empty GRanges object.
#' @param extend_summits Number of bases to extend on each side of the summit to
#'   create fixed-width peak ranges.
#' @param genome Genome build key used to choose chromosome sizes, blacklist resources, and external-tool parameters.
#' @param blacklist_GRanges GRanges object containing blacklist GRanges coordinates and metadata.
#' @return A GRanges object with coordinates and metadata columns expected by downstream ATAC helpers.
#' @keywords internal

get_peak_GRanges_w_fixed_width <- function(narrow_peak_path, extend_summits = 250, genome = "GRCh38", blacklist_GRanges) {
  # Set genome-specific parameters
  valid_chroms <- get_standard_chroms(genome)

  if (!fs::file_exists(narrow_peak_path) || fs::file_size(narrow_peak_path) == 0) {
    return(empty_peak_GRanges())
  }

  peak_GRanges <- utils::read.table(
    file = narrow_peak_path,
    col.names = c(
      "chr",
      "start",
      "end",
      "name",
      "score",
      "strand",
      "fold_change",
      "neg_log10pvalue_summit",
      "neg_log10qvalue_summit",
      "relative_summit_position"
    )
  ) %>%
    dplyr::mutate(
      start = start + relative_summit_position - extend_summits,
      end = start + 2 * extend_summits
    ) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE) %>%
    GenomeInfoDb::keepSeqlevels(intersect(GenomeInfoDb::seqlevels(.), valid_chroms), pruning.mode = "coarse")

  # Remove peaks overlapping with blacklisted regions
  blacklist_hits <- peak_GRanges %>% GenomicRanges::findOverlaps(blacklist_GRanges)
  peak_GRanges_wo_blacklist_hits <- peak_GRanges[-S4Vectors::queryHits(blacklist_hits)]

  return(peak_GRanges_wo_blacklist_hits)
}

#' Combine collapse GRanges ArchR
#'
#' Construct or transform genomic ranges with coordinates compatible with downstream ATAC helpers.
#'
#' @param GRanges_list GRanges object containing GRanges list coordinates and metadata.
#' @param by Metadata column used to prioritize overlapping GRanges while collapsing peaks.
#' @param decreasing Logical; when TRUE, higher values of `by` are preferred during overlap collapse.
#' @return A GRanges object with coordinates and metadata columns expected by downstream ATAC helpers.
#' @keywords internal

combine_collapse_GRanges_ArchR <- function(GRanges_list, by, decreasing = TRUE) {
  GRanges_list <- purrr::keep(GRanges_list, \(gr) length(gr) > 0)
  if (length(GRanges_list) == 0) {
    return(empty_peak_GRanges())
  }

  .clusterGRanges <- function(gr, by, decreasing) {
    gr <- sort(GenomeInfoDb::sortSeqlevels(gr))
    r <- GenomicRanges::reduce(gr, min.gapwidth = 0L, ignore.strand = TRUE)
    o <- GenomicRanges::findOverlaps(gr, r, ignore.strand = TRUE)
    GenomicRanges::mcols(gr)$cluster <- S4Vectors::subjectHits(o)
    gr <- gr[order(GenomicRanges::mcols(gr)[, by], decreasing = decreasing), ]
    gr <- gr[!duplicated(GenomicRanges::mcols(gr)$cluster), ]
    gr <- sort(GenomeInfoDb::sortSeqlevels(gr))
    GenomicRanges::mcols(gr)$cluster <- NULL
    return(gr)
  }

  gr <- purrr::reduce(.x = GRanges_list, .f = c)
  i <- 0
  grConverge <- gr
  while (length(grConverge) > 0) {
    i <- i + 1
    grSelect <- .clusterGRanges(grConverge, by, decreasing)
    grConverge <- IRanges::subsetByOverlaps(grConverge, grSelect, invert = TRUE, ignore.strand = TRUE) # blacklist selected gr

    if (i == 1) {
      grAll <- grSelect # if i=1 then set gr_all to clustered
    } else {
      grAll <- c(grAll, grSelect)
    }
  }
  grAll <- sort(GenomeInfoDb::sortSeqlevels(grAll))

  return(grAll)
}

#' Format peak GRanges
#'
#' Construct or transform genomic ranges with coordinates compatible with downstream ATAC helpers.
#'
#' @param GRanges GRanges object containing GRanges coordinates and metadata.
#' @param cluster_name Cluster label written into peak names and output filenames.
#' @return A GRanges object with coordinates and metadata columns expected by downstream ATAC helpers.
#' @keywords internal

format_peak_GRanges <- function(GRanges, cluster_name = NULL) {
  peak_df <- GRanges %>%
    GenomicRanges::as.data.frame()

  if (is.null(cluster_name)) {
    cluster <- if ("cluster" %in% names(peak_df)) {
      as.character(peak_df$cluster)
    } else if ("name" %in% names(peak_df)) {
      stringr::str_match(peak_df$name, "\\.([^\\.]+)\\.")[, 2]
    } else {
      NA_character_
    }
  } else {
    cluster <- as.character(cluster_name)
  }

  peak_df %>%
    dplyr::mutate(
      cluster = .env$cluster,
      summit_position = start + relative_summit_position,
      region_vec = paste(seqnames, start, end, sep = "-")
    ) %>%
    dplyr::select(-dplyr::any_of(c("name", "relative_summit_position"))) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE, ignore.strand = FALSE) %>%
    stats::setNames(.$region_vec)
}


#' Get ATAC QC metadata from BPCells
#'
#' Add BPCells-derived ATAC peak and blacklist QC metrics to cell metadata.
#'
#' @param metadata_df Cell metadata data frame. If `barcode_w_prefix` is absent,
#'   row names are promoted to that column before joining.
#' @param ATAC_peak_BPCells_matrix BPCells peak-by-cell ATAC matrix with peak names in rows and cell barcodes in columns.
#' @param ATAC_combined_BPCells_fragment_obj Combined BPCells fragment object for one aggregation, with prefixed cell names.
#' @param blacklist_GRanges GRanges object containing blacklist GRanges coordinates and metadata.
#' @param ATAC_peak_GRanges GRanges for ATAC peaks, aligned by peak name to the corresponding ATAC matrix rows.
#' @param genome Genome build key used to choose chromosome sizes, blacklist resources, and external-tool parameters.
#' @param peak_matrix_mode Fragment counting mode passed to `BPCells::peak_matrix()`
#'   for blacklist counts.
#' @return The input metadata as a tibble with `nCount_ATAC`,
#'   `atac_peak_counts_blacklist`, fraction/enrichment columns, and missing
#'   count values replaced by zero.
#' @keywords internal

get_ATAC_QC_metadata_from_BPCells <- function(
  metadata_df,
  ATAC_peak_BPCells_matrix,
  ATAC_combined_BPCells_fragment_obj,
  blacklist_GRanges,
  ATAC_peak_GRanges,
  genome,
  peak_matrix_mode
) {
  metadata_tibble <- tibble::as_tibble(metadata_df)
  if (!"barcode_w_prefix" %in% colnames(metadata_tibble)) {
    metadata_tibble <- tibble::rownames_to_column(as.data.frame(metadata_df), var = "barcode_w_prefix") |>
      tibble::as_tibble()
  }

  peak_counts <- BPCells::colSums(ATAC_peak_BPCells_matrix)

  blacklist_matrix <- BPCells::peak_matrix(
    fragments = ATAC_combined_BPCells_fragment_obj,
    ranges = blacklist_GRanges,
    mode = peak_matrix_mode
  )
  blacklist_counts <- BPCells::colSums(blacklist_matrix)

  qc_tibble <- tibble::tibble(
    barcode_w_prefix = colnames(ATAC_peak_BPCells_matrix),
    nCount_ATAC = as.numeric(peak_counts),
    atac_peak_counts_blacklist = as.numeric(blacklist_counts[colnames(ATAC_peak_BPCells_matrix)])
  ) |>
    dplyr::mutate(dplyr::across(c(nCount_ATAC, atac_peak_counts_blacklist), ~ tidyr::replace_na(.x, 0)))

  genome_peak_cov <- ATAC_peak_GRanges |>
    GenomicRanges::width() |>
    sum() |>
    magrittr::divide_by(get_effective_genome_size(genome))

  metadata_tibble |>
    dplyr::left_join(qc_tibble, by = "barcode_w_prefix") |>
    dplyr::mutate(
      dplyr::across(c(nCount_ATAC, atac_peak_counts_blacklist), ~ tidyr::replace_na(.x, 0)),
      atac_peak_counts_frac = nCount_ATAC / .data$atac_fragments,
      atac_peak_counts_blacklist_frac = atac_peak_counts_blacklist / nCount_ATAC,
      log10_nCount_ATAC = log10(nCount_ATAC),
      atac_peak_count_enrichment = atac_peak_counts_frac / genome_peak_cov
    )
}

#' Run ATAC LSI BPCells
#'
#' Compute TF-IDF transformed ATAC LSI embeddings with BPCells SVD.
#'
#' @param ATAC_peak_BPCells_matrix BPCells peak-by-cell ATAC matrix with peak names in rows and cell barcodes in columns.
#' @param n_components Number of singular vectors/components to compute.
#' @param scale_factor Multiplier applied before `log1p()` after column-depth
#'   and inverse-document-frequency scaling.
#' @param threads Number of threads passed to `BPCells::svds()`.
#' @return A list with cell embeddings (`cell_embeddings`), peak loadings
#'   (`feature_loadings`), and singular values (`singular_values`).
#' @keywords internal

run_ATAC_LSI_BPCells <- function(ATAC_peak_BPCells_matrix, n_components, scale_factor = 10000, threads = 1) {
  peak_matrix <- ATAC_peak_BPCells_matrix

  col_sums <- BPCells::colSums(peak_matrix)
  col_scaling <- ifelse(col_sums > 0, 1 / col_sums, 0)

  row_sums <- BPCells::rowSums(peak_matrix)
  idf <- ifelse(row_sums > 0, ncol(peak_matrix) / row_sums, 0)

  mat_lsi <- peak_matrix |>
    BPCells::multiply_cols(col_scaling) |>
    BPCells::multiply_rows(idf)
  mat_lsi <- log1p(scale_factor * mat_lsi)

  svd <- BPCells::svds(mat_lsi, k = n_components, threads = threads)

  cell_embeddings <- sweep(svd$v, 2, svd$d, FUN = "*")
  rownames(cell_embeddings) <- colnames(peak_matrix)
  colnames(cell_embeddings) <- paste0("LSI_", seq_len(ncol(cell_embeddings)))

  feature_loadings <- svd$u
  rownames(feature_loadings) <- rownames(peak_matrix)
  colnames(feature_loadings) <- paste0("LSI_", seq_len(ncol(feature_loadings)))

  list(
    cell_embeddings = cell_embeddings,
    feature_loadings = feature_loadings,
    singular_values = svd$d
  )
}

#' Run UMAP from embedding matrix
#'
#' Run cosine UMAP on selected embedding dimensions, preserving input cell row names in the output matrix.
#'
#' @param embedding_matrix Numeric matrix with cells/barcodes in rows and embedding dimensions in columns; row names are carried into downstream coordinates.
#' @param dims Integer dimension indices to use; combined with `dim_prefix` to select columns such as `PCA_1` or `LSI_2`.
#' @param n_neighbors UMAP neighbor count; clipped below the number of input cells where the helper does that internally.
#' @param min_dist UMAP minimum-distance parameter controlling how tightly local neighborhoods are packed.
#' @param seed Random seed passed to stochastic clustering, sampling, or embedding code for reproducibility.
#' @param dim_prefix Prefix used to translate `dims` into embedding column names, for example `PCA_`, `LSI_`, or `WNN_`.
#' @param col_prefix Prefix assigned to generated coordinate columns, for example `LSI_UMAP` gives `LSI_UMAP_1` and `LSI_UMAP_2`.
#' @param n_components Number of UMAP output dimensions to compute; use 2 for plotting panels and 3 for 3D widgets.
#' @return A numeric UMAP matrix with preserved cell row names and `<col_prefix>_<dimension>` columns.
#' @keywords internal

run_UMAP_from_embedding_matrix <- function(
  embedding_matrix,
  dims,
  n_neighbors,
  min_dist,
  seed = 1,
  dim_prefix = "LSI_",
  col_prefix = "LSI_UMAP",
  n_components = 2
) {
  umap_input <- select_embedding_dimensions(
    embedding_matrix = embedding_matrix,
    dims = dims,
    dim_prefix = dim_prefix
  )
  n_neighbors <- min(as.integer(n_neighbors), nrow(umap_input) - 1L)

  set.seed(seed)
  umap <- uwot::umap(
    X = umap_input,
    n_neighbors = n_neighbors,
    min_dist = min_dist,
    metric = "cosine",
    n_components = n_components
  )

  colnames(umap) <- paste0(col_prefix, "_", seq_len(ncol(umap)))
  rownames(umap) <- rownames(umap_input)
  umap
}

#' Run harmony on embedding matrix
#'
#' Apply Harmony correction to selected embedding dimensions.
#'
#' @param embedding_matrix Numeric matrix with cells/barcodes in rows and embedding dimensions in columns; row names are carried into downstream coordinates.
#' @param metadata_tibble Tibble with one row per cell or pseudobulk sample; must contain the barcode/grouping columns referenced by the helper arguments.
#' @param harmony_correction_metadata_col_names Metadata columns used as Harmony
#'   covariates. `NULL` or length zero returns `embedding_matrix` unchanged.
#' @param dims Integer dimension indices to use; combined with `dim_prefix` to select columns such as `PCA_1` or `LSI_2`.
#' @param cores Number of CPU cores requested for external tools or parallel work.
#' @param dim_prefix Prefix used to translate `dims` into embedding column names, for example `PCA_`, `LSI_`, or `WNN_`.
#' @return A corrected embedding matrix with the selected dimension columns. Cells
#'   with missing Harmony covariates are removed with a warning.
#' @keywords internal

run_harmony_on_embedding_matrix <- function(embedding_matrix, metadata_tibble, harmony_correction_metadata_col_names, dims, cores = 1, dim_prefix = "LSI_") {
  if (is.null(harmony_correction_metadata_col_names) || length(harmony_correction_metadata_col_names) == 0) {
    return(embedding_matrix)
  }

  metadata <- metadata_tibble |>
    dplyr::filter(.data$barcode_w_prefix %in% rownames(embedding_matrix)) |>
    dplyr::arrange(match(.data$barcode_w_prefix, rownames(embedding_matrix)))

  missing_covars <- setdiff(harmony_correction_metadata_col_names, colnames(metadata))
  if (length(missing_covars) > 0) {
    stop("Missing Harmony covariate column(s): ", paste(missing_covars, collapse = ", "))
  }
  barcodes_with_NAs <- metadata |>
    dplyr::filter(dplyr::if_any(dplyr::all_of(harmony_correction_metadata_col_names), ~ is.na(.x))) |>
    dplyr::pull(.data$barcode_w_prefix)
  if (length(barcodes_with_NAs) > 0) {
    keep_barcodes <- setdiff(rownames(embedding_matrix), barcodes_with_NAs)
    embedding_matrix <- embedding_matrix[keep_barcodes, , drop = FALSE]
    metadata <- metadata |>
      dplyr::filter(.data$barcode_w_prefix %in% keep_barcodes) |>
      dplyr::arrange(match(.data$barcode_w_prefix, rownames(embedding_matrix)))
    warning("Removed ", length(barcodes_with_NAs), " cells with NA values in Harmony correction covariates.")
  }

  vars_use <- harmony_correction_metadata_col_names
  if (length(vars_use) > 1) {
    # Harmony 2.0.2 can fail in moe_correct_ridge_cpp() with
    # "inv(): use of LAPACK must be enabled" for multiple covariates in the
    # pixi/conda build. Collapse them into one explicit interaction batch as a
    # pragmatic workaround. This is not identical to additive covariate
    # correction: each observed covariate combination is corrected as its own
    # batch.
    metadata <- metadata |>
      dplyr::mutate(
        harmony_interaction_batch = interaction(!!!rlang::syms(vars_use), drop = TRUE, sep = "__")
      )
    vars_use <- "harmony_interaction_batch"
  }

  harmony_input <- select_embedding_dimensions(
    embedding_matrix = embedding_matrix,
    dims = dims,
    dim_prefix = dim_prefix
  )
  harmonized_embeddings <- harmony::RunHarmony(
    data_mat = harmony_input,
    meta_data = metadata,
    vars_use = vars_use,
    max_iter = 25,
    lambda = 1,
    ncores = cores,
    verbose = TRUE
  )

  colnames(harmonized_embeddings) <- colnames(harmony_input)
  harmonized_embeddings
}

#' Cluster KNN SNN leiden
#'
#' Build an SNN graph from KNN output and cluster it with Leiden modularity.
#'
#' @param knn Nearest-neighbor result with index and distance components, usually from BPCells HNSW helpers.
#' @param resolution Leiden clustering resolution; higher values generally split clusters more finely.
#' @param seed Random seed passed to stochastic clustering, sampling, or embedding code for reproducibility.
#' @return A factor of Leiden cluster assignments in neighbor-graph vertex order.
#' @keywords internal

cluster_knn_snn_leiden <- function(knn, resolution, seed = 1) {
  # BPCells::cluster_graph_leiden() first converts a sparse SNN matrix with
  # igraph::graph_from_adjacency_matrix(), which can spike RAM on large datasets.
  snn <- BPCells::knn_to_snn_graph(knn, return_type = "list")
  snn_graph <- igraph::graph_from_data_frame(
    data.frame(
      from = snn$i + 1L,
      to = snn$j + 1L,
      weight = snn$weight
    ),
    directed = FALSE,
    vertices = data.frame(name = seq_len(snn$dim))
  )

  withr::with_seed(
    seed,
    igraph::cluster_leiden(
      snn_graph,
      weights = igraph::E(snn_graph)$weight,
      resolution = resolution,
      objective_function = "modularity"
    ) |>
      igraph::membership() |>
      as.factor()
  )
}

#' Filter clusters by min barcodes
#'
#' Drop cluster labels that contain too few barcodes.
#'
#' @param clusters Named vector or factor of cluster assignments. Names, when
#'   present, are preserved on the returned factor.
#' @param min_barcodes Minimum number of barcodes required for a cluster/group to be retained; smaller groups are dropped.
#' @return A mixed-sort ordered factor where clusters below `min_barcodes` are
#'   dropped; if the threshold is disabled, all labels are retained.
#' @keywords internal

filter_clusters_by_min_barcodes <- function(clusters, min_barcodes = 100) {
  if (is.null(min_barcodes) || is.na(min_barcodes) || min_barcodes <= 1) {
    return(get_mixsorted_factor(as.character(clusters)) |> purrr::set_names(names(clusters)))
  }

  cluster_sizes <- table(clusters)
  retained_clusters <- names(cluster_sizes)[cluster_sizes >= min_barcodes]
  retained_barcodes <- as.character(clusters) %in% retained_clusters

  if (!any(retained_barcodes)) {
    stop("No clusters have at least ", min_barcodes, " barcodes.")
  }

  n_removed_barcodes <- sum(!retained_barcodes)
  n_removed_clusters <- sum(cluster_sizes < min_barcodes)
  if (n_removed_barcodes > 0) {
    message(
      "Removing ", n_removed_barcodes, " barcodes from ", n_removed_clusters,
      " clusters with fewer than ", min_barcodes, " barcodes."
    )
  }

  retained_clusters <- clusters[retained_barcodes]
  get_mixsorted_factor(as.character(retained_clusters)) |> purrr::set_names(names(retained_clusters))
}

#' Cluster embedding matrix BPCells
#'
#' Cluster selected embedding dimensions with BPCells KNN/SNN and Leiden.
#'
#' @param embedding_matrix Numeric matrix with cells/barcodes in rows and embedding dimensions in columns; row names are carried into downstream coordinates.
#' @param dims Integer dimension indices to use; combined with `dim_prefix` to select columns such as `PCA_1` or `LSI_2`.
#' @param k Number of nearest neighbors to use for KNN/SNN construction.
#' @param resolution Leiden clustering resolution; higher values generally split clusters more finely.
#' @param threads Number of threads passed to BPCells, HNSW, or matrix-stat routines.
#' @param seed Random seed passed to stochastic clustering, sampling, or embedding code for reproducibility.
#' @param dim_prefix Prefix used to translate `dims` into embedding column names, for example `PCA_`, `LSI_`, or `WNN_`.
#' @param min_barcodes Minimum number of barcodes required for a cluster/group to be retained; smaller groups are dropped.
#' @return A named factor of cluster assignments for retained cells.
#' @keywords internal

cluster_embedding_matrix_BPCells <- function(embedding_matrix, dims, k, resolution, threads = 1, seed = 1, dim_prefix = "LSI_", min_barcodes = 100) {
  cluster_input <- select_embedding_dimensions(
    embedding_matrix = embedding_matrix,
    dims = dims,
    dim_prefix = dim_prefix
  )
  k <- min(as.integer(k), nrow(cluster_input) - 1L)
  clusters <- cluster_input |>
    BPCells::knn_hnsw(k = k, metric = "cosine", threads = threads, ef = 500) |>
    cluster_knn_snn_leiden(resolution = resolution, seed = seed)

  cluster_names <- rownames(cluster_input)
  names(clusters) <- cluster_names
  filter_clusters_by_min_barcodes(clusters, min_barcodes = min_barcodes)
}

plot_LSI_loadings_from_tibble <- function(LSI_loadings_tibble, dims, nfeatures = 50) {
  dims |>
    purrr::set_names(paste0("LSI_", dims)) |>
    purrr::map(\(dim_idx) {
      dim_col <- paste0("LSI_", dim_idx)
      LSI_loadings_tibble |>
        dplyr::select(peak, loading = dplyr::all_of(dim_col)) |>
        dplyr::slice_max(order_by = abs(.data$loading), n = nfeatures) |>
        dplyr::mutate(peak = forcats::fct_reorder(.data$peak, .data$loading)) |>
        ggplot2::ggplot(ggplot2::aes(x = loading, y = peak)) +
        ggplot2::geom_col() +
        ggplot2::labs(title = dim_col, x = "Loading", y = ggplot2::element_blank())
    })
}

plot_embedding_singular_values <- function(singular_values, dims) {
  plot_tibble <- tibble::tibble(
    dim = seq_along(singular_values),
    singular_value = singular_values
  ) |>
    dplyr::filter(.data$dim %in% dims)

  plot_tibble |>
    ggplot2::ggplot(ggplot2::aes(x = .data$dim, y = .data$singular_value)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::scale_x_continuous(breaks = plot_tibble$dim) +
    ggplot2::labs(x = "Dimension", y = "Singular value")
}

#' Plot embedding sdev
#'
#' Plot standard deviations of embedding dimensions before and after Harmony.
#'
#' @param embedding_matrix Numeric matrix with cells/barcodes in rows and embedding dimensions in columns; row names are carried into downstream coordinates.
#' @param dims Integer dimension indices to use; combined with `dim_prefix` to select columns such as `PCA_1` or `LSI_2`.
#' @param harmony_embedding_matrix Optional Harmony-corrected embedding matrix with the same row names and dimension naming convention as `embedding_matrix`.
#' @param dim_prefix Prefix used to translate `dims` into embedding column names, for example `PCA_`, `LSI_`, or `WNN_`.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for saving or composition.
#' @keywords internal

plot_embedding_sdev <- function(embedding_matrix, dims, harmony_embedding_matrix = NULL, dim_prefix = "LSI_") {
  get_sdev_tibble <- function(embedding_matrix, embedding_type) {
    dim_cols <- paste0(dim_prefix, dims)
    keep_dims <- dim_cols %in% colnames(embedding_matrix)
    if (!any(keep_dims)) {
      return(tibble::tibble(dim = integer(), sdev = numeric(), embedding_type = character()))
    }

    embedding_matrix[, dim_cols[keep_dims], drop = FALSE] |>
      tibble::as_tibble() |>
      dplyr::summarise(dplyr::across(dplyr::everything(), sd)) |>
      tidyr::pivot_longer(dplyr::everything(), names_to = "dimension", values_to = "sdev") |>
      dplyr::mutate(
        dim = dims[keep_dims],
        embedding_type = embedding_type
      ) |>
      dplyr::select(dim, sdev, embedding_type)
  }

  plot_tibble <- get_sdev_tibble(embedding_matrix, "Non-Harmony")

  if (!is.null(harmony_embedding_matrix)) {
    plot_tibble <- dplyr::bind_rows(plot_tibble, get_sdev_tibble(harmony_embedding_matrix, "Harmony"))
  }

  plot_tibble |>
    ggplot2::ggplot(ggplot2::aes(x = .data$dim, y = .data$sdev, color = .data$embedding_type)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::scale_x_continuous(breaks = sort(unique(plot_tibble$dim))) +
    ggplot2::scale_color_manual(values = c("Non-Harmony" = "#999999", "Harmony" = "#2166AC")) +
    ggplot2::labs(x = "Dimension", y = "Embedding SD", color = NULL) +
    ggplot2::theme(legend.position = "top")
}

empty_embedding_metadata_association_plot <- function(title, message) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = message) +
    ggplot2::labs(title = title) +
    ggplot2::theme_void()
}

#' Prepare embedding metadata association inputs
#'
#' Align embedding dimensions and metadata rows before association testing.
#'
#' @param embedding_matrix Numeric matrix with cells/barcodes in rows and embedding dimensions in columns; row names are carried into downstream coordinates.
#' @param metadata_tibble Tibble with one row per cell or pseudobulk sample; must contain the barcode/grouping columns referenced by the helper arguments.
#' @param dims Integer dimension indices to use; combined with `dim_prefix` to select columns such as `PCA_1` or `LSI_2`.
#' @param dim_prefix Prefix used to translate `dims` into embedding column names, for example `PCA_`, `LSI_`, or `WNN_`.
#' @return A list containing the subset embedding matrix, aligned metadata,
#'   retained dimension column names, and retained numeric dimension indices.
#' @keywords internal

prepare_embedding_metadata_association_inputs <- function(embedding_matrix, metadata_tibble, dims, dim_prefix) {
  dim_cols <- paste0(dim_prefix, dims)
  keep_dims <- dim_cols %in% colnames(embedding_matrix)
  dim_cols <- dim_cols[keep_dims]
  dims <- dims[keep_dims]
  metadata <- metadata_tibble |>
    dplyr::filter(.data$barcode_w_prefix %in% rownames(embedding_matrix)) |>
    dplyr::arrange(match(.data$barcode_w_prefix, rownames(embedding_matrix)))

  list(
    embedding_matrix = embedding_matrix[metadata$barcode_w_prefix, dim_cols, drop = FALSE],
    metadata = metadata,
    dim_cols = dim_cols,
    dims = dims
  )
}

continuous_embedding_r2 <- function(embedding_values, metadata_values) {
  if (!is.numeric(metadata_values)) {
    return(NA_real_)
  }
  keep <- is.finite(embedding_values) & is.finite(metadata_values)
  if (sum(keep) < 3 || stats::sd(embedding_values[keep]) == 0 || stats::sd(metadata_values[keep]) == 0) {
    return(NA_real_)
  }
  stats::cor(embedding_values[keep], metadata_values[keep])^2
}

categorical_embedding_r2 <- function(embedding_values, metadata_values) {
  keep <- is.finite(embedding_values) & !is.na(metadata_values)
  if (sum(keep) < 3 || stats::sd(embedding_values[keep]) == 0) {
    return(NA_real_)
  }
  groups <- factor(metadata_values[keep])
  if (nlevels(groups) < 2) {
    return(NA_real_)
  }
  group_means <- tapply(embedding_values[keep], groups, mean)
  group_n <- table(groups)[names(group_means)]
  total_ss <- sum((embedding_values[keep] - mean(embedding_values[keep]))^2)
  if (total_ss == 0) {
    return(NA_real_)
  }
  sum(as.numeric(group_n) * (group_means - mean(embedding_values[keep]))^2) / total_ss
}

#' Get embedding metadata association tibble
#'
#' Compute per-dimension association scores between embeddings and metadata columns.
#'
#' @param embedding_matrix Numeric matrix with cells/barcodes in rows and embedding dimensions in columns; row names are carried into downstream coordinates.
#' @param metadata Metadata data frame already aligned row-for-row with
#'   `embedding_matrix`.
#' @param dim_cols Character vector of embedding matrix columns to test.
#' @param dims Integer dimension indices to use; combined with `dim_prefix` to select columns such as `PCA_1` or `LSI_2`.
#' @param metadata_cols Character vector of metadata columns to test or plot.
#' @param variable_type Metadata type to summarize; expected values are `continuous` or `categorical`.
#' @return A tibble with `variable`, `dim`, and `metric`, where `metric` is
#'   squared Pearson correlation for continuous metadata and between-group
#'   variance fraction for categorical metadata.
#' @keywords internal

get_embedding_metadata_association_tibble <- function(embedding_matrix, metadata, dim_cols, dims, metadata_cols, variable_type = c("continuous", "categorical")) {
  variable_type <- match.arg(variable_type)
  metadata_cols <- intersect(metadata_cols, colnames(metadata))
  if (identical(variable_type, "continuous")) {
    metadata_cols <- metadata_cols[purrr::map_lgl(metadata_cols, \(metadata_col) is.numeric(metadata[[metadata_col]]))]
  }

  purrr::map_dfr(metadata_cols, \(metadata_col) {
    purrr::map2_dfr(dim_cols, dims, \(dim_col, dim) {
      metric <- switch(
        variable_type,
        continuous = continuous_embedding_r2(embedding_matrix[, dim_col], metadata[[metadata_col]]),
        categorical = categorical_embedding_r2(embedding_matrix[, dim_col], metadata[[metadata_col]])
      )
      tibble::tibble(variable = metadata_col, dim = dim, metric = metric)
    })
  })
}

#' Plot embedding metadata association barplot
#'
#' Plot association strength between one embedding matrix and metadata variables.
#'
#' @param embedding_matrix Numeric matrix with cells/barcodes in rows and embedding dimensions in columns; row names are carried into downstream coordinates.
#' @param harmony_embedding_matrix Optional Harmony-corrected embedding matrix with the same row names and dimension naming convention as `embedding_matrix`.
#' @param metadata_tibble Tibble with one row per cell or pseudobulk sample; must contain the barcode/grouping columns referenced by the helper arguments.
#' @param dims Integer dimension indices to use; combined with `dim_prefix` to select columns such as `PCA_1` or `LSI_2`.
#' @param metadata_cols Character vector of metadata columns to test or plot.
#' @param variable_type Metadata type to summarize; expected values are `continuous` or `categorical`.
#' @param dim_prefix Prefix used to translate `dims` into embedding column names, for example `PCA_`, `LSI_`, or `WNN_`.
#' @param title Plot title used for the assembled association panel.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for saving or composition.
#' @keywords internal

plot_embedding_metadata_association_barplot <- function(
  embedding_matrix,
  harmony_embedding_matrix = NULL,
  metadata_tibble,
  dims,
  metadata_cols,
  variable_type = c("continuous", "categorical"),
  dim_prefix = "LSI_",
  title = NULL
) {
  variable_type <- match.arg(variable_type)
  title <- title %||% variable_type
  inputs <- prepare_embedding_metadata_association_inputs(
    embedding_matrix = embedding_matrix,
    metadata_tibble = metadata_tibble,
    dims = dims,
    dim_prefix = dim_prefix
  )
  plot_tibble <- get_embedding_metadata_association_tibble(
    embedding_matrix = inputs$embedding_matrix,
    metadata = inputs$metadata,
    dim_cols = inputs$dim_cols,
    dims = inputs$dims,
    metadata_cols = metadata_cols,
    variable_type = variable_type
  ) |>
    dplyr::mutate(embedding_type = "Non-Harmony")

  if (!is.null(harmony_embedding_matrix)) {
    harmony_inputs <- prepare_embedding_metadata_association_inputs(
      embedding_matrix = harmony_embedding_matrix,
      metadata_tibble = metadata_tibble,
      dims = dims,
      dim_prefix = dim_prefix
    )
    harmony_plot_tibble <- get_embedding_metadata_association_tibble(
      embedding_matrix = harmony_inputs$embedding_matrix,
      metadata = harmony_inputs$metadata,
      dim_cols = harmony_inputs$dim_cols,
      dims = harmony_inputs$dims,
      metadata_cols = metadata_cols,
      variable_type = variable_type
    ) |>
      dplyr::mutate(embedding_type = "Harmony")
    plot_tibble <- dplyr::bind_rows(plot_tibble, harmony_plot_tibble)
  }

  if (nrow(plot_tibble) == 0 || all(is.na(plot_tibble$metric))) {
    return(empty_embedding_metadata_association_plot(title, "No usable metadata variables"))
  }

  plot_tibble <- plot_tibble |>
    dplyr::mutate(
      variable = factor(.data$variable, levels = unique(.data$variable)),
      embedding_type = factor(.data$embedding_type, levels = c("Non-Harmony", "Harmony"))
    )
  y_max <- max(plot_tibble$metric, na.rm = TRUE)
  y_upper <- if (is.finite(y_max) && y_max > 0) {
    min(1, max(pretty(c(0, y_max), n = 4)))
  } else {
    0.01
  }

  plot <- plot_tibble |>
    ggplot2::ggplot(ggplot2::aes(x = .data$dim, y = .data$metric, fill = .data$embedding_type)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.2, color = "grey60") +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.7, na.rm = TRUE) +
    ggplot2::facet_wrap(~variable, ncol = 1, strip.position = "right") +
    ggplot2::scale_x_continuous(breaks = inputs$dims) +
    ggplot2::scale_y_continuous(
      labels = scales::label_percent(accuracy = 1),
      limits = c(0, y_upper),
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::scale_fill_manual(values = c("Non-Harmony" = "#999999", "Harmony" = "#2166AC")) +
    ggplot2::labs(title = title, x = "Dimension", y = "Variance explained", fill = NULL) +
    ggplot2::theme(
      legend.position = "top",
      strip.text.y.right = ggplot2::element_text(angle = 0, hjust = 0)
    )
}

#' Plot embedding metadata association barplots
#'
#' Assemble technical and biological metadata-association barplots for embeddings.
#'
#' @param embedding_matrix Numeric matrix with cells/barcodes in rows and embedding dimensions in columns; row names are carried into downstream coordinates.
#' @param harmony_embedding_matrix Optional Harmony-corrected embedding matrix with the same row names and dimension naming convention as `embedding_matrix`.
#' @param metadata_tibble Tibble with one row per cell or pseudobulk sample; must contain the barcode/grouping columns referenced by the helper arguments.
#' @param dims Integer dimension indices to use; combined with `dim_prefix` to select columns such as `PCA_1` or `LSI_2`.
#' @param dim_prefix Prefix used to translate `dims` into embedding column names, for example `PCA_`, `LSI_`, or `WNN_`.
#' @param continuous_technical_cols Numeric technical metadata columns to test
#'   against embedding dimensions.
#' @param categorical_technical_cols Categorical technical metadata columns to
#'   test against embedding dimensions.
#' @param continuous_biological_cols Numeric biological metadata columns to test
#'   against embedding dimensions.
#' @param categorical_biological_cols Categorical biological metadata columns to
#'   test against embedding dimensions.
#' @return A ggplot, patchwork, or BPCells trackplot object ready for saving or composition.
#' @keywords internal

plot_embedding_metadata_association_barplots <- function(
  embedding_matrix,
  harmony_embedding_matrix = NULL,
  metadata_tibble,
  dims,
  dim_prefix,
  continuous_technical_cols = character(),
  categorical_technical_cols = character(),
  continuous_biological_cols = character(),
  categorical_biological_cols = character()
) {
  plot_association_barplot <- function(metadata_cols, variable_type, title) {
    plot_embedding_metadata_association_barplot(
      embedding_matrix = embedding_matrix,
      harmony_embedding_matrix = harmony_embedding_matrix,
      metadata_tibble = metadata_tibble,
      dims = dims,
      metadata_cols = metadata_cols,
      variable_type = variable_type,
      dim_prefix = dim_prefix,
      title = title
    )
  }

  list(
    continuous_technical = plot_association_barplot(
      metadata_cols = continuous_technical_cols,
      variable_type = "continuous",
      title = "Continuous technical variables"
    ),
    categorical_technical = plot_association_barplot(
      metadata_cols = categorical_technical_cols,
      variable_type = "categorical",
      title = "Categorical technical variables"
    ),
    continuous_biological = plot_association_barplot(
      metadata_cols = continuous_biological_cols,
      variable_type = "continuous",
      title = "Continuous biological variables"
    ),
    categorical_biological = plot_association_barplot(
      metadata_cols = categorical_biological_cols,
      variable_type = "categorical",
      title = "Categorical biological variables"
    )
  )
}
