args <- commandArgs(trailingOnly = TRUE)
cache <- args[[1]]
run <- args[[2]]
wells <- read.delim("cfg_GEM_wells.tsv", check.names = FALSE)
wells <- wells[wells$GEM_well_is_active, ]
config <- yaml::read_yaml("cfg_aggregations.yaml")
stopifnot(isTRUE(config$immune_human_2x$is_active), isTRUE(config$brain_mouse$is_active))
manifest <- read.delim("example_data/public_core_cellranger_count_manifest.tsv")
input_root <- file.path(cache, "inputs", unname(tools::md5sum("example_data/public_core_cellranger_count_manifest.tsv")))
writeLines(input_root, file.path(run, "input-root"))
manifest <- manifest[manifest$GEM_well_ID %in% wells$GEM_well_ID, ]
stopifnot(setequal(unique(manifest$GEM_well_ID), wells$GEM_well_ID))
write.table(manifest, file.path(run, "manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
for (i in seq_len(nrow(wells))) {
  source <- file.path(input_root, wells$GEM_well_ID[[i]])
  dir.create(source, recursive = TRUE, showWarnings = FALSE)
  stopifnot(file.symlink(source, wells$GEM_well_cellranger_arc_count_dir[[i]]))
}
