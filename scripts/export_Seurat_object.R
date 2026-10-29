# Export an aggregation's Seurat/Signac convenience object (multimodal_Seurat_object.8_multimodal_QC) with its own copy
# of the BPCells matrix folders it reads, so it loads without the pipeline's targets store. Its layers stay lazy BPCells
# operations, re-pointed at the copies; ATAC fragment files stay where the GEM wells keep them. Also writes the nucleus
# metadata as a TSV.
# Usage, from the repository root:
#   pixi run --use-environment-activation-cache Rscript scripts/export_Seurat_object.R <aggregation> <output folder>
args <- commandArgs(trailingOnly = TRUE)
aggregation <- args[[1]]
out_dir <- args[[2]]
matrices_dir <- file.path(out_dir, "BPCells_matrices")
fs::dir_create(matrices_dir)

object <- targets::tar_read_raw(paste0("multimodal_Seurat_object.8_multimodal_QC.", aggregation))
copies <- list()
copy_input <- function(input) {
  if (!inherits(input, "MatrixDir")) stop("Unexpected BPCells input of class ", class(input)[[1]])
  if (is.null(copies[[input@dir]])) {
    target <- file.path(matrices_dir, basename(input@dir))
    fs::dir_copy(input@dir, target, overwrite = TRUE)
    copies[[input@dir]] <<- target
  }
  BPCells::open_matrix_dir(copies[[input@dir]])
}
for (assay in SeuratObject::Assays(object)) {
  for (layer in SeuratObject::Layers(object[[assay]])) {
    matrix <- SeuratObject::LayerData(object, assay = assay, layer = layer)
    if (!inherits(matrix, "IterableMatrix")) next
    BPCells::all_matrix_inputs(matrix) <- lapply(BPCells::all_matrix_inputs(matrix), copy_input)
    SeuratObject::LayerData(object, assay = assay, layer = layer) <- matrix
  }
}
saveRDS(object, file.path(out_dir, paste0(aggregation, "_Seurat_object.rds")))
metadata <- tibble::rownames_to_column(object[[]], "cell")
readr::write_tsv(metadata, file.path(out_dir, paste0(aggregation, "_nucleus_metadata.tsv.gz")))
fragment_paths <- unlist(lapply(Signac::Fragments(object[["ATAC"]]), \(fragments) Signac::GetFragmentData(fragments, slot = "file.path")))
message(ncol(object), " nuclei; copied ", length(copies), " BPCells matrix folders; fragment files in ",
  paste(unique(dirname(dirname(dirname(fragment_paths)))), collapse = ", "))
