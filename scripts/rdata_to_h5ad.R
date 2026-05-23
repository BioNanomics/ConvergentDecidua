#!/usr/bin/env Rscript
# Convert a Seurat .RData(.gz) file to a 10X-style MTX directory + obs.csv.
#
# Usage:
#   Rscript scripts/rdata_to_h5ad.R <input.RData[.gz]> <output_dir> [object_name]
#
# The Python ingest layer (src/ingest/seurat_rdata.py) loads the resulting
# matrix.mtx + barcodes.tsv + features.tsv + obs.csv into AnnData.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: rdata_to_h5ad.R <input.RData[.gz]> <output_dir> [object_name]")
}
input <- args[1]
output_dir <- args[2]
obj_name <- if (length(args) >= 3) args[3] else NULL

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# load() handles .gz transparently when given a connection
con <- if (grepl("\\.gz$", input)) gzfile(input) else input
env <- new.env()
load(con, envir = env)
if (inherits(con, "connection")) close(con)

seurat_objs <- Filter(function(n) inherits(env[[n]], "Seurat"), ls(envir = env))
if (length(seurat_objs) == 0) {
  stop("No Seurat object found in ", input,
       ". Objects present: ", paste(ls(envir = env), collapse = ", "))
}
if (!is.null(obj_name)) {
  if (!obj_name %in% seurat_objs) {
    stop("Object '", obj_name, "' not in file. Available: ",
         paste(seurat_objs, collapse = ", "))
  }
  obj <- env[[obj_name]]
} else {
  obj <- env[[seurat_objs[1]]]
  message("Using Seurat object: ", seurat_objs[1])
}

assay <- DefaultAssay(obj)

# Seurat >= 5 uses `layer`; older versions use `slot`. Try both, prefer counts.
get_data <- function(obj, assay, key) {
  tryCatch(
    GetAssayData(obj, assay = assay, layer = key),
    error = function(e) GetAssayData(obj, assay = assay, slot = key)
  )
}
counts <- tryCatch(get_data(obj, assay, "counts"), error = function(e) NULL)
if (is.null(counts) || nrow(counts) == 0 || ncol(counts) == 0) {
  message("'counts' layer empty; falling back to 'data'.")
  counts <- get_data(obj, assay, "data")
}

writeMM(counts, file.path(output_dir, "matrix.mtx"))
writeLines(colnames(counts), file.path(output_dir, "barcodes.tsv"))
features <- data.frame(
  gene_id = rownames(counts),
  gene_symbol = rownames(counts),
  stringsAsFactors = FALSE
)
write.table(features, file.path(output_dir, "features.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

meta <- obj@meta.data
meta$barcode <- rownames(meta)
write.csv(meta, file.path(output_dir, "obs.csv"), row.names = FALSE)

message(sprintf("Wrote %d cells x %d features -> %s",
                ncol(counts), nrow(counts), output_dir))
