
# from obj_sketch
outdir = "paracocco/bonsai_input"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

counts = GetAssayData(obj_sketch,assay = "RNA",layer = "counts")
counts = as(counts, "dgCMatrix")

keep = Matrix::rowSums(counts > 0) > 0
counts = counts[keep, , drop = FALSE]

Matrix::writeMM(counts,
  file = file.path(outdir, "matrix.mtx")
)

writeLines(
  rownames(counts),
  file.path(outdir, "genes.tsv")
)

writeLines(
  colnames(counts),
  file.path(outdir, "cells.tsv")
)

# Metadata
metadata = obj_sketch@meta.data[colnames(counts), ,drop = FALSE]
metadata$cell_id = rownames(metadata)
metadata = metadata[
  , c(
    "cell_id",
    "cartig",
    "groups",
    "seurat_clusters",
    "subtype_cell"
  ),
  drop = FALSE
]
write.csv(
  metadata,
  file.path(outdir, "cell_metadata.csv"),
  row.names = FALSE
)

readLines(
  file.path(outdir, "matrix.mtx"),
  n = 3
)
