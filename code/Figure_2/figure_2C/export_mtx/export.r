# Figure/panel: 2C input
# Purpose: Export sparse counts, identifiers and metadata.
# Inputs (under data/ unless session objects): sketched_data_1.rds.
# Packages: Matrix, Seurat.
# Outputs: data/input/{matrix.mtx, genes.tsv, cells.tsv, cell_metadata.csv}.
# Prerequisites: code_sketch.r or verified historical sketch.
# Execution: repository root; see README and docs/SCRIPT_INVENTORY.md.


# Load the saved sketch; use the repository root as working directory.
library(Seurat)
library(Matrix)
obj_sketch = readRDS("data/sketched_data_1.rds")
outdir = "data/input"
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
