# Figure/panel: 2C input
# Purpose: 5,000-cell leverage-score sketch.
# Inputs (under data/ unless session objects): combined_data_1.rds.
# Packages: Seurat, SeuratObject.
# Outputs: sketched_data_1.rds.
# Prerequisites: None; overwrites sketched_data_1.rds if run.
# Execution: repository root; see README and docs/SCRIPT_INVENTORY.md.

library(Seurat)
library(SeuratObject)

set.seed(123)

rds_path = file.path("data", "combined_data_1.rds")
combined_data_integrated = readRDS(rds_path)

combined_data_integrated$RNA_snn_res.0.2 = factor(
  combined_data_integrated$RNA_snn_res.0.2,
  levels = as.character(0:10)
)

DefaultAssay(combined_data_integrated) <- "RNA"
obj = NormalizeData(combined_data_integrated,assay = "RNA",verbose = FALSE)
obj = FindVariableFeatures(
  obj,
  assay = "RNA",
  selection.method = "vst",
  nfeatures = 3000,
  verbose = FALSE
)

length(VariableFeatures(obj))

obj = SketchData(
  object = obj,
  assay = "RNA",
  ncells = 5000,
  method = "LeverageScore",
  sketched.assay = "sketch",
  seed = 123,
  verbose = TRUE
)

sketch_cells = colnames(obj[["sketch"]])

obj_sketch = subset(
  x = combined_data_integrated,
  cells = sketch_cells
)

DefaultAssay(obj_sketch) = "RNA"

dim(obj_sketch)
table(obj_sketch$groups)
table(obj_sketch$cartig, obj_sketch$groups)

saveRDS(obj_sketch, file = file.path("data", "sketched_data_1.rds"))