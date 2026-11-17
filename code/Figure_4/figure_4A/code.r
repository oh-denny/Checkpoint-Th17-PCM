# Figure/panel: 4A
# Purpose: Myeloid–T-cell communication population preparation.
# Inputs (under data/ unless session objects): combined_data_1.rds; t_cells_CD4.rds.
# Packages: Seurat.
# Outputs: DIMT_sub.tiff.
# Prerequisites: None; review hard-coded cluster selections before regeneration.
# Execution: repository root; see README and docs/SCRIPT_INVENTORY.md.


######
library(Seurat)
rds_path = file.path("data", "combined_data_1.rds")
combined_data_integrated = readRDS(rds_path)
t_cells_CD4 = readRDS("./data/t_cells_CD4.rds")

Idents(combined_data_integrated) = "subtype_cell"
head(combined_data_integrated)

obj_comm = subset(combined_data_integrated,idents = c("Macrophage","Neutrophils"))

Idents(t_cells_CD4) = "seurat_clusters"
head(t_cells_CD4)

obj_t = subset(t_cells_CD4,idents = c("1"))
head(obj_comm)
head(obj_t)

DefaultAssay(obj_comm) = "RNA"
DefaultAssay(obj_t) = "RNA"

obj_comm$comm_celltype = as.character(obj_comm$subtype_cell)
obj_t$comm_celltype = "Treg and Th17"

obj_cellchat = merge(
  x = obj_comm,
  y = obj_t,
  project = "Th17_Myeloid"
)

DefaultAssay(obj_cellchat) = "RNA"
obj_cellchat = NormalizeData(
  obj_cellchat,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
obj_cellchat = FindVariableFeatures(
  obj_cellchat,
  selection.method = "vst",
  nfeatures = 2000
)
obj_cellchat = ScaleData(
  obj_cellchat,
  vars.to.regress = "percent.mt"
)
obj_cellchat = RunPCA(
  obj_cellchat,
  npcs = 30
)
obj_cellchat = RunUMAP(
  obj_cellchat,
  dims = 1:10
)
obj_cellchat = FindNeighbors(
  obj_cellchat,
  dims = 1:10
)
obj_cellchat = FindClusters(
  obj_cellchat,
  resolution = 0.2
)

DimPlot(
  obj_cellchat,
  group.by = "comm_celltype",
  reduction = "umap",
  label = TRUE
)

obj_cellchat_clean = subset(
  obj_cellchat,
  cells = colnames(obj_cellchat)[
    (
      obj_cellchat$comm_celltype == "Treg and Th17" &
      obj_cellchat$seurat_clusters == "1"
    ) |
    (
      obj_cellchat$comm_celltype == "Neutrophils" &
      obj_cellchat$seurat_clusters %in% c("0", "2", "6")
    ) |
    (
      obj_cellchat$comm_celltype == "Macrophage" &
      obj_cellchat$seurat_clusters %in% c("4", "5", "7")
    )
  ]
)

DefaultAssay(obj_cellchat_clean) = "RNA"
obj_cellchat_clean = NormalizeData(obj_cellchat_clean)
obj_cellchat_clean = FindVariableFeatures(
  obj_cellchat_clean,
  nfeatures = 2000
)
obj_cellchat_clean = ScaleData(
  obj_cellchat_clean,
  vars.to.regress = "percent.mt"
)
obj_cellchat_clean = RunPCA(
  obj_cellchat_clean,
  npcs = 30
)
obj_cellchat_clean = RunUMAP(
  obj_cellchat_clean,
  dims = 1:10
)
obj_cellchat_clean = FindNeighbors(
  obj_cellchat_clean,
  dims = 1:10
)
obj_cellchat_clean = FindClusters(
  obj_cellchat_clean,
  resolution = 0.2
)

Idents(obj_cellchat_clean) = "comm_celltype"
cols_comm = c(
  "Neutrophils" = "#E53637",
  "Macrophage"  = "#396544",
  "Treg and Th17" = "#22707E"
)

tiff("DIMT_sub.tiff",units = 'in', width = 7.5, height = 5.5, res = 300)
DimPlot(
  obj_cellchat_clean,
  reduction = "umap",
  group.by = "comm_celltype",
  label = F,
  repel = F,
  cols = cols_comm
) +
  theme_classic()
dev.off()