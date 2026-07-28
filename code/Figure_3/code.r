
library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)

rds_path = file.path("data", "combined_data_1.rds")
combined_data_integrated = readRDS(rds_path)

output_dir = file.path("results", "figure1_audit")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)



t_cells = subset(
  combined_data_integrated,
  subset = cluster_id %in% c("1", "3", "10")
)

DefaultAssay(t_cells) <- "RNA"

t_cells = t_cells |>
  NormalizeData(verbose = FALSE) |>
  FindVariableFeatures(
    selection.method = "vst",
    nfeatures = 2000,
    verbose = FALSE
  ) |>
  ScaleData(
    vars.to.regress = "percent.mt",
    verbose = FALSE
  ) |>
  RunPCA(
    npcs = 30,
    verbose = FALSE
  ) |>
  RunUMAP(
    dims = 1:15,
    n.neighbors = 30,
    min.dist = 0.3,
    verbose = FALSE
  ) |>
  FindNeighbors(
    dims = 1:15,
    verbose = FALSE
  ) |>
  FindClusters(
    resolution = 0.2,
    verbose = FALSE
  )


DimPlot(
  t_cells,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.6
) +
  NoLegend() +
  theme_classic() +
  labs(
    x = "UMAP1",
    y = "UMAP2"
  ) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

FeaturePlot(t_cells, features = c("Cd8a"))




t_cells_CD4 = subset(
  t_cells,
  subset = !cluster_id %in% c("1")
)


DefaultAssay(t_cells_CD4) <- "RNA"

t_cells_CD4 = t_cells_CD4 |>
  NormalizeData(verbose = FALSE) |>
  FindVariableFeatures(
    selection.method = "vst",
    nfeatures = 2000,
    verbose = FALSE
  ) |>
  ScaleData(
    vars.to.regress = "percent.mt",
    verbose = FALSE
  ) |>
  RunPCA(
    npcs = 30,
    verbose = FALSE
  ) |>
  RunUMAP(
    dims = 1:20,
    n.neighbors = 30,
    min.dist = 0.3,
    verbose = FALSE
  ) |>
  FindNeighbors(
    dims = 1:20,
    verbose = FALSE
  ) |>
  FindClusters(
    resolution = 0.2,
    verbose = FALSE
  )


DimPlot(
  t_cells_CD4,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.6
) +
  NoLegend() +
  theme_classic() +
  labs(
    x = "UMAP1",
    y = "UMAP2"
  ) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

FeaturePlot(t_cells_CD4, features = c("Foxp3", "Il17a", "Rora"))
