
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


cols <- c(
  "0" = "#8B2A3E",  # vinho
  "1" = "#22707E",  # azul petróleo
  "2" = "#244847",  # verde petróleo escuro
  "3" = "#5F7574",  # cinza azulado
  "4" = "#5E6A30",  # verde oliva
  "5" = "#A5927B"   # bege acinzentado
)


DimPlot(
  t_cells_CD4,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = F,
  repel = F,
  pt.size = 1,
  cols = cols
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



genes <- c(

# Naive
"Tcf7",
"Lef1",
"Il7r",

# Proliferação
"Mki67",
"Top2a",

# Treg
"Foxp3",
"Il2ra",
"Ctla4",
"Trac",

# Th17
"Rorc",
"Il17a",
"Rora",

# Ativação
"Icos",
"Cd44",
"Cd69",

# Hipóxia
"Hif1a"
)


feature_plot = FeaturePlot(
  object = t_cells_CD4,
  features = genes,
  reduction = "umap",
  pt.size = 1.5,
  min.cutoff = "q01",
  max.cutoff = "q95",
  order = TRUE,
  cols = c("grey85", "#B2182B"),
  ncol = 4,
  keep.scale = "feature"
) &
  NoAxes() &
  NoLegend() &
  theme(
    plot.title = element_text(
      color = "#B2182B",
      size = 18,
      face = "bold",
      hjust = 0.5
    )
  )

pdf("feature_plot.pdf", width = 8, height = 7)
feature_plot
dev.off()
