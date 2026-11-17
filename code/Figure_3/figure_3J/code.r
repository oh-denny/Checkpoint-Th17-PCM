# Figure/panel: 3J
# Purpose: Th17 module score split by experimental group.
# Inputs (under data/ unless session objects): t_cells_CD4.rds.
# Packages: Seurat, ggplot2.
# Outputs: feature_plot.pdf.
# Prerequisites: None; historical CD4 object required.
# Execution: repository root; see README and docs/SCRIPT_INVENTORY.md.



library(ggplot2)
library(Seurat)

t_cells_CD4 = readRDS("./data/t_cells_CD4.rds")
Th17 = c("Rorc","Maf","Ccr6","Il17a","Il17f","Batf","Rora","Il23r","Ahr","Irf4")

Th17_present = intersect(
  Th17,
  rownames(t_cells_CD4)
)

Th17_missing = setdiff(
  Th17,
  rownames(t_cells_CD4)
)

Th17_present
Th17_missing
DefaultAssay(t_cells_CD4) = "RNA"

t_cells_CD4 = AddModuleScore(
  object = t_cells_CD4,
  features = list(Th17_present),
  name = "Th17_Score",
  assay = "RNA",
  ctrl = 100,
  seed = 1234,
  search = FALSE
)

t_cells_CD4@meta.data$groups_plot <- factor(
  t_cells_CD4@meta.data$groups,
  levels = c(
    "naive",
    "pbs",
    "igg",
    "antipd1",
    "antictla4"
  ),
  labels = c(
    "Naive",
    "PBS",
    "IgG",
    "anti-PD-1",
    "anti-CTLA-4"
  )
)

feature_plot = FeaturePlot(
  object = t_cells_CD4,
  features = "Th17_Score1",
  reduction = "umap",
  pt.size = 2,
  min.cutoff = "q01",
  max.cutoff = "q95",
  order = TRUE,
  cols = c("grey85", "#B2182B"),
  ncol = 1,
  split.by = 'groups_plot',
  keep.scale = "feature"
) &
  NoAxes() &
  NoLegend() &
  theme(
    plot.title = element_text(
      color = "#B2182B",
      size = 20,
      face = "bold",
      hjust = 0.5
    )
  )

pdf("feature_plot.pdf", width = 15, height = 3)
print(feature_plot)
dev.off()