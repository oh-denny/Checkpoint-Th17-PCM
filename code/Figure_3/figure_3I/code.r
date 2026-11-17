# Figure/panel: 3I
# Purpose: CD4 marker feature plots.
# Inputs (under data/ unless session objects): t_cells_CD4.rds.
# Packages: Seurat, ggplot2.
# Outputs: feature_plot.pdf.
# Prerequisites: None; historical CD4 object required.
# Execution: repository root; see README and docs/SCRIPT_INVENTORY.md.


library(ggplot2)
library(Seurat)

t_cells_CD4 = readRDS("./data/t_cells_CD4.rds")

genes = c(
"Tcf7",
"Lef1",
"Il7r",
"Mki67",
"Top2a",
"Foxp3",
"Il2ra",
"Ctla4",
"Trac",
"Rorc",
"Il17a",
"Rora",
"Icos",
"Cd44",
"Cd69",
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
print(feature_plot)
dev.off()