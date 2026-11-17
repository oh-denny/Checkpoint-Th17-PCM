# Figure/panel: 1C
# Purpose: Annotated lung-cell UMAP.
# Inputs (under data/ unless session objects): combined_data_1.rds.
# Packages: Seurat, dplyr, ggplot2, purrr, ragg, readr, stringr, tidyr.
# Outputs: Figure1B_UMAP_600dpi.tiff.
# Prerequisites: None; annotations and UMAP must already exist.
# Execution: repository root; see README and docs/SCRIPT_INVENTORY.md.



library(ggplot2)
library(purrr)
library(stringr)
library(readr)
library(Seurat)
library(dplyr)
library(tidyr)

rds_path = file.path("data", "combined_data_1.rds")
combined_data_integrated = readRDS(rds_path)

combined_data_integrated$RNA_snn_res.0.2 = factor(
  as.character(combined_data_integrated$RNA_snn_res.0.2),
  levels = as.character(0:10)
)

combined_data_integrated$cluster_celltype = paste0(
  "Cluster ", combined_data_integrated$RNA_snn_res.0.2, " - ", combined_data_integrated$subtype_cell
)

ordered_levels = unique(combined_data_integrated$cluster_celltype[order(as.numeric(as.character(combined_data_integrated$RNA_snn_res.0.2)))])

combined_data_integrated$cluster_celltype = factor(
  combined_data_integrated$cluster_celltype, 
  levels = ordered_levels
)

Idents(combined_data_integrated) = "cluster_celltype"

exact_colors = c(
  "#F06C45", # 0
  "#296DB0", # 1
  "#E93E3F", # 2
  "#88D3D3", # 3
  "#2C8C61", # 4
  "#155C42", # 5
  "#41704D", # 6
  "#0A6E75", # 7
  "#84C868", # 8
  "#B15929", # 9
  "#A99099"  # 10
)

names(exact_colors) <- ordered_levels

p_umap = DimPlot(
  object = combined_data_integrated,
  reduction = "umap",
  group.by = "cluster_celltype",
  cols = exact_colors,
  pt.size = .15,
  raster = FALSE
) +
  labs(
    color = "Clusters"
  ) +
  guides(
    color = guide_legend(
      nrow = 3,
      byrow = FALSE,
      override.aes = list(size = 3)
    )
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(
      size = 10,
      face = "plain",
      hjust = 0.5
    ),
    legend.text = element_text(size = 9),
    legend.key.width = unit(0.4, "cm"),
    legend.key.height = unit(0.4, "cm"),
    legend.background = element_blank()
  )

ggsave(
  filename = "Figure1B_UMAP_600dpi.tiff",
  plot = p_umap,
  device = ragg::agg_tiff,
  width = 6,
  height = 5,
  units = "in",
  res = 600,
  compression = "lzw",
  background = "white"
)
