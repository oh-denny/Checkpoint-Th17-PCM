

library(ggplot2)
library(purrr)
library(stringr)
library(readr)
library(Seurat)
library(dplyr)
library(tidyr)

rds_path = file.path("data", "combined_data_1.rds")
combined_data_integrated = readRDS(rds_path)

cluster_colors = c(
  "0"  = "#F06C45",
  "1"  = "#296DB0",
  "2"  = "#E93E3F",
  "3"  = "#88D3D3",
  "4"  = "#2C8C61",
  "5"  = "#155C42",
  "6"  = "#41704D",
  "7"  = "#0A6E75",
  "8"  = "#84C868",
  "9"  = "#B15929",
  "10" = "#A99099"
)

combined_data_integrated$RNA_snn_res.0.2 = factor(
  combined_data_integrated$RNA_snn_res.0.2,
  levels = as.character(0:10)
)

p_umap = DimPlot(
  object = combined_data_integrated,
  reduction = "umap",
  group.by = "RNA_snn_res.0.2",
  cols = cluster_colors,
  pt.size = 0.15,
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
    legend.position = c(0.72, 0.18),
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
