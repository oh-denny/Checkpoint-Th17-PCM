#### dotplot
library(Seurat)
library(ggplot2)
library(dplyr)

rds_path = file.path("data", "combined_data_1.rds")
combined_data_integrated = readRDS(rds_path)

obj = combined_data_integrated
DefaultAssay(obj) = "RNA"

markers = list(
  "T cells" = c(
    "Cd3e", "Cd4", "Cd8a", "Cd8b1",
    "Tcf7", "Il7r", "Mki67", "Il17a", "Foxp3"
  ),
  "Macrophages" = c(
    "Adgre1", "Cd68", "Csf1r", "Mrc1", "Lgals3"
  ),
  "Neutrophils" = c(
    "S100a8", "S100a9", "Ly6g", "Csf3r", "Cxcr2"
  )
)

markers = lapply(markers, function(x) intersect(x, rownames(obj)))

p_dotplot = DotPlot(
  object = obj,
  features = markers,
  group.by = "subtype_cell",
  assay = "RNA",
  dot.scale = 10,   
  scale = TRUE
) +
  RotatedAxis() +
  scale_color_gradientn(
    colours = c(
      "#F5F1EE", 
      "#EBCDD4", 
      "#F28CB4", 
      "#D94F8A", 
      "#7A0177"  
    ),
    name = "Average Expression"
  ) +
  scale_size_continuous(
    range = c(1.2, 10), 
    name = "Percent Expressed"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 11,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 12,
      color = "black"
    ),
    axis.title = element_blank(),
    strip.text = element_blank(),
    strip.background = element_blank(),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

p_dotplot
ggsave(
  "Figure2_DotPlot_BonsaiColors.tiff",
  plot = p_dotplot,
  device = ragg::agg_tiff,
  width = 9,
  height = 4,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)




