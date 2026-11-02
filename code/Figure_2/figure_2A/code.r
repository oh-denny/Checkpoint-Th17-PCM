

library(ggplot2)
library(purrr)
library(stringr)
library(readr)
library(Seurat)
library(dplyr)
library(tidyr)
library(scales) 


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

markers = FindAllMarkers(combined_data_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

genes = c(
  "Cd74", "H2-Aa", "Igkc", "Cd79a", "H2-Ab1", "Trbc2", "Tcf7", "Trac", 
  "Gm10275", "Rpl12", "S100a9", "Cxcl2", "S100a8", "Il1b", "Acod1", 
  "S100a4", "Maf", "Thy1", "Icos", "Hba-a1", "Hba-a2", "Hbb-bs", 
  "Hbb-bt", "Alas2", "Gzma", "Ccl5", "Nkg7", "Fcer1g", "Ccl3", "Lyz2", 
  "Ctss", "Psap", "Ctsb", "Lgals3", "Tspan7", "Ifitm3", "Cd36", "Epas1", 
  "Ramp2", "Sftpc", "Cbr2", "Sftpa1", "Cxcl15", "Sftpd", "Chil3", 
  "Atp6v0d2", "Lpl", "Ear2", "Plet1", "Stmn1", "Top2a", "Mki67", "Pclaf", 
  "Birc5", "Cd3e", "Cd4", "Sell", "Il7r", "Cd8a", "Cd8b1", "Cst3", 
  "Cd19", "Nfkbia", "Gzmb", "Klrk1", "Ly6c1", "Itgax", "Mpeg1"
)

genes_sig = markers[markers$gene %in% genes,]

Idents(combined_data_integrated) = "RNA_snn_res.0.2"
heat = DoHeatmap(combined_data_integrated, features = genes_sig$gene, slot = "data") + 
  scale_fill_gradientn(colors = c("white", "#d39c83", "#5c1d24"),
  limits = c(0, 6),                  
    oob = squish)

ggsave(
  filename = "Figure2A_UMAP_600dpi.tiff",
  plot = heat,
  device = ragg::agg_tiff,
  width = 6,
  height = 5,
  units = "in",
  res = 600,
  compression = "lzw",
  background = "white"
)


