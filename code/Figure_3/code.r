

#Supplementary figura

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





# supp
genes <- c(
# Naive
"Cd3e",
"Cd3d",
"Cd4",

# Proliferação
"Cd8a",
"Cd8b1",

# Treg
"Il7r",
"Il2ra",
"Il2rb",
"Thy1",

# Th17
"Eomes",
"Ifng",
"Xcl1"
)


feature_plot = FeaturePlot(
  object = t_cells,
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

pdf("feature_plot_supp.pdf", width = 8, height = 7)
feature_plot
dev.off()
