

library(Seurat)
library(CellChat)
library(dplyr)


cellchat_list = readRDS("./data/cellchat_list.rds")

names(cellchat_list)

cellchat_list$pbs
cellchat_list$antipd1
cellchat_list$antictla4


levels(cellchat_list$pbs@idents)


pbs_treg_neut = subsetCommunication(
  cellchat_list$pbs,
  sources.use = "Treg and Th17",
  targets.use = "Neutrophils"
)

pd1_treg_neut = subsetCommunication(
  cellchat_list$antipd1,
  sources.use = "Treg and Th17",
  targets.use = "Neutrophils"
)

ctla4_treg_neut = subsetCommunication(
  cellchat_list$antictla4,
  sources.use = "Treg and Th17",
  targets.use = "Neutrophils"
)



pbs_treg_neut <- pbs_treg_neut |>
  mutate(group = "PBS")

pd1_treg_neut <- pd1_treg_neut |>
  mutate(group = "anti-PD-1")

ctla4_treg_neut <- ctla4_treg_neut |>
  mutate(group = "anti-CTLA-4")

treg_neut_all <- bind_rows(
  pbs_treg_neut,
  pd1_treg_neut,
  ctla4_treg_neut
)



treg_neut_all |>
  dplyr::select(
    group,
    ligand,
    receptor,
    pathway_name,
    prob,
    pval
  ) |>
  arrange(group, desc(prob))


netVisual_bubble(
  cellchat_list$pbs,
  sources.use = "Treg and Th17",
  targets.use = "Neutrophils",
  remove.isolate = FALSE
)

netVisual_bubble(
  cellchat_list$antipd1,
  sources.use = "Treg and Th17",
  targets.use = "Neutrophils",
  remove.isolate = FALSE
)

netVisual_bubble(
  cellchat_list$antictla4,
  sources.use = "Treg and Th17",
  targets.use = "Neutrophils",
  remove.isolate = FALSE
)


sort(unique(treg_neut_all$pathway_name))
table(treg_neut_all$pathway_name)



neut_treg_pbs <- subsetCommunication(
  cellchat_list$pbs,
  sources.use = "Neutrophils",
  targets.use = "Treg and Th17"
)

neut_treg_pd1 <- subsetCommunication(
  cellchat_list$antipd1,
  sources.use = "Neutrophils",
  targets.use = "Treg and Th17"
)

neut_treg_ctl <- subsetCommunication(
  cellchat_list$antictla4,
  sources.use = "Neutrophils",
  targets.use = "Treg and Th17"
)


neut_treg_all <- bind_rows(
  neut_treg_pbs |> mutate(group = "PBS"),
  neut_treg_pd1 |> mutate(group = "anti-PD-1"),
  neut_treg_ctl |> mutate(group = "anti-CTLA-4")
)

neut_treg_all |>
  select(
    group,
    ligand,
    receptor,
    pathway_name,
    prob,
    pval
  ) |>
  arrange(group, desc(prob))



treg_mac_pbs <- subsetCommunication(
  cellchat_list$pbs,
  sources.use = "Treg and Th17",
  targets.use = "Macrophage"
)

treg_mac_pd1 <- subsetCommunication(
  cellchat_list$antipd1,
  sources.use = "Treg and Th17",
  targets.use = "Macrophage"
)

treg_mac_ctl <- subsetCommunication(
  cellchat_list$antictla4,
  sources.use = "Treg and Th17",
  targets.use = "Macrophage"
)


treg_mac_all <- bind_rows(
  treg_mac_pbs |> mutate(group = "PBS"),
  treg_mac_pd1 |> mutate(group = "anti-PD-1"),
  treg_mac_ctl |> mutate(group = "anti-CTLA-4")
)

treg_mac_all |>
  select(
    group,
    ligand,
    receptor,
    pathway_name,
    prob,
    pval
  ) |>
  arrange(group, desc(prob))



mac_treg_pbs <- subsetCommunication(
  cellchat_list$pbs,
  sources.use = "Macrophage",
  targets.use = "Treg and Th17"
)

mac_treg_pd1 <- subsetCommunication(
  cellchat_list$antipd1,
  sources.use = "Macrophage",
  targets.use = "Treg and Th17"
)

mac_treg_ctl <- subsetCommunication(
  cellchat_list$antictla4,
  sources.use = "Macrophage",
  targets.use = "Treg and Th17"
)


mac_treg_all <- bind_rows(
  mac_treg_pbs |> mutate(group = "PBS"),
  mac_treg_pd1 |> mutate(group = "anti-PD-1"),
  mac_treg_ctl |> mutate(group = "anti-CTLA-4")
)

mac_treg_all |>
  select(
    group,
    ligand,
    receptor,
    pathway_name,
    prob,
    pval
  ) |>
  arrange(group, desc(prob))






cellchat_summary <- data.frame(
  Interaction = c(
    "Treg/Th17 → Neutrophils",
    "Neutrophils → Treg/Th17"
  ),

  PBS = c(
    cellchat_list$pbs@net$weight[
      "Treg and Th17",
      "Neutrophils"
    ],
    cellchat_list$pbs@net$weight[
      "Neutrophils",
      "Treg and Th17"
    ]
  ),

  anti_PD1 = c(
    cellchat_list$antipd1@net$weight[
      "Treg and Th17",
      "Neutrophils"
    ],
    cellchat_list$antipd1@net$weight[
      "Neutrophils",
      "Treg and Th17"
    ]
  ),

  anti_CTLA4 = c(
    cellchat_list$antictla4@net$weight[
      "Treg and Th17",
      "Neutrophils"
    ],
    cellchat_list$antictla4@net$weight[
      "Neutrophils",
      "Treg and Th17"
    ]
  )
)

cellchat_summary


library(dplyr)
library(tidyr)
library(ggplot2)

cellchat_long <- cellchat_summary |>
  pivot_longer(
    cols = c("PBS", "anti_PD1", "anti_CTLA4"),
    names_to = "group",
    values_to = "weight"
  ) |>
  mutate(
    group = factor(
      group,
      levels = c(
        "PBS",
        "anti_PD1",
        "anti_CTLA4"
      ),
      labels = c(
        "PBS",
        "anti-PD-1",
        "anti-CTLA-4"
      )
    )
  )

interaction_colors <- c(
  "Treg/Th17 → Neutrophils" = "#D95F8D",
  "Neutrophils → Treg/Th17" = "#355C3A"
)

p_neut <- ggplot(
  cellchat_long,
  aes(
    x = group,
    y = weight,
    group = Interaction,
    color = Interaction
  )
) +
  geom_line(
    linewidth = 1
  ) +
  geom_point(
    size = 4
  ) +
  scale_color_manual(
    values = interaction_colors
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = NULL,
    y = "Interaction strength",
    color = NULL
  ) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5
    )
  )

p_neut






cellchat_summary_mac <- data.frame(

  Interaction = c(
    "Treg/Th17 → Macrophage",
    "Macrophage → Treg/Th17"
  ),

  PBS = c(
    cellchat_list$pbs@net$weight[
      "Treg and Th17",
      "Macrophage"
    ],
    cellchat_list$pbs@net$weight[
      "Macrophage",
      "Treg and Th17"
    ]
  ),

  anti_PD1 = c(
    cellchat_list$antipd1@net$weight[
      "Treg and Th17",
      "Macrophage"
    ],
    cellchat_list$antipd1@net$weight[
      "Macrophage",
      "Treg and Th17"
    ]
  ),

  anti_CTLA4 = c(
    cellchat_list$antictla4@net$weight[
      "Treg and Th17",
      "Macrophage"
    ],
    cellchat_list$antictla4@net$weight[
      "Macrophage",
      "Treg and Th17"
    ]
  )
)

cellchat_summary_mac



cellchat_mac_long <- cellchat_summary_mac |>
  pivot_longer(
    cols = c(
      "PBS",
      "anti_PD1",
      "anti_CTLA4"
    ),
    names_to = "group",
    values_to = "weight"
  ) |>
  mutate(
    group = factor(
      group,
      levels = c(
        "PBS",
        "anti_PD1",
        "anti_CTLA4"
      ),
      labels = c(
        "PBS",
        "anti-PD-1",
        "anti-CTLA-4"
      )
    )
  )

mac_colors <- c(
  "Treg/Th17 → Macrophage" = "#D95F8D",
  "Macrophage → Treg/Th17" = "#41704D"
)

p_mac <- ggplot(
  cellchat_mac_long,
  aes(
    x = group,
    y = weight,
    group = Interaction,
    color = Interaction
  )
) +
  geom_line(
    linewidth = 1
  ) +
  geom_point(
    size = 4
  ) +
  scale_color_manual(
    values = mac_colors
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = NULL,
    y = "Interaction strength",
    color = NULL
  ) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5
    )
  )

p_mac

cellchat_combined <- bind_rows(

  cellchat_long |>
    mutate(
      Myeloid = "Neutrophils"
    ),

  cellchat_mac_long |>
    mutate(
      Myeloid = "Macrophage"
    )
)


ggplot(
  cellchat_combined,
  aes(
    x = group,
    y = weight,
    group = Interaction,
    color = Interaction
  )
) +
  geom_line(
    linewidth = 1
  ) +
  geom_point(
    size = 3.5
  ) +
  facet_wrap(
    ~ Myeloid,
    scales = "free_y"
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = NULL,
    y = "Interaction strength",
    color = NULL
  ) +
  theme(
    legend.position = "top",
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold"
    )
  )












sapply(
  cellchat_list[c("pbs", "antipd1", "antictla4")],
  function(x) {
    sum(x@net$weight["Neutrophils", ])
  }
)


sapply(
  cellchat_list[c("pbs", "antipd1", "antictla4")],
  function(x) {
    sum(x@net$weight[, "Treg and Th17"])
  }
)


neut_outgoing <- bind_rows(
  subsetCommunication(
    cellchat_list$pbs,
    sources.use = "Neutrophils"
  ) |> mutate(group = "PBS"),

  subsetCommunication(
    cellchat_list$antipd1,
    sources.use = "Neutrophils"
  ) |> mutate(group = "anti-PD-1"),

  subsetCommunication(
    cellchat_list$antictla4,
    sources.use = "Neutrophils"
  ) |> mutate(group = "anti-CTLA-4")
)

neut_outgoing |>
  group_by(group, target) |>
  summarise(
    weight = sum(prob),
    n_interactions = n(),
    .groups = "drop"
  ) |>
  arrange(group, desc(weight)) |> data.frame()






######

t_cells_CD4
Idents(combined_data_integrated) <- "subtype_cell"
head(combined_data_integrated)

obj_comm <- subset(
  combined_data_integrated,
  idents = c(
    "Macrophage",
    "Neutrophils"
  )
)

Idents(t_cells_CD4) <- "seurat_clusters"
head(t_cells_CD4)

obj_t <- subset(
  t_cells_CD4,
  idents = c(
    "1"
  )
)


head(obj_comm)
head(obj_t)



library(Seurat)

# garantir RNA como assay principal
DefaultAssay(obj_comm) <- "RNA"
DefaultAssay(obj_t) <- "RNA"

# colocar uma identidade biológica explícita
obj_comm$comm_celltype <- as.character(obj_comm$subtype_cell)
obj_t$comm_celltype <- "Treg and Th17"


obj_cellchat <- merge(
  x = obj_comm,
  y = obj_t,
  project = "Th17_Myeloid"
)
anyDuplicated(colnames(obj_cellchat))


table(obj_cellchat$comm_celltype)
table(obj_cellchat$groups)
table(obj_cellchat$comm_celltype, obj_cellchat$groups)


DefaultAssay(obj_cellchat) <- "RNA"

obj_cellchat <- NormalizeData(
  obj_cellchat,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

obj_cellchat <- FindVariableFeatures(
  obj_cellchat,
  selection.method = "vst",
  nfeatures = 2000
)

obj_cellchat <- ScaleData(
  obj_cellchat,
  vars.to.regress = "percent.mt"
)

obj_cellchat <- RunPCA(
  obj_cellchat,
  npcs = 30
)

obj_cellchat <- RunUMAP(
  obj_cellchat,
  dims = 1:10
)

obj_cellchat <- FindNeighbors(
  obj_cellchat,
  dims = 1:10
)

obj_cellchat <- FindClusters(
  obj_cellchat,
  resolution = 0.2
)


DimPlot(
  obj_cellchat,
  group.by = "comm_celltype",
  reduction = "umap",
  label = TRUE
)


FeaturePlot(
  obj_cellchat,
  features = c(
    "Cd3d", "Cd3e", "Trac",
    "Cd4", "Il7r",
    "S100a8", "S100a9", "Ly6g"
  ),
  ncol = 3
)


table(
  obj_cellchat$seurat_clusters,
  obj_cellchat$comm_celltype
)




obj_cellchat_clean <- subset(
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

table(
  obj_cellchat_clean$seurat_clusters,
  obj_cellchat_clean$comm_celltype
)

table(obj_cellchat_clean$comm_celltype)


DefaultAssay(obj_cellchat_clean) <- "RNA"

obj_cellchat_clean <- NormalizeData(obj_cellchat_clean)

obj_cellchat_clean <- FindVariableFeatures(
  obj_cellchat_clean,
  nfeatures = 2000
)

obj_cellchat_clean <- ScaleData(
  obj_cellchat_clean,
  vars.to.regress = "percent.mt"
)

obj_cellchat_clean <- RunPCA(
  obj_cellchat_clean,
  npcs = 30
)

obj_cellchat_clean <- RunUMAP(
  obj_cellchat_clean,
  dims = 1:10
)

obj_cellchat_clean <- FindNeighbors(
  obj_cellchat_clean,
  dims = 1:10
)

obj_cellchat_clean <- FindClusters(
  obj_cellchat_clean,
  resolution = 0.2
)

Idents(obj_cellchat_clean) <- "comm_celltype"


cols_comm <- c(
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






# objeto já limpo
obj_plot <- obj_cellchat_clean

# dataframe do UMAP
umap_df <- Embeddings(obj_plot, reduction = "umap") |>
  as.data.frame()

colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")

umap_df$celltype <- obj_plot$comm_celltype[
  rownames(umap_df)
]

# centroides
centroids <- umap_df |>
  group_by(celltype) |>
  summarise(
    x = median(UMAP_1),
    y = median(UMAP_2),
    .groups = "drop"
  )

centroids


cc <- cellchat_list$antictla4

edges <- data.frame(
  source = c(
    "Macrophage",
    "Treg and Th17",
    "Neutrophils",
    "Treg and Th17"
  ),
  target = c(
    "Treg and Th17",
    "Macrophage",
    "Treg and Th17",
    "Neutrophils"
  ),
  weight = c(
    cc@net$weight["Macrophage", "Treg and Th17"],
    cc@net$weight["Treg and Th17", "Macrophage"],
    cc@net$weight["Neutrophils", "Treg and Th17"],
    cc@net$weight["Treg and Th17", "Neutrophils"]
  )
)

edges



edges_plot <- edges |>
  left_join(
    centroids,
    by = c("source" = "celltype")
  ) |>
  rename(
    x = x,
    y = y
  ) |>
  left_join(
    centroids,
    by = c("target" = "celltype"),
    suffix = c("", "_end")
  )

edges_plot


cols_comm <- c(
  "Macrophage"    = "#396544",
  "Neutrophils"   = "#E53637",
  "Treg and Th17" = "#2C7FB8"
)

p_comm <- ggplot(
  umap_df,
  aes(
    x = UMAP_1,
    y = UMAP_2
  )
) +

  geom_point(
    aes(color = celltype),
    size = 0.65,
    alpha = 0.75
  ) +

  geom_curve(
    data = edges_plot,
    aes(
      x = x,
      y = y,
      xend = x_end,
      yend = y_end,
      linewidth = weight
    ),
    inherit.aes = FALSE,
    curvature = 0.28,
    color = "#355C9A",
    alpha = 0.95,
    arrow = arrow(
      type = "closed",
      length = unit(0.22, "cm")
    ),
    lineend = "round"
  ) +

  scale_color_manual(
    values = cols_comm
  ) +

  scale_linewidth_continuous(
    range = c(0.4, 3.5),
    name = "CellChat\ninteraction strength"
  ) +

  theme_void(base_size = 13) +

  theme(
    legend.position = "right"
  ) +

  labs(
    color = NULL
  )

p_comm







Idents(obj_cellchat_clean) <- "comm_celltype"

table(obj_cellchat_clean$comm_celltype)
table(obj_cellchat_clean$groups)


obj_pbs <- subset(
  obj_cellchat_clean,
  subset = groups == "pbs"
)

obj_pd1 <- subset(
  obj_cellchat_clean,
  subset = groups == "antipd1"
)

obj_ctl <- subset(
  obj_cellchat_clean,
  subset = groups == "antictla4"
)


library(CellChat)

run_cellchat <- function(obj) {

  DefaultAssay(obj) <- "RNA"

  data.input <- GetAssayData(
    obj,
    assay = "RNA",
    layer = "data"
  )

  meta <- obj@meta.data

  cellchat <- createCellChat(
    object = data.input,
    meta = meta,
    group.by = "comm_celltype"
  )

  cellchat@DB <- CellChatDB.mouse

  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)

  cellchat <- computeCommunProb(
    cellchat,
    type = "triMean"
  )

  cellchat <- filterCommunication(
    cellchat,
    min.cells = 10
  )

  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)

  return(cellchat)
}


cc_pbs_clean <- run_cellchat(obj_pbs)
cc_pd1_clean <- run_cellchat(obj_pd1)
cc_ctl_clean <- run_cellchat(obj_ctl)



subsetCommunication(
  cc_pbs_clean,
  sources.use = "Neutrophils",
  targets.use = "Treg and Th17"
)

subsetCommunication(
  cc_pd1_clean,
  sources.use = "Neutrophils",
  targets.use = "Treg and Th17"
)

subsetCommunication(
  cc_ctl_clean,
  sources.use = "Neutrophils",
  targets.use = "Treg and Th17"
)


subsetCommunication(
  cc_pbs_clean,
  sources.use = "Treg and Th17",
  targets.use = "Neutrophils"
)




library(CellChat)

cols_comm <- c(
  "Macrophage"    = "#396544",
  "Neutrophils"   = "#E53637",
  "Treg and Th17" = "#2C7FB8"
)


max_weight <- max(
  cc_pbs_clean@net$weight,
  cc_pd1_clean@net$weight,
  cc_ctl_clean@net$weight,
  na.rm = TRUE
)

max_weight


par(mfrow = c(1, 3))


tiff("communication_pbs.tiff",units = 'in', width = 7.5, height = 5.5, res = 300)
netVisual_circle(
  cc_pbs_clean@net$weight,
  vertex.weight = as.numeric(table(cc_pbs_clean@idents)),
  weight.scale = TRUE,
  edge.weight.max = max_weight,
  vertex.label.cex = 0,
  label.edge = FALSE,
  vertex.label.color = "transparent",
  color.use = cols_comm,
  title.name = "PBS"
)
dev.off()

tiff("communication_pd1.tiff",units = 'in', width = 7.5, height = 5.5, res = 300)
netVisual_circle(
  cc_pd1_clean@net$weight,
  vertex.weight = as.numeric(table(cc_pd1_clean@idents)),
  weight.scale = TRUE,
  edge.weight.max = max_weight,
  vertex.label.cex = 0,
  label.edge = FALSE,
  vertex.label.color = "transparent",
  color.use = cols_comm,
  title.name = "anti-PD-1"
)
dev.off()

tiff("communication_ctla4.tiff",units = 'in', width = 7.5, height = 5.5, res = 300)
netVisual_circle(
  cc_ctl_clean@net$weight,
  vertex.weight = as.numeric(table(cc_ctl_clean@idents)),
  weight.scale = TRUE,
  edge.weight.max = max_weight,
  vertex.label.cex = 0,
  label.edge = FALSE,
  vertex.label.color = "transparent",
  color.use = cols_comm,
  title.name = "anti-CTLA-4"
)
dev.off()

par(mfrow = c(1, 1))


get_comm_weight <- function(cc, source, target) {

  mat <- cc@net$weight

  if (
    source %in% rownames(mat) &&
    target %in% colnames(mat)
  ) {
    return(mat[source, target])
  }

  return(0)
}


cellchat_summary_clean <- data.frame(

  Interaction = c(
    "Macrophage → Treg/Th17",
    "Treg/Th17 → Macrophage",
    "Neutrophils → Treg/Th17",
    "Treg/Th17 → Neutrophils"
  ),

  PBS = c(
    get_comm_weight(
      cc_pbs_clean,
      "Macrophage",
      "Treg and Th17"
    ),

    get_comm_weight(
      cc_pbs_clean,
      "Treg and Th17",
      "Macrophage"
    ),

    get_comm_weight(
      cc_pbs_clean,
      "Neutrophils",
      "Treg and Th17"
    ),

    get_comm_weight(
      cc_pbs_clean,
      "Treg and Th17",
      "Neutrophils"
    )
  ),

  anti_PD1 = c(
    get_comm_weight(
      cc_pd1_clean,
      "Macrophage",
      "Treg and Th17"
    ),

    get_comm_weight(
      cc_pd1_clean,
      "Treg and Th17",
      "Macrophage"
    ),

    get_comm_weight(
      cc_pd1_clean,
      "Neutrophils",
      "Treg and Th17"
    ),

    get_comm_weight(
      cc_pd1_clean,
      "Treg and Th17",
      "Neutrophils"
    )
  ),

  anti_CTLA4 = c(
    get_comm_weight(
      cc_ctl_clean,
      "Macrophage",
      "Treg and Th17"
    ),

    get_comm_weight(
      cc_ctl_clean,
      "Treg and Th17",
      "Macrophage"
    ),

    get_comm_weight(
      cc_ctl_clean,
      "Neutrophils",
      "Treg and Th17"
    ),

    get_comm_weight(
      cc_ctl_clean,
      "Treg and Th17",
      "Neutrophils"
    )
  )
)

cellchat_summary_clean


library(tidyr)
library(dplyr)

cellchat_long <- cellchat_summary_clean |>
  pivot_longer(
    cols = c(
      PBS,
      anti_PD1,
      anti_CTLA4
    ),
    names_to = "group",
    values_to = "weight"
  ) |>
  mutate(
    group = factor(
      group,
      levels = c(
        "PBS",
        "anti_PD1",
        "anti_CTLA4"
      ),
      labels = c(
        "PBS",
        "anti-PD-1",
        "anti-CTLA-4"
      )
    ),

    Myeloid = case_when(
      grepl("Macrophage", Interaction) ~ "Macrophage",
      grepl("Neutrophils", Interaction) ~ "Neutrophils"
    )
  )


library(ggplot2)

interaction_colors <- c(
  "Macrophage → Treg/Th17"   = "#396544",
  "Treg/Th17 → Macrophage"   = "#73A681",
  "Neutrophils → Treg/Th17"  = "#E53637",
  "Treg/Th17 → Neutrophils"  = "#E98989"
)

p_line <- ggplot(
  cellchat_long,
  aes(
    x = group,
    y = weight,
    group = Interaction,
    color = Interaction
  )
) +

  geom_line(
    linewidth = 1
  ) +

  geom_point(
    size = 3.5
  ) +

  facet_wrap(
    ~ Myeloid,
    scales = "free_y"
  ) +

  scale_color_manual(
    values = interaction_colors
  ) +

  labs(
    x = NULL,
    y = "CellChat interaction strength",
    color = NULL
  ) +

  theme_classic(base_size = 14) +

  theme(
    legend.position = "top",
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    )
  )


tiff("communication_panel.tiff",units = 'in', width = 6, height = 4.5, res = 300)
p_line
dev.off()


tiff("communication_heatmap_pbs.tiff",units = 'in', width = 6, height = 4.5, res = 300)
netVisual_heatmap(
  cc_pbs_clean,
  measure = "weight",
  color.heatmap = "Reds",
  title.name = "PBS"
)
dev.off()

tiff("communication_heatmap_pd1.tiff",units = 'in', width = 6, height = 4.5, res = 300)
netVisual_heatmap(
  cc_pd1_clean,
  measure = "weight",
  color.heatmap = "Reds",
  title.name = "anti-PD-1"
)
dev.off()

tiff("communication_heatmap_ctl.tiff",units = 'in', width = 6, height = 4.5, res = 300)
netVisual_heatmap(
  cc_ctl_clean,
  measure = "weight",
  color.heatmap = "Reds",
  title.name = "anti-CTLA-4"
)
dev.off()



mat_pbs <- cc_pbs_clean@net$weight
mat_pd1 <- cc_pd1_clean@net$weight
mat_ctl <- cc_ctl_clean@net$weight


cell_order <- c(
  "Macrophage",
  "Neutrophils",
  "Treg and Th17"
)

mat_pbs <- mat_pbs[cell_order, cell_order]
mat_pd1 <- mat_pd1[cell_order, cell_order]
mat_ctl <- mat_ctl[cell_order, cell_order]


max_comm <- max(
  mat_pbs,
  mat_pd1,
  mat_ctl,
  na.rm = TRUE
)


library(ComplexHeatmap)
library(circlize)
library(grid)

col_fun <- colorRamp2(
  c(0, .2, .4),
  c("white", "#F4A6A6", "#B2182B")
)


ht_pbs <- Heatmap(
  mat_pbs,
  name = "Interaction\nstrength",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_title = "Sender",
  column_title = "PBS",
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10)
)


ht_pd1 <- Heatmap(
  mat_pd1,
  name = "Interaction\nstrength",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_title = NULL,
  column_title = "anti-PD-1",
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10)
)


ht_ctl <- Heatmap(
  mat_ctl,
  name = "Interaction\nstrength",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_title = NULL,
  column_title = "anti-CTLA-4",
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10)
)

tiff("communication_heatmap.tiff",units = 'in', width = 8, height = 4, res = 300)
draw(
  ht_pbs + ht_pd1 + ht_ctl,
  merge_legends = TRUE
)
dev.off()




























Idents(obj_cellchat_clean) <- "comm_celltype"

neut <- subset(
  obj_cellchat_clean,
  idents = "Neutrophils"
)

Idents(neut) <- "groups"

de_neut_pd1 <- FindMarkers(
  neut,
  ident.1 = "antipd1",
  ident.2 = "pbs",
  assay = "RNA",
  min.pct = 0.10,
  logfc.threshold = 0.25
)

de_neut_ctl <- FindMarkers(
  neut,
  ident.1 = "antictla4",
  ident.2 = "pbs",
  assay = "RNA",
  min.pct = 0.10,
  logfc.threshold = 0.25
)







library(clusterProfiler)
library(org.Mm.eg.db)

genes_pd1_neut <- rownames(
  de_neut_pd1[
    de_neut_pd1$p_val < 0.05 &
    de_neut_pd1$avg_log2FC > 0,
  ]
)

genes_ctl_neut <- rownames(
  de_neut_ctl[
    de_neut_ctl$p_val < 0.05 &
    de_neut_ctl$avg_log2FC > 0,
  ]
)

go_pd1_neut <- enrichGO(
  gene = genes_pd1_neut,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP"
)

go_ctl_neut <- enrichGO(
  gene = genes_ctl_neut,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP"
)


grep(
  "chemotaxis|migration|chemokine|neutrophil|granulocyte",
  go_pd1_neut@result$Description,
  ignore.case = TRUE,
  value = TRUE
)

grep(
  "chemotaxis|migration|chemokine|neutrophil|granulocyte",
  go_ctl_neut@result$Description,
  ignore.case = TRUE,
  value = TRUE
)


pathways_focus <- c(
  "neutrophil chemotaxis",
  "neutrophil migration",
  "granulocyte chemotaxis",
  "granulocyte migration",
  "leukocyte chemotaxis",
  "leukocyte migration involved in inflammatory response",
  "response to chemokine",
  "chemokine-mediated signaling pathway",
  "regulation of neutrophil chemotaxis",
  "regulation of neutrophil migration"
)


library(dplyr)
library(tidyr)
library(ggplot2)

go_pd1_df <- go_pd1_neut@result |>
  as.data.frame() |>
  filter(Description %in% pathways_focus) |>
  mutate(Treatment = "anti-PD-1")

go_ctl_df <- go_ctl_neut@result |>
  as.data.frame() |>
  filter(Description %in% pathways_focus) |>
  mutate(Treatment = "anti-CTLA-4")

go_neut_focus <- bind_rows(
  go_pd1_df,
  go_ctl_df
)


go_neut_focus <- go_neut_focus |>
  mutate(
    GeneRatio_num = sapply(
      GeneRatio,
      function(x) {
        z <- strsplit(x, "/")[[1]]
        as.numeric(z[1]) / as.numeric(z[2])
      }
    ),

    minus_log10_padj = -log10(p.adjust)
  )



p_go_neut <- ggplot(
  go_neut_focus,
  aes(
    x = Treatment,
    y = Description
  )
) +
  geom_point(
    aes(
      size = GeneRatio_num,
      color = minus_log10_padj
    )
  ) +
  scale_color_gradient(
    low = "#F4B6B6",
    high = "#B2182B"
  ) +
  scale_size(
    range = c(3, 10)
  ) +
  labs(
    x = NULL,
    y = NULL,
    size = "Gene ratio",
    color = expression(-log[10]("adjusted p-value"))
  ) +
  theme_classic(base_size = 18) +
  theme(
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 18
    )
  )

tiff("enrichment.tiff",units = 'in', width = 12, height = 6, res = 300)
p_go_neut
dev.off()







go_alluvial <- go_neut_focus |>
  dplyr::select(
    Treatment,
    Description,
    geneID,
    p.adjust
  ) |>
  separate_rows(
    geneID,
    sep = "/"
  )

colnames(go_alluvial)[colnames(go_alluvial) == "Description"] = "Pathway"
colnames(go_alluvial)[colnames(go_alluvial) == "geneID"] = "Gene"

gene_frequency <- go_alluvial |>
  count(Gene, sort = TRUE)

head(gene_frequency, 30)


top_genes <- gene_frequency |>
  slice_head(n = 20) |>
  pull(Gene)

go_alluvial_plot <- go_alluvial |>
  filter(Gene %in% top_genes)


go_alluvial_plot <- go_alluvial_plot |>
  distinct(
    Treatment,
    Pathway,
    Gene
  ) |>
  mutate(
    Weight = 1
  )



library(ggalluvial)

group_colors = c(
  "anti-PD-1"   = "#28331A", 
  "anti-CTLA-4" = "#9BC5F3"  
) 


p_alluvial <- ggplot(
  go_alluvial_plot,
  aes(
    axis1 = Pathway,
    axis2 = Gene,
    y = Weight
  )
) +
  geom_alluvium(
    aes(fill = Treatment),
    width = 0.08,
    alpha = 0.65,
    knot.pos = 0.45
  ) +
  geom_stratum(
    width = 0.10,
    fill = "white",
    color = "black",
    linewidth = 0.6
  ) +
  theme_minimal(base_size = 16) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  scale_x_discrete(limits = c("Description", "Gene"), expand = c(.1, .1)) +
  labs(x = "", y = "") +
  theme(
    legend.position = 'none'
  ) +
  scale_fill_manual(
    values = group_colors
  )



p_alluvial



pb = ggplot_build(p_alluvial)
strata_data = pb$data[[2]] 


tiff("enrichment_alluvial.tiff",units = 'in', width = 9, height = 5, res = 300)
p_alluvial + 
  geom_text_repel(
    data = strata_data |> filter(x == 1),
    aes(x = x, y = y, label = stratum),
    inherit.aes = FALSE,
    force = 1,
    hjust = 1,
    nudge_x = -1.4,
    min.segment.length = 1,
    segment.curvature = 2,
    box.padding = unit(.02, "lines"),
    direction = "y",
    size = 4,
    segment.size = .3,
    color = 'black'
  ) +
  geom_text_repel(
    data = strata_data |> filter(x == 2),
    aes(x = x, y = y, label = stratum),
    inherit.aes = FALSE,
    force = 1,
    nudge_x = .4,
    color = 'black',
    #nudge_y = 2,
    min.segment.length = 1,
    segment.curvature = 2,
    direction = "y",
    box.padding = unit(.02, "lines"),
    size = 4,
    segment.size = .3,
    max.overlaps = getOption("ggrepel.max.overlaps", default = 1000)
  ) 
dev.off()


