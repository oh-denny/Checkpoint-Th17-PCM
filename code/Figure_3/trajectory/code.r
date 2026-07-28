
library(Seurat)
library(monocle)
library(dplyr)
library(ggplot2)

source("functions/patch/patch.r")
patch_monocle2()

seurat_rds <- t_cells_CD4

SC_obj_sub = t_cells_CD4

DefaultAssay(SC_obj_sub) <- "RNA"
Idents(SC_obj_sub) <- "seurat_clusters"

counts_matrix = GetAssayData(
  SC_obj_sub,
  assay = "RNA",
  layer = "counts"
)

cell_metadata = SC_obj_sub@meta.data[
  colnames(counts_matrix),
  ,
  drop = FALSE
]

pd = new(
  "AnnotatedDataFrame",
  data = cell_metadata
)

gene_metadata = data.frame(
  gene_short_name = rownames(counts_matrix),
  row.names = rownames(counts_matrix),
  stringsAsFactors = FALSE
)

fd = new(
  "AnnotatedDataFrame",
  data = gene_metadata
)

cds = newCellDataSet(
  counts_matrix,
  phenoData = pd,
  featureData = fd,
  lowerDetectionLimit = 1,
  expressionFamily = VGAM::negbinomial.size()
)

patch_monocle2()

markers <- FindAllMarkers(
  SC_obj_sub,
  assay = "RNA",
  slot = "data",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  recorrect_umi = FALSE
)

markers_filtered = markers |>
  filter(
    p_val_adj < 0.01,
    avg_log2FC > 0.5,
    pct.1 >= 0.25
  ) |>
  filter(
    !grepl("^mt-", gene, ignore.case = TRUE),
    !grepl("^Rpl|^Rps", gene),
    !grepl("^Hba|^Hbb", gene),
    !grepl("^Igh|^Igk|^Igl", gene),
    !grepl("^Gm[0-9]", gene)
  )

ordering_genes_markers = markers_filtered |>
  group_by(cluster) |>
  slice_max(
    order_by = avg_log2FC,
    n = 50,
    with_ties = FALSE
  ) |>
  pull(gene) |>
  unique()

ordering_genes_markers <- intersect(
  ordering_genes_markers,
  rownames(cds)
)
#ordering_genes_markers <- intersect(ordering_genes_markers, rownames(cds))

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- setOrderingFilter(cds, ordering_genes_markers)

set.seed(1234)

cds <- reduceDimension(
  cds,
  max_components = 2,
  method = "DDRTree"
)

patch_monocle2()
cds <- monocle::orderCells(cds)

saveRDS(cds, "data/seu_obj/monocle2_cds_ordered.rds")

table(
  pData(cds)$State,
  pData(cds)$seurat_clusters
)


cols <- c(
  "0" = "#A61C4B",
  "1" = "#0F6C84",
  "2" = "#2C6B6B",
  "3" = "#AEB8B7",
  "4" = "#6C7A3A"
)

pData(cds)$seurat_clusters <- factor(
  pData(cds)$seurat_clusters,
  levels = names(cols)
)

p_cluster <- monocle::plot_cell_trajectory(
  cds,
  color_by = "seurat_clusters"
) +
  scale_color_manual(values = cols) +
  theme_classic(base_size = 14) +
  labs(color = "Cluster")

p_cluster

library(viridis)

p_pseudotime <- monocle::plot_cell_trajectory(
  cds,
  color_by = "Pseudotime"
) +
  scale_color_viridis_c(option = "magma", direction = 1) +
  theme_classic(base_size = 14) +
  labs(color = "Pseudotime")

p_pseudotime



p_group <- monocle::plot_cell_trajectory(
  cds,
  color_by = "groups"
) +
  #scale_color_manual(values = c("Healthy" = "#1F5C5C", "CD" = "#E18487")) +
  theme_classic(base_size = 14) +
  labs(color = "groups")

p_group


p_state <- monocle::plot_cell_trajectory(
  cds,
  color_by = "State"
) +
  #scale_color_manual(values = c("Healthy" = "#1F5C5C", "CD" = "#E18487")) +
  theme_classic(base_size = 14) +
  labs(color = "State")

p_state




monocle::plot_cell_trajectory(cds, color_by = "seurat_clusters")
monocle::plot_cell_trajectory(cds, color_by = "Pseudotime")
monocle::plot_cell_trajectory(cds, color_by = "State")

genes_to_plot <- intersect(c("HIF1A", "FOXP3"), rownames(cds))

if ("HIF1A" %in% genes_to_plot) {
  monocle::plot_genes_in_pseudotime(cds["HIF1A", ], color_by = "State")
  monocle::plot_genes_jitter(cds["HIF1A", ], grouping = "State")
}

if ("FOXP3" %in% genes_to_plot) {
  monocle::plot_genes_in_pseudotime(cds["FOXP3", ], color_by = "State")
}

if ("RORgt_Treg_Score1" %in% colnames(SC_obj_sub@meta.data)) {
  pData(cds)$RORgt_Treg_Score1 <- SC_obj_sub$RORgt_Treg_Score1[colnames(cds)]
  monocle::plot_cell_trajectory(cds, color_by = "RORgt_Treg_Score1")
}



tab_state_group <- table(
  State = pData(cds)$State,
  Group = pData(cds)$groups
)

prop_state_within_group <- prop.table(
  tab_state_group,
  margin = 2
)

round(prop_state_within_group * 100, 2)


library(dplyr)
library(ggplot2)

state_df <- pData(cds) |>
  as.data.frame() |>
  count(groups, State, name = "n") |>
  group_by(groups) |>
  mutate(
    proportion = n / sum(n),
    percent = proportion * 100
  ) |>
  ungroup()

ggplot(
  state_df,
  aes(
    x = groups,
    y = percent,
    fill = factor(State)
  )
) +
  geom_col(position = "fill") +
  scale_y_continuous(
    labels = scales::percent
  ) +
  labs(
    x = NULL,
    y = "Proportion of cells",
    fill = "Monocle State"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )




SC_obj_sub$MonocleState <- factor(
    pData(cds)$State[
        match(colnames(SC_obj_sub), rownames(pData(cds)))
    ]
)


tab <- table(
    Cluster = SC_obj_sub$seurat_clusters,
    State = SC_obj_sub$MonocleState
)

tab
round(
    prop.table(tab, margin = 2) * 100,
    1
)


Idents(SC_obj_sub) <- "MonocleState"


genes <- c(

# Naive
"Tcf7",
"Lef1",
"Ccr7",
"Sell",
"Il7r",

# Proliferação
"Mki67",
"Top2a",
"Birc5",
"Cenpf",

# Treg
"Foxp3",
"Il2ra",
"Ctla4",
"Tigit",

# Th17
"Rorc",
"Il17a",
"Il17f",
"Ccr6",
"Rora",

# Ativação
"Icos",
"Cd44",
"Cd69",

# Th1
"Ifng",
"Tbx21",

# Hipóxia
"Hif1a"
)

genes <- intersect(
    genes,
    rownames(SC_obj_sub)
)

DotPlot(
    SC_obj_sub,
    features = genes
) +
RotatedAxis()


FeaturePlot(
    SC_obj_sub,
    features=c(
        "Tcf7",
        "Foxp3",
        "Rorc",
        "Mki67"
    ),
    order=TRUE
)
