




cds_pbs = readRDS("./data/cds_condition_PBS_naive.rds")
cds_igg = readRDS("./data/cds_condition_IGG_pbs.rds")
cds_pd1 = readRDS("./data/cds_condition_PD1_pbs.rds")
cds_ctl = readRDS("./data/cds_condition_CTLA4_pbs.rds")


tab_pd1 <- table(
  State = pData(cds_pd1)$State,
  Group = pData(cds_pd1)$groups
)

round(
  prop.table(tab_pd1, margin = 2) * 100,
  2
)

tab_ctla4 <- table(
  State = pData(cds_ctl)$State,
  Group = pData(cds_ctl)$groups
)

round(
  prop.table(tab_ctla4, margin = 2) * 100,
  2
)

state_pd1 <- setNames(
  as.character(pData(cds_pd1)$State),
  rownames(pData(cds_pd1))
)

t_cells_CD4$PD1_State <- state_pd1[colnames(t_cells_CD4)]

state_ctla4 <- setNames(
  as.character(pData(cds_ctl)$State),
  rownames(pData(cds_ctl))
)

t_cells_CD4$CTLA4_State <- state_ctla4[colnames(t_cells_CD4)]

head(t_cells_CD4)


Idents(t_cells_CD4) <- "PD1_State"

markers_pd1_state <- FindAllMarkers(
  t_cells_CD4,
  assay = "RNA",
  layer = "data",
  only.pos = TRUE,
  min.pct = 0.15,
  logfc.threshold = 0.25
)

markers_pd1_state <- markers_pd1_state |>
  filter(p_val < 0.05) |>
  arrange(desc(avg_log2FC))


Idents(t_cells_CD4) <- "CTLA4_State"

markers_ctla4_state <- FindAllMarkers(
  t_cells_CD4,
  assay = "RNA",
  layer = "data",
  only.pos = TRUE,
  min.pct = 0.15,
  logfc.threshold = 0.25
)

markers_ctla4_state <- markers_ctla4_state |>
  filter(p_val < 0.05) |>
  arrange(desc(avg_log2FC))

genes_pd1_state <- markers_pd1_state |>
  slice_head(n = 100) |>
  pull(gene)

genes_ctla4_state <- markers_ctla4_state |>
  slice_head(n = 100) |>
  pull(gene)




markers_pd1_state_bio <-
    markers_pd1_state |>
    filter(
        !grepl("^Gm",gene),
        !grepl("-ps",gene),
        !grepl("^Rpl",gene),
        !grepl("^Rps",gene)
    )

lista_pd1 <- split(
  markers_pd1_state_bio$gene,
  markers_pd1_state_bio$cluster
)

markers_ctla4_state_bio <-
    markers_ctla4_state |>
    filter(
        !grepl("^Gm",gene),
        !grepl("-ps",gene),
        !grepl("^Rpl",gene),
        !grepl("^Rps",gene)
    )

lista_ctl <- split(
  markers_ctla4_state_bio$gene,
  markers_ctla4_state_bio$cluster
)

library(clusterProfiler)
library(org.Mm.eg.db)
bp_res = compareCluster(
  geneCluster = lista_pd1,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP"
)
dotplot(bp_res)

bp_ctl = compareCluster(
  geneCluster = lista_ctl,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP"
)
dotplot(bp_ctl)


res_enrch = bp_res@compareClusterResult |> data.frame()
res_enrch1 = bp_ctl@compareClusterResult |> data.frame()

class(res_enrch1)

res_enrch$group = 'Anti-PD-1'
res_enrch1$group = 'Anti-CTLA-4'

res = rbind(res_enrch, res_enrch1)

openxlsx::write.xlsx(res, "~/Desktop/res.xlsx")


### alguma figura X

Th17 <- c(
    "Rorc",
    "Maf",
    "Ccr6",
    "Il17a",
    "Il17f",
    "Batf",
    "Rora",
    "Il23r",
    "Ahr",
    "Irf4"
    )


Th17_present <- intersect(
  Th17,
  rownames(t_cells_CD4)
)

Th17_missing <- setdiff(
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

FeaturePlot(
  t_cells_CD4,
  features = "Th17_Score1",
  order = TRUE,
  split.by = 'groups'
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
feature_plot
dev.off()





#### sankey plot

library(dplyr)
library(ggplot2)
library(stringr)
library(forcats)

go_ctl <- as.data.frame(bp_ctl)

head(go_ctl)
colnames(go_ctl)


ratio_to_numeric <- function(x) {
  vapply(
    strsplit(as.character(x), "/"),
    function(z) as.numeric(z[1]) / as.numeric(z[2]),
    numeric(1)
  )
}

go_ctl <- go_ctl |>
  mutate(
    State = paste0("State ", Cluster),
    GeneRatio_numeric = ratio_to_numeric(GeneRatio),
    Description = str_wrap(Description, width = 45)
  )

go_ctl_selected <- go_ctl |>
  filter(
    p.adjust < 0.05,
    str_detect(
      str_to_lower(Description),
      paste0(
        "t cell activation|",
        "t cell differentiation|",
        "lymphocyte activation|",
        "lymphocyte differentiation|",
        "regulatory t cell|",
        "interleukin-17|",
        "il-17|",
        "th17|",
        "nf-kappa|",
        "nf-kappab|",
        "tumor necrosis factor|",
        "inflammatory response|",
        "chemotaxis|",
        "neutrophil|",
        "granulocyte|",
        "macrophage|",
        "monocyte|",
        "myeloid|",
        "antigen processing|",
        "antigen presentation"
      )
    )
  )


go_ctl_plot <- go_ctl_selected |>
  group_by(State) |>
  slice_min(
    order_by = p.adjust,
    n = 6,
    with_ties = FALSE
  ) |>
  ungroup()




p_go_ctl <- ggplot(
  go_ctl_plot,
  aes(
    x = State,
    y = fct_reorder(Description, GeneRatio_numeric),
    size = GeneRatio_numeric,
    color = p.adjust
  )
) +
  geom_point() +
  scale_color_gradient(
    low = "#B2182B",
    high = "#2166AC",
    trans = "reverse"
  ) +
  scale_size(
    range = c(3, 10)
  ) +
  theme_classic(base_size = 13) +
  labs(
    x = "Anti-CTLA-4 trajectory state",
    y = NULL,
    size = "Gene ratio",
    color = "Adjusted p-value"
  ) +
  theme(
    axis.text.y = element_text(size = 10),
    legend.position = "right"
  )

p_go_ctl





pd1_states <- setNames(
  as.character(pData(cds_pd1)$State),
  rownames(pData(cds_pd1))
)

t_cells_CD4$PD1_State <- pd1_states[colnames(t_cells_CD4)]


ctl_states <- setNames(
  as.character(pData(cds_ctl)$State),
  rownames(pData(cds_ctl))
)

t_cells_CD4$CTLA4_State <- ctl_states[colnames(t_cells_CD4)]


t_cells_CD4$Trajectory_State <- NA_character_

# Anti-PD-1: ajuste os números conforme sua orientação real
t_cells_CD4$Trajectory_State[
  t_cells_CD4$groups == "antipd1" &
  t_cells_CD4$PD1_State == "1"
] <- "anti-PD-1 | Early"

t_cells_CD4$Trajectory_State[
  t_cells_CD4$groups == "antipd1" &
  t_cells_CD4$PD1_State == "3"
] <- "anti-PD-1 | Intermediate"

t_cells_CD4$Trajectory_State[
  t_cells_CD4$groups == "antipd1" &
  t_cells_CD4$PD1_State == "5"
] <- "anti-PD-1 | Late"

# Anti-CTLA-4
t_cells_CD4$Trajectory_State[
  t_cells_CD4$groups == "antictla4" &
  t_cells_CD4$CTLA4_State == "1"
] <- "anti-CTLA-4 | Early"

t_cells_CD4$Trajectory_State[
  t_cells_CD4$groups == "antictla4" &
  t_cells_CD4$CTLA4_State == "2"
] <- "anti-CTLA-4 | Intermediate"

t_cells_CD4$Trajectory_State[
  t_cells_CD4$groups == "antictla4" &
  t_cells_CD4$CTLA4_State == "3"
] <- "anti-CTLA-4 | Late"


table(
  t_cells_CD4$Trajectory_State,
  useNA = "ifany"
)




trajectory_cells <- subset(
  t_cells_CD4,
  subset = !is.na(Trajectory_State)
)

trajectory_cells$Trajectory_State <- factor(
  trajectory_cells$Trajectory_State,
  levels = c(
    "anti-PD-1 | Early",
    "anti-PD-1 | Intermediate",
    "anti-PD-1 | Late",
    "anti-CTLA-4 | Early",
    "anti-CTLA-4 | Intermediate",
    "anti-CTLA-4 | Late"
  )
)

genes_trajectory <- list(
  "Less differentiated" = c(
    "Tcf7", "Lef1", "Sell", "Klf2", "S1pr1"
  ),
  "T-cell activation" = c(
    "Cd28", "Icos", "Cd44", "Cd69", "Dusp10", "Ccl5"
  ),
  "Th17" = c(
    "Rorc", "Maf", "Ccr6", "Il17a",
    "Il17f", "Rora", "Hif1a"
  ),
  "Treg" = c(
    "Foxp3", "Il2ra", "Ctla4", "Tigit", "Ikzf2"
  ),
  "Proteostasis" = c(
    "Hspa1a", "Hspa1b", "Hsp90aa1", "Dnajb1", "Hspe1"
  ),
  "Antigen presentation" = c(
    "H2-Aa", "H2-Eb1", "Cd74", "Cd79a"
  )
)

genes_trajectory <- lapply(
  genes_trajectory,
  intersect,
  y = rownames(trajectory_cells)
)





p_dot <- DotPlot(
  trajectory_cells,
  features = genes_trajectory,
  group.by = "Trajectory_State",
  assay = "RNA",
  dot.scale = 10,
  scale = TRUE
) +
  scale_size(
    range = c(0.7, 10),
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100)
  ) +
  scale_color_gradient2(
    low = "grey90",
    mid = "white",
    high = "#B2182B",
    midpoint = 0
  ) +
  theme_classic(base_size = 18) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(
      face = "italic"
    )
  )

p_dot



p_dot <- DotPlot(
  trajectory_cells,
  features = genes_trajectory,
  group.by = "Trajectory_State",
  assay = "RNA",
  dot.scale = 10,
  scale = TRUE
) +
  scale_size(
    range = c(0.7, 10),
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100)
  ) +
  scale_color_gradient2(
    low = "grey90",
    mid = "white",
    high = "#B2182B",
    midpoint = 0
  ) +
  labs(
    x = NULL,
    y = NULL,
    color = "Average expression",
    size = "Percent expressed"
  ) +
  theme_classic(base_size = 11)

































library(Seurat)
library(monocle)
library(dplyr)
library(ggplot2)

# ============================================================
# 1. Função para extrair Early, Intermediate e Late
# ============================================================

extract_trajectory_stage <- function(
    cds,
    group_name,
    group_label,
    reverse_pseudotime = FALSE
) {

  metadata <- pData(cds) |>
    as.data.frame() |>
    tibble::rownames_to_column("cell") |>
    dplyr::filter(
      groups == group_name,
      !is.na(Pseudotime)
    )

  if (nrow(metadata) == 0) {
    stop(
      "Nenhuma célula encontrada para o grupo '",
      group_name,
      "' no CDS informado."
    )
  }

  if (reverse_pseudotime) {
    metadata$Pseudotime_use <-
      max(metadata$Pseudotime, na.rm = TRUE) -
      metadata$Pseudotime
  } else {
    metadata$Pseudotime_use <- metadata$Pseudotime
  }

  # Dividir o pseudotempo em três partes com números
  # aproximadamente semelhantes de células
  metadata <- metadata |>
    dplyr::mutate(
      Trajectory_stage = dplyr::ntile(
        Pseudotime_use,
        3
      ),

      Trajectory_stage = dplyr::case_when(
        Trajectory_stage == 1 ~ "Early",
        Trajectory_stage == 2 ~ "Intermediate",
        Trajectory_stage == 3 ~ "Late"
      ),

      Trajectory_State = paste(
        group_label,
        Trajectory_stage,
        sep = " | "
      )
    )

  metadata |>
    dplyr::select(
      cell,
      Trajectory_State,
      Trajectory_stage,
      Pseudotime = Pseudotime_use
    )
}




states_pbs <- extract_trajectory_stage(
  cds = cds_condition_PBS_naive,
  group_name = "pbs",
  group_label = "PBS",
  reverse_pseudotime = FALSE
)

states_igg <- extract_trajectory_stage(
  cds = cds_condition_IGG_pbs,
  group_name = "igg",
  group_label = "IgG",
  reverse_pseudotime = FALSE
)

states_pd1 <- extract_trajectory_stage(
  cds = cds_condition_PD1_pbs,
  group_name = "antipd1",
  group_label = "anti-PD-1",
  reverse_pseudotime = FALSE
)

states_ctla4 <- extract_trajectory_stage(
  cds = cds_condition_CTLA4_pbs,
  group_name = "antictla4",
  group_label = "anti-CTLA-4",
  reverse_pseudotime = FALSE
)

trajectory_metadata <- bind_rows(
  states_pbs,
  states_igg,
  states_pd1,
  states_ctla4
)




# ============================================================
# Transferir os estados para o Seurat de forma robusta
# ============================================================

# Garantir que cada célula apareça apenas uma vez
duplicated_cells <- trajectory_metadata$cell[
  duplicated(trajectory_metadata$cell)
]

if (length(duplicated_cells) > 0) {
  stop(
    "Existem células duplicadas em trajectory_metadata: ",
    paste(head(unique(duplicated_cells), 10), collapse = ", ")
  )
}

# Ordem das células no objeto Seurat
seurat_cells <- colnames(t_cells_CD4)

# Localizar cada célula do Seurat em trajectory_metadata
cell_match <- match(
  seurat_cells,
  trajectory_metadata$cell
)

# Criar metadado com exatamente a mesma ordem do Seurat
trajectory_metadata_seurat <- data.frame(

  Trajectory_State = trajectory_metadata$Trajectory_State[
    cell_match
  ],

  Trajectory_stage = trajectory_metadata$Trajectory_stage[
    cell_match
  ],

  Trajectory_Pseudotime = trajectory_metadata$Pseudotime[
    cell_match
  ],

  row.names = seurat_cells,
  stringsAsFactors = FALSE
)

# Remover colunas antigas apenas quando existirem
columns_to_remove <- intersect(
  c(
    "Trajectory_State",
    "Trajectory_stage",
    "Trajectory_Pseudotime"
  ),
  colnames(t_cells_CD4@meta.data)
)

if (length(columns_to_remove) > 0) {
  t_cells_CD4@meta.data[
    ,
    columns_to_remove
  ] <- NULL
}

# Adicionar o novo metadado
t_cells_CD4 <- AddMetaData(
  object = t_cells_CD4,
  metadata = trajectory_metadata_seurat
)




trajectory_state_levels <- c(
  "PBS | Early",
  "PBS | Intermediate",
  "PBS | Late",

  "IgG | Early",
  "IgG | Intermediate",
  "IgG | Late",

  "anti-PD-1 | Early",
  "anti-PD-1 | Intermediate",
  "anti-PD-1 | Late",

  "anti-CTLA-4 | Early",
  "anti-CTLA-4 | Intermediate",
  "anti-CTLA-4 | Late"
)

cells_trajectory <- rownames(
  t_cells_CD4@meta.data
)[
  !is.na(t_cells_CD4$Trajectory_State)
]


t_cells_CD4$Trajectory_State <- factor(
  t_cells_CD4$Trajectory_State,
  levels = trajectory_state_levels
)

trajectory_cells <- subset(
  t_cells_CD4,
  cells = cells_trajectory
)


table(
  trajectory_cells$Trajectory_State,
  useNA = "ifany"
)




genes_trajectory <- list(

  "Less differentiated" = c(
    "Lef1",
    "Klf3",
    "Sell",
    "Foxp1",
    "S1pr1",
    "Klf2",
    "Ccr7",
    "Tob1",
    "Tcf7",
    "Tob2",
    "Zfp36l2",
    "Btg2"
  ),

  "T-cell activation" = c(
    "Tnfrsf1b",
    "Icos",
    "Cyld",
    "Nr4a1",
    "Itk",
    "Tnfaip3",
    "Tnfsf8",
    "Cd28",
    "Zap70",
    "Lat"
  ),

  "Migration" = c(
    "Rgs1",
    "Cxcr4",
    "Ccr2",
    "Spn",
    "Adam8",
    "Itgb1",
    "Itga4",
    "S1pr4",
    "Rac2",
    "Ripor2",
    "Cxcr6",
    "S100a4",
    "Ccl5"
  ),

  "Th17" = c(
    "Il17a",
    "Maf",
    "Hif1a",
    "Bhlhe40"
  ),

  "Treg" = c(
    "Ctla4",
    "Il2ra",
    "Tmem176a",
    "Tmem176b",
    "Tgfbr2"
  )
)



all_dotplot_genes <- unlist(
  genes_trajectory,
  use.names = FALSE
)

duplicated_dotplot_genes <- unique(
  all_dotplot_genes[
    duplicated(all_dotplot_genes)
  ]
)

duplicated_dotplot_genes



gene_frequency <- table(all_dotplot_genes)

gene_frequency[
  gene_frequency > 1
]


genes_trajectory_present <- lapply(
  genes_trajectory,
  function(genes) {
    intersect(
      genes,
      rownames(trajectory_cells)
    )
  }
)



genes_present_flat <- unlist(
  genes_trajectory_present,
  use.names = FALSE
)

if (anyDuplicated(genes_present_flat) > 0) {

  duplicated_genes <- unique(
    genes_present_flat[
      duplicated(genes_present_flat)
    ]
  )

  stop(
    "Genes repetidos entre classificações: ",
    paste(duplicated_genes, collapse = ", ")
  )
}

trajectory_cells$Trajectory_Group <- sub(
  " \\| .*",
  "",
  as.character(trajectory_cells$Trajectory_State)
)

trajectory_cells$Trajectory_Stage <- sub(
  ".*\\| ",
  "",
  as.character(trajectory_cells$Trajectory_State)
)

trajectory_cells$Trajectory_Group <- factor(
  trajectory_cells$Trajectory_Group,
  levels = c(
    "PBS",
    "IgG",
    "anti-PD-1",
    "anti-CTLA-4"
  )
)

trajectory_cells$Trajectory_Stage <- factor(
  trajectory_cells$Trajectory_Stage,
  levels = c(
    "Early",
    "Intermediate",
    "Late"
  )
)


DefaultAssay(trajectory_cells) <- "RNA"

p_dot <- DotPlot(
  object = trajectory_cells,
  features = genes_trajectory_present,
  group.by = "Trajectory_State",
  assay = "RNA",
  dot.scale = 9,
  scale = TRUE
)


p_dot$data <- p_dot$data |>
  dplyr::filter(
    avg.exp.scaled > .5
  )


p_dot$data <- p_dot$data |>
  dplyr::mutate(
    Trajectory_Group = sub(
      " \\| .*",
      "",
      as.character(id)
    ),
    Trajectory_Stage = sub(
      ".*\\| ",
      "",
      as.character(id)
    ),
    Trajectory_Group = factor(
      Trajectory_Group,
      levels = c(
        "PBS",
        "IgG",
        "anti-PD-1",
        "anti-CTLA-4"
      )
    ),
    Trajectory_Stage = factor(
      Trajectory_Stage,
      levels = c(
        "Late",
        "Intermediate",
        "Early"
      )
    )
  )



p_dot <- p_dot +
  aes(y = Trajectory_Stage) +
  facet_grid(
    rows = vars(Trajectory_Group),
    cols = vars(feature.groups),
    scales = "free",
    space = "free",
    switch = "y"
  ) +
  scale_size(
    range = c(0.5, 13),
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100)
  ) +
  scale_color_gradient(
    low = "white",
    high = "#B2182B",
    limits = c(1, 2),
    breaks = c(1, 1.5, 2),
    oob = scales::squish,
    name = "Average Expression"
  ) +
  labs(
    x = NULL,
    y = NULL,
    color = "Average Expression",
    size = "Percent Expressed"
  ) +
  theme_classic(base_size = 15) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(
      face = "italic",
      color = "black"
    ),
    strip.background = element_rect(
      fill = "white",
      color = "black"
    ),
    strip.text.x = element_text(
      face = "bold"
    ),
    strip.text.y.left = element_text(
      face = "bold",
      angle = 0
    ),
    strip.placement = "outside",
    panel.spacing.y = grid::unit(
      0.8,
      "lines"
    ),
    legend.position = "right"
  )

p_dot


ggsave(
  filename = "./p_dot.tiff",
  plot = p_dot,
  device = ragg::agg_tiff,
  width = 20,
  height = 4,
  units = "in",
  res = 600,
  compression = "lzw",
  background = "white"
)
