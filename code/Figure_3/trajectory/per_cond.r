
source("functions/patch/patch.r")
patch_monocle2()

group_colors = c(
  "naive"     = "#F06A8A", 
  "pbs"       = "#6F7F5A",  
  "igg"       = "#4C2A78", 
  "antipd1"   = "#28331A", 
  "antictla4" = "#9BC5F3"  
) 

seurat_rds = t_cells_CD4

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


SC_condition <- t_cells_CD4

DefaultAssay(SC_condition) <- "RNA"
Idents(SC_condition) <- "groups"

table(Idents(SC_condition))

de_pbs_naive <- FindMarkers(
  SC_condition,
  ident.1 = "pbs",
  ident.2 = "naive",
  assay = "RNA",
  layer = "data",
  min.pct = 0.10,
  logfc.threshold = 0.10,
  only.pos = FALSE
)

de_pd1_pbs <- FindMarkers(
  SC_condition,
  ident.1 = "antipd1",
  ident.2 = "pbs",
  assay = "RNA",
  layer = "data",
  min.pct = 0.10,
  logfc.threshold = 0.10,
  only.pos = FALSE
)

de_ctla4_pbs <- FindMarkers(
  SC_condition,
  ident.1 = "antictla4",
  ident.2 = "pbs",
  assay = "RNA",
  layer = "data",
  min.pct = 0.10,
  logfc.threshold = 0.10,
  only.pos = FALSE
)

de_igg_pbs <- FindMarkers(
  SC_condition,
  ident.1 = "igg",
  ident.2 = "pbs",
  assay = "RNA",
  layer = "data",
  min.pct = 0.10,
  logfc.threshold = 0.10,
  only.pos = FALSE
)


select_de_genes <- function(
    de_table,
    n_up = 50,
    n_down = 50,
    padj_cutoff = 0.05,
    logfc_cutoff = 0.25
) {

  de_df <- de_table |>
    as.data.frame() |>
    tibble::rownames_to_column("gene")

  fc_column <- if ("avg_log2FC" %in% colnames(de_df)) {
    "avg_log2FC"
  } else if ("avg_logFC" %in% colnames(de_df)) {
    "avg_logFC"
  } else {
    stop("Não foi encontrada a coluna avg_log2FC ou avg_logFC.")
  }

  # Filtrar genes significativos
  de_filtered <- de_df |>
    filter(
      !is.na(p_val),
      p_val < padj_cutoff,
      abs(.data[[fc_column]]) > logfc_cutoff
    ) |>
    filter(
      !grepl(
        "^mt-|^Rpl|^Rps|^Hba|^Hbb|^Igh|^Igk|^Igl",
        gene,
        ignore.case = TRUE
      )
    )

  # Genes aumentados
  genes_up <- de_filtered |>
    filter(.data[[fc_column]] > 0) |>
    arrange(desc(.data[[fc_column]])) |>
    slice_head(n = n_up) |>
    pull(gene)

  # Genes reduzidos
  genes_down <- de_filtered |>
    filter(.data[[fc_column]] < 0) |>
    arrange(.data[[fc_column]]) |>
    slice_head(n = n_down) |>
    pull(gene)

  unique(c(genes_up, genes_down))
}


genes_infection <- select_de_genes(
  de_pbs_naive,
  n_up = 50,
  n_down = 50
)

genes_pd1 <- select_de_genes(
  de_pd1_pbs,
  n_up = 40,
  n_down = 40
)

genes_ctla4 <- select_de_genes(
  de_ctla4_pbs,
  n_up = 40,
  n_down = 40
)

genes_igg <- select_de_genes(
  de_igg_pbs,
  n_up = 40,
  n_down = 40
)


ordering_genes_condition <- unique(c(
  genes_infection,
  genes_pd1,
  genes_ctla4,
  genes_igg
))

ordering_genes_condition <- intersect(
  genes_ctla4,
  rownames(cds)
)

length(ordering_genes_condition)


cds_condition <- cds

cds_condition <- setOrderingFilter(
  cds_condition,
  ordering_genes_condition
)


set.seed(1234)

cds_condition <- reduceDimension(
  cds_condition,
  max_components = 2,
  method = "DDRTree"
)

cds_condition <- orderCells(cds_condition)

pData(cds_condition)$groups <- factor(
  pData(cds_condition)$groups,
  levels = names(group_colors)
)

saveRDS(cds_condition, "./data/cds_condition_PBS_naive.rds")
saveRDS(cds_condition, "./data/cds_condition_IGG_pbs.rds")
saveRDS(cds_condition, "./data/cds_condition_PD1_pbs.rds")
saveRDS(cds_condition, "./data/cds_condition_CTLA4_pbs.rds")


p_pseudotime <- monocle::plot_cell_trajectory(
  cds_condition,
  color_by = "Pseudotime"
) +
  scale_color_viridis_c(option = "magma", direction = 1) +
  theme_classic(base_size = 14) +
  labs(color = "Pseudotime")

p_pseudotime

ggsave(
  filename = "igg_pbs.tiff",
  plot = p_pseudotime,
  device = ragg::agg_tiff,
  width = 6,
  height = 5,
  units = "in",
  res = 600,
  compression = "lzw",
  background = "white"
)


p_state <- monocle::plot_cell_trajectory(
  cds_condition,
  color_by = "State"
) +
  #scale_color_manual(values = c("Healthy" = "#1F5C5C", "CD" = "#E18487")) +
  theme_classic(base_size = 14) +
  labs(color = "State")

p_state


ggsave(
  filename = "igg_pbs_p_state.tiff",
  plot = p_state,
  device = ragg::agg_tiff,
  width = 6,
  height = 5,
  units = "in",
  res = 600,
  compression = "lzw",
  background = "white"
)















###### HEATMAP

genes_all <- unique(c(
  genes_infection,
  genes_igg,
  genes_pd1,
  genes_ctla4
))

genes_all <- intersect(
  genes_all,
  rownames(t_cells_CD4)
)

genes_heatmap <- genes_all[
  !grepl(
    "^Gm[0-9]|-ps[0-9]*$|^Rpl|^Rps|^mt-",
    genes_all,
    ignore.case = TRUE
  )
]

length(genes_all)
length(genes_heatmap)


genes_heatmap <- setdiff(
  genes_heatmap,
  c(
    "Gapdh", "Actb", "Tmsb4x", "Hprt",
    "Pgk1", "Ogt", "Uba1", "Hspa8"
  )
)


t_cells_CD4$groups_plot <- factor(
  t_cells_CD4$groups,
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


avg_expr <- AggregateExpression(
  object = t_cells_CD4,
  assays = "RNA",
  features = genes_heatmap,
  group.by = "groups_plot",
  return.seurat = FALSE
)$RNA


avg_expr <- AverageExpression(
  object = t_cells_CD4,
  assays = "RNA",
  features = genes_heatmap,
  group.by = "groups_plot",
  slot = "data",
  return.seurat = FALSE
)$RNA




library(Matrix)

rna_data <- GetAssayData(
  object = t_cells_CD4,
  assay = "RNA",
  layer = "data"
)

genes_heatmap <- intersect(
  genes_heatmap,
  rownames(rna_data)
)

groups <- droplevels(t_cells_CD4$groups_plot)

avg_expr <- sapply(
  levels(groups),
  function(grp) {
    cells_group <- colnames(t_cells_CD4)[groups == grp]

    Matrix::rowMeans(
      rna_data[
        genes_heatmap,
        cells_group,
        drop = FALSE
      ]
    )
  }
)

rownames(avg_expr) <- genes_heatmap
colnames(avg_expr) <- levels(groups)


avg_expr_scaled <- t(scale(t(avg_expr)))

avg_expr_scaled <- avg_expr_scaled[
  apply(
    avg_expr_scaled,
    1,
    function(x) all(is.finite(x))
  ),
  ,
  drop = FALSE
]


library(ComplexHeatmap)
library(circlize)
library(grid)


col_fun <- colorRamp2(
    c(-1,0,1),
    c("#0D47A1","white","#CF4173")
)


group_cols <- c(
    "Naive"      = "#F06A8A",
    "PBS"        = "#6F7F5A",
    "IgG"        = "#4C2A78",
    "anti-PD-1"  = "#28331A",
    "anti-CTLA-4"= "#9BC5F3"
)

ha <- HeatmapAnnotation(
    Group = colnames(avg_expr_scaled),
    col = list(
        Group = group_cols
    ),
    show_annotation_name = FALSE
)

ht <- Heatmap(
    avg_expr_scaled,
    name = "Z-score",
    col = col_fun,
    #column_split = colnames(avg_expr_scaled),
    top_annotation = ha,
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(
        fontsize = 10,
        fontface = "bold"
    ),
    column_title = NULL,
    row_title = NULL,
    heatmap_legend_param = list(
        title = "Expression"
    )
)

tiff("heatmap.tiff",
  width = 6,
  height = 5,
  units = "in",
  res = 600, )
draw(ht)
dev.off()


ggsave(
  filename = "heatmap.tiff",
  plot = cc,
  device = ragg::agg_tiff,
  width = 6,
  height = 5,
  units = "in",
  res = 600,
  compression = "lzw",
  background = "white"
)








#### pSEUDO BULK


library(Seurat)
library(DESeq2)
library(dplyr)
library(tibble)

DefaultAssay(t_cells_CD4) <- "RNA"

pb_global <- AggregateExpression(
  t_cells_CD4,
  assays = "RNA",
  group.by = c("cartig", "groups"),
  return.seurat = FALSE,
  slot = "counts"
)

pb_counts <- pb_global$RNA



t_cells_CD4$pb_id <- paste(
  t_cells_CD4$cartig,
  t_cells_CD4$groups,
  sep = "__"
)

pb_global <- AggregateExpression(
  t_cells_CD4,
  assays = "RNA",
  group.by = "pb_id",
  return.seurat = FALSE,
  slot = "counts"
)

pb_counts <- pb_global$RNA

pb_counts[1:3,1:3]

pb_metadata <- t_cells_CD4@meta.data |>
  dplyr::select(pb_id, cartig, groups) |>
  dplyr::distinct() |>
  dplyr::mutate(
    pb_id_counts = gsub("_", "-", pb_id)
  )

pb_metadata <- pb_metadata[
  match(colnames(pb_counts), pb_metadata$pb_id_counts),
  ,
  drop = FALSE
]

rownames(pb_metadata) <- pb_metadata$pb_id_counts





keep_genes <- rowSums(pb_counts >= 10) >= 2

pb_counts_filtered <- pb_counts[
  keep_genes,
  ,
  drop = FALSE
]


dds <- DESeqDataSetFromMatrix(
  countData = round(pb_counts_filtered),
  colData = pb_metadata,
  design = ~ groups
)

dds$groups <- factor(
  dds$groups,
  levels = c(
    "naive",
    "pbs",
    "igg",
    "antipd1",
    "antictla4"
  )
)

dds <- DESeq(dds)


res_pbs_naive <- results(
  dds,
  contrast = c(
    "groups",
    "pbs",
    "naive"
  ),
  alpha = 0.05
) |>
  as.data.frame() |>
  rownames_to_column("gene") |>
  arrange(padj)


res_igg_pbs <- results(
  dds,
  contrast = c(
    "groups",
    "igg",
    "pbs"
  ),
  alpha = 0.05
) |>
  as.data.frame() |>
  rownames_to_column("gene") |>
  arrange(padj)



res_pd1_pbs <- results(
  dds,
  contrast = c(
    "groups",
    "antipd1",
    "pbs"
  ),
  alpha = 0.05
) |>
  as.data.frame() |>
  rownames_to_column("gene") |>
  arrange(padj)


res_ctla4_pbs <- results(
  dds,
  contrast = c(
    "groups",
    "antictla4",
    "pbs"
  ),
  alpha = 0.05
) |>
  as.data.frame() |>
  rownames_to_column("gene") |>
  arrange(padj)

head(res_ctla4_pbs)
