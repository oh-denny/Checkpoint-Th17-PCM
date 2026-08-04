
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


cds_condition_pbs = readRDS("./data/cds_condition_PBS_naive.rds")
cds_condition_igg = readRDS("./data/cds_condition_IGG_pbs.rds")
cds_condition_pd1 = readRDS("./data/cds_condition_PD1_pbs.rds")
cds_condition_ctl = readRDS("./data/cds_condition_CTLA4_pbs.rds")



p_pseudotime <- monocle::plot_cell_trajectory(
  cds_condition_ctl,
  color_by = "Pseudotime"
) +
  scale_color_viridis_c(option = "magma", direction = 1) +
  theme_classic(base_size = 20) +
  labs(color = "Pseudotime") +
  theme(
    legend.position = 'bottom'
  )

p_pseudotime

ggsave(
  filename = "ctl_pbs.tiff",
  plot = p_pseudotime,
  device = ragg::agg_tiff,
  width = 6,
  height = 7,
  units = "in",
  res = 600,
  compression = "lzw",
  background = "white"
)


p_state <- monocle::plot_cell_trajectory(
  cds_condition_ctl,
  color_by = "State"
) +
  #scale_color_manual(values = c("Healthy" = "#1F5C5C", "CD" = "#E18487")) +
  theme_classic(base_size = 20) +
  labs(color = "State") +
   theme(
    legend.position = 'bottom'
  )

p_state


ggsave(
  filename = "ctl_pbs_p_state.tiff",
  plot = p_state,
  device = ragg::agg_tiff,
  width = 6,
  height = 8,
  units = "in",
  res = 600,
  compression = "lzw",
  background = "white"
)















### enriquecimento

library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(dplyr)
library(ggplot2)

gene_lists <- list(
  Infection   = unique(genes_infection),
  IgG         = unique(genes_igg),
  `anti-PD-1` = unique(genes_pd1),
  `anti-CTLA-4` = unique(genes_ctla4)
)

all_symbols <- unique(unlist(gene_lists))



cc_go <- compareCluster(
  geneCluster = gene_lists,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  minGSSize = 10,
  maxGSSize = 500,
  readable = TRUE
)

dotplot(
  cc_go,
  showCategory = 8,
  includeAll = TRUE
) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )



openxlsx::write.xlsx(cc_go@compareClusterResult, "~/Desktop/cc_go.xlsx")



pathways_infection <- c(
  "T cell differentiation",
  "regulation of T cell activation",
  "chemokine-mediated signaling pathway",
  "T cell migration",
  "regulation of inflammatory response"
)

pathways_igg <- c(
  "protein folding",
  "regulation of leukocyte differentiation",
  "cellular response to interleukin-2",
  "lymphocyte migration"
)

pathways_pd1 <- c(
  "protein folding",
  "T cell differentiation",
  "regulation of T cell activation",
  "glycolytic process",
  "leukocyte migration"
)

pathways_ctla4 <- c(
  "regulation of T cell activation",
  "regulation of leukocyte differentiation",
  "leukocyte cell-cell adhesion",
  "T cell migration",
  "response to chemokine",
  "glycolytic process",
  "glutamine metabolic process",
  "protein folding"
)
selected_pathways <- list(
  Infection = pathways_infection,
  IgG = pathways_igg,
  `anti-PD-1` = pathways_pd1,
  `anti-CTLA-4` = pathways_ctla4
)

selected_pathways <- list(
  Infection = pathways_infection,
  IgG = pathways_igg,
  `anti-PD-1` = pathways_pd1,
  `anti-CTLA-4` = pathways_ctla4
)



library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(forcats)
library(circlize)
library(ComplexHeatmap)

cc_df <- as.data.frame(cc_go)

unique(cc_df$Cluster)




selected_pathways <- list(

  Infection = c(
    "T cell differentiation",
    "regulation of T cell activation",
    "chemokine-mediated signaling pathway",
    "T cell migration",
    "regulation of inflammatory response"
  ),

  IgG = c(
    "protein folding",
    "regulation of leukocyte differentiation",
    "cellular response to interleukin-2",
    "lymphocyte migration"
  ),

  `anti-PD-1` = c(
    "protein folding",
    "T cell differentiation",
    "regulation of T cell activation",
    "glycolytic process",
    "leukocyte migration"
  ),

  `anti-CTLA-4` = c(
    "regulation of T cell activation",
    "regulation of leukocyte differentiation",
    "leukocyte cell-cell adhesion",
    "T cell migration",
    "response to chemokine",
    "glycolytic process",
    "glutamine metabolic process",
    "protein folding"
  )
)


selected_table <- bind_rows(
  lapply(names(selected_pathways), function(contrast_name) {

    pathways_contrast <- selected_pathways[[contrast_name]]

    cc_df |>
      filter(Cluster == contrast_name) |>
      filter(
        str_to_lower(Description) %in%
          str_to_lower(pathways_contrast)
      )
  })
)


selected_table <- selected_table |>
  mutate(
    GeneRatio_numeric = sapply(
      strsplit(GeneRatio, "/"),
      function(x) as.numeric(x[1]) / as.numeric(x[2])
    ),
    minus_log10_padj = -log10(p.adjust),
    Description = str_to_sentence(Description)
  )


selected_table <- selected_table |>
  arrange(Cluster, GeneRatio_numeric) |>
  mutate(
    pathway_plot = paste(Description, Cluster, sep = "___")
  )



p_dotplot <- ggplot(
  selected_table,
  aes(
    x = Cluster,
    y = reorder(Description, Cluster)
  )
) +
  geom_point(
    aes(
      size = Count,
      color = minus_log10_padj
    ),
    alpha = 0.9
  ) +

  theme_classic(base_size = 13) +
  theme(
    axis.text.y = element_text(
      size = 10,
      color = "black"
    ),
    axis.text.x = element_text(
      color = "black"
    ),
    legend.position = "right",
    panel.spacing = grid::unit(1.2, "lines")
  )

p_dotplot




plot_pathway_gene_circos <- function(
    pathway_gene_df,
    contrast_name,
    max_genes_per_pathway = 15,
    filename = NULL
) {

  df_contrast <- pathway_gene_df |>
    filter(Contrast == contrast_name) |>
    group_by(Pathway) |>
    slice_head(n = max_genes_per_pathway) |>
    ungroup() |>
    distinct(Pathway, Gene)

  if (nrow(df_contrast) == 0) {
    stop(
      paste(
        "Nenhuma via ou gene encontrado para:",
        contrast_name
      )
    )
  }

  pathways <- unique(df_contrast$Pathway)
  genes <- unique(df_contrast$Gene)

  pathway_colors <- setNames(
    grDevices::hcl.colors(
      length(pathways),
      palette = "Dark 3"
    ),
    pathways
  )

  gene_colors <- setNames(
    rep("#D9D9D9", length(genes)),
    genes
  )

  grid_colors <- c(
    pathway_colors,
    gene_colors
  )

  if (!is.null(filename)) {
    pdf(
      filename,
      width = 12,
      height = 12
    )
  }

  circos.clear()

  circos.par(
    start.degree = 90,
    gap.degree = 2,
    track.margin = c(0.01, 0.01),
    cell.padding = c(0, 0, 0, 0)
  )

  chordDiagram(
    x = df_contrast |>
      dplyr::select(Pathway, Gene),
    grid.col = grid_colors,
    transparency = 0.35,
    directional = 0,
    annotationTrack = "grid",
    preAllocateTracks = 1,
    link.sort = TRUE,
    link.decreasing = FALSE
  )

  circos.trackPlotRegion(
    track.index = 1,
    bg.border = NA,
    panel.fun = function(x, y) {

      sector_name <- CELL_META$sector.index

      circos.text(
        x = CELL_META$xcenter,
        y = CELL_META$ylim[1] + 0.1,
        labels = sector_name,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5),
        cex = ifelse(
          sector_name %in% pathways,
          0.7,
          0.6
        ),
        font = ifelse(
          sector_name %in% pathways,
          2,
          1
        )
      )
    }
  )

  title(
    main = contrast_name,
    cex.main = 1.6,
    font.main = 2
  )

  circos.clear()

  if (!is.null(filename)) {
    dev.off()
  }

  invisible(df_contrast)
}

pathway_gene_table <- selected_table |>
  dplyr::select(
    Cluster,
    Description,
    geneID
  ) |>
  separate_rows(
    geneID,
    sep = "/"
  ) |>
  dplyr::rename(
    Contrast = Cluster,
    Pathway = Description,
    Gene = geneID
  ) |>
  distinct()

circos_infection <- plot_pathway_gene_circos(
  pathway_gene_df = pathway_gene_table,
  contrast_name = "Infection",
  max_genes_per_pathway = 15,
  filename = "Circos_Infection.pdf"
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































library(clusterProfiler)
library(dplyr)
library(tidyr)
library(stringr)
library(forcats)
library(ggplot2)
library(scales)
library(tibble)

# Converter resultado para data.frame
go_table <- as.data.frame(cc_go)

if (nrow(go_table) == 0) {
  stop(
    "Nenhum termo GO foi identificado pelo compareCluster."
  )
}

required_go_columns <- c(
  "Cluster",
  "ID",
  "Description",
  "GeneRatio",
  "p.adjust",
  "geneID",
  "Count"
)

missing_go_columns <- setdiff(
  required_go_columns,
  colnames(go_table)
)

if (length(missing_go_columns) > 0) {
  stop(
    "As seguintes colunas não foram encontradas no resultado: ",
    paste(missing_go_columns, collapse = ", ")
  )
}


# ============================================================
# 4. FUNÇÃO PARA CONVERTER GENE RATIO
# ============================================================

ratio_to_numeric <- function(x) {

  x <- as.character(x)

  vapply(
    strsplit(x, "/", fixed = TRUE),
    function(z) {

      if (length(z) != 2) {
        return(NA_real_)
      }

      numerator <- suppressWarnings(as.numeric(z[1]))
      denominator <- suppressWarnings(as.numeric(z[2]))

      if (
        is.na(numerator) ||
        is.na(denominator) ||
        denominator == 0
      ) {
        return(NA_real_)
      }

      numerator / denominator
    },
    numeric(1)
  )
}


# ============================================================
# 5. PREPARAR TABELA GO
# ============================================================

go_table <- go_table |>
  dplyr::mutate(

    Cluster = as.character(.data$Cluster),
    Description = as.character(.data$Description),

    Contrast = dplyr::recode(
      .data$Cluster,
      "Infection" = "PBS vs Naive",
      "IgG" = "IgG vs PBS",
      "anti-PD-1" = "Anti-PD-1 vs PBS",
      "anti-CTLA-4" = "Anti-CTLA-4 vs PBS",
      .default = .data$Cluster
    ),

    description_lower = stringr::str_to_lower(
      .data$Description
    ),

    GeneRatio_numeric = ratio_to_numeric(
      .data$GeneRatio
    ),

    minus_log10_padj = -log10(
      pmax(
        as.numeric(.data$p.adjust),
        .Machine$double.xmin
      )
    )
  )


# ============================================================
# 6. CLASSIFICAR TERMOS GO NOS CINCO PROGRAMAS
#
# A ordem do case_when é importante.
# Primeiro entram os termos mais específicos.
# ============================================================

go_classified <- go_table |>
  dplyr::mutate(

    Program = dplyr::case_when(

      # ------------------------------------------------------
      # TREG-ASSOCIATED
      # Termos mais específicos devem vir primeiro
      # ------------------------------------------------------

      stringr::str_detect(
        .data$description_lower,
        "regulatory t cell differentiation"
      ) ~ "Treg-associated",

      stringr::str_detect(
        .data$description_lower,
        "regulation of regulatory t cell"
      ) ~ "Treg-associated",

      stringr::str_detect(
        .data$description_lower,
        "immune tolerance"
      ) ~ "Treg-associated",

      stringr::str_detect(
        .data$description_lower,
        "negative regulation of immune response"
      ) ~ "Treg-associated",

      # ------------------------------------------------------
      # TH17-ASSOCIATED
      # ------------------------------------------------------

      stringr::str_detect(
        .data$description_lower,
        "th17|interleukin-17|il-17"
      ) ~ "Th17-associated",

      .data$description_lower %in% c(
        "regulation of inflammatory response",
        "positive regulation of inflammatory response",
        "inflammatory response"
      ) ~ "Th17-associated",

      # ------------------------------------------------------
      # T-CELL ACTIVATION
      # ------------------------------------------------------

      stringr::str_detect(
        .data$description_lower,
        "regulation of t cell activation"
      ) ~ "T-cell activation",

      stringr::str_detect(
        .data$description_lower,
        "positive regulation of t cell activation"
      ) ~ "T-cell activation",

      stringr::str_detect(
        .data$description_lower,
        "negative regulation of t cell activation"
      ) ~ "T-cell activation",

      stringr::str_detect(
        .data$description_lower,
        "alpha-beta t cell activation"
      ) ~ "T-cell activation",

      stringr::str_detect(
        .data$description_lower,
        "t cell receptor signaling"
      ) ~ "T-cell activation",

      stringr::str_detect(
        .data$description_lower,
        "antigen receptor-mediated signaling"
      ) ~ "T-cell activation",

      stringr::str_detect(
        .data$description_lower,
        "lymphocyte activation"
      ) ~ "T-cell activation",

      # ------------------------------------------------------
      # MIGRATION AND ADHESION
      # ------------------------------------------------------

      stringr::str_detect(
        .data$description_lower,
        "t cell migration"
      ) ~ "Migration and adhesion",

      stringr::str_detect(
        .data$description_lower,
        "lymphocyte migration"
      ) ~ "Migration and adhesion",

      stringr::str_detect(
        .data$description_lower,
        "leukocyte migration"
      ) ~ "Migration and adhesion",

      stringr::str_detect(
        .data$description_lower,
        "chemokine-mediated signaling"
      ) ~ "Migration and adhesion",

      stringr::str_detect(
        .data$description_lower,
        "chemotaxis"
      ) ~ "Migration and adhesion",

      stringr::str_detect(
        .data$description_lower,
        "lamellipodium"
      ) ~ "Migration and adhesion",

      stringr::str_detect(
        .data$description_lower,
        "cell adhesion"
      ) ~ "Migration and adhesion",

      stringr::str_detect(
        .data$description_lower,
        "integrin-mediated"
      ) ~ "Migration and adhesion",

      stringr::str_detect(
        .data$description_lower,
        "actin filament organization"
      ) ~ "Migration and adhesion",

      # ------------------------------------------------------
      # LESS DIFFERENTIATED
      # Termos gerais de diferenciação ficam por último
      # ------------------------------------------------------

      .data$description_lower %in% c(
        "t cell differentiation",
        "lymphocyte differentiation"
      ) ~ "Less differentiated",

      stringr::str_detect(
        .data$description_lower,
        "^regulation of t cell differentiation$"
      ) ~ "Less differentiated",

      stringr::str_detect(
        .data$description_lower,
        "^regulation of lymphocyte differentiation$"
      ) ~ "Less differentiated",

      TRUE ~ NA_character_
    )
  )







# ============================================================
# 7. CONFERIR A CLASSIFICAÇÃO
# ============================================================

program_count <- go_classified |>
  dplyr::filter(!is.na(.data$Program)) |>
  dplyr::group_by(.data$Program) |>
  dplyr::summarise(
    N_terms = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::arrange(
    dplyr::desc(.data$N_terms)
  )

print(program_count)


classified_terms_table <- go_classified |>
  dplyr::filter(
    !is.na(.data$Program)
  ) |>
  dplyr::select(
    dplyr::all_of(
      c(
        "Contrast",
        "Program",
        "Description",
        "Count",
        "GeneRatio",
        "p.adjust",
        "geneID"
      )
    )
  ) |>
  dplyr::arrange(
    .data$Program,
    .data$Contrast,
    .data$p.adjust
  )

print(
  classified_terms_table,
  n = nrow(classified_terms_table)
)


# ============================================================
# 8. DEFINIR ORDEM DOS PROGRAMAS E CONTRASTES
# ============================================================

program_levels <- c(
  "Less differentiated",
  "T-cell activation",
  "Migration and adhesion",
  "Th17-associated",
  "Treg-associated"
)

contrast_levels <- c(
  "PBS vs Naive",
  "IgG vs PBS",
  "Anti-PD-1 vs PBS",
  "Anti-CTLA-4 vs PBS"
)


# ============================================================
# 9. SELECIONAR OS MELHORES TERMOS
#
# Atualmente:
# - p.adjust < 0.05
# - até 2 termos por programa e contraste
# ============================================================

n_terms_per_program <- 2
padj_cutoff <- 0.05

go_selected <- go_classified |>
  dplyr::filter(
    !is.na(.data$Program),
    !is.na(.data$p.adjust),
    .data$p.adjust < padj_cutoff
  ) |>
  dplyr::group_by(
    .data$Contrast,
    .data$Program
  ) |>
  dplyr::arrange(
    .data$p.adjust,
    dplyr::desc(.data$Count),
    .by_group = TRUE
  ) |>
  dplyr::slice_head(
    n = n_terms_per_program
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(

    Program = factor(
      .data$Program,
      levels = program_levels
    ),

    Contrast = factor(
      .data$Contrast,
      levels = contrast_levels
    ),

    Description_wrapped = stringr::str_wrap(
      .data$Description,
      width = 42
    ),

    minus_log10_padj = -log10(
      pmax(
        as.numeric(.data$p.adjust),
        .Machine$double.xmin
      )
    )
  )


if (nrow(go_selected) == 0) {
  stop(
    "Nenhum termo permaneceu depois da classificação e do filtro."
  )
}


# Conferir tabela final
selected_terms_table <- go_selected |>
  dplyr::select(
    dplyr::all_of(
      c(
        "Contrast",
        "Program",
        "Description",
        "Count",
        "GeneRatio",
        "p.adjust",
        "geneID"
      )
    )
  ) |>
  dplyr::arrange(
    .data$Program,
    .data$Contrast,
    .data$p.adjust
  )

print(
  selected_terms_table,
  n = nrow(selected_terms_table)
)


# Salvar tabela
utils::write.csv(
  selected_terms_table,
  file = "selected_GO_terms_trajectory_programs.csv",
  row.names = FALSE
)


# ============================================================
# 10. DOTPLOT PRINCIPAL
#
# Eixo X = contraste
# Eixo Y = termo GO
# Tamanho = número de genes
# Cor = -log10(p.adjust)
# ============================================================

p_go_dotplot <- ggplot2::ggplot(
  go_selected,
  ggplot2::aes(
    x = .data$Contrast,
    y = .data$Description_wrapped
  )
) +
  ggplot2::geom_point(
    ggplot2::aes(
      size = .data$Count,
      color = .data$minus_log10_padj
    ),
    alpha = 0.90
  ) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(.data$Program),
    scales = "free_y",
    space = "free_y",
    switch = "y",
    drop = TRUE
  ) +
  ggplot2::scale_size_continuous(
    name = "Gene count",
    range = c(2.5, 8)
  ) +
  ggplot2::scale_color_viridis_c(
    name = expression(-log[10]("adjusted p-value")),
    option = "magma",
    direction = 1
  ) +
  ggplot2::labs(
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_bw(
    base_size = 12
  ) +
  ggplot2::theme(

    panel.grid.major = ggplot2::element_line(
      linewidth = 0.25,
      color = "grey88"
    ),

    panel.grid.minor = ggplot2::element_blank(),

    strip.background = ggplot2::element_rect(
      fill = "grey95",
      color = "grey75"
    ),

    strip.placement = "outside",

    strip.text.y.left = ggplot2::element_text(
      angle = 0,
      face = "bold",
      size = 10
    ),

    axis.text.x = ggplot2::element_text(
      angle = 35,
      hjust = 1,
      vjust = 1,
      face = "bold",
      color = "black"
    ),

    axis.text.y = ggplot2::element_text(
      size = 9,
      color = "black"
    ),

    legend.position = "right",

    panel.spacing.y = grid::unit(
      0.6,
      "lines"
    ),

    plot.margin = ggplot2::margin(
      t = 6,
      r = 10,
      b = 6,
      l = 6
    )
  )

# ============================================================
# 11. DOTPLOT HORIZONTAL E COMPACTO
#
# Mais adequado para espaços horizontais na figura
# ============================================================

p_go_dotplot_horizontal <- ggplot2::ggplot(
  go_selected,
  ggplot2::aes(
    x = .data$Description_wrapped,
    y = .data$Contrast
  )
) +
  ggplot2::geom_point(
    ggplot2::aes(
      size = .data$Count,
      color = .data$minus_log10_padj
    ),
    alpha = 0.90
  ) +
  ggplot2::facet_wrap(
    facets = ggplot2::vars(.data$Program),
    scales = "free_x",
    nrow = 1,
    drop = TRUE
  ) +
  ggplot2::scale_size_continuous(
    name = "Gene count",
    range = c(2.5, 7)
  ) +
  ggplot2::scale_color_viridis_c(
    name = expression(-log[10]("adjusted p-value")),
    option = "magma",
    direction = 1
  ) +
  ggplot2::labs(
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_bw(
    base_size = 11
  ) +
  ggplot2::theme(

    strip.background = ggplot2::element_blank(),

    strip.text = ggplot2::element_text(
      face = "bold",
      size = 10
    ),

    axis.text.x = ggplot2::element_text(
      angle = 55,
      hjust = 1,
      vjust = 1,
      size = 8,
      color = "black"
    ),

    axis.text.y = ggplot2::element_text(
      face = "bold",
      size = 9,
      color = "black"
    ),

    panel.grid.minor = ggplot2::element_blank(),

    panel.spacing.x = grid::unit(
      1,
      "lines"
    ),

    legend.position = "right"
  )

print(p_go_dotplot_horizontal)














# ============================================================
# HEATMAP RESUMIDO DE PROGRAMAS BIOLÓGICOS
#
# Agrupa termos GO semanticamente semelhantes.
#
# Linhas:
#   Less differentiated
#   T-cell activation
#   Migration and adhesion
#   Th17-associated
#   Treg-associated
#
# Colunas:
#   PBS vs Naive
#   IgG vs PBS
#   Anti-PD-1 vs PBS
#   Anti-CTLA-4 vs PBS
#
# Cor:
#   melhor -log10(p.adjust) dentro do programa
#
# Número:
#   quantidade de genes únicos associados ao programa
# ============================================================


# ============================================================
# 0. PACOTES
# ============================================================

library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(tibble)


# ============================================================
# 1. EXTRAIR RESULTADO DO compareCluster
# ============================================================

go_table <- as.data.frame(cc_go)

if (nrow(go_table) == 0) {
  stop("O objeto cc_go não contém termos enriquecidos.")
}


# ============================================================
# 2. PREPARAR TABELA
# ============================================================

go_table <- go_table |>
  dplyr::mutate(

    Cluster = as.character(.data$Cluster),

    Description = as.character(
      .data$Description
    ),

    description_lower = stringr::str_to_lower(
      .data$Description
    ),

    Contrast = dplyr::recode(
      .data$Cluster,

      "Infection" = "PBS vs Naive",
      "IgG" = "IgG vs PBS",
      "anti-PD-1" = "Anti-PD-1 vs PBS",
      "anti-CTLA-4" = "Anti-CTLA-4 vs PBS",

      .default = .data$Cluster
    ),

    p.adjust = as.numeric(
      .data$p.adjust
    )
  )


# ============================================================
# 3. CLASSIFICAR CADA TERMO EM UM PROGRAMA AMPLO
#
# Importante:
# "regulation of inflammatory response" e
# "positive regulation of inflammatory response"
# serão agrupados no mesmo programa.
# ============================================================

go_classified <- go_table |>
  dplyr::mutate(

    Program = dplyr::case_when(

      # ------------------------------------------------------
      # TREG-ASSOCIATED
      # Termos específicos primeiro
      # ------------------------------------------------------

      stringr::str_detect(
        .data$description_lower,
        paste(
          c(
            "regulatory t cell",
            "regulatory t-cell",
            "immune tolerance",
            "tolerance induction",
            "negative regulation of immune response",
            "negative regulation of t cell activation",
            "negative regulation of lymphocyte activation",
            "suppression of immune response"
          ),
          collapse = "|"
        )
      ) ~ "Treg-associated",

      # ------------------------------------------------------
      # TH17 / INFLAMMATORY
      #
      # Todas as versões de inflammatory response entram juntas:
      # - inflammatory response
      # - regulation of inflammatory response
      # - positive regulation of inflammatory response
      # ------------------------------------------------------

      stringr::str_detect(
        .data$description_lower,
        paste(
          c(
            "th17",
            "interleukin-17",
            "il-17",
            "inflammatory response",
            "inflammatory cytokine",
            "response to interleukin-17",
            "cellular response to interleukin-17"
          ),
          collapse = "|"
        )
      ) ~ "Th17-associated",

      # ------------------------------------------------------
      # T-CELL ACTIVATION
      # ------------------------------------------------------

      stringr::str_detect(
        .data$description_lower,
        paste(
          c(
            "t cell activation",
            "t-cell activation",
            "alpha-beta t cell activation",
            "lymphocyte activation",
            "antigen receptor-mediated signaling",
            "t cell receptor signaling",
            "t-cell receptor signaling",
            "immune response-activating",
            "costimulation",
            "co-stimulation"
          ),
          collapse = "|"
        )
      ) ~ "T-cell activation",

      # ------------------------------------------------------
      # MIGRATION AND ADHESION
      # ------------------------------------------------------

      stringr::str_detect(
        .data$description_lower,
        paste(
          c(
            "t cell migration",
            "t-cell migration",
            "lymphocyte migration",
            "leukocyte migration",
            "cell migration",
            "chemotaxis",
            "chemokine-mediated signaling",
            "cell-cell adhesion",
            "cell adhesion",
            "integrin-mediated",
            "lamellipodium",
            "cell motility",
            "actin filament organization",
            "actin cytoskeleton organization"
          ),
          collapse = "|"
        )
      ) ~ "Migration and adhesion",

      # ------------------------------------------------------
      # LESS DIFFERENTIATED
      #
      # Termos gerais de diferenciação ficam por último.
      # Isso impede que "regulatory T cell differentiation"
      # seja erroneamente colocado aqui.
      # ------------------------------------------------------

      stringr::str_detect(
        .data$description_lower,
        paste(
          c(
            "^t cell differentiation$",
            "^lymphocyte differentiation$",
            "^regulation of t cell differentiation$",
            "^regulation of lymphocyte differentiation$",
            "^positive regulation of t cell differentiation$",
            "^negative regulation of t cell differentiation$"
          ),
          collapse = "|"
        )
      ) ~ "Less differentiated",

      TRUE ~ NA_character_
    )
  )


# ============================================================
# 4. MANTER SOMENTE TERMOS CLASSIFICADOS E SIGNIFICATIVOS
# ============================================================

go_classified <- go_classified |>
  dplyr::filter(
    !is.na(.data$Program),
    !is.na(.data$p.adjust),
    .data$p.adjust < 0.05
  )


if (nrow(go_classified) == 0) {
  stop(
    paste(
      "Nenhum termo permaneceu após a classificação.",
      "Revise os padrões usados no case_when()."
    )
  )
}


# ============================================================
# 5. FUNÇÃO PARA EXTRAIR GENES ÚNICOS
# ============================================================

extract_unique_genes <- function(gene_strings) {

  gene_strings <- gene_strings[
    !is.na(gene_strings) &
      nzchar(gene_strings)
  ]

  if (length(gene_strings) == 0) {
    return(character(0))
  }

  genes <- unlist(
    strsplit(
      gene_strings,
      split = "/",
      fixed = TRUE
    ),
    use.names = FALSE
  )

  genes <- genes[
    !is.na(genes) &
      nzchar(genes)
  ]

  unique(genes)
}


# ============================================================
# 6. RESUMIR CADA PROGRAMA EM CADA CONTRASTE
#
# best_padj:
#   menor p.adjust entre os termos do programa
#
# Best_term:
#   termo mais significativo dentro do programa
#
# N_unique_genes:
#   união dos genes de todos os termos daquele programa
#
# N_GO_terms:
#   quantidade de termos agrupados
# ============================================================

program_summary <- go_classified |>
  dplyr::group_by(
    .data$Contrast,
    .data$Program
  ) |>
  dplyr::summarise(

    best_padj = min(
      .data$p.adjust,
      na.rm = TRUE
    ),

    Best_term = .data$Description[
      which.min(.data$p.adjust)
    ][1],

    N_GO_terms = dplyr::n_distinct(
      .data$ID
    ),

    Genes = paste(
      extract_unique_genes(
        .data$geneID
      ),
      collapse = "/"
    ),

    N_unique_genes = length(
      extract_unique_genes(
        .data$geneID
      )
    ),

    .groups = "drop"
  ) |>
  dplyr::mutate(

    minus_log10_padj = -log10(
      pmax(
        .data$best_padj,
        .Machine$double.xmin
      )
    )
  )


# ============================================================
# 7. DEFINIR ORDEM
# ============================================================

program_levels <- c(
  "Less differentiated",
  "T-cell activation",
  "Migration and adhesion",
  "Th17-associated",
  "Treg-associated"
)

contrast_levels <- c(
  "PBS vs Naive",
  "IgG vs PBS",
  "Anti-PD-1 vs PBS",
  "Anti-CTLA-4 vs PBS"
)


# ============================================================
# 8. COMPLETAR COMBINAÇÕES AUSENTES
#
# Combinações sem termo significativo receberão:
# - score = 0
# - genes = 0
# ============================================================

program_summary_complete <- tidyr::expand_grid(

  Program = program_levels,

  Contrast = contrast_levels

) |>
  dplyr::left_join(

    program_summary,

    by = c(
      "Program",
      "Contrast"
    )

  ) |>
  dplyr::mutate(

    best_padj = tidyr::replace_na(
      .data$best_padj,
      1
    ),

    minus_log10_padj = tidyr::replace_na(
      .data$minus_log10_padj,
      0
    ),

    N_unique_genes = tidyr::replace_na(
      .data$N_unique_genes,
      0L
    ),

    N_GO_terms = tidyr::replace_na(
      .data$N_GO_terms,
      0L
    ),

    Best_term = tidyr::replace_na(
      .data$Best_term,
      "Not enriched"
    ),

    Genes = tidyr::replace_na(
      .data$Genes,
      ""
    ),

    Program = factor(
      .data$Program,
      levels = rev(program_levels)
    ),

    Contrast = factor(
      .data$Contrast,
      levels = contrast_levels
    )
  )


# ============================================================
# 9. INSPECIONAR RESULTADO
# ============================================================

print(
  program_summary_complete |>
    dplyr::select(
      dplyr::all_of(
        c(
          "Contrast",
          "Program",
          "Best_term",
          "N_GO_terms",
          "N_unique_genes",
          "best_padj"
        )
      )
    ),
  n = Inf
)


# ============================================================
# 10. HEATMAP RESUMIDO
#
# Agora há somente cinco linhas.
#
# Cor:
#   melhor significância dentro do programa
#
# Número:
#   genes únicos associados ao programa
# ============================================================

p_program_heatmap <- ggplot2::ggplot(

  program_summary_complete,

  ggplot2::aes(
    x = .data$Contrast,
    y = .data$Program,
    fill = .data$minus_log10_padj
  )

) +

  ggplot2::geom_tile(
    color = "white",
    linewidth = 1
  ) +

  ggplot2::geom_text(

    ggplot2::aes(
      label = dplyr::if_else(
        .data$N_unique_genes > 0,
        as.character(.data$N_unique_genes),
        ""
      )
    ),

    size = 4,
    fontface = "bold"
  ) +

  ggplot2::scale_fill_viridis_c(

    name = expression(
      -log[10]("adjusted p-value")
    ),

    option = "magma",

    limits = c(
      0,
      max(
        program_summary_complete$minus_log10_padj,
        na.rm = TRUE
      )
    )
  ) +

  ggplot2::labs(
    x = NULL,
    y = NULL
  ) +

  ggplot2::theme_classic(
    base_size = 13
  ) +

  ggplot2::theme(

    axis.text.x = ggplot2::element_text(
      angle = 40,
      hjust = 1,
      vjust = 1,
      face = "bold",
      color = "black"
    ),

    axis.text.y = ggplot2::element_text(
      face = "bold",
      color = "black"
    ),

    axis.ticks = ggplot2::element_blank(),

    legend.position = "right",

    plot.margin = ggplot2::margin(
      8,
      8,
      8,
      8
    )
  )

print(p_program_heatmap)


# ============================================================
# 11. SALVAR HEATMAP
# ============================================================

ggplot2::ggsave(

  filename = "GO_programs_summary_heatmap.tiff",

  plot = p_program_heatmap,

  device = ragg::agg_tiff,

  width = 7.5,
  height = 4.8,
  units = "in",

  res = 600,

  compression = "lzw",

  background = "white"
)


ggplot2::ggsave(

  filename = "GO_programs_summary_heatmap.pdf",

  plot = p_program_heatmap,

  width = 7.5,
  height = 4.8,
  units = "in",

  device = grDevices::cairo_pdf
)


# ============================================================
# 12. DOTPLOT RESUMIDO
#
# Também mostra apenas os cinco programas.
#
# Tamanho:
#   número de genes únicos
#
# Cor:
#   melhor significância
# ============================================================

p_program_dotplot <- ggplot2::ggplot(

  program_summary_complete |>
    dplyr::filter(
      .data$N_unique_genes > 0
    ),

  ggplot2::aes(
    x = .data$Contrast,
    y = .data$Program
  )

) +

  ggplot2::geom_point(

    ggplot2::aes(
      size = .data$N_unique_genes,
      color = .data$minus_log10_padj
    ),

    alpha = 0.9
  ) +

  ggplot2::scale_size_continuous(

    name = "Unique genes",

    range = c(
      3,
      11
    )
  ) +

  ggplot2::scale_color_viridis_c(

    name = expression(
      -log[10]("adjusted p-value")
    ),

    option = "magma"
  ) +

  ggplot2::labs(
    x = NULL,
    y = NULL
  ) +

  ggplot2::theme_classic(
    base_size = 13
  ) +

  ggplot2::theme(

    axis.text.x = ggplot2::element_text(
      angle = 40,
      hjust = 1,
      vjust = 1,
      face = "bold",
      color = "black"
    ),

    axis.text.y = ggplot2::element_text(
      face = "bold",
      color = "black"
    ),

    axis.line = ggplot2::element_blank(),

    axis.ticks = ggplot2::element_blank(),

    legend.position = "right"
  )

print(p_program_dotplot)


# ============================================================
# 13. SALVAR DOTPLOT
# ============================================================

ggplot2::ggsave(

  filename = "GO_programs_summary_dotplot.tiff",

  plot = p_program_dotplot,

  device = ragg::agg_tiff,

  width = 7.5,
  height = 4.8,
  units = "in",

  res = 600,

  compression = "lzw",

  background = "white"
)


# ============================================================
# 14. SALVAR TABELAS DE SUPORTE
# ============================================================

utils::write.csv(

  go_classified,

  file = "GO_terms_classified_by_program.csv",

  row.names = FALSE
)


utils::write.csv(

  program_summary_complete,

  file = "GO_programs_summary.csv",

  row.names = FALSE
)


saveRDS(

  program_summary_complete,

  file = "GO_programs_summary.rds"
)























# Criar uma coluna indicando se Il17a está presente no termo GO
go_table <- go_table |>
  dplyr::mutate(
    has_Il17a = stringr::str_detect(
      paste0("/", .data$geneID, "/"),
      stringr::fixed("/Il17a/")
    )
  )



go_classified <- go_table |>
  dplyr::mutate(

    Program = dplyr::case_when(

      # ------------------------------------------------------
      # TREG-ASSOCIATED
      # ------------------------------------------------------

      stringr::str_detect(
        .data$description_lower,
        "regulatory t cell|regulatory t-cell|immune tolerance"
      ) ~ "Treg-associated",

      # ------------------------------------------------------
      # TH17-ASSOCIATED
      #
      # Exige:
      # 1. termo relacionado à inflamação/IL-17
      # 2. Il17a entre os genes responsáveis pelo enriquecimento
      # ------------------------------------------------------

      (
        stringr::str_detect(
          .data$description_lower,
          paste(
            c(
              "th17",
              "interleukin-17",
              "il-17",
              "inflammatory response",
              "leukocyte differentiation",
              "helper t cell differentiation"
            ),
            collapse = "|"
          )
        ) &
          .data$has_Il17a
      ) ~ "Th17-associated",

      # ------------------------------------------------------
      # T-CELL ACTIVATION
      # ------------------------------------------------------

      stringr::str_detect(
        .data$description_lower,
        paste(
          c(
            "t cell activation",
            "t-cell activation",
            "lymphocyte activation",
            "antigen receptor-mediated signaling",
            "t cell receptor signaling",
            "t-cell receptor signaling"
          ),
          collapse = "|"
        )
      ) ~ "T-cell activation",

      # ------------------------------------------------------
      # MIGRATION AND ADHESION
      # ------------------------------------------------------

      stringr::str_detect(
        .data$description_lower,
        paste(
          c(
            "t cell migration",
            "t-cell migration",
            "lymphocyte migration",
            "leukocyte migration",
            "chemokine-mediated signaling",
            "cell-cell adhesion",
            "cell adhesion",
            "lamellipodium",
            "chemotaxis"
          ),
          collapse = "|"
        )
      ) ~ "Migration and adhesion",

      # ------------------------------------------------------
      # LESS DIFFERENTIATED
      # ------------------------------------------------------

      .data$description_lower %in% c(
        "t cell differentiation",
        "lymphocyte differentiation",
        "regulation of t cell differentiation",
        "regulation of lymphocyte differentiation"
      ) ~ "Less differentiated",

      TRUE ~ NA_character_
    )
  )


go_classified |>
  dplyr::filter(
    .data$Program == "Th17-associated"
  ) |>
  dplyr::select(
    dplyr::all_of(
      c(
        "Contrast",
        "Description",
        "geneID",
        "Count",
        "p.adjust"
      )
    )
  ) |>
  dplyr::arrange(
    .data$Contrast,
    .data$p.adjust
  ) 



th17_signature <- c(
  "Il17a",
  "Maf",
  "Hif1a",
  "Bhlhe40"
)

count_signature_genes <- function(
    gene_string,
    signature
) {

  if (
    is.na(gene_string) ||
    !nzchar(gene_string)
  ) {
    return(0L)
  }

  term_genes <- unique(
    strsplit(
      gene_string,
      split = "/",
      fixed = TRUE
    )[[1]]
  )

  length(
    intersect(
      term_genes,
      signature
    )
  )
}

get_signature_genes <- function(
    gene_string,
    signature
) {

  if (
    is.na(gene_string) ||
    !nzchar(gene_string)
  ) {
    return("")
  }

  term_genes <- unique(
    strsplit(
      gene_string,
      split = "/",
      fixed = TRUE
    )[[1]]
  )

  paste(
    intersect(
      term_genes,
      signature
    ),
    collapse = "/"
  )
}

go_table <- go_table |>
  dplyr::mutate(

    N_Th17_genes = vapply(
      .data$geneID,
      count_signature_genes,
      integer(1),
      signature = th17_signature
    ),

    Th17_genes = vapply(
      .data$geneID,
      get_signature_genes,
      character(1),
      signature = th17_signature
    ),

    has_Il17a = stringr::str_detect(
      paste0("/", .data$geneID, "/"),
      stringr::fixed("/Il17a/")
    )
  )







signatures <- list(

  Less_differentiated = c(
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

  T_cell_activation = c(
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

  Migration_and_adhesion = c(
    "Rgs1",
    "Cxcr4",
    "Ccr2",
    "Ccr7",
    "S1pr1",
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

  Th17_associated = c(
    "Il17a",
    "Maf",
    "Hif1a",
    "Bhlhe40"
  ),

  Treg_associated = c(
    "Ctla4",
    "Il2ra",
    "Tmem176a",
    "Tmem176b",
    "Tgfbr2"
  )
)
















split_gene_ids <- function(gene_string) {

  if (
    length(gene_string) == 0 ||
    is.na(gene_string) ||
    !nzchar(gene_string)
  ) {
    return(character(0))
  }

  genes <- unlist(
    strsplit(
      gene_string,
      split = "/",
      fixed = TRUE
    ),
    use.names = FALSE
  )

  genes <- trimws(genes)

  unique(
    genes[
      !is.na(genes) &
        nzchar(genes)
    ]
  )
}


count_signature_genes <- function(
    gene_string,
    signature
) {

  term_genes <- split_gene_ids(gene_string)

  length(
    intersect(
      term_genes,
      signature
    )
  )
}


get_signature_genes <- function(
    gene_string,
    signature
) {

  term_genes <- split_gene_ids(gene_string)

  matched_genes <- intersect(
    term_genes,
    signature
  )

  paste(
    matched_genes,
    collapse = "/"
  )
}


contains_gene <- function(
    gene_string,
    gene
) {

  gene %in% split_gene_ids(gene_string)
}


# ============================================================
# 4. PREPARAR TABELA GO
# ============================================================

go_table <- go_table |>
  dplyr::mutate(

    Cluster = as.character(.data$Cluster),

    Description = as.character(
      .data$Description
    ),

    geneID = as.character(
      .data$geneID
    ),

    p.adjust = as.numeric(
      .data$p.adjust
    ),

    description_lower = stringr::str_to_lower(
      .data$Description
    ),

    Contrast = dplyr::recode(
      .data$Cluster,

      "Infection" = "PBS vs Naive",
      "IgG" = "IgG vs PBS",
      "anti-PD-1" = "Anti-PD-1 vs PBS",
      "anti-CTLA-4" = "Anti-CTLA-4 vs PBS",

      .default = .data$Cluster
    )
  )


# ============================================================
# 5. CONTAR GENES DE CADA ASSINATURA EM CADA TERMO GO
# ============================================================

go_table <- go_table |>
  dplyr::mutate(

    N_Less_genes = vapply(
      .data$geneID,
      count_signature_genes,
      integer(1),
      signature = signatures$Less_differentiated
    ),

    Less_genes = vapply(
      .data$geneID,
      get_signature_genes,
      character(1),
      signature = signatures$Less_differentiated
    ),

    N_Activation_genes = vapply(
      .data$geneID,
      count_signature_genes,
      integer(1),
      signature = signatures$T_cell_activation
    ),

    Activation_genes = vapply(
      .data$geneID,
      get_signature_genes,
      character(1),
      signature = signatures$T_cell_activation
    ),

    N_Migration_genes = vapply(
      .data$geneID,
      count_signature_genes,
      integer(1),
      signature = signatures$Migration_and_adhesion
    ),

    Migration_genes = vapply(
      .data$geneID,
      get_signature_genes,
      character(1),
      signature = signatures$Migration_and_adhesion
    ),

    N_Th17_genes = vapply(
      .data$geneID,
      count_signature_genes,
      integer(1),
      signature = signatures$Th17_associated
    ),

    Th17_genes = vapply(
      .data$geneID,
      get_signature_genes,
      character(1),
      signature = signatures$Th17_associated
    ),

    N_Treg_genes = vapply(
      .data$geneID,
      count_signature_genes,
      integer(1),
      signature = signatures$Treg_associated
    ),

    Treg_genes = vapply(
      .data$geneID,
      get_signature_genes,
      character(1),
      signature = signatures$Treg_associated
    ),

    has_Il17a = vapply(
      .data$geneID,
      contains_gene,
      logical(1),
      gene = "Il17a"
    )
  )



go_classified <- go_table |>
  dplyr::mutate(

    Program = dplyr::case_when(

      # ------------------------------------------------------
      # TH17-ASSOCIATED
      #
      # Exige:
      # - termo compatível com Th17/inflamação
      # - presença obrigatória de Il17a
      # ------------------------------------------------------

      (
        stringr::str_detect(
          .data$description_lower,
          paste(
            c(
              "th17",
              "interleukin-17",
              "il-17",
              "inflammatory response",
              "helper t cell differentiation",
              "leukocyte differentiation"
            ),
            collapse = "|"
          )
        ) &
          .data$has_Il17a
      ) ~ "Th17-associated",

      # ------------------------------------------------------
      # TREG-ASSOCIATED
      #
      # Exige:
      # - termo relacionado a Treg/regulação imune
      # - pelo menos 1 gene da assinatura Treg
      # ------------------------------------------------------

      (
        stringr::str_detect(
          .data$description_lower,
          paste(
            c(
              "regulatory t cell",
              "regulatory t-cell",
              "immune tolerance",
              "negative regulation of immune response",
              "negative regulation of t cell activation",
              "regulation of immune response"
            ),
            collapse = "|"
          )
        ) &
          .data$N_Treg_genes >= 1
      ) ~ "Treg-associated",

      # ------------------------------------------------------
      # T-CELL ACTIVATION
      #
      # Exige:
      # - termo de ativação/TCR
      # - pelo menos 1 gene da assinatura
      # ------------------------------------------------------

      (
        stringr::str_detect(
          .data$description_lower,
          paste(
            c(
              "t cell activation",
              "t-cell activation",
              "lymphocyte activation",
              "alpha-beta t cell activation",
              "antigen receptor-mediated signaling",
              "t cell receptor signaling",
              "t-cell receptor signaling",
              "costimulation",
              "co-stimulation"
            ),
            collapse = "|"
          )
        ) &
          .data$N_Activation_genes >= 1
      ) ~ "T-cell activation",

      # ------------------------------------------------------
      # MIGRATION AND ADHESION
      #
      # Exige:
      # - termo de migração/adesão/citoesqueleto
      # - pelo menos 1 gene da assinatura
      # ------------------------------------------------------

      (
        stringr::str_detect(
          .data$description_lower,
          paste(
            c(
              "t cell migration",
              "t-cell migration",
              "lymphocyte migration",
              "leukocyte migration",
              "cell migration",
              "chemotaxis",
              "chemokine-mediated signaling",
              "cell adhesion",
              "cell-cell adhesion",
              "integrin-mediated",
              "lamellipodium",
              "cell motility",
              "actin filament",
              "actin cytoskeleton"
            ),
            collapse = "|"
          )
        ) &
          .data$N_Migration_genes >= 1
      ) ~ "Migration and adhesion",

      # ------------------------------------------------------
      # LESS DIFFERENTIATED
      #
      # Exige:
      # - termo de diferenciação
      # - pelo menos 1 gene da assinatura
      # ------------------------------------------------------

      (
        stringr::str_detect(
          .data$description_lower,
          paste(
            c(
              "^t cell differentiation$",
              "^lymphocyte differentiation$",
              "^regulation of t cell differentiation$",
              "^regulation of lymphocyte differentiation$"
            ),
            collapse = "|"
          )
        ) &
          .data$N_Less_genes >= 1
      ) ~ "Less differentiated",

      TRUE ~ NA_character_
    ),

    Signature_genes = dplyr::case_when(

      .data$Program == "Less differentiated" ~
        .data$Less_genes,

      .data$Program == "T-cell activation" ~
        .data$Activation_genes,

      .data$Program == "Migration and adhesion" ~
        .data$Migration_genes,

      .data$Program == "Th17-associated" ~
        .data$Th17_genes,

      .data$Program == "Treg-associated" ~
        .data$Treg_genes,

      TRUE ~ NA_character_
    ),

    N_signature_genes = dplyr::case_when(

      .data$Program == "Less differentiated" ~
        .data$N_Less_genes,

      .data$Program == "T-cell activation" ~
        .data$N_Activation_genes,

      .data$Program == "Migration and adhesion" ~
        .data$N_Migration_genes,

      .data$Program == "Th17-associated" ~
        .data$N_Th17_genes,

      .data$Program == "Treg-associated" ~
        .data$N_Treg_genes,

      TRUE ~ NA_integer_
    )
  )







go_classified_filtered <- go_classified |>
  dplyr::filter(
    !is.na(.data$Program),
    !is.na(.data$p.adjust),
    .data$p.adjust < 0.05,
    .data$N_signature_genes > 0
  )


if (nrow(go_classified_filtered) == 0) {
  stop(
    paste(
      "Nenhum termo GO permaneceu.",
      "Verifique as assinaturas ou relaxe os critérios."
    )
  )
}





classified_check <- go_classified_filtered |>
  dplyr::select(
    dplyr::all_of(
      c(
        "Contrast",
        "Program",
        "Description",
        "Signature_genes",
        "N_signature_genes",
        "p.adjust",
        "geneID"
      )
    )
  ) |>
  dplyr::arrange(
    .data$Program,
    .data$Contrast,
    .data$p.adjust
  )






program_summary <- go_classified_filtered |>
  dplyr::group_by(
    .data$Contrast,
    .data$Program
  ) |>
  dplyr::arrange(
    .data$p.adjust,
    dplyr::desc(.data$N_signature_genes),
    .by_group = TRUE
  ) |>
  dplyr::slice_head(n = 1) |>
  dplyr::ungroup() |>
  dplyr::transmute(

    Contrast = .data$Contrast,
    Program = .data$Program,

    Best_term = .data$Description,

    Signature_genes = .data$Signature_genes,

    N_signature_genes = .data$N_signature_genes,

    best_padj = .data$p.adjust,

    minus_log10_padj = -log10(
      pmax(
        .data$p.adjust,
        .Machine$double.xmin
      )
    )
  )


# ============================================================
# 10. ORDENS
# ============================================================

program_levels <- c(
  "Less differentiated",
  "T-cell activation",
  "Migration and adhesion",
  "Th17-associated",
  "Treg-associated"
)

contrast_levels <- c(
  "PBS vs Naive",
  "IgG vs PBS",
  "Anti-PD-1 vs PBS",
  "Anti-CTLA-4 vs PBS"
)


# ============================================================
# 11. COMPLETAR COMBINAÇÕES AUSENTES
# ============================================================

program_summary_complete <- tidyr::expand_grid(

  Program = program_levels,

  Contrast = contrast_levels

) |>
  dplyr::left_join(

    program_summary,

    by = c(
      "Program",
      "Contrast"
    )

  ) |>
  dplyr::mutate(

    Best_term = tidyr::replace_na(
      .data$Best_term,
      "Not enriched"
    ),

    Signature_genes = tidyr::replace_na(
      .data$Signature_genes,
      ""
    ),

    N_signature_genes = tidyr::replace_na(
      .data$N_signature_genes,
      0L
    ),

    best_padj = tidyr::replace_na(
      .data$best_padj,
      1
    ),

    minus_log10_padj = tidyr::replace_na(
      .data$minus_log10_padj,
      0
    ),

    Program = factor(
      .data$Program,
      levels = rev(program_levels)
    ),

    Contrast = factor(
      .data$Contrast,
      levels = contrast_levels
    )
  )


# ============================================================
# 12. HEATMAP
#
# Cor:
#   -log10(p.adjust) do melhor termo válido
#
# Texto:
#   genes da assinatura que sustentam o programa
# ============================================================


program_summary_complete = program_summary_complete |> subset( Program != 'Treg-associated')

p_program_heatmap <- ggplot2::ggplot(

  program_summary_complete,

  ggplot2::aes(
    x = .data$Contrast,
    y = .data$Program,
    fill = .data$minus_log10_padj
  )

) +

  ggplot2::geom_tile(
    color = "white",
    linewidth = 1
  ) +

  ggplot2::scale_fill_viridis_c(

    name = expression(
      -log[10]("adjusted p-value")
    ),

    option = "magma",

    limits = c(
      0,
      max(
        program_summary_complete$minus_log10_padj,
        na.rm = TRUE
      )
    )
  ) +

  ggplot2::labs(
    x = NULL,
    y = NULL
  ) +

  ggplot2::theme_test(
    base_size = 18
  ) +

  ggplot2::theme(

    axis.text.x = ggplot2::element_text(
      angle = 40,
      hjust = 1,
      vjust = 1,
      color = "black"
    ),

    axis.text.y = ggplot2::element_blank(),

    axis.ticks = ggplot2::element_blank(),

    legend.position = "right",

    plot.margin = ggplot2::margin(
      8,
      8,
      8,
      8
    )
  )

print(p_program_heatmap)


# ============================================================
# 13. SALVAR
# ============================================================

ggplot2::ggsave(
  filename = "GO_programs_gene_supported_heatmap.tiff",
  plot = p_program_heatmap,
  device = ragg::agg_tiff,
  width = 6.5,
  height = 5.2,
  units = "in",
  res = 600,
  compression = "lzw",
  background = "white"
)









program_levels <- c(
  "Less differentiated",
  "T-cell activation",
  "Migration and adhesion",
  "Th17-associated",
  "Treg-associated"
)

contrast_levels <- c(
  "PBS vs Naive",
  "IgG vs PBS",
  "Anti-PD-1 vs PBS",
  "Anti-CTLA-4 vs PBS"
)


alluvial_data <- program_summary_complete |>
  dplyr::mutate(

    Program = factor(
      as.character(Program),
      levels = program_levels
    ),

    Contrast = factor(
      as.character(Contrast),
      levels = contrast_levels
    ),

    Enriched = (
      Best_term != "Not enriched" &
      best_padj < 0.05 &
      minus_log10_padj > 0
    ),

    raw_weight = dplyr::if_else(
      Enriched,
      minus_log10_padj,
      0
    )
  ) |>
  dplyr::filter(
    Enriched
  )



alluvial_data <- alluvial_data |>
  dplyr::group_by(Program) |>
  dplyr::mutate(

    program_total = sum(
      raw_weight,
      na.rm = TRUE
    ),

    Flow_weight = dplyr::if_else(
      program_total > 0,
      raw_weight / program_total,
      0
    )
  ) |>
  dplyr::ungroup()



alluvial_data |>
  dplyr::group_by(Program) |>
  dplyr::summarise(
    Total_node_width = sum(Flow_weight),
    .groups = "drop"
  )



program_colors <- c(
  "Less differentiated"   = "#D99B7F",
  "T-cell activation"     = "#A56F63",
  "Migration and adhesion" = "#464858",
  "Th17-associated"       = "#618764",
  "Treg-associated"       = "#9CB080"
)

contrast_colors <- c(
  "PBS vs Naive"        = "#6F7F5A",
  "IgG vs PBS"          = "#4C2A78",
  "Anti-PD-1 vs PBS"    = "#28331A",
  "Anti-CTLA-4 vs PBS"  = "#9BC5F3"
)


gap_factor <- 0.80

alluvial_data <- alluvial_data |>
  dplyr::mutate(
    flow_weight_plot = Flow_weight * gap_factor
  )

alluvial_data$Program = factor(alluvial_data$Program, levels = c("Th17-associated","Migration and adhesion", "T-cell activation" ,"Less differentiated"))


p_program_alluvial <- ggplot(
  alluvial_data,
  aes(axis1 = Best_term,axis2 = Program,y = Flow_weight)) +
  geom_alluvium(
    aes(fill = Program),
    width = 0.10,
    alpha = 0.75,
    knot.pos = 0.45,
    color = NA
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
  scale_fill_manual(values = program_colors)



pb = ggplot_build(p_program_alluvial)
strata_data = pb$data[[2]] 

p_program_alluvial = p_program_alluvial + 
  geom_text_repel(
    data = strata_data |> filter(x == 1),
    aes(x = x, y = y, label = stratum),
    inherit.aes = FALSE,
    force = 1,
    hjust = 1,
    nudge_x = -1,
    min.segment.length = 0,
    segment.curvature = 0,
    box.padding = unit(.02, "lines"),
    direction = "y",
    size = 0,
    segment.size = 0,
    color = 'white'
  ) +
  geom_text_repel(
    data = strata_data |> filter(x == 2),
    aes(x = x, y = y, label = stratum),
    inherit.aes = FALSE,
    force = 1,
    nudge_x = .6,
    color = 'white',
    #nudge_y = 2,
    min.segment.length = 0,
    segment.curvature = 0,
    direction = "y",
    box.padding = unit(.02, "lines"),
    size = 0,
    segment.size = 00,
    max.overlaps = getOption("ggrepel.max.overlaps", default = 1000)
  ) 


tiff("alluvial.tiff",units = 'in', width = 6, height = 5, res = 300)
p_program_alluvial
dev.off()



