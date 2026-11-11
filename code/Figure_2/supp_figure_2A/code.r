

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


remotes::install_github("neurorestore/Augur")




##### tripartite

library(Seurat)
library(Matrix)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(patchwork)

# ============================================================
# 1. PREPARAR O OBJETO
# ============================================================

obj <- combined_data_integrated
DefaultAssay(obj) <- "RNA"

obj$cell_class <- case_when(
  obj$subtype_cell == "Macrophage" ~ "Macrophages",
  obj$subtype_cell == "Neutrophils" ~ "Neutrophils",
  grepl("^T cells", obj$subtype_cell) ~ "T cells",
  TRUE ~ NA_character_
)

expr <- GetAssayData(
  obj,
  assay = "RNA",
  layer = "data"
)

meta <- obj@meta.data
meta$cell_id <- rownames(meta)

# ============================================================
# 2. DEFINIR COMPARAÇÕES E CORES
# ============================================================

comparisons <- list(
  "Naive / PBS / anti-PD-1" = c(
    "naive", "pbs", "antipd1"
  ),
  "Naive / PBS / anti-CTLA-4" = c(
    "naive", "pbs", "antictla4"
  )
)

celltypes <- c(
  "Macrophages",
  "Neutrophils",
  "T cells"
)

group_labels <- c(
  naive = "Naive",
  pbs = "PBS",
  antipd1 = "anti-PD-1",
  antictla4 = "anti-CTLA-4"
)

group_colors <- c(
  "Naive" = "#F15B83",
  "PBS" = "#788A64",
  "anti-PD-1" = "#26361E",
  "anti-CTLA-4" = "#99C4EF"
)

# ============================================================
# 3. CALCULAR EXPRESSÃO RELATIVA POR GRUPO
# ============================================================

get_ternary_data <- function(
    celltype,
    groups_use,
    comparison_name,
    min_expression = 0.1
) {

  md <- meta %>%
    filter(
      cell_class == celltype,
      groups %in% groups_use
    )

  if (nrow(md) == 0) {
    stop("Nenhuma célula encontrada para: ", celltype)
  }

  # Calcular expressão média por animal
  sample_ids <- unique(md$cartig)

  sample_expr <- lapply(sample_ids, function(s) {

    cells <- md$cell_id[md$cartig == s]

    Matrix::rowMeans(
      expr[, cells, drop = FALSE]
    )
  })

  sample_mat <- do.call(cbind, sample_expr)
  colnames(sample_mat) <- sample_ids

  sample_groups <- vapply(
    sample_ids,
    function(s) {
      unique(as.character(
        md$groups[md$cartig == s]
      ))[1]
    },
    character(1)
  )

  # Média entre animais de cada condição
  group_mat <- sapply(
    groups_use,
    function(g) {

      sample_cols <- sample_groups == g

      if (!any(sample_cols)) {
        stop("Grupo ausente: ", g)
      }

      rowMeans(
        sample_mat[, sample_cols, drop = FALSE]
      )
    }
  )

  # Filtrar genes com expressão muito baixa
  keep <- apply(
    group_mat,
    1,
    max
  ) >= min_expression

  group_mat <- group_mat[
    keep,
    ,
    drop = FALSE
  ]

  # Converter para proporções relativas
  totals <- rowSums(group_mat)

  keep_total <- is.finite(totals) & totals > 0

  group_mat <- group_mat[
    keep_total,
    ,
    drop = FALSE
  ]

  relative <- sweep(
    group_mat,
    1,
    rowSums(group_mat),
    "/"
  )

  result <- data.frame(
    gene = rownames(relative),
    Left = relative[, groups_use[1]],
    Right = relative[, groups_use[2]],
    Top = relative[, groups_use[3]],
    cell_type = celltype,
    Comparison = comparison_name,
    stringsAsFactors = FALSE
  )

  # Identificar condição dominante
  labels_use <- unname(group_labels[groups_use])

  proportions <- as.matrix(
    result[, c("Left", "Right", "Top")]
  )

  result$dominant <- labels_use[
    max.col(proportions, ties.method = "first")
  ]

  result$max_prop <- apply(
    proportions, 1, max
  )

  result$second_prop <- apply(
    proportions,
    1,
    function(x) sort(x, decreasing = TRUE)[2]
  )

  result$specificity <- (
    result$max_prop - result$second_prop
  )

  # Coordenadas do triângulo
  result$x <- result$Right + 0.5 * result$Top

  result$y <- sqrt(3) / 2 * result$Top

  result
}

# ============================================================
# 4. EXECUTAR AS DUAS COMPARAÇÕES
# ============================================================

ternary_list <- list()

for (comparison_name in names(comparisons)) {

  groups_use <- comparisons[[comparison_name]]

  for (celltype in celltypes) {

    message(
      "Processing: ",
      comparison_name,
      " | ",
      celltype
    )

    result <- get_ternary_data(
      celltype = celltype,
      groups_use = groups_use,
      comparison_name = comparison_name
    )

    key <- paste(
      comparison_name,
      celltype,
      sep = "_"
    )

    ternary_list[[key]] <- result
  }
}

ternary_df <- bind_rows(ternary_list)

# Conferir quantidade de genes
ternary_df %>%
  count(Comparison, cell_type)

# ============================================================
# 5. SELECIONAR GENES MAIS ESPECÍFICOS
# ============================================================

top_vertex_genes <- ternary_df %>%
  group_by(
    Comparison,
    cell_type,
    dominant
  ) %>%
  slice_max(
    order_by = specificity,
    n = 8,
    with_ties = FALSE
  ) %>%
  ungroup()

# Tabela auxiliar
top_vertex_table <- top_vertex_genes %>%
  select(
    Comparison,
    cell_type,
    dominant,
    gene,
    Left,
    Right,
    Top,
    max_prop,
    specificity
  ) %>%
  arrange(
    Comparison,
    cell_type,
    dominant,
    desc(specificity)
  )

print(top_vertex_table, n = Inf)

# ============================================================
# 6. FUNÇÃO PARA PLOTAR O TERNARY
# ============================================================

triangle <- data.frame(
  x = c(0, 1, 0.5, 0),
  y = c(0, 0, sqrt(3) / 2, 0)
)

plot_ternary <- function(
    comparison_name,
    show_labels = TRUE,
    n_labels = 8,
    show_vertex_labels = FALSE
) {

  plot_df <- ternary_df %>%
    filter(Comparison == comparison_name)

  label_df <- top_vertex_genes %>%
    filter(Comparison == comparison_name) %>%
    group_by(cell_type, dominant) %>%
    slice_head(n = n_labels) %>%
    ungroup()

  groups_use <- comparisons[[comparison_name]]

  left_label <- group_labels[groups_use[1]]
  right_label <- group_labels[groups_use[2]]
  top_label <- group_labels[groups_use[3]]

  p <- ggplot() +
    geom_path(
      data = triangle,
      aes(x, y),
      color = "grey40",
      linewidth = 0.5
    ) +
    geom_point(
      data = plot_df,
      aes(
        x = x,
        y = y,
        color = dominant
      ),
      size = 0.55,
      alpha = 0.6
    ) +
    scale_color_manual(
      values = group_colors,
      name = "Dominant condition"
    ) +
    facet_wrap(
      ~cell_type,
      nrow = 1
    ) +
    coord_equal(
      clip = "off",
      xlim = c(-0.12, 1.12),
      ylim = c(-0.12, 1.02)
    ) +
    theme_void(base_size = 12) +
    theme(
      strip.text = element_text(
        face = "bold",
        size = 13
      ),
      legend.position = "bottom",
      plot.margin = margin(15, 20, 15, 20)
    )

  if (show_vertex_labels) {
    p <- p +
      annotate("text", x = 0, y = -0.055, label = left_label, size = 3.5) +
      annotate("text", x = 1, y = -0.055, label = right_label, size = 3.5) +
      annotate("text", x = 0.5, y = 0.92, label = top_label, size = 3.5)
  }

  if (show_labels) {
    p <- p +
      geom_text_repel(
        data = label_df,
        aes(
          x = x,
          y = y,
          label = gene,
        ),
        color = "black",
        size = 4,
        box.padding = 0.25,
        point.padding = 0.15,
        segment.color = "grey55",
        max.overlaps = Inf,
        seed = 123,
        show.legend = FALSE
      )
  }

  p
}

# ============================================================
# 7. GERAR OS GRÁFICOS
# ============================================================

# Naive / PBS / anti-PD-1
p_pd1 <- plot_ternary(
  "Naive / PBS / anti-PD-1",
  show_labels = TRUE,
  n_labels = 10,
  show_vertex_labels = FALSE
)

print(p_pd1)

# Naive / PBS / anti-CTLA-4
p_ctla4 <- plot_ternary(
  "Naive / PBS / anti-CTLA-4",
  show_labels = TRUE,
  n_labels = 10,
  show_vertex_labels = FALSE
)

print(p_ctla4)


# Salvar cada comparação separadamente
ggsave(
  filename = "./Ternary_Naive_PBS_PD1.tiff",
  plot = p_pd1,
  device = ragg::agg_tiff,
  width = 15,
  height = 5,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggsave(
  filename = "./Ternary_Naive_PBS_CTLA4.tiff",
  plot = p_ctla4,
  device = ragg::agg_tiff,
  width = 15,
  height = 5,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)










library(Seurat)
library(Matrix)
library(dplyr)
library(tidyr)
library(ggplot2)

make_facet_dotplot_df <- function(
    obj,
    genes_use,
    groups_use,
    celltype_col = "subtype_cell",
    group_col = "groups",
    assay = "RNA",
    layer = "data"
) {

  DefaultAssay(obj) <- assay

  genes_use <- intersect(genes_use, rownames(obj))

  expr <- GetAssayData(
    obj,
    assay = assay,
    layer = layer
  )[genes_use, , drop = FALSE]

  meta <- obj@meta.data %>%
    mutate(cell_id = colnames(obj)) %>%
    filter(.data[[group_col]] %in% groups_use)

  expr <- expr[, meta$cell_id, drop = FALSE]

  res_list <- lapply(seq_along(genes_use), function(i) {
    g <- genes_use[i]
    vals <- expr[g, ]

    df <- meta %>%
      mutate(expr = as.numeric(vals)) %>%
      group_by(.data[[group_col]], .data[[celltype_col]]) %>%
      summarise(
        avg.exp = mean(expr, na.rm = TRUE),
        pct.exp = mean(expr > 0, na.rm = TRUE) * 100,
        .groups = "drop"
      ) %>%
      mutate(gene = g)

    df
  })

  plot_df <- bind_rows(res_list)

  # escalar por gene, como o DotPlot do Seurat faz
  plot_df <- plot_df %>%
    group_by(gene) %>%
    mutate(
      avg.exp.scaled = as.numeric(scale(avg.exp))
    ) %>%
    ungroup()

  plot_df[[group_col]] <- factor(plot_df[[group_col]], levels = groups_use)

  plot_df
}

plot_facet_dotplot <- function(
    plot_df,
    genes_order = NULL,
    group_col = "groups",
    celltype_col = "subtype_cell",
    title = NULL
) {

  if (!is.null(genes_order)) {
    plot_df$gene <- factor(plot_df$gene, levels = genes_order)
  }

  ggplot(
    plot_df,
    aes(
      x = gene,
      y = .data[[celltype_col]],
      size = pct.exp,
      color = avg.exp.scaled
    )
  ) +
    geom_point() +
    facet_wrap(
      as.formula(paste("~", group_col)),
      nrow = 1
    ) +
    scale_color_gradient2(
      low = "#EAEAEA",
      mid = "#F3B4C4",
      high = "#BD3361",
      midpoint = 0,
      name = "Average Expression"
    ) +
    scale_size(
      range = c(0, 6),
      name = "Percent Expressed"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 8
      ),
      axis.text.y = element_text(size = 10),
      strip.text = element_text(face = "bold", size = 12),
      axis.title = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5)
    ) +
    labs(title = title)
}

genes_pd1 <- top_vertex_genes %>%
  filter(Comparison == "Naive / PBS / anti-PD-1") %>%
  group_by(cell_type, dominant) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  pull(gene) %>%
  unique()

genes_pd1
plot_df_pd1 <- make_facet_dotplot_df(
  obj = combined_data_integrated,
  genes_use = genes_pd1,
  groups_use = c("naive", "pbs", "antipd1"),
  celltype_col = "subtype_cell",
  group_col = "groups",
  assay = "RNA",
  layer = "data"
)



p_pd1 <- plot_facet_dotplot(
  plot_df = plot_df_pd1,
  genes_order = genes_pd1,
  group_col = "groups",
  celltype_col = "subtype_cell",
  title = "Naive / PBS / anti-PD-1"
)

print(p_pd1)




genes_ctla4 <- top_vertex_genes %>%
  filter(Comparison == "Naive / PBS / anti-CTLA-4") %>%
  group_by(cell_type, dominant) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  pull(gene) %>%
  unique()

genes_ctla4













#########


library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

obj <- combined_data_integrated
DefaultAssay(obj) <- "RNA"

obj$cell_class <- dplyr::case_when(
  obj$subtype_cell == "Macrophage" ~ "Macrophages",
  obj$subtype_cell == "Neutrophils" ~ "Neutrophils",

  obj$subtype_cell %in% c(
    "T cells CD4+/CD3+",
    "T cells CD3+"
  ) ~ "T cells",

  TRUE ~ NA_character_
)

# Conferir as populações selecionadas
table(obj$cell_class, useNA = "ifany")
table(obj$subtype_cell[obj$cell_class == "T cells"])
# ====================================================
# 2. ASSINATURAS DO DOTPLOT
# ====================================================

signatures <- list(
  "T cells" = c(
    "Cd3e", "Cd4", "Cd8a", "Cd8b1",
    "Tcf7", "Il7r", "Mki67",
    "Il17a", "Foxp3"
  ),
  "Macrophages" = c(
    "Adgre1", "Cd68", "Csf1r",
    "Mrc1", "Lgals3"
  ),
  "Neutrophils" = c(
    "S100a8", "S100a9", "Ly6g",
    "Csf3r", "Cxcr2"
  )
)

# ====================================================
# 3. CALCULAR SCORES POR POPULAÇÃO
# ====================================================

set.seed(123)

score_list <- list()

for (celltype in names(signatures)) {

  message("Processing: ", celltype)

  cells_use <- colnames(obj)[
    !is.na(obj$cell_class) &
    obj$cell_class == celltype
  ]

  obj_sub <- subset(
    obj,
    cells = cells_use
  )

  genes_use <- intersect(
    signatures[[celltype]],
    rownames(obj_sub)
  )

  print(genes_use)

  if (length(genes_use) < 3) {
    stop("Poucos genes para: ", celltype)
  }

  obj_sub <- AddModuleScore(
    object = obj_sub,
    features = list(genes_use),
    assay = "RNA",
    name = "Signature_score",
    seed = 123
  )

  # Média por animal
  df_sub <- obj_sub@meta.data %>%
    group_by(cartig, groups) %>%
    summarise(
      Score = mean(Signature_score1, na.rm = TRUE),
      n_cells = n(),
      .groups = "drop"
    ) %>%
    mutate(
      Signature = celltype
    )

  score_list[[celltype]] <- df_sub
}

df_long <- bind_rows(score_list)

# ====================================================
# 4. ORGANIZAR OS GRUPOS
# ====================================================

df_long <- df_long %>%
  mutate(
    groups = factor(
      groups,
      levels = c(
        "naive", "pbs", "igg",
        "antipd1", "antictla4"
      ),
      labels = c(
        "Naive", "PBS", "IgG",
        "anti-PD-1", "anti-CTLA-4"
      )
    ),
    Signature = factor(
      Signature,
      levels = c(
        "T cells",
        "Macrophages",
        "Neutrophils"
      )
    )
  )

# ====================================================
# 5. MÉDIA E SEM POR GRUPO
# ====================================================

df_summary <- df_long %>%
  group_by(Signature, groups) %>%
  summarise(
    mean = mean(Score, na.rm = TRUE),
    sd = sd(Score, na.rm = TRUE),
    n = sum(!is.na(Score)),
    sem = sd / sqrt(n),
    .groups = "drop"
  )

print(df_summary, n = Inf)

# ====================================================
# 6. CORES DO MANUSCRITO
# ====================================================

group_colors <- c(
  "Naive"       = "#F15B83",
  "PBS"         = "#788A64",
  "IgG"         = "#573080",
  "anti-PD-1"   = "#26361E",
  "anti-CTLA-4" = "#99C4EF"
)

# ====================================================
# 7. GRÁFICO — MESMO ESTILO ORIGINAL
# ====================================================

p_scores <- ggplot(
  df_long,
  aes(x = groups, y = Score)
) +

  geom_errorbar(
    data = df_summary,
    aes(
      x = groups,
      ymin = mean - sem,
      ymax = mean + sem
    ),
    inherit.aes = FALSE,
    width = 0.13,
    linewidth = 0.5,
    color = "black"
  ) +

  geom_crossbar(
    data = df_summary,
    aes(
      x = groups,
      y = mean,
      ymin = mean,
      ymax = mean
    ),
    inherit.aes = FALSE,
    width = 0.32,
    linewidth = 0.55,
    color = "black"
  ) +

  geom_jitter(
    aes(color = groups),
    width = 0.09,
    height = 0,
    size = 2.5,
    alpha = 0.9
  ) +

  scale_color_manual(
    values = group_colors
  ) +

  facet_wrap(
    ~ Signature,
    scales = "free_y",
    nrow = 1
  ) +

  labs(
    x = NULL,
    y = "Mean module score"
  ) +

  theme_classic(base_size = 12) +

  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 9,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 9,
      color = "black"
    ),
    axis.title.y = element_text(
      size = 11,
      color = "black"
    ),
    axis.line = element_line(
      linewidth = 0.5,
      color = "black"
    ),
    panel.spacing = grid::unit(0.8, "cm")
  )

print(p_scores)
