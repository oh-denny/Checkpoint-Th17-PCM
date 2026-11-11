




library(dplyr)
library(tidyr)
library(ggplot2)
library(ggalluvial)

# Assinaturas utilizadas na Figura 2B
signatures_fig2 <- list(
  "Less differentiated" = c(
    "Tcf7", "Lef1", "Sell", "Ccr7", "Klf2"
  ),
  "T-cell activation" = c(
    "Cd44", "Cd69", "Icos", "Tnfrsf4", "Tnfrsf9"
  ),
  "Proliferation" = c(
    "Mki67", "Top2a", "Pcna", "Mcm5", "Stmn1"
  ),
  "Th17" = c(
    "Rorc", "Rora", "Il17a", "Il17f",
    "Il23r", "Ccr6", "Maf"
  ),
  "Treg" = c(
    "Foxp3", "Il2ra", "Ctla4",
    "Ikzf2", "Tigit", "Il10"
  )
)

program_order <- names(signatures_fig2)

program_colors <- c(
  "Less differentiated" = "#7AA5C4",
  "T-cell activation" = "#B879A5",
  "Proliferation" = "#80A477",
  "Th17" = "#DE7968",
  "Treg" = "#8162AD"
)

# Converter assinaturas em tabela
alluvial_df <- tibble(
  Program = rep(
    names(signatures_fig2),
    lengths(signatures_fig2)
  ),
  Gene = unlist(
    signatures_fig2,
    use.names = FALSE
  ),
  Weight = 1
) %>%
  mutate(
    Program = factor(
      Program,
      levels = program_order
    ),
    Gene = factor(
      Gene,
      levels = rev(unique(Gene))
    )
  )

# Conferir os genes
print(alluvial_df, n = Inf)

# Construir o alluvial
p_alluvial <- ggplot(
  alluvial_df,
  aes(
    axis1 = Gene,
    axis2 = Program,
    y = Weight
  )
) +
  geom_alluvium(
    aes(fill = Program),
    width = 0.13,
    alpha = 0.7,
    knot.pos = 0.45
  ) +
  geom_stratum(
    width = 0.13,
    fill = "grey95",
    color = "grey40",
    linewidth = 0.2
  ) +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    size = 2.6
  ) +
  scale_x_discrete(
    limits = c("Genes", "Transcriptional programs"),
    expand = c(0.20, 0.20)
  ) +
  scale_fill_manual(
    values = program_colors
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      color = "black",
      face = "bold",
      size = 11
    ),
    plot.margin = margin(10, 35, 10, 35)
  )

print(p_alluvial)






library(Seurat)
library(dplyr)
library(ggplot2)

obj <- combined_data_integrated
DefaultAssay(obj) <- "RNA"

# Conferir os nomes das populações
table(obj$subtype_cell)

# Selecionar populações T
t_cells <- colnames(obj)[
  grepl("T cells", obj$subtype_cell, ignore.case = TRUE)
]

obj_T <- subset(obj, cells = t_cells)

table(obj_T$groups)
table(obj_T$subtype_cell)


signatures <- list(
  Th17 = c(
    "Rorc", "Rora", "Il17a", "Il17f",
    "Il23r", "Ccr6", "Maf"
  ),
  Treg = c(
    "Foxp3", "Il2ra", "Ctla4",
    "Ikzf2", "Tigit", "Il10"
  ),
  Activation = c(
    "Cd44", "Cd69", "Icos",
    "Tnfrsf4", "Tnfrsf9"
  ),
  Naive = c(
    "Tcf7", "Lef1", "Sell",
    "Ccr7", "Klf2"
  ),
  Proliferation = c(
    "Mki67", "Top2a", "Pcna",
    "Mcm5", "Stmn1"
  )
)

# Manter genes presentes na matriz
signatures <- lapply(
  signatures,
  function(x) intersect(x, rownames(obj_T))
)

sapply(signatures, length)


set.seed(123)

for (nm in names(signatures)) {
  if (length(signatures[[nm]]) < 3) next

  obj_T <- AddModuleScore(
    object = obj_T,
    features = list(signatures[[nm]]),
    assay = "RNA",
    name = paste0(nm, "_score"),
    seed = 123
  )
}

score_cols <- grep(
  "_score1$",
  colnames(obj_T@meta.data),
  value = TRUE
)

df <- obj_T@meta.data %>%
  group_by(cartig, groups) %>%
  summarise(
    across(
      all_of(score_cols),
      ~ mean(.x, na.rm = TRUE)
    ),
    n_cells = n(),
    .groups = "drop"
  )

df_long <- tidyr::pivot_longer(
  df,
  cols = all_of(score_cols),
  names_to = "Signature",
  values_to = "Score"
)

library(ggplot2)
library(dplyr)

# Organizar os grupos
df_long <- df_long %>%
  mutate(
    groups = factor(
      groups,
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
    ),
    Signature = factor(
      Signature,
      levels = c(
        "Naive_score1",
        "Activation_score1",
        "Proliferation_score1",
        "Th17_score1",
        "Treg_score1"
      ),
      labels = c(
        "Less differentiated",
        "T-cell activation",
        "Proliferation",
        "Th17",
        "Treg"
      )
    )
  )

# Calcular média e SEM por grupo
df_summary <- df_long %>%
  group_by(Signature, groups) %>%
  summarise(
    mean = mean(Score, na.rm = TRUE),
    sd = sd(Score, na.rm = TRUE),
    n = sum(!is.na(Score)),
    sem = sd / sqrt(n),
    .groups = "drop"
  )

# Paleta utilizada no manuscrito
group_colors <- c(
  "Naive"       = "#F15B83",
  "PBS"         = "#788A64",
  "IgG"         = "#573080",
  "anti-PD-1"   = "#26361E",
  "anti-CTLA-4" = "#99C4EF"
)

# Gráfico
p_scores <- ggplot(
  df_long,
  aes(x = groups, y = Score)
) +

  # Barras de erro SEM
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

  # Média de cada grupo
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

  # Pontos dos animais
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














library(dplyr)
library(tidyr)
library(ggplot2)

# Garantir que os grupos estejam padronizados
df_long <- df_long %>%
  mutate(
    groups = as.character(groups),
    groups = recode(
      groups,
      "naive" = "Naive",
      "pbs" = "PBS",
      "igg" = "IgG",
      "antipd1" = "anti-PD-1",
      "antictla4" = "anti-CTLA-4"
    )
  )

# Confirmar número de animais por grupo
df_long %>%
  distinct(Signature, groups, cartig) %>%
  count(Signature, groups)


comparisons <- list(
  c("IgG", "anti-PD-1"),
  c("IgG", "anti-CTLA-4"),
  c("PBS", "anti-PD-1"),
  c("PBS", "anti-CTLA-4"),
  c("Naive", "PBS")
)

# Função para teste t de Welch
run_test <- function(data, comparison) {

  g1 <- comparison[1]
  g2 <- comparison[2]

  x <- data$Score[data$groups == g1]
  y <- data$Score[data$groups == g2]

  x <- x[is.finite(x)]
  y <- y[is.finite(y)]

  if (length(x) < 2 || length(y) < 2) {
    return(tibble(
      group1 = g1, group2 = g2,
      n1 = length(x), n2 = length(y),
      mean_diff = NA_real_,
      p = NA_real_
    ))
  }

  test <- tryCatch(
    t.test(y, x, var.equal = FALSE),
    error = function(e) NULL
  )

  tibble(
    group1 = g1,
    group2 = g2,
    n1 = length(x),
    n2 = length(y),
    mean_diff = mean(y) - mean(x),
    p = if (is.null(test)) NA_real_ else test$p.value
  )
}

# Executar para cada assinatura
results_stats <- df_long %>%
  group_by(Signature) %>%
  group_modify(~ {
    bind_rows(
      lapply(
        comparisons,
        function(comp) run_test(.x, comp)
      )
    )
  }) %>%
  ungroup() %>%
  mutate(
    p_adj = p.adjust(p, method = "BH"),
    significance = case_when(
      is.na(p_adj) ~ "NA",
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(Signature, p_adj)

print(results_stats, n = Inf)
results_stats %>%
  filter(!is.na(p), p < 0.05) %>%
  select(
    Signature,
    group1,
    group2,
    mean_diff,
    p,
    p_adj,
    significance
  )





library(scales)

# Posições das anotações
plot_ranges <- df_long %>%
  group_by(Signature) %>%
  summarise(
    ymax = max(Score, na.rm = TRUE),
    ymin = min(Score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    span = pmax(ymax - ymin, 0.01)
  )

annotations <- results_stats %>%
  filter(
    group1 == "IgG",
    group2 %in% c("anti-PD-1", "anti-CTLA-4")
  ) %>%
  left_join(plot_ranges, by = "Signature") %>%
  mutate(
    x = ifelse(group2 == "anti-PD-1", 4, 5),
    y = ymax + ifelse(
      group2 == "anti-PD-1",
      0.15 * span,
      0.30 * span
    ),
    label = paste0(
      "FDR = ",
      ifelse(
        is.na(p_adj),
        "NA",
        format.pval(p_adj, digits = 2, eps = 0.001)
      )
    )
  )

p_scores_stats <- p_scores +
  geom_text(
    data = annotations,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 2.5,
    hjust = 0.5,
    color = "black"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.4))
  )

print(p_scores_stats)







# Instalar somente se necessário
if (!requireNamespace("rstatix", quietly = TRUE)) {
  install.packages("rstatix")
}

library(rstatix)
library(dplyr)
library(tidyr)

# Organizar os grupos
df_long <- df_long %>%
  mutate(
    groups = factor(
      as.character(groups),
      levels = c(
        "Naive",
        "PBS",
        "IgG",
        "anti-PD-1",
        "anti-CTLA-4"
      )
    )
  )

# Conferir réplicas biológicas
df_long %>%
  group_by(Signature, groups) %>%
  summarise(
    n_animals = n_distinct(cartig),
    .groups = "drop"
  )


df_long <- df_long %>%
  mutate(
    groups = as.character(groups),
    groups = recode(
      groups,
      "naive" = "Naive",
      "pbs" = "PBS",
      "igg" = "IgG",
      "antipd1" = "anti-PD-1",
      "antictla4" = "anti-CTLA-4"
    ),
    groups = factor(
      groups,
      levels = c(
        "Naive", "PBS", "IgG",
        "anti-PD-1", "anti-CTLA-4"
      )
    )
  )

kw_results <- df_long %>%
  group_by(Signature) %>%
  kruskal_test(Score ~ groups) %>%
  ungroup() %>%
  mutate(
    p_adj_BH = p.adjust(p, method = "BH")
  )

print(kw_results, n = Inf)
dunn_results <- df_long %>%
  group_by(Signature) %>%
  dunn_test(
    Score ~ groups,
    p.adjust.method = "none"
  ) %>%
  ungroup() %>%
  mutate(
    p.adj = p.adjust(p, method = "BH"),
    significance = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01 ~ "**",
      p.adj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(Signature, p.adj)

library(dplyr)

dunn_results_raw <- dunn_results %>%
  mutate(
    significance_raw = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  select(
    Signature,
    group1,
    group2,
    statistic,
    p,
    significance_raw
  ) %>%
  arrange(Signature, p)

print(dunn_results_raw, n = Inf)


library(ggpubr)
library(dplyr)

# Anotações com p não ajustado
stat_plot <- dunn_results %>%
  filter(p < 0.05) %>%
  mutate(
    p_label = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*"
    )
  )

# Posições verticais por assinatura
y_positions <- df_long %>%
  group_by(Signature) %>%
  summarise(
    ymax = max(Score, na.rm = TRUE),
    ymin = min(Score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    span = pmax(ymax - ymin, 0.01)
  )

stat_plot <- stat_plot %>%
  left_join(y_positions, by = "Signature") %>%
  group_by(Signature) %>%
  arrange(p, .by_group = TRUE) %>%
  mutate(
    y.position = ymax + row_number() * 0.15 * span
  ) %>%
  ungroup()

p_scores_dunn <- p_scores +
  stat_pvalue_manual(
    stat_plot,
    label = "p_label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    bracket.size = 0.4,
    size = 6,
    inherit.aes = FALSE
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.45))
  )

print(p_scores_dunn)


ggsave(
  filename = "Figure2_Tcell_ModuleScores_Dunn.tiff",
  plot = p_scores_dunn,
  device = "tiff",
  width = 10,
  height = 3,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)







