# Figure/panel: 2F
# Purpose: Animal-averaged T-cell module scores; Kruskal–Wallis and Dunn tests.
# Inputs (under data/ unless session objects): combined_data_1.rds.
# Packages: Seurat, dplyr, ggplot2, ggpubr, grid, rstatix, tidyr.
# Outputs: Figure2_Tcell_ModuleScores_Dunn.tiff.
# Prerequisites: None; processed inputs required.
# Execution: repository root; see README and docs/SCRIPT_INVENTORY.md.


library(Seurat)
library(dplyr)
library(ggplot2)
library(rstatix)
library(tidyr)
library(ggpubr)

rds_path = file.path("data", "combined_data_1.rds")
obj = readRDS(rds_path)
DefaultAssay(obj) = "RNA"

table(obj$subtype_cell)

t_cells = colnames(obj)[
  grepl("T cells", obj$subtype_cell, ignore.case = TRUE)
]
obj_T = subset(obj, cells = t_cells)

table(obj_T$groups)
table(obj_T$subtype_cell)

signatures = list(
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

signatures = lapply(
  signatures,
  function(x) intersect(x, rownames(obj_T))
)

sapply(signatures, length)
set.seed(123)

for (nm in names(signatures)) {
  if (length(signatures[[nm]]) < 3) next

  obj_T = AddModuleScore(
    object = obj_T,
    features = list(signatures[[nm]]),
    assay = "RNA",
    name = paste0(nm, "_score"),
    seed = 123
  )
}

score_cols = grep("_score1$", colnames(obj_T@meta.data),value = TRUE)
df = obj_T@meta.data %>%
  group_by(cartig, groups) %>%
  summarise(
    across(
      all_of(score_cols),
      ~ mean(.x, na.rm = TRUE)
    ),
    n_cells = n(),
    .groups = "drop"
  )

df_long = tidyr::pivot_longer(
  df,
  cols = all_of(score_cols),
  names_to = "Signature",
  values_to = "Score"
)

df_long = df_long %>%
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


df_summary = df_long %>%
  group_by(Signature, groups) %>%
  summarise(
    mean = mean(Score, na.rm = TRUE),
    sd = sd(Score, na.rm = TRUE),
    n = sum(!is.na(Score)),
    sem = sd / sqrt(n),
    .groups = "drop"
  )

group_colors = c(
  "Naive"       = "#F15B83",
  "PBS"         = "#788A64",
  "IgG"         = "#573080",
  "anti-PD-1"   = "#26361E",
  "anti-CTLA-4" = "#99C4EF"
)

p_scores = ggplot(
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


df_long = df_long %>%
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

kw_results = df_long %>%
  group_by(Signature) %>%
  kruskal_test(Score ~ groups) %>%
  ungroup() %>%
  mutate(
    p_adj_BH = p.adjust(p, method = "BH")
  )

print(kw_results, n = Inf)
dunn_results = df_long %>%
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

dunn_results_raw = dunn_results %>%
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


stat_plot = dunn_results %>%
  filter(p < 0.05) %>%
  mutate(
    p_label = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*"
    )
  )

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

