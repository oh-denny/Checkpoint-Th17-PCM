


install.packages("DirichletReg")
install.packages("dplyr")
install.packages("tidyr")

library(DirichletReg)
library(dplyr)
library(tidyr)




meta <- combined_data_integrated@meta.data

meta$cluster <- meta$RNA_snn_res.0.2 

counts_df <- meta %>%
  group_by(cartig, groups, cluster) %>%
  summarise(count = n(), .groups = "drop")

wide_counts <- counts_df %>%
  pivot_wider(names_from = cluster, values_from = count, values_fill = 0)

cluster_cols <- setdiff(colnames(wide_counts), c("cartig", "groups"))

prop_matrix <- as.matrix(wide_counts[, cluster_cols])
prop_matrix <- prop_matrix / rowSums(prop_matrix)

# Criar objeto de proporção formatado para DirichletReg
Y <- DR_data(prop_matrix)

# Criar dataframe final com as variáveis independentes
dirichlet_data <- data.frame(
  groups = factor(wide_counts$groups),
  cartig = wide_counts$cartig
)
dirichlet_data$Y <- Y


# Ajustando a Regressão de Dirichlet
# Assumindo 'naive' ou 'PBS' como referência (ajuste o relevel conforme seu grupo controle principal)
dirichlet_data$groups <- relevel(dirichlet_data$groups, ref = "naive")

model <- DirichReg(Y ~ groups, data = dirichlet_data)

# Resumo do modelo com P-valores para cada cluster vs controle
summary(model)





library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(scales)

# Grupos na ordem usada pelo modelo
levels(dirichlet_data$groups)

# Clusters/componentes na mesma ordem da matriz composicional
cluster_names <- colnames(prop_matrix)

# Termos do modelo na mesma ordem do summary()
model_terms <- c(
  "(Intercept)",
  paste0(
    "groups",
    levels(dirichlet_data$groups)[-1]
  )
)

# Conferência
length(cluster_names)
length(model_terms)
length(model$coefficients)

stopifnot(
  length(model$coefficients) ==
    length(cluster_names) * length(model_terms)
)

dirichlet_results <- tidyr::expand_grid(
  cluster = cluster_names,
  term = model_terms
) |>
  mutate(
    estimate = unname(model$coefficients),
    std_error = unname(model$se),
    z_value = estimate / std_error,
    p_value = 2 * pnorm(
      abs(z_value),
      lower.tail = FALSE
    ),
    conf_low = estimate - 1.96 * std_error,
    conf_high = estimate + 1.96 * std_error
  )

head(dirichlet_results)


dirichlet_contrasts <- dirichlet_results |>
  filter(term != "(Intercept)") |>
  mutate(
    group = sub("^groups", "", term),

    group = recode(
      group,
      "antictla4" = "anti-CTLA-4",
      "antipd1"   = "anti-PD-1",
      "igg"       = "IgG",
      "pbs"       = "PBS"
    ),

    cluster = factor(
      cluster,
      levels = as.character(0:10)
    ),

    group = factor(
      group,
      levels = c(
        "PBS",
        "IgG",
        "anti-PD-1",
        "anti-CTLA-4"
      )
    )
  )


dirichlet_contrasts <- dirichlet_contrasts |>
  mutate(
    padj_global = p.adjust(
      p_value,
      method = "BH"
    ),
    significance = case_when(
      padj_global < 0.0001 ~ "****",
      padj_global < 0.001  ~ "***",
      padj_global < 0.01   ~ "**",
      padj_global < 0.05   ~ "*",
      TRUE                 ~ ""
    )
  )


dirichlet_contrasts |>
  filter(padj_global < 0.05) |>
  arrange(cluster, padj_global)



heatmap_dirichlet <- ggplot(
  dirichlet_contrasts,
  aes(
    x = group,
    y = cluster,
    fill = estimate
  )
) +
  geom_tile(
    colour = "black",
    linewidth = 0.2
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    name = expression(beta)
  ) +
  labs(
    x = NULL,
    y = "Cluster"
  ) +
  coord_fixed(ratio = 0.75) +
  theme_classic(18) +
  theme(
    axis.text.x = element_blank(),
    legend.position = "none"
  )


heatmap_dirichlet
pdf("heatmap_dirichlet.pdf", width = 5, height = 4)
heatmap_dirichlet
dev.off()




forest_dirichlet <- ggplot(
  dirichlet_contrasts,
  aes(
    x = estimate,
    y = fct_rev(cluster),
    colour = group
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_errorbarh(
    aes(
      xmin = conf_low,
      xmax = conf_high
    ),
    height = 0.16,
    position = position_dodge(
      width = 0.65
    ),
    linewidth = 0.55
  ) +
  geom_point(
    aes(
      shape = padj_global < 0.05
    ),
    position = position_dodge(
      width = 0.65
    ),
    size = 2.8
  ) +
  scale_shape_manual(
    values = c(
      `FALSE` = 1,
      `TRUE` = 16
    ),
    labels = c(
      `FALSE` = "Not significant",
      `TRUE` = "BH-adjusted p < 0.05"
    ),
    name = NULL
  ) +
  scale_colour_manual(
    values = c(
      "PBS"         = "#6F7F5A",
      "IgG"         = "#4C2A78",
      "anti-PD-1"   = "#6AAED6",
      "anti-CTLA-4" = "#5A86B8"
    ),
    name = NULL
  ) +
  theme_classic(18) +
  theme(
    axis.text.x = element_blank(),
    legend.position = "none"
  )











# Referência = PBS
dirichlet_data_pbs <- dirichlet_data

dirichlet_data_pbs$groups <- relevel(
  dirichlet_data_pbs$groups,
  ref = "pbs"
)

model_pbs <- DirichReg(
  Y ~ groups,
  data = dirichlet_data_pbs
)

sum_mod_pbs <- summary(model_pbs)

sum_mod_pbs



library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

cluster_names <- colnames(prop_matrix)

pbs_terms <- c(
  "(Intercept)",
  paste0(
    "groups",
    levels(dirichlet_data_pbs$groups)[-1]
  )
)

# Verifique a ordem
levels(dirichlet_data_pbs$groups)
pbs_terms


dirichlet_results_pbs <- tidyr::expand_grid(
  cluster = cluster_names,
  term = pbs_terms
) |>
  mutate(
    estimate = unname(model_pbs$coefficients),
    std_error = unname(model_pbs$se),
    z_value = estimate / std_error,
    p_value = 2 * pnorm(
      abs(z_value),
      lower.tail = FALSE
    ),
    conf_low = estimate - 1.96 * std_error,
    conf_high = estimate + 1.96 * std_error
  )


dirichlet_contrasts_pbs <- dirichlet_results_pbs |>
  filter(term != "(Intercept)") |>
  mutate(
    group = sub("^groups", "", term),

    group = recode(
      group,
      "naive"     = "Naive",
      "igg"       = "IgG",
      "antipd1"   = "anti-PD-1",
      "antictla4" = "anti-CTLA-4"
    ),

    cluster = factor(
      cluster,
      levels = as.character(0:10)
    ),

    group = factor(
      group,
      levels = c(
        "Naive",
        "IgG",
        "anti-PD-1",
        "anti-CTLA-4"
      )
    )
  )


dirichlet_contrasts_pbs <- dirichlet_contrasts_pbs |>
  mutate(
    padj_global = p.adjust(
      p_value,
      method = "BH"
    ),
    significance = case_when(
      padj_global < 0.0001 ~ "****",
      padj_global < 0.001  ~ "***",
      padj_global < 0.01   ~ "**",
      padj_global < 0.05   ~ "*",
      TRUE                 ~ ""
    )
  )



group_colors_pbs <- c(
  "Naive"       = "#F06A8A",
  "IgG"         = "#4C2A78",
  "anti-PD-1"   = "#6AAED6",
  "anti-CTLA-4" = "#5A86B8"
)

forest_pbs <- ggplot(
  dirichlet_contrasts_pbs,
  aes(
    x = estimate,
    y = fct_rev(cluster),
    colour = group
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_errorbarh(
    aes(
      xmin = conf_low,
      xmax = conf_high
    ),
    height = 0.16,
    position = position_dodge(width = 0.65),
    linewidth = 0.55
  ) +
  geom_point(
    aes(
      shape = padj_global < 0.05
    ),
    position = position_dodge(width = 0.65),
    size = 2.8
  ) +
  scale_colour_manual(
    values = group_colors_pbs,
    name = NULL
  ) +
  scale_shape_manual(
    values = c(
      `FALSE` = 1,
      `TRUE` = 16
    ),
    labels = c(
      `FALSE` = "Not significant",
      `TRUE` = "BH-adjusted p < 0.05"
    ),
    name = NULL
  ) +
  labs(
    x = "Dirichlet regression coefficient vs PBS (95% CI)",
    y = "Cluster",
    title = "Differential cellular composition",
    subtitle = "Positive coefficients indicate relative enrichment; negative coefficients indicate relative depletion"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

forest_pbs



heatmap_pbs <- ggplot(
  dirichlet_contrasts_pbs,
  aes(
    x = group,
    y = cluster,
    fill = estimate
  )
) +
  geom_tile(
    colour = "black",
    linewidth = 0.2
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    name = expression(beta)
  ) +
  labs(
    x = NULL,
    y = "Cluster"
  ) +
  coord_fixed(ratio = 0.75) +
  theme_classic(18) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

heatmap_pbs




dirichlet_contrasts_pbs$comp = "pbs"
dirichlet_contrasts$comp = "naive"

all_contrast = rbind(dirichlet_contrasts, dirichlet_contrasts_pbs)


ggplot(
  all_contrast,
  aes(x = group,y = cluster,fill = estimate)
) +
  geom_tile(colour = "black",linewidth = 0.2
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    name = expression(beta)
  ) +
  facet_grid(
    ~comp, switch = "x"
  ) +
  labs(
    x = NULL,
    y = "Cluster"
  ) +
  coord_fixed(ratio = 0.75) +
  theme_classic(18) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )


figure1C_group <- meta |>
  dplyr::count(
    groups,
    RNA_snn_res.0.2,
    name = "cell_count"
  ) |>
  dplyr::group_by(groups) |>
  dplyr::mutate(
    total_cells_group = sum(cell_count),
    proportion = cell_count / total_cells_group,
    percentage = proportion * 100
  ) |>
  dplyr::ungroup() |>
  dplyr::rename(
    group = groups,
    cluster = RNA_snn_res.0.2
  )


figure1C_group$group = factor(figure1C_group$group, levels = c("naive", "pbs", "igg", "antipd1", "antictla4"))

p_figure1C <- ggplot(
  figure1C_group,
  aes(
    x = group,
    y = percentage,
    fill = group
  )
) +
  geom_col(
    width = 0.75,
    colour = "white",
    linewidth = 0.15
  ) +
  geom_line(
    aes(group = 1),
    colour = "black",
    linewidth = 0.65
  ) +
  geom_point(
    colour = "black",
    size = 1.7
  ) +
  facet_wrap(
    ~ cluster,
    scales = "free_y",
    nrow = 2
  ) +
  scale_fill_manual(
    values = group_colors,
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = "Cell proportion (%)",
    fill = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    strip.background = element_rect(
      fill = "white",
      colour = "black"
    )
  )

p_figure1C
