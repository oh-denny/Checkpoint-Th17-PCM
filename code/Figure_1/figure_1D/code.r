# Figure/panel: 1D
# Purpose: Animal-level composition and Dirichlet regression.
# Inputs (under data/ unless session objects): combined_data_1.rds.
# Packages: DirichletReg, dplyr, forcats, ggplot2, ragg, scales, tidyr.
# Outputs: heatmap_dirichlet.tiff.
# Prerequisites: None; processed inputs required.
# Execution: repository root; see README and docs/SCRIPT_INVENTORY.md.

library(DirichletReg)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(scales)


rds_path = file.path("data", "combined_data_1.rds")
combined_data_integrated = readRDS(rds_path)
meta = combined_data_integrated@meta.data

meta$cluster = meta$RNA_snn_res.0.2 

counts_df = meta %>%
  group_by(cartig, groups, cluster) %>%
  summarise(count = n(), .groups = "drop")

wide_counts = counts_df %>%
  pivot_wider(names_from = cluster, values_from = count, values_fill = 0)

cluster_cols = setdiff(colnames(wide_counts), c("cartig", "groups"))
prop_matrix = as.matrix(wide_counts[, cluster_cols])
prop_matrix = prop_matrix / rowSums(prop_matrix)

Y = DR_data(prop_matrix)

dirichlet_data = data.frame(
  groups = factor(wide_counts$groups),
  cartig = wide_counts$cartig
)
dirichlet_data$Y = Y

dirichlet_data$groups = relevel(dirichlet_data$groups, ref = "naive")

model = DirichReg(Y ~ groups, data = dirichlet_data)
summary(model)

cluster_names = colnames(prop_matrix)
model_terms = c(
  "(Intercept)",
  paste0(
    "groups",
    levels(dirichlet_data$groups)[-1]
  )
)

dirichlet_results = tidyr::expand_grid(
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

dirichlet_contrasts = dirichlet_results |>
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

dirichlet_contrasts = dirichlet_contrasts |>
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

heatmap_dirichlet = ggplot(
  dirichlet_contrasts,
  aes(x = group,y = cluster,fill = estimate)) +
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
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.position = "none"
  )


ggsave(
  filename = "heatmap_dirichlet.tiff",
  plot = heatmap_dirichlet,
  device = ragg::agg_tiff,
  width = 6,
  height = 5,
  units = "in",
  res = 600,
  compression = "lzw",
  background = "white"
)

