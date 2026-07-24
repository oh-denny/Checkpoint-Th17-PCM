library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)


group_colors = c(
  "naive"     = "#F06A8A", 
  "pbs"       = "#6F7F5A",  
  "igg"       = "#4C2A78", 
  "antipd1"   = "#28331A", 
  "antictla4" = "#9BC5F3"  
)


rds_path = file.path("data", "combined_data_1.rds")
combined_data_integrated = readRDS(rds_path)

output_dir = file.path("results", "figure1_audit")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

meta = combined_data_integrated@meta.data |>
  tibble::rownames_to_column("cell_barcode")

meta |>
  distinct(groups, cartig) |>
  arrange(groups, cartig)


cluster_levels = levels(meta$RNA_snn_res.0.2)

barplot_animal = meta |>
  count(
    groups,
    cartig,
    RNA_snn_res.0.2,
    name = "total"
  ) |>
  complete(
    nesting(groups, cartig),
    RNA_snn_res.0.2 = cluster_levels,
    fill = list(total = 0)
  ) |>
  mutate(
    groups = factor(
      groups,
      levels = c(
        "naive",
        "pbs",
        "igg",
        "antipd1",
        "antictla4"
      )
    )
  )
barplot_animal$groups = factor(
  barplot_animal$groups,
  levels = c("naive", "pbs", "igg", "antipd1", "antictla4")
)

barplot_animal$RNA_snn_res.0.2 = factor(
  barplot_animal$RNA_snn_res.0.2,
  levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
)

summary_effect = barplot_animal |>
  group_by(groups, RNA_snn_res.0.2) |>
  summarise(
    mean = mean(total),
    sd = sd(total),
    median = median(total),
    .groups = "drop"
  )

naive_mean = summary_effect |>
  filter(groups == "naive") |>
  select(
    RNA_snn_res.0.2,
    naive_mean = mean
  )

summary_effect = summary_effect |>
  left_join(
    naive_mean,
    by = "RNA_snn_res.0.2"
  ) |>
  mutate(
    FoldChange = mean / naive_mean,
    Log2FC = log2(FoldChange)
  )

summary_effect

fc_labels = summary_effect |>
  filter(groups != "naive") |>
  mutate(
    label = sprintf("%.2fx", FoldChange)
  )

summary_effect = summary_effect |>
  mutate(
    groups = factor(
      groups,
      levels = c("naive","pbs", "igg", "antipd1", "antictla4"),
      labels = c("Naive","PBS", "IgG", "anti-PD-1", "anti-CTLA-4")
    )
  )

group_colors = c(
  "Naive"     = "#F06A8A", 
  "PBS"       = "#6F7F5A",  
  "IgG"       = "#4C2A78", 
  "anti-PD-1"   = "#28331A", 
  "anti-CTLA-4" = "#9BC5F3"  
) 

barplot_logFC = ggplot(summary_effect, aes(groups, Log2FC, fill = groups)) +
    geom_col() +
    facet_grid(~RNA_snn_res.0.2, switch = 'y') +
    scale_fill_manual(
            values = group_colors,
            drop = FALSE
        ) +
    theme_classic(18) +
    theme(
        legend.position = 'bottom',
        axis.text.x = element_text(angle = 45, hjust = 1)
    )

ggsave(
  filename = file.path(
    output_dir,
    "barplot_logFC.tiff"
  ),
  plot = barplot_logFC,
  device = ragg::agg_tiff,
  width = 20,
  height = 5.5,
  units = "in",
  res = 600,
  compression = "lzw",
  background = "white"
)

###
library(tidyverse)
library(Seurat)
library(ggpubr)

combined_data_integrated$cluster_id = combined_data_integrated$RNA_snn_res.0.2
combined_data_integrated$original_animal = mapping_exp1$original_animal[
  match(combined_data_integrated$cartig, mapping_exp1$cartig)
]
table(combined_data_integrated$cartig, combined_data_integrated$original_animal, useNA = "ifany")
cell_proportions = combined_data_integrated@meta.data %>%
  filter(!is.na(original_animal)) %>%
  group_by(original_animal, groups, cluster_id) %>%
  summarise(count = n(), .groups = "drop_last") %>%
  mutate(
    total_cells = sum(count),
    percentage = (count / total_cells) * 100
  ) %>%
  ungroup()

data_for_cor = cell_proportions %>%
  left_join(cfu_exp1, by = c("original_animal", "groups"))

pro_cell = ggplot(data_for_cor, 
       aes(x = percentage, y = log10_cfu_g)) +
  geom_point(aes(color = groups), size = 5) +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = 'gray80') +
  facet_wrap(~cluster_id, scales = "free_x", nrow = 1) +
  theme_classic(18) +
  labs(
    x = "Cell proportion (%)",
    y = expression(Log[10] * " CFU/g"),
  ) +
    scale_color_manual(
    values = group_colors,
    drop = FALSE
  ) +
  theme(
    legend.position = 'botom'
  )

ggsave(
  filename = file.path(
    output_dir,
    "pro_cell.tiff"
  ),
  plot = pro_cell,
  device = ragg::agg_tiff,
  width = 20,
  height = 3,
  units = "in",
  res = 600,
  compression = "lzw",
  background = "white"
)
