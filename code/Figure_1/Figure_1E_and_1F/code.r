
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(emmeans)
library(forcats)
library(ggrepel)
library(tibble)

source("./functions/functions.r")

cfu_exp1 = read.csv("data/CFU/cfu_exp1.csv")
rds_path = file.path("data", "combined_data_1.rds")
combined_data_integrated <- readRDS(rds_path)

output_dir = file.path("results", "figure1_audit")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


meta = combined_data_integrated@meta.data |>
  rownames_to_column(
    "cell_barcode"
  ) |>
  mutate(
    groups = factor(
      groups,
      levels = group_levels
    ),

    RNA_snn_res.0.2 = factor(
      RNA_snn_res.0.2,
      levels = sort(
        unique(
          as.numeric(
            as.character(
              RNA_snn_res.0.2
            )
          )
        )
      )
    )
  )

sample_information = meta |>
  distinct(
    groups,
    cartig
  ) |>
  arrange(
    groups,
    cartig
  )

cluster_levels = levels(meta$RNA_snn_res.0.2)
animal_cluster = meta |>
  count(
    groups,
    cartig,
    RNA_snn_res.0.2,
    name = "cluster_cells"
  ) |>
  complete(
    nesting(
      groups,
      cartig
    ),
    RNA_snn_res.0.2 = cluster_levels,
    fill = list(
      cluster_cells = 0
    )
  ) |>
  group_by(
    groups,
    cartig
  ) |>
  mutate(
    total_cells = sum(
      cluster_cells
    ),

    proportion =
      cluster_cells /
      total_cells,

    percentage =
      proportion * 100
  ) |>
  ungroup()


proportion_check = animal_cluster |>
  group_by(
    groups,
    cartig
  ) |>
  summarise(
    sum_proportion = sum(
      proportion
    ),
    .groups = "drop"
  )

summary_effect = animal_cluster |>
  group_by(
    groups,
    RNA_snn_res.0.2
  ) |>
  summarise(
    mean_proportion = mean(
      proportion
    ),
    sd_proportion = sd(
      proportion
    ),
    median_proportion = median(
      proportion
    ),
    n_animals = n_distinct(
      cartig
    ),
    .groups = "drop"
  ) |>
  group_by(
    RNA_snn_res.0.2
  ) |>
  mutate(
    naive_mean = mean_proportion[
      groups == "naive"
    ],
    Log2FC = log2(
      (
        mean_proportion +
          pseudocount
      ) /
        (
          naive_mean +
            pseudocount
        )
    )
  ) |>
  ungroup()

cluster_wide = summary_effect |>
  filter(
    groups %in% c(
      "pbs",
      "igg",
      "antipd1",
      "antictla4"
    )
  ) |>
  select(
    RNA_snn_res.0.2,
    groups,
    Log2FC
  ) |>
  pivot_wider(
    names_from = groups,
    values_from = Log2FC
  )


rank_ctla4 = create_prioritization_ranking(
  cluster_wide = cluster_wide,
  treatment = "antictla4",
  comparison_groups = c(
    "pbs",
    "igg",
    "antipd1"
  )
)

rank_pd1 = create_prioritization_ranking(
  cluster_wide = cluster_wide,
  treatment = "antipd1",
  comparison_groups = c(
    "pbs",
    "igg",
    "antictla4"
  )
)

plot_ctla4 = plot_prioritization_ranking(
  ranking_data = rank_ctla4,
  treatment_label = "anti-CTLA-4",
  cluster_colors = cluster_colors,
  show_legend = TRUE
)

plot_pd1 = plot_prioritization_ranking(
  ranking_data = rank_pd1,
  treatment_label = "anti-PD-1",
  cluster_colors = cluster_colors,
  show_legend = TRUE
)

plot_ctla4
plot_pd1


ggsave(
  filename = file.path(
    output_dir,
    "cellular_prioritization_antictla4.tiff"
  ),
  plot = plot_ctla4,
  device = ragg::agg_tiff,
  width = 6,
  height = 6.5,
  units = "in",
  res = 600,
  compression = "lzw",
  background = "white"
)

ggsave(
  filename = file.path(
    output_dir,
    "cellular_prioritization_antipd1.tiff"
  ),
  plot = plot_pd1,
  device = ragg::agg_tiff,
  width = 6,
  height = 6.5,
  units = "in",
  res = 600,
  compression = "lzw",
  background = "white"
)
