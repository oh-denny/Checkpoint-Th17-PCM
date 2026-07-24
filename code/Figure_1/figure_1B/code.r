
library(emmeans)
library(ggplot2)
library(purrr)
library(stringr)
library(readr)
library(Seurat)
library(dplyr)
library(tidyr)

rds_path = file.path("data", "combined_data_1.rds")
combined_data_integrated = readRDS(rds_path)

output_dir = file.path("results", "figure1_audit")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

meta = combined_data_integrated@meta.data |>
  tibble::rownames_to_column("cell_barcode")

mapping_exp1 = read.csv("data/CFU/mapping_exp1.csv")
cfu_exp1 = read.csv("data/CFU/cfu_exp1.csv")

meta = combined_data_integrated@meta.data |>
  tibble::rownames_to_column("cell_barcode") |>
  mutate(
    cluster = as.character(RNA_snn_res.0.2),
    cartig = tolower(cartig),
    groups = tolower(groups)
  )

animal_map_scrna = meta |>
  distinct(cartig, groups) |>
  arrange(groups, cartig)

animal_map_scrna
animal_map_scrna |>
  count(groups)

cluster_prop_animal = meta |>
  count(
    cartig,
    groups,
    cluster,
    name = "cluster_cells"
  ) |>
  complete(
    nesting(cartig, groups),
    cluster,
    fill = list(cluster_cells = 0)
  ) |>
  group_by(cartig, groups) |>
  mutate(
    total_cells_animal = sum(cluster_cells),
    proportion = cluster_cells / total_cells_animal,
    percentage = 100 * proportion
  ) |>
  ungroup()


cfu_scrna_exp1 = mapping_exp1 |>
  left_join(
    cfu_exp1,
    by = c("original_animal", "groups")
  )

cfu_scrna_exp1

cluster_cfu = cluster_prop_animal |>
  filter(groups != "naive") |>
  left_join(
    cfu_scrna_exp1,
    by = c("cartig", "groups")
  )

cluster_cfu |>
  distinct(
    cartig,
    groups,
    original_animal,
    log10_cfu_g
  ) |>
  arrange(groups, cartig)

fit = lm(log10_cfu_g ~ groups,data = cluster_cfu)
emmeans::emmeans(fit, pairwise ~ groups)

group_colors = c(
  "PBS"       = "#6F7F5A",  
  "IgG"       = "#4C2A78", 
  "anti-PD-1"   = "#28331A", 
  "anti-CTLA-4" = "#9BC5F3"  
)

emm = emmeans(fit, ~ groups)
emm_df = as.data.frame(emm)
emm_df$groups = factor(
  emm_df$groups,
  levels = c("pbs", "igg", "antipd1", "antictla4"),
  labels = c("PBS", "IgG", "anti-PD-1", "anti-CTLA-4")
)

forestplot = ggplot(
  emm_df, aes(x = groups,y = emmean, color = groups)) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = lower.CL,ymax = upper.CL),
    width = 0.2,
    linewidth = 0.9) +
  labs(x = NULL,y = expression(Log[10] * " CFU/g")) +
  theme_classic(22) +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    ),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  ) +
  scale_color_manual(values = group_colors)


ggsave(
  filename = "forestplot.tiff",
  plot = forestplot,
  device = ragg::agg_tiff,
  width = 6,
  height = 5,
  units = "in",
  res = 600,
  compression = "lzw",
  background = "white"
)
