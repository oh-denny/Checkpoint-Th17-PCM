


ligands_5B <- c(
  "Il1b", "Il1a", "Tnf",
  "Cxcl16", "Mif",
  "Tgfb1", "Lgals3",
  "Spp1"
)


receptors_5B <- c(
  "Il1r1", "Il1rap",
  "Tnfrsf1a", "Tnfrsf1b",
  "Cxcr6",
  "Cd74",
  "Tgfbr1", "Tgfbr2",
  "Itgb1",
  "Cxcr4",
  "Il15ra", "Il2rb", "Il2rg"
)


library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

obj <- obj_nichenet

DefaultAssay(obj) <- "RNA"
Idents(obj) <- "cell_nichenet"

obj$cell_nichenet <- as.character(obj$cell_nichenet)
obj$cell_nichenet[obj$cell_nichenet == "Macrophage"] <- "Macrophages"
obj$cell_nichenet <- factor(
  obj$cell_nichenet,
  levels = c("Macrophages", "Neutrophils", "Th17/Treg")
)

obj$cell_condition <- paste(
  obj$cell_nichenet,
  obj$groups,
  sep = "_"
)

obj$cell_condition <- factor(
  obj$cell_condition,
  levels = c(
    "Macrophages_antictla4",
    "Macrophages_antipd1",
    "Macrophages_pbs",
    "Macrophages_naive",
    "Neutrophils_antictla4",
    "Neutrophils_antipd1",
    "Neutrophils_pbs",
    "Neutrophils_naive",
    "Th17/Treg_antictla4",
    "Th17/Treg_antipd1",
    "Th17/Treg_pbs",
    "Th17/Treg_naive"
  )
)

ligands_5B <- intersect(ligands_5B, rownames(obj))
receptors_5B <- intersect(receptors_5B, rownames(obj))

ligands_5B
receptors_5B


obj_ligands <- subset(
  obj,
  subset = cell_nichenet %in% c("Macrophages", "Neutrophils")
)

obj_ligands$cell_condition <- factor(
  obj_ligands$cell_condition,
  levels = c(
    "Macrophages_antictla4",
    "Macrophages_antipd1",
    "Macrophages_pbs",
    "Macrophages_naive",
    "Neutrophils_antictla4",
    "Neutrophils_antipd1",
    "Neutrophils_pbs",
    "Neutrophils_naive"
  )
)


obj_receptors <- subset(
  obj,
  subset = cell_nichenet == "Th17/Treg"
)

obj_receptors$cell_condition <- factor(
  obj_receptors$cell_condition,
  levels = c(
    "Th17/Treg_antictla4",
    "Th17/Treg_antipd1",
    "Th17/Treg_pbs",
    "Th17/Treg_naive"
  )
)


p_ligands <- DotPlot(
  obj_ligands,
  features = ligands_5B,
  group.by = "cell_condition",
  scale = FALSE
) +
  RotatedAxis() +
  labs(
    x = NULL,
    y = "Ligands",
    title = "Figure 5B",
    subtitle = "Expression of prioritized ligands in sender populations"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 10
    ),
    axis.text.y = element_text(
      face = "italic",
      size = 12
    ),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )





















library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

ligands_5B <- c("Il1b", "Il1a", "Tnf", "Cxcl16", "Mif", "Tgfb1", "Lgals3", "Spp1", "Slpi", "Il15", "Il15ra")
receptors_5B <- c("Il1r1", "Il1rap", "Tnfrsf1a", "Tnfrsf1b", "Cxcr6", "Cd74", "Tgfbr1", "Tgfbr2", "Itgb1", "Cxcr4", "Il15ra", "Il2rb", "Il2rg")

obj <- obj_nichenet
DefaultAssay(obj) <- "RNA"
obj$cell_nichenet <- as.character(obj$cell_nichenet)
obj$cell_nichenet[obj$cell_nichenet == "Macrophage"] <- "Macrophages"
groups_use <- c("naive", "pbs", "antipd1", "antictla4")
obj <- subset(obj, subset = groups %in% groups_use)

obj$cell_nichenet <- factor(obj$cell_nichenet, levels = c("Macrophages", "Neutrophils", "Th17/Treg"))

Idents(obj) <- "cell_nichenet"

obj$cell_condition <- paste(obj$cell_nichenet, obj$groups, sep = "_")

obj$cell_condition <- factor(obj$cell_condition, levels = c("Macrophages_antictla4", "Macrophages_antipd1", "Macrophages_pbs", "Macrophages_naive", "Neutrophils_antictla4", "Neutrophils_antipd1", "Neutrophils_pbs", "Neutrophils_naive", "Th17/Treg_antictla4", "Th17/Treg_antipd1", "Th17/Treg_pbs", "Th17/Treg_naive"))

table(obj$cell_condition, useNA = "always")

stopifnot(!anyNA(obj$cell_condition))

ligands_5B <- intersect(ligands_5B, rownames(obj[["RNA"]]))

receptors_5B <- intersect(receptors_5B, rownames(obj[["RNA"]]))

ligands_5B

receptors_5B

# Sender populations

obj_ligands <- subset(obj, subset = cell_nichenet %in% c("Macrophages", "Neutrophils"))

obj_ligands$cell_condition <- droplevels(obj_ligands$cell_condition)

Idents(obj_ligands) <- "cell_condition"

table(obj_ligands$cell_condition, useNA = "always")

# Receiver population

obj_receptors <- subset(obj, subset = cell_nichenet == "Th17/Treg")

obj_receptors$cell_condition <- droplevels(obj_receptors$cell_condition)

Idents(obj_receptors) <- "cell_condition"

table(obj_receptors$cell_condition, useNA = "always")

# Ligand expression

p_ligands <- DotPlot(obj_ligands, features = ligands_5B, assay = "RNA", group.by = "cell_condition", scale = FALSE, scale.by = "size", dot.scale = 7) +
scale_color_gradient(low = "grey90", high = "#B2182B") +
labs(x = NULL, y = NULL, title = "Sender populations", subtitle = "Expression of prioritized ligands", color = "Average expression", size = "Percent expressed") +
theme_classic(base_size = 14) +
theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 11), axis.text.y = element_text(size = 12, colour = "black"), plot.title = element_text(face = "bold"), legend.position = "right")

print(p_ligands)

# Receptor expression

p_receptors <- DotPlot(obj_receptors, features = receptors_5B, assay = "RNA", group.by = "cell_condition", scale = FALSE, scale.by = "size", dot.scale = 7) +
scale_color_gradient(low = "grey90", high = "#B2182B") +
labs(x = NULL, y = NULL, title = "Receiver population", subtitle = "Expression of corresponding receptors in Th17/Treg cells", color = "Average expression", size = "Percent expressed") +
theme_classic(base_size = 14) +
theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 11), axis.text.y = element_text(size = 12, colour = "black"), plot.title = element_text(face = "bold"), legend.position = "right")

print(p_receptors)

# Figure 5B

p_5B <- p_ligands / p_receptors + plot_layout(heights = c(2, 1))

print(p_5B)

# Save Figure 5B

dir.create("results/nichenet", recursive = TRUE, showWarnings = FALSE)

ggsave("results/nichenet/Figure5B.pdf", plot = p_5B, width = 14, height = 10)

ggsave("results/nichenet/Figure5B.png", plot = p_5B, width = 14, height = 10, dpi = 600)



max_expression <- max(c(p_ligands$data$avg.exp.scaled, p_receptors$data$avg.exp.scaled), na.rm = TRUE)

p_ligands <- p_ligands +
scale_size_area(max_size = 8, limits = c(0, 100), breaks = c(10, 25, 50, 75, 100), name = "Percent expressed") +
scale_color_gradient(low = "grey90", high = "#B2182B", limits = c(0, max_expression), name = "Average expression") +
geom_hline(yintercept = 4.5, colour = "grey60", linewidth = 0.4, linetype = "dashed")

p_receptors <- p_receptors +
scale_size_area(max_size = 8, limits = c(0, 100), breaks = c(10, 25, 50, 75, 100), name = "Percent expressed") +
scale_color_gradient(low = "grey90", high = "#B2182B", limits = c(0, max_expression), name = "Average expression")

p_5B <- (p_ligands / p_receptors) + plot_layout(heights = c(2, 1), guides = "collect")

p_5B



expr_il15 <- FetchData(obj_receptors, vars = c("Il15ra", "Il2rb", "Il2rg", "Cxcr6", "Il17a", "Foxp3", "groups"))

il15_summary <- expr_il15 %>%
group_by(groups) %>%
summarise(n_cells = n(), pct_Il15ra = mean(Il15ra > 0) * 100, pct_Il2rb = mean(Il2rb > 0) * 100, pct_Il2rg = mean(Il2rg > 0) * 100, pct_double_positive = mean(Il2rb > 0 & Il2rg > 0) * 100, pct_triple_positive = mean(Il15ra > 0 & Il2rb > 0 & Il2rg > 0) * 100, .groups = "drop")

il15_summary
