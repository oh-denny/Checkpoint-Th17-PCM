# Figure/panel: 3D–G companion
# Purpose: Bootstrap relative effects for selected CD4 genes.
# Inputs (under data/ unless session objects): t_cells_CD4.rds.
# Packages: inherited through prerequisite/helper; see inventory.
# Outputs: CD4_relative_effects.tiff.
# Prerequisites: functions/rel.r; REffect installed in advance.
# Execution: repository root; see README and docs/SCRIPT_INVENTORY.md.



if (!requireNamespace("REffect", quietly = TRUE)) {
  stop("Install REffect from oh-denny/relative_effect before running this script.")
}
source("./functions/rel.r")
t_cells_CD4 = readRDS("./data/t_cells_CD4.rds")

genes_use = c(
  "Sell",
  "Foxp3",
  "Pdcd1",
  "Il10",
  "Cxcr6",
  "Ccr2",
  "Il17a",
  "Ctla4",
  "Cd44"
)

DefaultAssay(t_cells_CD4) = "RNA"
genes_present = intersect(genes_use,rownames(t_cells_CD4[["RNA"]]))
df_expression = FetchData(
  object = t_cells_CD4,
  vars = c(genes_present, "groups"),
  layer = "data"
)
head(df_expression)
table(df_expression$groups)

comparisons = list(
  Infection = list(
    reference = "naive",
    treatment = "pbs",
    colors = c("#E88EAE", "#778973")
  ),
  IgG = list(
    reference = "pbs",
    treatment = "igg",
    colors = c("#778973", "#8755B1")
  ),
  PD1 = list(
    reference = "pbs",
    treatment = "antipd1",
    colors = c("#778973", "#34482A")
  ),
  CTLA4 = list(
    reference = "pbs",
    treatment = "antictla4",
    colors = c("#778973", "#91BFF0")
  )
)

set.seed(123)

relative_results = lapply(
  comparisons,
  function(comp) {

    run_relative_effect(
      data = df_expression,
      reference = comp$reference,
      treatment = comp$treatment,
      genes = genes_present,
      colors = comp$colors,
      n_boot = 1000
    )

  }
)

relative_theme = theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(
      color = "black"
    ),
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.grid = element_blank()
  )

p_infection = relative_results$Infection$plot
p_igg       = relative_results$IgG$plot
p_pd1       = relative_results$PD1$plot
p_ctla4     = relative_results$CTLA4$plot

p_infection = p_infection + relative_theme
p_igg       = p_igg + relative_theme
p_pd1       = p_pd1 + relative_theme
p_ctla4     = p_ctla4 + relative_theme


p_final = (
  p_infection /
  p_igg /
  p_pd1 /
  p_ctla4
)

p_final

ggsave(
  filename = "CD4_relative_effects.tiff",
  plot = p_final,
  width = 4,
  height = 12,
  units = "in",
  dpi = 600,
  compression = "lzw"
)
