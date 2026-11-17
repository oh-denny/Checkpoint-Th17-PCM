# Figure/panel: 4E
# Purpose: CellChat interaction-strength heatmaps.
# Inputs (under data/ unless session objects): session cc_pbs_clean, cc_pd1_clean, cc_ctl_clean.
# Packages: ComplexHeatmap, circlize, grid.
# Outputs: communication_heatmap.tiff.
# Prerequisites: Run figure_4A then figure_4B-D.
# Execution: repository root; see README and docs/SCRIPT_INVENTORY.md.

if (!exists("cc_pbs_clean", inherits = TRUE)) stop("Missing cc_pbs_clean. Run Figure 4A then 4B-D in this R session first.")
if (!exists("cc_pd1_clean", inherits = TRUE)) stop("Missing cc_pd1_clean. Run Figure 4A then 4B-D in this R session first.")
if (!exists("cc_ctl_clean", inherits = TRUE)) stop("Missing cc_ctl_clean. Run Figure 4A then 4B-D in this R session first.")



library(ComplexHeatmap)
library(circlize)
library(grid)

# run Figure 4B-D
mat_pbs = cc_pbs_clean@net$weight
mat_pd1 = cc_pd1_clean@net$weight
mat_ctl = cc_ctl_clean@net$weight
cell_order = c(
  "Macrophage",
  "Neutrophils",
  "Treg and Th17"
)
mat_pbs = mat_pbs[cell_order, cell_order]
mat_pd1 = mat_pd1[cell_order, cell_order]
mat_ctl = mat_ctl[cell_order, cell_order]
max_comm = max(
  mat_pbs,
  mat_pd1,
  mat_ctl,
  na.rm = TRUE
)
col_fun = colorRamp2(
  c(0, .2, .4),
  c("white", "#F4A6A6", "#B2182B")
)
ht_pbs = Heatmap(
  mat_pbs,
  name = "Interaction\nstrength",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_title = "Sender",
  column_title = "PBS",
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10)
)
ht_pd1 = Heatmap(
  mat_pd1,
  name = "Interaction\nstrength",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_title = NULL,
  column_title = "anti-PD-1",
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10)
)
ht_ctl = Heatmap(
  mat_ctl,
  name = "Interaction\nstrength",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_title = NULL,
  column_title = "anti-CTLA-4",
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10)
)

tiff("communication_heatmap.tiff",units = 'in', width = 8, height = 4, res = 300)
draw(
  ht_pbs + ht_pd1 + ht_ctl,
  merge_legends = TRUE
)
dev.off()

