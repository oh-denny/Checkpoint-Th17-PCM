# Figure/panel: 4B–D
# Purpose: Group-specific CellChat networks.
# Inputs (under data/ unless session objects): session obj_cellchat_clean.
# Packages: CellChat, Seurat.
# Outputs: communication_pbs.tiff; communication_pd1.tiff; communication_ctla4.tiff.
# Prerequisites: Run figure_4A; functions/cellchat.r.
# Execution: repository root; see README and docs/SCRIPT_INVENTORY.md.

library(Seurat)
if (!exists("obj_cellchat_clean", inherits = TRUE)) stop("Missing obj_cellchat_clean. Run Figure 4A in this R session first.")

# run figure 4_A
library(CellChat)
source("./functions/cellchat.r")

Idents(obj_cellchat_clean) = "comm_celltype"
table(obj_cellchat_clean$comm_celltype)
table(obj_cellchat_clean$groups)

obj_pbs = subset(obj_cellchat_clean,subset = groups == "pbs")
obj_pd1 = subset(obj_cellchat_clean,subset = groups == "antipd1")
obj_ctl = subset(obj_cellchat_clean,subset = groups == "antictla4")

cc_pbs_clean = run_cellchat(obj_pbs)
cc_pd1_clean = run_cellchat(obj_pd1)
cc_ctl_clean = run_cellchat(obj_ctl)

subsetCommunication(cc_pbs_clean,sources.use = "Neutrophils",targets.use = "Treg and Th17")
subsetCommunication(cc_pd1_clean,sources.use = "Neutrophils",targets.use = "Treg and Th17")
subsetCommunication(cc_ctl_clean,sources.use = "Neutrophils",targets.use = "Treg and Th17")
subsetCommunication(cc_pbs_clean,sources.use = "Treg and Th17",targets.use = "Neutrophils")

cols_comm = c(
  "Macrophage"    = "#396544",
  "Neutrophils"   = "#E53637",
  "Treg and Th17" = "#2C7FB8"
)

max_weight = max(
  cc_pbs_clean@net$weight,
  cc_pd1_clean@net$weight,
  cc_ctl_clean@net$weight,
  na.rm = TRUE
)

par(mfrow = c(1, 3))
tiff("communication_pbs.tiff",units = 'in', width = 7.5, height = 5.5, res = 300)
netVisual_circle(
  cc_pbs_clean@net$weight,
  vertex.weight = as.numeric(table(cc_pbs_clean@idents)),
  weight.scale = TRUE,
  edge.weight.max = max_weight,
  vertex.label.cex = 0,
  label.edge = FALSE,
  vertex.label.color = "transparent",
  color.use = cols_comm,
  title.name = "PBS"
)
dev.off()

tiff("communication_pd1.tiff",units = 'in', width = 7.5, height = 5.5, res = 300)
netVisual_circle(
  cc_pd1_clean@net$weight,
  vertex.weight = as.numeric(table(cc_pd1_clean@idents)),
  weight.scale = TRUE,
  edge.weight.max = max_weight,
  vertex.label.cex = 0,
  label.edge = FALSE,
  vertex.label.color = "transparent",
  color.use = cols_comm,
  title.name = "anti-PD-1"
)
dev.off()

tiff("communication_ctla4.tiff",units = 'in', width = 7.5, height = 5.5, res = 300)
netVisual_circle(
  cc_ctl_clean@net$weight,
  vertex.weight = as.numeric(table(cc_ctl_clean@idents)),
  weight.scale = TRUE,
  edge.weight.max = max_weight,
  vertex.label.cex = 0,
  label.edge = FALSE,
  vertex.label.color = "transparent",
  color.use = cols_comm,
  title.name = "anti-CTLA-4"
)
dev.off()

par(mfrow = c(1, 1))
