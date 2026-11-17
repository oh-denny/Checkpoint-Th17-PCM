# Figure/panel: 3C
# Purpose: DE ordering genes and DDRTree trajectories.
# Inputs: data/combined_data_1.rds; original CD3+ and CD4+/CD3+ annotations.
# Packages: Matrix, Seurat, VGAM, monocle, openxlsx, ragg.
# Outputs: results/trajectory/ (four ordering lists, four contrast-specific CDS objects, plots and workbook).
# Prerequisites: functions/functions.r, trajectory_input.r and patch/patch.r; fresh output directory.
# Execution: repository root; see README and docs/SCRIPT_INVENTORY.md.

source("functions/functions.r")
source("functions/patch/patch.r")
source("functions/trajectory_input.r")
patch_monocle2()

library(Seurat)
library(monocle)

group_colors = c(
  "naive"     = "#F06A8A", 
  "pbs"       = "#6F7F5A",  
  "igg"       = "#4C2A78", 
  "antipd1"   = "#28331A", 
  "antictla4" = "#9BC5F3"  
) 

SC_condition = load_trajectory_population()
output_dir = getOption("paracocco.trajectory_output_dir", "results/trajectory")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

DefaultAssay(SC_condition) = "RNA"
counts_matrix = GetAssayData(
  SC_condition,
  assay = "RNA",
  layer = "counts"
)

cell_metadata = SC_condition@meta.data[
  colnames(counts_matrix),
  ,
  drop = FALSE
]

pd = new(
  "AnnotatedDataFrame",
  data = cell_metadata
)

gene_metadata = data.frame(
  gene_short_name = rownames(counts_matrix),
  row.names = rownames(counts_matrix),
  stringsAsFactors = FALSE
)

fd = new(
  "AnnotatedDataFrame",
  data = gene_metadata
)

cds = newCellDataSet(
  counts_matrix,
  phenoData = pd,
  featureData = fd,
  lowerDetectionLimit = 1,
  expressionFamily = VGAM::negbinomial.size()
)

patch_monocle2()

DefaultAssay(SC_condition) = "RNA"
Idents(SC_condition) = "groups"
table(Idents(SC_condition))

de_pbs_naive = FindMarkers(
  SC_condition,
  ident.1 = "pbs",
  ident.2 = "naive",
  assay = "RNA",
  layer = "data",
  min.pct = 0.10,
  logfc.threshold = 0.10,
  only.pos = FALSE
)

de_pd1_pbs = FindMarkers(
  SC_condition,
  ident.1 = "antipd1",
  ident.2 = "pbs",
  assay = "RNA",
  layer = "data",
  min.pct = 0.10,
  logfc.threshold = 0.10,
  only.pos = FALSE
)

de_ctla4_pbs = FindMarkers(
  SC_condition,
  ident.1 = "antictla4",
  ident.2 = "pbs",
  assay = "RNA",
  layer = "data",
  min.pct = 0.10,
  logfc.threshold = 0.10,
  only.pos = FALSE
)

de_igg_pbs = FindMarkers(
  SC_condition,
  ident.1 = "igg",
  ident.2 = "pbs",
  assay = "RNA",
  layer = "data",
  min.pct = 0.10,
  logfc.threshold = 0.10,
  only.pos = FALSE
)

genes_infection = select_de_genes(
  de_pbs_naive,
  n_up = 50,
  n_down = 40
)

genes_pd1 = select_de_genes(
  de_pd1_pbs,
  n_up = 40,
  n_down = 40
)

genes_ctla4 = select_de_genes(
  de_ctla4_pbs,
  n_up = 40,
  n_down = 40
)

genes_igg = select_de_genes(
  de_igg_pbs,
  n_up = 40,
  n_down = 40
)

saveRDS(genes_infection, file.path(output_dir, "genes_infection.rds"))
saveRDS(genes_pd1, file.path(output_dir, "genes_pd1.rds"))
saveRDS(genes_ctla4, file.path(output_dir, "genes_ctla4.rds"))
saveRDS(genes_igg, file.path(output_dir, "genes_igg.rds"))

ordering_lists = list(
  PBS_naive = genes_infection,
  IGG_pbs = genes_igg,
  PD1_pbs = genes_pd1,
  CTLA4_pbs = genes_ctla4
)
for (contrast in names(ordering_lists)) {
  ordering_genes_condition = intersect(ordering_lists[[contrast]], rownames(cds))
  set.seed(1234)
  cds_condition = cds
  cds_condition = estimateSizeFactors(cds_condition)
  cds_condition = estimateDispersions(cds_condition)
  cds_condition = setOrderingFilter(cds_condition, ordering_genes_condition)
  cds_condition = reduceDimension(
    cds_condition, max_components = 2, method = "DDRTree", ncenter = 30
  )
  cds_condition = orderCells(cds_condition)
  pData(cds_condition)$groups = factor(
    pData(cds_condition)$groups, levels = names(group_colors)
  )
  saveRDS(cds_condition, file.path(output_dir, paste0("cds_condition_", contrast, ".rds")))
}


cds_condition_pbs = readRDS(file.path(output_dir, "cds_condition_PBS_naive.rds"))
cds_condition_igg = readRDS(file.path(output_dir, "cds_condition_IGG_pbs.rds"))
cds_condition_pd1 = readRDS(file.path(output_dir, "cds_condition_PD1_pbs.rds"))
cds_condition_ctl = readRDS(file.path(output_dir, "cds_condition_CTLA4_pbs.rds"))


p_pseudotime = monocle::plot_cell_trajectory(
  cds_condition_ctl,
  color_by = "Pseudotime"
) +
  scale_color_viridis_c(option = "magma", direction = 1) +
  theme_classic(base_size = 20) +
  labs(color = "Pseudotime") +
  theme(
    legend.position = 'bottom'
  )

p_pseudotime

ggsave(
  filename = file.path(output_dir, "ctl_pbs.tiff"),
  plot = p_pseudotime,
  device = ragg::agg_tiff,
  width = 6,
  height = 7,
  units = "in",
  res = 600,
  compression = "lzw",
  background = "white"
)

p_state = monocle::plot_cell_trajectory(
  cds_condition_ctl,
  color_by = "State"
) +
  theme_classic(base_size = 20) +
  labs(color = "State") +
   theme(
    legend.position = 'bottom'
  )

p_state

ggsave(
  filename = file.path(output_dir, "ctl_pbs_p_state.tiff"),
  plot = p_state,
  device = ragg::agg_tiff,
  width = 6,
  height = 8,
  units = "in",
  res = 600,
  compression = "lzw",
  background = "white"
)



# Supp table
# Tables are exported with openxlsx below; no writexl calls are used.

max_len <- max(
  length(genes_infection),
  length(genes_pd1),
  length(genes_ctla4),
  length(genes_igg)
)

genes_table <- data.frame(
  Infection = c(
    genes_infection,
    rep(NA, max_len - length(genes_infection))
  ),
  
  anti_PD1 = c(
    genes_pd1,
    rep(NA, max_len - length(genes_pd1))
  ),
  
  anti_CTLA4 = c(
    genes_ctla4,
    rep(NA, max_len - length(genes_ctla4))
  ),
  
  IgG = c(
    genes_igg,
    rep(NA, max_len - length(genes_igg))
  )
)

openxlsx::write.xlsx(
  genes_table,
  file.path(output_dir, "genes_trajectory.xlsx")
)
