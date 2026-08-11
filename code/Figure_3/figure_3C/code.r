source("functions/functions.r")
source("functions/patch/patch.r")
patch_monocle2()

group_colors = c(
  "naive"     = "#F06A8A", 
  "pbs"       = "#6F7F5A",  
  "igg"       = "#4C2A78", 
  "antipd1"   = "#28331A", 
  "antictla4" = "#9BC5F3"  
) 

SC_obj_sub = t_cells_CD4
DefaultAssay(SC_obj_sub) = "RNA"
Idents(SC_obj_sub) = "seurat_clusters"

counts_matrix = GetAssayData(
  SC_obj_sub,
  assay = "RNA",
  layer = "counts"
)

cell_metadata = SC_obj_sub@meta.data[
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
SC_condition = t_cells_CD4

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
  n_down = 50
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

ordering_genes_condition = intersect(
  genes_infection,
  rownames(cds)
) # 1
 
ordering_genes_condition = intersect(
  genes_igg,
  rownames(cds)
) # 2

ordering_genes_condition = intersect(
  genes_pd1,
  rownames(cds)
) # 3

ordering_genes_condition = intersect(
  genes_ctla4,
  rownames(cds)
) # 4

cds_condition = cds
cds_condition = setOrderingFilter(
  cds_condition,
  ordering_genes_condition
)

set.seed(1234)

cds_condition = reduceDimension(
  cds_condition,
  max_components = 2,
  method = "DDRTree"
)
cds_condition = orderCells(cds_condition)
pData(cds_condition)$groups = factor(
  pData(cds_condition)$groups,
  levels = names(group_colors)
)

saveRDS(cds_condition, "./data/cds_condition_PBS_naive.rds")
saveRDS(cds_condition, "./data/cds_condition_IGG_pbs.rds")
saveRDS(cds_condition, "./data/cds_condition_PD1_pbs.rds")
saveRDS(cds_condition, "./data/cds_condition_CTLA4_pbs.rds")


cds_condition_pbs = readRDS("./data/cds_condition_PBS_naive.rds")
cds_condition_igg = readRDS("./data/cds_condition_IGG_pbs.rds")
cds_condition_pd1 = readRDS("./data/cds_condition_PD1_pbs.rds")
cds_condition_ctl = readRDS("./data/cds_condition_CTLA4_pbs.rds")


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
  filename = "ctl_pbs.tiff",
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
  filename = "ctl_pbs_p_state.tiff",
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
library(writexl)

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
  "~/Desktop/genes_trajectory.xlsx"
)
