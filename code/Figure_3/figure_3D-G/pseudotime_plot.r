# Figure/panel: 3D–G
# Purpose: Fixed mean-expression module scores along pseudotime.
# Inputs: four trajectory CDS objects and ordering lists; data/combined_data_1.rds for original annotation-based scores.
# Packages: patchwork, ragg.
# Outputs: results/trajectory_scores/ (plot and four scored CDS objects).
# Prerequisites: Historical trajectories; functions/trajectory.r; provenance review.
# Execution: repository root; see README and docs/SCRIPT_INVENTORY.md.

### Fixed module scores
source("./functions/trajectory.r")
library(patchwork)
trajectory_data_dir = getOption("paracocco.trajectory_data_dir", "data")
score_output_dir = getOption("paracocco.trajectory_score_output_dir", "results/trajectory_scores")
dir.create(score_output_dir, recursive = TRUE, showWarnings = FALSE)

cds_condition_PBS_naive = readRDS(file.path(trajectory_data_dir, "cds_condition_PBS_naive.rds"))
cds_condition_IGG_pbs = readRDS(file.path(trajectory_data_dir, "cds_condition_IGG_pbs.rds"))
cds_condition_PD1_pbs = readRDS(file.path(trajectory_data_dir, "cds_condition_PD1_pbs.rds"))
cds_condition_CTLA4_pbs = readRDS(file.path(trajectory_data_dir, "cds_condition_CTLA4_pbs.rds"))

genes_infection = readRDS(file.path(trajectory_data_dir, "genes_infection.rds"))
genes_pd1 = readRDS(file.path(trajectory_data_dir, "genes_pd1.rds"))
genes_ctla4 = readRDS(file.path(trajectory_data_dir, "genes_ctla4.rds"))
genes_igg = readRDS(file.path(trajectory_data_dir, "genes_igg.rds"))

source("functions/trajectory_input.r")
trajectory_cells = load_trajectory_population()

gene_lists = list(
  Infection   = unique(genes_infection),
  IgG         = unique(genes_igg),
  `anti-PD-1` = unique(genes_pd1),
  `anti-CTLA-4` = unique(genes_ctla4)
)

group_colors = c(
  "naive"     = "#F06A8A", 
  "pbs"       = "#6F7F5A",  
  "igg"       = "#4C2A78", 
  "antipd1"   = "#28331A", 
  "antictla4" = "#9BC5F3"  
) 


module_signatures = lapply(module_signatures, unique)
all_contrast_genes = unique(unlist(gene_lists, use.names = FALSE))
not_from_contrasts = lapply(module_signatures,function(genes) setdiff(genes, all_contrast_genes))


module_signatures_present = lapply(
  module_signatures,
  function(genes) intersect(genes, rownames(trajectory_cells))
)

signature_summary = tibble(
  Module = names(module_signatures_present),
  N_genes = lengths(module_signatures_present),
  Genes = vapply(
    module_signatures_present,
    paste,
    collapse = ", ",
    FUN.VALUE = character(1)
  )
)

fixed_module_scores = calculate_fixed_module_scores(
  seurat_object = trajectory_cells,
  modules = module_signatures_present,
  assay = "RNA",
  layer = "data"
)

fixed_score_columns = c("Less_differentiated_score","T_cell_activation_score","Migration_score","Th17_score","Treg_score")
print(colnames(fixed_module_scores))

cds_list = list(
  Infection = cds_condition_PBS_naive,
  IgG = cds_condition_IGG_pbs,
  PD1 = cds_condition_PD1_pbs,
  CTLA4 = cds_condition_CTLA4_pbs
)

for (cds_item in cds_list) {
  if (!setequal(colnames(cds_item), colnames(trajectory_cells))) {
    stop("Trajectory cell identities do not match the original annotation-based population.")
  }
}
cds_list = lapply(cds_list, remove_old_score_columns)
cds_list = lapply(
  cds_list,
  transfer_fixed_scores,
  score_metadata = fixed_module_scores,
  score_columns = fixed_score_columns
)

cds_condition_PBS_naive = cds_list$Infection
cds_condition_IGG_pbs = cds_list$IgG
cds_condition_PD1_pbs = cds_list$PD1
cds_condition_CTLA4_pbs = cds_list$CTLA4

score_check = lapply(
  cds_list,
  function(cds) {
    grep("_score$", colnames(pData(cds)), value = TRUE)
  }
)

score_check


scores_infection = prepare_fixed_scores(
  cds = cds_condition_PBS_naive,
  groups_keep = c("naive", "pbs")
)

scores_igg = prepare_fixed_scores(
  cds = cds_condition_IGG_pbs,
  groups_keep = c("pbs", "igg")
)

scores_pd1 = prepare_fixed_scores(
  cds = cds_condition_PD1_pbs,
  groups_keep = c("pbs", "antipd1")
)

scores_ctla4 = prepare_fixed_scores(
  cds = cds_condition_CTLA4_pbs,
  groups_keep = c("pbs", "antictla4")
)

plot_module_check = lapply(
  list(
    Infection = scores_infection,
    IgG = scores_igg,
    PD1 = scores_pd1,
    CTLA4 = scores_ctla4
  ),
  function(x) levels(x$Module)
)

plot_module_check

p_infection = plot_fixed_module_scores(
  plot_df = scores_infection,
  plot_title = "PBS vs Naive",
  group_colors = group_colors
)

p_igg = plot_fixed_module_scores(
  plot_df = scores_igg,
  plot_title = "IgG vs PBS",
  group_colors = group_colors
)

p_pd1 = plot_fixed_module_scores(
  plot_df = scores_pd1,
  plot_title = "Anti-PD-1 vs PBS",
  group_colors = group_colors
)

p_ctla4 = plot_fixed_module_scores(
  plot_df = scores_ctla4,
  plot_title = "Anti-CTLA-4 vs PBS",
  group_colors = group_colors
)


p_all_fixed_scores = (
  p_infection /
    p_igg /
    p_pd1 /
    p_ctla4
) +
  plot_layout(
    ncol = 1,
    guides = "collect",
    heights = rep(1, 4)
  ) +
  theme(
    legend.position = "top"
  )

print(p_all_fixed_scores)

ggsave(
  filename = file.path(score_output_dir, "fixed_module_scores_pseudotime.tiff"),
  plot = p_all_fixed_scores,
  device = ragg::agg_tiff,
  width = 14,
  height = 13,
  units = "in",
  res = 600,
  compression = "lzw",
  background = "white"
)

saveRDS(cds_condition_PBS_naive,file.path(score_output_dir, "cds_condition_PBS_naive_fixed_scores.rds"))
saveRDS(cds_condition_IGG_pbs,file.path(score_output_dir, "cds_condition_IGG_pbs_fixed_scores.rds"))
saveRDS(cds_condition_PD1_pbs,file.path(score_output_dir, "cds_condition_PD1_pbs_fixed_scores.rds"))
saveRDS(cds_condition_CTLA4_pbs,file.path(score_output_dir, "cds_condition_CTLA4_pbs_fixed_scores.rds"))