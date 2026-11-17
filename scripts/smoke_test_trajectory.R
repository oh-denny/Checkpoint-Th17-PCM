# Read-only verification of the original trajectory population and score transfer.
suppressPackageStartupMessages({library(Seurat);library(monocle)})
source('functions/trajectory_input.r')
source('functions/trajectory.r')
trajectory_cells <- load_trajectory_population()
stopifnot(ncol(trajectory_cells) == 1908L)
stopifnot(sum(trajectory_cells$subtype_cell == 'T cells CD3+') == 123L,
          sum(trajectory_cells$subtype_cell == 'T cells CD4+/CD3+') == 1785L)
cat('PASS: original annotation-selected population contains 1,908 cells.\n')
scores <- calculate_fixed_module_scores(trajectory_cells, module_signatures)
stopifnot(nrow(scores) == 1908L, !anyNA(scores))
for (contrast in c('PBS_naive','IGG_pbs','PD1_pbs','CTLA4_pbs')) {
  cds <- readRDS(file.path('data',paste0('cds_condition_',contrast,'.rds')))
  stopifnot(setequal(colnames(cds),colnames(trajectory_cells)))
  cds <- transfer_fixed_scores(remove_old_score_columns(cds), scores, colnames(scores))
  stopifnot(!anyNA(pData(cds)[,colnames(scores),drop=FALSE]))
  cat('PASS:',contrast,'exact cell identity and complete fixed-score transfer.\n')
}
characterization <- readRDS('data/t_cells_CD4.rds')
stopifnot(ncol(characterization) == 6219L)
cat('PASS: the 6,219-cell CD4 characterization object is a distinct input.\n')
cat('No trajectory fitting or scientific output writes performed.\n')
