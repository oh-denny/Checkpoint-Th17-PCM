
load_trajectory_population = function(path = 'data/combined_data_1.rds') {
  if (!file.exists(path)) stop('Missing annotated integrated Seurat input: ', path)
  obj <- readRDS(path)
  if (!'subtype_cell' %in% colnames(obj[[]])) {
    stop('The trajectory input must retain the original subtype_cell annotations.')
  }
  cells <- colnames(obj)[obj$subtype_cell %in% c('T cells CD3+', 'T cells CD4+/CD3+')]
  if (!length(cells)) stop('No cells match the original trajectory annotations.')
  obj <- subset(obj, cells = cells)
  Seurat::DefaultAssay(obj) <- 'RNA'
  obj
}
