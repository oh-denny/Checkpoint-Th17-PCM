
library(CellChat)

run_cellchat <- function(obj) {

  DefaultAssay(obj) <- "RNA"

  data.input <- GetAssayData(
    obj,
    assay = "RNA",
    layer = "data"
  )

  meta <- obj@meta.data

  cellchat <- createCellChat(
    object = data.input,
    meta = meta,
    group.by = "comm_celltype"
  )

  cellchat@DB <- CellChatDB.mouse

  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)

  cellchat <- computeCommunProb(
    cellchat,
    type = "triMean"
  )

  cellchat <- filterCommunication(
    cellchat,
    min.cells = 10
  )

  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)

  return(cellchat)
}



get_comm_weight <- function(cc, source, target) {

  mat <- cc@net$weight

  if (
    source %in% rownames(mat) &&
    target %in% colnames(mat)
  ) {
    return(mat[source, target])
  }

  return(0)
}