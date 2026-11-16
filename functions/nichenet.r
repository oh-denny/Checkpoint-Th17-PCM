

library(purrr)
library(stringr)
library(readr)
library(Seurat)
library(dplyr)
library(tidyr)
library(scales)
library(Matrix)
library(ggplot2)
library(nichenetr)
library(edgeR)

library(OmnipathR)
library(dplyr)
library(VennDiagram)
library(readr)


get_expressed_genes_custom <- function(
    celltype,
    condition,
    pct = expression_pct
) {

  cells <- rownames(meta)[
    meta$cell_nichenet == celltype &
    meta$groups == condition
  ]

  if (length(cells) == 0) {
    return(character(0))
  }

  mat <- expr[, cells, drop = FALSE]

  detected <- Matrix::rowMeans(mat > 0)

  names(detected)[detected >= pct]
}


get_expressed_lr <- function(
    sender,
    condition,
    receiver = "Th17/Treg",
    pct = expression_pct
) {

  expressed_sender <- get_expressed_genes_custom(
    sender,
    condition,
    pct
  )

  expressed_receiver <- get_expressed_genes_custom(
    receiver,
    condition,
    pct
  )

  expressed_ligands <- intersect(
    expressed_sender,
    ligands
  )

  expressed_receptors <- intersect(
    expressed_receiver,
    receptors
  )

  lr_filtered <- lr_Network_Nichenet %>%
    filter(
      from %in% expressed_ligands,
      to %in% expressed_receptors,
      from %in% colnames(ligand_target_matrix)
    )

  # Additional criterion for IL-15 signaling

  lr_filtered <- lr_filtered %>%
    filter(
      from != "Il15" |
      (
        all(c("Il2rb", "Il2rg") %in% expressed_receiver) &
        to %in% c("Il2rb", "Il2rg", "Il15ra")
      )
    )

  lr_filtered %>%
    distinct(from, to)
}



get_target_genes <- function(
    deg,
    top_n = 100,
    min_logFC = 0.25
) {

  deg %>%
    filter(
      log2FoldChange > min_logFC,
      Gene %in% rownames(ligand_target_matrix)
    ) %>%
    arrange(desc(log2FoldChange)) %>%
    slice_head(n = top_n) %>%
    pull(Gene) %>%
    unique()
}




run_pseudobulk_edger <- function(
    pb,
    condition,
    reference
) {

  # Select samples
  md <- pb$metadata %>%
    filter(
      groups %in% c(condition, reference)
    )

  md$groups <- factor(
    md$groups,
    levels = c(reference, condition)
  )

  # Check biological replicates
  if (any(table(md$groups) < 2)) {

    stop(
      paste(
        "Insufficient biological replicates:",
        condition,
        "vs",
        reference
      )
    )
  }

  # Select pseudobulk counts
  mat <- pb$counts[
    ,
    rownames(md),
    drop = FALSE
  ]

  # Check count validity
  if (any(!is.finite(mat))) {
    stop("Non-finite values detected.")
  }

  if (any(mat < 0)) {
    stop("Negative values detected.")
  }

  if (any(colSums(mat) == 0)) {
    stop("Samples with zero library size detected.")
  }

  # Report fractional counts
  if (any(abs(mat - round(mat)) > 1e-8)) {

    message(
      "Fractional counts detected. ",
      "Preserving original values for edgeR."
    )
  }

  # Design matrix
  design <- model.matrix(
    ~ groups,
    data = md
  )

  # DGEList
  y <- edgeR::DGEList(
    counts = mat,
    group = md$groups
  )

  # Filter lowly expressed genes
  keep <- edgeR::filterByExpr(
    y,
    design = design
  )

  y <- y[
    keep,
    ,
    keep.lib.sizes = FALSE
  ]

  # TMM normalization
  y <- edgeR::calcNormFactors(y)

  # Estimate dispersion
  y <- edgeR::estimateDisp(
    y,
    design,
    robust = TRUE
  )

  # Quasi-likelihood model
  fit <- edgeR::glmQLFit(
    y,
    design,
    robust = TRUE
  )

  qlf <- edgeR::glmQLFTest(
    fit,
    coef = 2
  )

  # Differential expression results
  results <- edgeR::topTags(
    qlf,
    n = Inf
  )$table %>%
    as.data.frame() %>%
    tibble::rownames_to_column(
      var = "Gene"
    ) %>%
    dplyr::rename(
      log2FoldChange = logFC,
      padj = FDR,
      pvalue = PValue
    )

  return(results)
}






run_nichenet <- function(
    sender,
    comparison,
    potential_ligands,
    geneset,
    background,
    top_ligands = 20,
    top_targets = 250
) {

  message(
    "Running: ",
    sender,
    " | ",
    comparison
  )

  background <- intersect(
    background,
    rownames(ligand_target_matrix)
  )

  geneset <- intersect(
    geneset,
    background
  )

  potential_ligands <- intersect(
    potential_ligands,
    colnames(ligand_target_matrix)
  )

  if (
    length(geneset) < 10 ||
    length(geneset) >= length(background) ||
    length(potential_ligands) == 0
  ) {
    warning("Insufficient genes or ligands: ", comparison)
    return(NULL)
  }

  # Predict ligand activities

  ligand_activities <- predict_ligand_activities(
    geneset = geneset,
    background_expressed_genes = background,
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = potential_ligands
  ) %>%
    arrange(desc(aupr_corrected)) %>%
    mutate(rank = row_number())

  # Select prioritized ligands

  best_ligands <- ligand_activities %>%
    filter(aupr_corrected > 0) %>%
    slice_head(n = top_ligands) %>%
    pull(test_ligand)

  if (length(best_ligands) == 0) {
    return(NULL)
  }

  # Predict ligand-target links

  ligand_target_links <- map_dfr(

    best_ligands,

    function(ligand) {

      get_weighted_ligand_target_links(
        ligand = ligand,
        geneset = geneset,
        ligand_target_matrix = ligand_target_matrix,
        n = top_targets
      )
    }
  ) %>%
    drop_na() %>%
    filter(weight > 0) %>%
    mutate(
      sender = sender,
      receiver = "Th17/Treg",
      comparison = comparison
    ) %>%
    left_join(
      ligand_activities %>%
        select(test_ligand, aupr_corrected, rank),
      by = c("ligand" = "test_ligand")
    )

  return(
    list(
      ligand_activity = ligand_activities,
      ligand_target = ligand_target_links
    )
  )
}

