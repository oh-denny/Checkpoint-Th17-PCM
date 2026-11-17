
suppressPackageStartupMessages({
  library(Seurat)
  library(monocle)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
})

module_levels = c(
  "Less_differentiated",
  "T_cell_activation",
  "Migration",
  "Th17",
  "Treg"
)

module_labels = c(
  "Less differentiated",
  "T-cell activation",
  "Migration",
  "Th17",
  "Treg"
)

module_signatures <- list(

  Th17 = c(
    "Il17a",
    "Maf",
    "Hif1a",
    "Bhlhe40"
  ),

  Treg = c(
    "Ctla4",
    "Il2ra",
    "Tmem176a",
    "Tmem176b",
    "Tgfbr2"
  ),

  T_cell_activation = c(
    "Tnfrsf1b",
    "Icos",
    "Cyld",
    "Nr4a1",
    "Itk",
    "Tnfaip3",
    "Tnfsf8",
    "Cd28",
    "Zap70",
    "Lat"
  ),

  Migration = c(
    "Rgs1",
    "Cxcr4",
    "Ccr2",
    "Ccr7",
    "S1pr1",
    "Spn",
    "Adam8",
    "Itgb1",
    "Itga4",
    "S1pr4",
    "Rac2",
    "Ripor2",
    "Cxcr6",
    "S100a4",
    "Ccl5"
  ),

  Less_differentiated = c(
    "Lef1",
    "Klf3",
    "Sell",
    "Foxp1",
    "S1pr1",
    "Klf2",
    "Ccr7",
    "Tob1",
    "Tcf7",
    "Tob2",
    "Zfp36l2",
    "Btg2"
  )
)



# ============================================================
# 3. Calculate the SAME five scores once for all cells
# ============================================================

calculate_fixed_module_scores <- function(
    seurat_object,
    modules,
    assay = "RNA",
    layer = "data"
) {

  assay_names <- names(seurat_object@assays)

  if (!(assay %in% assay_names)) {
    stop(
      "Assay '", assay, "' não encontrado. Assays disponíveis: ",
      paste(assay_names, collapse = ", ")
    )
  }

  if (!is.list(modules)) {
    stop("'modules' precisa ser uma lista nomeada.")
  }

  if (is.null(names(modules)) || any(names(modules) == "")) {
    stop("'modules' precisa ter nomes para todas as assinaturas.")
  }

  expression_matrix <- SeuratObject::LayerData(
    object = seurat_object,
    assay = assay,
    layer = layer
  )

  if (is.null(expression_matrix)) {
    stop(
      "A layer '", layer,
      "' não foi encontrada no assay '", assay, "'."
    )
  }

  if (nrow(expression_matrix) == 0 || ncol(expression_matrix) == 0) {
    stop("A matriz de expressão selecionada está vazia.")
  }

  available_genes <- rownames(expression_matrix)

  if (is.null(available_genes)) {
    stop("A matriz de expressão não possui nomes de genes.")
  }

  scores <- data.frame(
    row.names = colnames(expression_matrix)
  )

  for (module_name in names(modules)) {

    module_genes <- unique(
      as.character(
        unlist(
          modules[[module_name]],
          use.names = FALSE
        )
      )
    )

    module_genes <- module_genes[
      !is.na(module_genes) &
        nzchar(module_genes)
    ]

    genes_use <- module_genes[
      module_genes %in% available_genes
    ]

    message(
      module_name,
      ": ",
      length(genes_use),
      " genes utilizados — ",
      paste(genes_use, collapse = ", ")
    )

    if (length(genes_use) == 0) {
      stop(
        "Nenhum gene da assinatura '",
        module_name,
        "' foi encontrado na matriz."
      )
    }

    score_name <- paste0(
      module_name,
      "_score"
    )

    scores[[score_name]] <- Matrix::colMeans(
      expression_matrix[
        genes_use,
        ,
        drop = FALSE
      ]
    )
  }

  if (anyNA(scores)) {
    stop("Foram produzidos valores NA no cálculo dos scores.")
  }

  scores
}




# ============================================================
# 4. Remove old scores and transfer only the fixed scores
# ============================================================

remove_old_score_columns <- function(cds) {

  old_score_columns <- grep(
    "_score$",
    colnames(pData(cds)),
    value = TRUE
  )

  if (length(old_score_columns) > 0) {
    for (score_name in old_score_columns) {
      pData(cds)[[score_name]] <- NULL
    }
  }

  cds
}

transfer_fixed_scores <- function(
    cds,
    score_metadata,
    score_columns
) {

  cds_cells <- rownames(pData(cds))
  score_cells <- rownames(score_metadata)

  missing_cells <- setdiff(cds_cells, score_cells)

  if (length(missing_cells) > 0) {
    warning(
      length(missing_cells),
      " CDS cells were not found in fixed_module_scores. ",
      "Missing cells will receive NA."
    )
  }

  missing_score_columns <- setdiff(
    score_columns,
    colnames(score_metadata)
  )

  if (length(missing_score_columns) > 0) {
    stop(
      "Missing score columns: ",
      paste(missing_score_columns, collapse = ", ")
    )
  }

  # Match cells by barcode
  cell_index <- match(cds_cells, score_cells)

  if (all(is.na(cell_index))) {
    stop("No matching cells found.")
  }

  # Transfer scores preserving CDS cell order
  for (score_name in score_columns) {

    pData(cds)[[score_name]] <- unname(
      score_metadata[[score_name]][cell_index]
    )

  }

  message(
    sum(!is.na(cell_index)),
    "/",
    length(cds_cells),
    " cells successfully matched."
  )

  return(cds)
}




prepare_fixed_scores <- function(
    cds,
    groups_keep,
    score_columns = fixed_score_columns,
    reverse_pseudotime = FALSE
) {

  metadata_df <- pData(cds) |>
    as.data.frame()

  required_metadata <- c(
    "groups",
    "Pseudotime",
    score_columns
  )

  missing_metadata <- setdiff(
    required_metadata,
    colnames(metadata_df)
  )

  if (length(missing_metadata) > 0) {
    stop(
      "Missing columns in pData(cds): ",
      paste(missing_metadata, collapse = ", ")
    )
  }

  metadata_df <- metadata_df |>
    rownames_to_column("cell") |>
    filter(
      groups %in% groups_keep,
      !is.na(Pseudotime)
    )

  if (nrow(metadata_df) == 0) {
    stop(
      "No cells remained for groups: ",
      paste(groups_keep, collapse = ", ")
    )
  }

  missing_groups <- setdiff(
    groups_keep,
    unique(as.character(metadata_df$groups))
  )

  if (length(missing_groups) > 0) {
    stop(
      "Groups absent after filtering: ",
      paste(missing_groups, collapse = ", ")
    )
  }

  pseudotime_range <- range(
    metadata_df$Pseudotime,
    na.rm = TRUE
  )

  if (!all(is.finite(pseudotime_range)) || diff(pseudotime_range) == 0) {
    stop("Pseudotime is invalid or has no variation.")
  }

  metadata_df <- metadata_df |>
    mutate(
      Pseudotime_scaled =
        (Pseudotime - pseudotime_range[1]) /
        diff(pseudotime_range)
    )

  if (reverse_pseudotime) {
    metadata_df <- metadata_df |>
      mutate(
        Pseudotime_scaled = 1 - Pseudotime_scaled
      )
  }

  plot_df <- metadata_df |>
    dplyr::select(
      cell,
      groups,
      any_of("State"),
      Pseudotime,
      Pseudotime_scaled,
      all_of(score_columns)
    ) |>
    pivot_longer(
      cols = all_of(score_columns),
      names_to = "Module",
      values_to = "Score"
    ) |>
    mutate(
      Module = sub("_score$", "", Module),
      Module = factor(
        Module,
        levels = module_levels,
        labels = module_labels
      ),
      groups = factor(
        groups,
        levels = groups_keep
      )
    ) |>
    filter(
      !is.na(Module),
      is.finite(Score),
      is.finite(Pseudotime_scaled)
    )

  expected_modules <- factor(
    module_labels,
    levels = module_labels
  )

  modules_present <- unique(plot_df$Module)

  if (!all(expected_modules %in% modules_present)) {
    stop(
      "Not all fixed modules are present after data preparation. Present: ",
      paste(as.character(modules_present), collapse = ", ")
    )
  }

  plot_df
}





# ============================================================
# 6. Plotting function
# ============================================================

plot_fixed_module_scores <- function(
    plot_df,
    plot_title,
    group_colors,
    ncol = 5,
    show_legend = TRUE
) {

  groups_present <- levels(droplevels(plot_df$groups))
  missing_colors <- setdiff(groups_present, names(group_colors))

  if (length(missing_colors) > 0) {
    stop(
      "Missing colors for group(s): ",
      paste(missing_colors, collapse = ", ")
    )
  }

  p <- ggplot(
    plot_df,
    aes(
      x = Pseudotime_scaled,
      y = Score,
      color = groups,
      fill = groups
    )
  ) +
    geom_smooth(
      method = "gam",
      formula = y ~ s(x, bs = "cs"),
      se = TRUE,
      linewidth = 0.9,
      alpha = 0.15,
      na.rm = TRUE
    ) +
    facet_wrap(
      ~Module,
      scales = "free_y",
      ncol = ncol,
      drop = FALSE
    ) +
    scale_color_manual(
      values = group_colors[groups_present],
      breaks = groups_present,
      drop = TRUE
    ) +
    scale_fill_manual(
      values = group_colors[groups_present],
      breaks = groups_present,
      drop = TRUE
    ) +
    theme_classic(base_size = 16) +
    labs(
      title = plot_title,
      x = "Scaled pseudotime",
      y = "Mean expression score",
      color = NULL,
      fill = NULL
    ) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0
      ),
      strip.background = element_blank(),
      strip.text = element_text(
        face = "bold",
        color = "black",
        size = 14
      ),
      axis.text = element_text(
        color = "black"
      ),
      panel.spacing = grid::unit(0.8, "lines"),
      legend.position = if (show_legend) "top" else "none"
    )

  p
}
