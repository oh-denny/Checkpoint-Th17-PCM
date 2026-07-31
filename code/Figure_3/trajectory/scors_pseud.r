



### modeule score
library(Seurat)
library(monocle)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

# ============================================================
# 1. Ler os CDS já criados
# ============================================================

cds_condition_PBS_naive <- readRDS(
  "./data/cds_condition_PBS_naive.rds"
)

cds_condition_IGG_pbs <- readRDS(
  "./data/cds_condition_IGG_pbs.rds"
)

cds_condition_PD1_pbs <- readRDS(
  "./data/cds_condition_PD1_pbs.rds"
)

cds_condition_CTLA4_pbs <- readRDS(
  "./data/cds_condition_CTLA4_pbs.rds"
)


# ============================================================
# 2. Definir assinaturas
# ============================================================

module_signatures <- list(

  Th17 = c(
    "Il17a",
    "Il17f",
    "Rorc",
    "Rora",
    "Maf",
    "Ccr6",
    "Il23r",
    "Il1r1",
    "Il21",
    "Il22",
    "Hif1a",
    "Bhlhe40"
  ),

  Treg = c(
    "Foxp3",
    "Il2ra",
    "Ctla4",
    "Tigit",
    "Ikzf2",
    "Il10",
    "Entpd1",
    "Tnfrsf18"
  ),

  Activation = c(
    "Cd69",
    "Icos",
    "Cd28",
    "Zap70",
    "Lat",
    "Tnfrsf4",
    "Tnfrsf9",
    "Nr4a1",
    "Batf"
  ),

  Migration = c(
    "Cxcr6",
    "Ccr2",
    "Ccr6",
    "Cxcr4",
    "Ccl5",
    "Itga4",
    "Itgb1",
    "S1pr1",
    "Rgs1"
  ),

  Less_differentiated = c(
    "Tcf7",
    "Lef1",
    "Sell",
    "Ccr7",
    "Klf2",
    "Klf3",
    "S1pr1",
    "Tob1",
    "Foxp1"
  )
)


# ============================================================
# 3. Calcular os module scores no Seurat
# ============================================================

SC_scores <- t_cells_CD4
DefaultAssay(SC_scores) <- "RNA"

module_signatures_present <- lapply(
  module_signatures,
  function(genes) {
    intersect(genes, rownames(SC_scores))
  }
)

print(module_signatures_present)
print(sapply(module_signatures_present, length))


for (module_name in names(module_signatures_present)) {

  genes_present <- module_signatures_present[[module_name]]

  if (length(genes_present) < 3) {

    warning(
      module_name,
      " possui menos de três genes presentes e não será calculado."
    )

    next
  }

  temporary_prefix <- paste0(module_name, "_tmp")

  SC_scores <- AddModuleScore(
    object = SC_scores,
    features = list(genes_present),
    assay = "RNA",
    name = temporary_prefix
  )

  generated_column <- paste0(temporary_prefix, "1")
  final_column <- paste0(module_name, "_score")

  SC_scores@meta.data[[final_column]] <-
    SC_scores@meta.data[[generated_column]]

  SC_scores@meta.data[[generated_column]] <- NULL
}


score_columns <- grep(
  "_score$",
  colnames(SC_scores@meta.data),
  value = TRUE
)

score_columns


transfer_scores_to_cds <- function(
    cds,
    seurat_object,
    score_columns = c(
      "Th17_score",
      "Treg_score",
      "Activation_score",
      "Migration_score",
      "Less_differentiated_score"
    )
) {

  cds_cells <- rownames(pData(cds))
  seurat_cells <- rownames(seurat_object@meta.data)

  missing_cells <- setdiff(
    cds_cells,
    seurat_cells
  )

  if (length(missing_cells) > 0) {

    stop(
      paste0(
        length(missing_cells),
        " células do CDS não foram encontradas no objeto Seurat.\n",
        "Exemplos: ",
        paste(head(missing_cells), collapse = ", ")
      )
    )
  }

  available_scores <- intersect(
    score_columns,
    colnames(seurat_object@meta.data)
  )

  if (length(available_scores) == 0) {

    stop(
      "Nenhuma coluna de module score foi encontrada no objeto Seurat."
    )
  }

  score_metadata <- seurat_object@meta.data[
    cds_cells,
    available_scores,
    drop = FALSE
  ]

  for (score_name in available_scores) {

    pData(cds)[[score_name]] <-
      score_metadata[[score_name]]
  }

  cds
}

cds_condition_PBS_naive <- transfer_scores_to_cds(
  cds = cds_condition_PBS_naive,
  seurat_object = SC_scores
)

cds_condition_IGG_pbs <- transfer_scores_to_cds(
  cds = cds_condition_IGG_pbs,
  seurat_object = SC_scores
)

cds_condition_PD1_pbs <- transfer_scores_to_cds(
  cds = cds_condition_PD1_pbs,
  seurat_object = SC_scores
)

cds_condition_CTLA4_pbs <- transfer_scores_to_cds(
  cds = cds_condition_CTLA4_pbs,
  seurat_object = SC_scores
)

prepare_pseudotime_scores <- function(cds) {

  score_columns <- c(
    "Th17_score",
    "Treg_score",
    "Activation_score",
    "Migration_score",
    "Less_differentiated_score"
  )

  available_scores <- intersect(
    score_columns,
    colnames(pData(cds))
  )

  if (length(available_scores) == 0) {

    stop(
      paste0(
        "Nenhum module score foi encontrado no pData(cds).\n",
        "Colunas disponíveis: ",
        paste(colnames(pData(cds)), collapse = ", ")
      )
    )
  }

  metadata_df <- pData(cds) |>
    as.data.frame() |>
    rownames_to_column("cell") |>
    filter(
      !is.na(Pseudotime)
    )

  # State pode não existir em algum CDS
  columns_to_keep <- c(
    "cell",
    "groups",
    "Pseudotime",
    available_scores
  )

  if ("State" %in% colnames(metadata_df)) {
    columns_to_keep <- c(
      "cell",
      "groups",
      "State",
      "Pseudotime",
      available_scores
    )
  }

  metadata_df |>
    dplyr::select(
      all_of(columns_to_keep)
    ) |>
    pivot_longer(
      cols = all_of(available_scores),
      names_to = "Module",
      values_to = "Score"
    ) |>
    mutate(
      Module = sub("_score$", "", Module),

      Module = factor(
        Module,
        levels = c(
          "Less_differentiated",
          "Activation",
          "Migration",
          "Treg",
          "Th17"
        ),

        labels = c(
          "Less differentiated",
          "T-cell activation",
          "Migration",
          "Treg",
          "Th17"
        )
      )
    )
}

scores_ctla4 <- prepare_pseudotime_scores(
  cds_condition_CTLA4_pbs
)

head(scores_ctla4)

table(
  scores_ctla4$Module,
  scores_ctla4$groups
)


p_scores_ctla4 <- ggplot(
  scores_ctla4,
  aes(
    x = Pseudotime,
    y = Score,
    color = groups,
    fill = groups
  )
) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, bs = "cs"),
    se = TRUE,
    linewidth = 1
  ) +
  facet_wrap(
    ~Module,
    scales = "free_y",
    ncol = 2
  ) +
  scale_color_manual(
    values = group_colors
  ) +
  scale_fill_manual(
    values = group_colors
  ) +
  theme_classic(base_size = 13) +
  labs(
    title = "Anti-CTLA-4 vs PBS",
    x = "Pseudotime",
    y = "Module score",
    color = NULL,
    fill = NULL
  ) +
  theme(
    strip.background = element_rect(
      fill = "grey95",
      color = "grey40"
    ),
    strip.text = element_text(
      face = "bold"
    ),
    legend.position = "top"
  )

p_scores_ctla4




plot_module_scores_pseudotime <- function(
    cds,
    plot_title,
    group_colors,
    ncol = 2
) {

  plot_df <- prepare_pseudotime_scores(cds)

  ggplot(
    plot_df,
    aes(
      x = Pseudotime,
      y = Score,
      color = groups,
      fill = groups
    )
  ) +
    geom_smooth(
      method = "gam",
      formula = y ~ s(x, bs = "cs"),
      se = TRUE,
      linewidth = 1
    ) +
    facet_wrap(
      ~Module,
      scales = "free_y",
      ncol = ncol
    ) +
    scale_color_manual(
      values = group_colors
    ) +
    scale_fill_manual(
      values = group_colors
    ) +
    theme_classic(base_size = 13) +
    labs(
      title = plot_title,
      x = "Pseudotime",
      y = "Module score",
      color = NULL,
      fill = NULL
    ) +
    theme(
      strip.background = element_rect(
        fill = "grey95",
        color = "grey40"
      ),
      strip.text = element_text(
        face = "bold"
      ),
      legend.position = "top"
    )
}


p_modules_infection <- plot_module_scores_pseudotime(
  cds = cds_condition_PBS_naive,
  plot_title = "PBS vs Naive",
  group_colors = group_colors
)

p_modules_igg <- plot_module_scores_pseudotime(
  cds = cds_condition_IGG_pbs,
  plot_title = "IgG vs PBS",
  group_colors = group_colors
)

p_modules_pd1 <- plot_module_scores_pseudotime(
  cds = cds_condition_PD1_pbs,
  plot_title = "Anti-PD-1 vs PBS",
  group_colors = group_colors
)

p_modules_ctla4 <- plot_module_scores_pseudotime(
  cds = cds_condition_CTLA4_pbs,
  plot_title = "Anti-CTLA-4 vs PBS",
  group_colors = group_colors
)









prepare_pseudotime_scores <- function(
    cds,
    groups_keep,
    reverse_pseudotime = FALSE
) {

  score_columns <- c(
    "Th17_score",
    "Treg_score",
    "Activation_score",
    "Migration_score",
    "Less_differentiated_score"
  )

  available_scores <- intersect(
    score_columns,
    colnames(pData(cds))
  )

  if (length(available_scores) == 0) {
    stop("Nenhum module score foi encontrado no pData(cds).")
  }

  metadata_df <- pData(cds) |>
    as.data.frame() |>
    tibble::rownames_to_column("cell") |>
    dplyr::filter(
      !is.na(Pseudotime),
      groups %in% groups_keep
    )

  # Padronizar o pseudotempo entre 0 e 1
  pt_min <- min(metadata_df$Pseudotime, na.rm = TRUE)
  pt_max <- max(metadata_df$Pseudotime, na.rm = TRUE)

  metadata_df$Pseudotime_scaled <- (
    metadata_df$Pseudotime - pt_min
  ) / (
    pt_max - pt_min
  )

  if (reverse_pseudotime) {
    metadata_df$Pseudotime_scaled <-
      1 - metadata_df$Pseudotime_scaled
  }

  metadata_df |>
    dplyr::select(
      cell,
      groups,
      dplyr::any_of("State"),
      Pseudotime,
      Pseudotime_scaled,
      dplyr::all_of(available_scores)
    ) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(available_scores),
      names_to = "Module",
      values_to = "Score"
    ) |>
    dplyr::mutate(
      Module = sub("_score$", "", Module),

      Module = factor(
        Module,
        levels = c(
          "Less_differentiated",
          "Activation",
          "Migration",
          "Treg",
          "Th17"
        ),
        labels = c(
          "Less differentiated",
          "T-cell activation",
          "Migration",
          "Treg",
          "Th17"
        )
      ),

      groups = droplevels(factor(groups))
    )
}


scores_infection <- prepare_pseudotime_scores(
  cds = cds_condition_PBS_naive,
  groups_keep = c("naive", "pbs")
)

scores_igg <- prepare_pseudotime_scores(
  cds = cds_condition_IGG_pbs,
  groups_keep = c("pbs", "igg")
)

scores_pd1 <- prepare_pseudotime_scores(
  cds = cds_condition_PD1_pbs,
  groups_keep = c("pbs", "antipd1")
)

scores_ctla4 <- prepare_pseudotime_scores(
  cds = cds_condition_CTLA4_pbs,
  groups_keep = c("pbs", "antictla4")
)



plot_module_scores <- function(
    plot_df,
    plot_title,
    group_colors,
    ncol = 2
) {

  groups_present <- unique(as.character(plot_df$groups))

  colors_present <- group_colors[
    groups_present
  ]

  ggplot(
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
      se = FALSE,
      linewidth = 1
    ) +
    facet_wrap(
      ~Module,
      scales = "free_y",
      ncol = 5
    ) +
    scale_color_manual(
      values = colors_present,
      drop = TRUE
    ) +
    scale_fill_manual(
      values = colors_present,
      drop = TRUE
    ) +
    theme_classic(base_size = 13) +
    labs(
      title = plot_title,
      x = "Scaled pseudotime",
      y = "Module score",
      color = NULL,
      fill = NULL
    ) +
    theme(
      strip.background = element_rect(
        fill = "grey95",
        color = "grey40"
      ),
      strip.text = element_text(
        face = "bold"
      ),
      legend.position = "top"
    )
}





p_modules_infection <- plot_module_scores(
  plot_df = scores_infection,
  plot_title = "PBS vs Naive",
  group_colors = group_colors
)

p_modules_igg <- plot_module_scores(
  plot_df = scores_igg,
  plot_title = "IgG vs PBS",
  group_colors = group_colors
)

p_modules_pd1 <- plot_module_scores(
  plot_df = scores_pd1,
  plot_title = "Anti-PD-1 vs PBS",
  group_colors = group_colors
)

p_modules_ctla4 <- plot_module_scores(
  plot_df = scores_ctla4,
  plot_title = "Anti-CTLA-4 vs PBS",
  group_colors = group_colors
)


library(patchwork)

p_all_modules <- p_modules_infection /
  p_modules_igg /
  p_modules_pd1 /
  p_modules_ctla4 +
  plot_layout(
    ncol = 1,
    guides = "collect"
  ) &
  theme(
    legend.position = "top"
  )

p_all_modules



























# ============================================================
# FIXED MODULE SCORES ALONG PSEUDOTIME
# Same five signatures are used in every contrast
# Signatures were assembled only from genes present in gene_lists
# ============================================================

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

# ============================================================
# 1. Read CDS objects already created
# ============================================================

cds_condition_PBS_naive <- readRDS(
  "./data/cds_condition_PBS_naive.rds"
)

cds_condition_IGG_pbs <- readRDS(
  "./data/cds_condition_IGG_pbs.rds"
)

cds_condition_PD1_pbs <- readRDS(
  "./data/cds_condition_PD1_pbs.rds"
)

cds_condition_CTLA4_pbs <- readRDS(
  "./data/cds_condition_CTLA4_pbs.rds"
)

# Required objects expected in the current R environment:
#   t_cells_CD4
#   gene_lists
#   group_colors

required_objects <- c("t_cells_CD4", "gene_lists", "group_colors")
missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), inherits = TRUE)
]

if (length(missing_objects) > 0) {
  stop(
    "Missing required objects: ",
    paste(missing_objects, collapse = ", ")
  )
}

# ============================================================
# 2. Fixed signatures used for ALL four contrasts
# ============================================================

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

module_signatures <- lapply(module_signatures, unique)

# Confirm that every signature gene occurred in at least one contrast.
all_contrast_genes <- unique(unlist(gene_lists, use.names = FALSE))
not_from_contrasts <- lapply(
  module_signatures,
  function(genes) setdiff(genes, all_contrast_genes)
)

if (any(lengths(not_from_contrasts) > 0)) {
  problem_text <- paste(
    vapply(
      names(not_from_contrasts),
      function(module_name) {
        missing_genes <- not_from_contrasts[[module_name]]
        if (length(missing_genes) == 0) return("")
        paste0(module_name, ": ", paste(missing_genes, collapse = ", "))
      },
      character(1)
    ),
    collapse = "\n"
  )

  stop(
    "The following signature genes are absent from gene_lists:\n",
    problem_text
  )
}

# Keep only genes available in the Seurat object.
module_signatures_present <- lapply(
  module_signatures,
  function(genes) intersect(genes, rownames(t_cells_CD4))
)

if (any(lengths(module_signatures_present) == 0)) {
  empty_modules <- names(module_signatures_present)[
    lengths(module_signatures_present) == 0
  ]

  stop(
    "No genes found in t_cells_CD4 for module(s): ",
    paste(empty_modules, collapse = ", ")
  )
}

signature_summary <- tibble(
  Module = names(module_signatures_present),
  N_genes = lengths(module_signatures_present),
  Genes = vapply(
    module_signatures_present,
    paste,
    collapse = ", ",
    FUN.VALUE = character(1)
  )
)

print(signature_summary)

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

fixed_module_scores <- calculate_fixed_module_scores(
  seurat_object = t_cells_CD4,
  modules = module_signatures_present,
  assay = "RNA",
  layer = "data"
)

fixed_score_columns <- c(
  "Less_differentiated_score",
  "T_cell_activation_score",
  "Migration_score",
  "Th17_score",
  "Treg_score"
)

if (!identical(sort(colnames(fixed_module_scores)), sort(fixed_score_columns))) {
  stop(
    "Unexpected score columns. Found: ",
    paste(colnames(fixed_module_scores), collapse = ", ")
  )
}

print(colnames(fixed_module_scores))

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
    stop(
      length(missing_cells),
      " CDS cells were not found in fixed_module_scores. Examples: ",
      paste(head(missing_cells, 5), collapse = ", ")
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

  score_metadata_use <- score_metadata[
    cds_cells,
    score_columns,
    drop = FALSE
  ]

  for (score_name in score_columns) {
    pData(cds)[[score_name]] <- score_metadata_use[[score_name]]
  }

  cds
}

cds_list <- list(
  Infection = cds_condition_PBS_naive,
  IgG = cds_condition_IGG_pbs,
  PD1 = cds_condition_PD1_pbs,
  CTLA4 = cds_condition_CTLA4_pbs
)

cds_list <- lapply(cds_list, remove_old_score_columns)

cds_list <- lapply(
  cds_list,
  transfer_fixed_scores,
  score_metadata = fixed_module_scores,
  score_columns = fixed_score_columns
)

cds_condition_PBS_naive <- cds_list$Infection
cds_condition_IGG_pbs <- cds_list$IgG
cds_condition_PD1_pbs <- cds_list$PD1
cds_condition_CTLA4_pbs <- cds_list$CTLA4

# Strict validation: every CDS must contain exactly the same five scores.
score_check <- lapply(
  cds_list,
  function(cds) {
    grep("_score$", colnames(pData(cds)), value = TRUE)
  }
)

print(score_check)

if (!all(vapply(
  score_check,
  function(x) identical(sort(x), sort(fixed_score_columns)),
  logical(1)
))) {
  stop("The CDS objects do not contain exactly the same fixed score columns.")
}

# ============================================================
# 5. Prepare plotting data
# ============================================================

module_levels <- c(
  "Less_differentiated",
  "T_cell_activation",
  "Migration",
  "Th17",
  "Treg"
)

module_labels <- c(
  "Less differentiated",
  "T-cell activation",
  "Migration",
  "Th17",
  "Treg"
)

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

scores_infection <- prepare_fixed_scores(
  cds = cds_condition_PBS_naive,
  groups_keep = c("naive", "pbs")
)

scores_igg <- prepare_fixed_scores(
  cds = cds_condition_IGG_pbs,
  groups_keep = c("pbs", "igg")
)

scores_pd1 <- prepare_fixed_scores(
  cds = cds_condition_PD1_pbs,
  groups_keep = c("pbs", "antipd1")
)

scores_ctla4 <- prepare_fixed_scores(
  cds = cds_condition_CTLA4_pbs,
  groups_keep = c("pbs", "antictla4")
)

# Verify that all four data frames contain exactly the same five modules.
plot_module_check <- lapply(
  list(
    Infection = scores_infection,
    IgG = scores_igg,
    PD1 = scores_pd1,
    CTLA4 = scores_ctla4
  ),
  function(x) levels(x$Module)
)

print(plot_module_check)

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

p_infection <- plot_fixed_module_scores(
  plot_df = scores_infection,
  plot_title = "PBS vs Naive",
  group_colors = group_colors
)

p_igg <- plot_fixed_module_scores(
  plot_df = scores_igg,
  plot_title = "IgG vs PBS",
  group_colors = group_colors
)

p_pd1 <- plot_fixed_module_scores(
  plot_df = scores_pd1,
  plot_title = "Anti-PD-1 vs PBS",
  group_colors = group_colors
)

p_ctla4 <- plot_fixed_module_scores(
  plot_df = scores_ctla4,
  plot_title = "Anti-CTLA-4 vs PBS",
  group_colors = group_colors
)

# ============================================================
# 7. Combine plots
# ============================================================

p_all_fixed_scores <- (
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

# ============================================================
# 8. Save outputs
# ============================================================

if (!dir.exists("./results")) {
  dir.create("./results", recursive = TRUE)
}


ggsave(
  filename = "./fixed_module_scores_pseudotime.tiff",
  plot = p_all_fixed_scores,
  device = ragg::agg_tiff,
  width = 14,
  height = 13,
  units = "in",
  res = 600,
  compression = "lzw",
  background = "white"
)

# Save CDS objects containing only the five fixed score columns.
saveRDS(
  cds_condition_PBS_naive,
  "./data/cds_condition_PBS_naive_fixed_scores.rds"
)

saveRDS(
  cds_condition_IGG_pbs,
  "./data/cds_condition_IGG_pbs_fixed_scores.rds"
)

saveRDS(
  cds_condition_PD1_pbs,
  "./data/cds_condition_PD1_pbs_fixed_scores.rds"
)

saveRDS(
  cds_condition_CTLA4_pbs,
  "./data/cds_condition_CTLA4_pbs_fixed_scores.rds"
)