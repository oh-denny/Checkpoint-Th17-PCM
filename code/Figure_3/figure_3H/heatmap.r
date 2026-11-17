# Figure/panel: 3H
# Purpose: GO enrichment and functional-program classification.
# Inputs (under data/ unless session objects): four cds_condition_*.rds.
# Packages: clusterProfiler, dplyr, enrichplot, ggplot2, monocle, org.Mm.eg.db, ragg, stringr, tibble, tidyr.
# Outputs: GO_programs_gene_supported_heatmap.tiff.
# Prerequisites: Historical trajectories; functions/classify_enrich.r.
# Execution: repository root; see README and docs/SCRIPT_INVENTORY.md.



library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(tibble)
library(enrichplot)
library(monocle)

source("./functions/classify_enrich.r")
trajectory_data_dir = getOption("paracocco.trajectory_data_dir", "data")

cds_condition_pbs = readRDS(file.path(trajectory_data_dir, "cds_condition_PBS_naive.rds"))
cds_condition_igg = readRDS(file.path(trajectory_data_dir, "cds_condition_IGG_pbs.rds"))
cds_condition_pd1 = readRDS(file.path(trajectory_data_dir, "cds_condition_PD1_pbs.rds"))
cds_condition_ctl = readRDS(file.path(trajectory_data_dir, "cds_condition_CTLA4_pbs.rds"))


genes_infection = get_ordering_genes(cds_condition_pbs)
genes_igg       = get_ordering_genes(cds_condition_igg)
genes_pd1       = get_ordering_genes(cds_condition_pd1)
genes_ctla4     = get_ordering_genes(cds_condition_ctl)

gene_lists = list(
  Infection   = unique(genes_infection),
  IgG         = unique(genes_igg),
  `anti-PD-1` = unique(genes_pd1),
  `anti-CTLA-4` = unique(genes_ctla4)
)

lengths(gene_lists)

all_symbols = unique(unlist(gene_lists))

bp_go = compareCluster(
  geneCluster = gene_lists,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  minGSSize = 10,
  maxGSSize = 500,
  readable = TRUE
)

go_table = as.data.frame(bp_go)
th17_signature = c("Il17a","Maf","Hif1a","Bhlhe40")

go_table = go_table |>
  dplyr::mutate(
    Cluster = as.character(.data$Cluster),
    Description = as.character(
      .data$Description
    ),
    geneID = as.character(
      .data$geneID
    ),
    p.adjust = as.numeric(
      .data$p.adjust
    ),
    description_lower = stringr::str_to_lower(
      .data$Description
    ),
    Contrast = dplyr::recode(
      .data$Cluster,
      "Infection" = "PBS vs Naive",
      "IgG" = "IgG vs PBS",
      "anti-PD-1" = "Anti-PD-1 vs PBS",
      "anti-CTLA-4" = "Anti-CTLA-4 vs PBS",
      .default = .data$Cluster
    )
  )

go_table = go_table |>
  dplyr::mutate(
    N_Less_genes = vapply(
      .data$geneID,
      count_signature_genes,
      integer(1),
      signature = signatures$Less_differentiated
    ),
    Less_genes = vapply(
      .data$geneID,
      get_signature_genes,
      character(1),
      signature = signatures$Less_differentiated
    ),
    N_Activation_genes = vapply(
      .data$geneID,
      count_signature_genes,
      integer(1),
      signature = signatures$T_cell_activation
    ),
    Activation_genes = vapply(
      .data$geneID,
      get_signature_genes,
      character(1),
      signature = signatures$T_cell_activation
    ),
    N_Migration_genes = vapply(
      .data$geneID,
      count_signature_genes,
      integer(1),
      signature = signatures$Migration_and_adhesion
    ),
    Migration_genes = vapply(
      .data$geneID,
      get_signature_genes,
      character(1),
      signature = signatures$Migration_and_adhesion
    ),
    N_Th17_genes = vapply(
      .data$geneID,
      count_signature_genes,
      integer(1),
      signature = signatures$Th17_associated
    ),
    Th17_genes = vapply(
      .data$geneID,
      get_signature_genes,
      character(1),
      signature = signatures$Th17_associated
    ),
    N_Treg_genes = vapply(
      .data$geneID,
      count_signature_genes,
      integer(1),
      signature = signatures$Treg_associated
    ),
    Treg_genes = vapply(
      .data$geneID,
      get_signature_genes,
      character(1),
      signature = signatures$Treg_associated
    ),
    has_Il17a = vapply(
      .data$geneID,
      contains_gene,
      logical(1),
      gene = "Il17a"
    )
  )

go_classified = go_table |>
  dplyr::mutate(
    Program = dplyr::case_when(
      (
        stringr::str_detect(
          .data$description_lower,
          paste(
            c(
              "th17",
              "interleukin-17",
              "il-17",
              "inflammatory response",
              "helper t cell differentiation",
              "leukocyte differentiation"
            ),
            collapse = "|"
          )
        ) &
          .data$has_Il17a
      ) ~ "Th17-associated",
      (
        stringr::str_detect(
          .data$description_lower,
          paste(
            c(
              "regulatory t cell",
              "regulatory t-cell",
              "immune tolerance",
              "negative regulation of immune response",
              "negative regulation of t cell activation",
              "regulation of immune response"
            ),
            collapse = "|"
          )
        ) &
          .data$N_Treg_genes >= 1
      ) ~ "Treg-associated",
      (
        stringr::str_detect(
          .data$description_lower,
          paste(
            c(
              "t cell activation",
              "t-cell activation",
              "lymphocyte activation",
              "alpha-beta t cell activation",
              "antigen receptor-mediated signaling",
              "t cell receptor signaling",
              "t-cell receptor signaling",
              "costimulation",
              "co-stimulation"
            ),
            collapse = "|"
          )
        ) &
          .data$N_Activation_genes >= 1
      ) ~ "T-cell activation",
      (
        stringr::str_detect(
          .data$description_lower,
          paste(
            c(
              "t cell migration",
              "t-cell migration",
              "lymphocyte migration",
              "leukocyte migration",
              "cell migration",
              "chemotaxis",
              "chemokine-mediated signaling",
              "cell adhesion",
              "cell-cell adhesion",
              "integrin-mediated",
              "lamellipodium",
              "cell motility",
              "actin filament",
              "actin cytoskeleton"
            ),
            collapse = "|"
          )
        ) &
          .data$N_Migration_genes >= 1
      ) ~ "Migration and adhesion",
      (
        stringr::str_detect(
          .data$description_lower,
          paste(
            c(
              "^t cell differentiation$",
              "^lymphocyte differentiation$",
              "^regulation of t cell differentiation$",
              "^regulation of lymphocyte differentiation$"
            ),
            collapse = "|"
          )
        ) &
          .data$N_Less_genes >= 1
      ) ~ "Less differentiated",

      TRUE ~ NA_character_
    ),
    Signature_genes = dplyr::case_when(
      .data$Program == "Less differentiated" ~
        .data$Less_genes,
      .data$Program == "T-cell activation" ~
        .data$Activation_genes,
      .data$Program == "Migration and adhesion" ~
        .data$Migration_genes,
      .data$Program == "Th17-associated" ~
        .data$Th17_genes,
      .data$Program == "Treg-associated" ~
        .data$Treg_genes,
      TRUE ~ NA_character_
    ),
    N_signature_genes = dplyr::case_when(
      .data$Program == "Less differentiated" ~
        .data$N_Less_genes,
      .data$Program == "T-cell activation" ~
        .data$N_Activation_genes,
      .data$Program == "Migration and adhesion" ~
        .data$N_Migration_genes,
      .data$Program == "Th17-associated" ~
        .data$N_Th17_genes,
      .data$Program == "Treg-associated" ~
        .data$N_Treg_genes,
      TRUE ~ NA_integer_
    )
  )


go_classified_filtered = go_classified |>
  dplyr::filter(
    !is.na(.data$Program),
    !is.na(.data$p.adjust),
    .data$p.adjust < 0.05,
    .data$N_signature_genes > 0
  )

classified_check = go_classified_filtered |>
  dplyr::select(
    dplyr::all_of(
      c(
        "Contrast",
        "Program",
        "Description",
        "Signature_genes",
        "N_signature_genes",
        "p.adjust",
        "geneID"
      )
    )
  ) |>
  dplyr::arrange(
    .data$Program,
    .data$Contrast,
    .data$p.adjust
  )

program_summary = go_classified_filtered |>
  dplyr::group_by(
    .data$Contrast,
    .data$Program
  ) |>
  dplyr::arrange(
    .data$p.adjust,
    dplyr::desc(.data$N_signature_genes),
    .by_group = TRUE
  ) |>
  dplyr::slice_head(n = 1) |>
  dplyr::ungroup() |>
  dplyr::transmute(
    Contrast = .data$Contrast,
    Program = .data$Program,
    Best_term = .data$Description,
    Signature_genes = .data$Signature_genes,
    N_signature_genes = .data$N_signature_genes,
    best_padj = .data$p.adjust,
    minus_log10_padj = -log10(
      pmax(
        .data$p.adjust,
        .Machine$double.xmin
      )
    )
  )

program_levels = c("Less differentiated","T-cell activation","Migration and adhesion","Th17-associated","Treg-associated")
contrast_levels = c("PBS vs Naive","IgG vs PBS","Anti-PD-1 vs PBS","Anti-CTLA-4 vs PBS")

program_summary_complete = tidyr::expand_grid(
  Program = program_levels,
  Contrast = contrast_levels
) |>
  dplyr::left_join(
    program_summary,
    by = c(
      "Program",
      "Contrast"
    )
  ) |>
  dplyr::mutate(
    Best_term = tidyr::replace_na(
      .data$Best_term,
      "Not enriched"
    ),
    Signature_genes = tidyr::replace_na(
      .data$Signature_genes,
      ""
    ),
    N_signature_genes = tidyr::replace_na(
      .data$N_signature_genes,
      0L
    ),
    best_padj = tidyr::replace_na(
      .data$best_padj,
      1
    ),
    minus_log10_padj = tidyr::replace_na(
      .data$minus_log10_padj,
      0
    ),
    Program = factor(
      .data$Program,
      levels = rev(program_levels)
    ),
    Contrast = factor(
      .data$Contrast,
      levels = contrast_levels
    )
  )


program_summary_complete = program_summary_complete |> subset( Program != 'Treg-associated')

p_program_heatmap = ggplot2::ggplot(
  program_summary_complete,
  ggplot2::aes(
    x = .data$Contrast,
    y = .data$Program,
    fill = .data$minus_log10_padj
  )
) +
  ggplot2::geom_tile(
    color = "white",
    linewidth = 1
  ) +
  ggplot2::scale_fill_viridis_c(
    name = expression(
      -log[10]("adjusted p-value")
    ),
    option = "magma",
    limits = c(
      0,
      max(
        program_summary_complete$minus_log10_padj,
        na.rm = TRUE
      )
    )
  ) +
  ggplot2::labs(
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_test(
    base_size = 18
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 40,
      hjust = 1,
      vjust = 1,
      color = "black"
    ),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    legend.position = "right",
    plot.margin = ggplot2::margin(
      8,
      8,
      8,
      8
    )
  )

p_program_heatmap

ggplot2::ggsave(
  filename = "GO_programs_gene_supported_heatmap.tiff",
  plot = p_program_heatmap,
  device = ragg::agg_tiff,
  width = 6.5,
  height = 5.2,
  units = "in",
  res = 600,
  compression = "lzw",
  background = "white"
)


