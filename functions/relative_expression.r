
library(Seurat)
library(Matrix)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(patchwork)

get_ternary_data <- function(
    celltype,
    groups_use,
    comparison_name,
    min_expression = 0.1
) {

  md <- meta %>%
    filter(
      cell_class == celltype,
      groups %in% groups_use
    )

  if (nrow(md) == 0) {
    stop("Nenhuma célula encontrada para: ", celltype)
  }

  sample_ids <- unique(md$cartig)

  sample_expr <- lapply(sample_ids, function(s) {

    cells <- md$cell_id[md$cartig == s]

    Matrix::rowMeans(
      expr[, cells, drop = FALSE]
    )
  })

  sample_mat <- do.call(cbind, sample_expr)
  colnames(sample_mat) <- sample_ids

  sample_groups <- vapply(
    sample_ids,
    function(s) {
      unique(as.character(
        md$groups[md$cartig == s]
      ))[1]
    },
    character(1)
  )

  group_mat <- sapply(
    groups_use,
    function(g) {

      sample_cols <- sample_groups == g

      if (!any(sample_cols)) {
        stop("Grupo ausente: ", g)
      }

      rowMeans(
        sample_mat[, sample_cols, drop = FALSE]
      )
    }
  )

  keep <- apply(
    group_mat,
    1,
    max
  ) >= min_expression

  group_mat <- group_mat[
    keep,
    ,
    drop = FALSE
  ]

  totals <- rowSums(group_mat)

  keep_total <- is.finite(totals) & totals > 0

  group_mat <- group_mat[
    keep_total,
    ,
    drop = FALSE
  ]

  relative <- sweep(
    group_mat,
    1,
    rowSums(group_mat),
    "/"
  )

  result <- data.frame(
    gene = rownames(relative),
    Left = relative[, groups_use[1]],
    Right = relative[, groups_use[2]],
    Top = relative[, groups_use[3]],
    cell_type = celltype,
    Comparison = comparison_name,
    stringsAsFactors = FALSE
  )

  labels_use <- unname(group_labels[groups_use])

  proportions <- as.matrix(
    result[, c("Left", "Right", "Top")]
  )

  result$dominant <- labels_use[
    max.col(proportions, ties.method = "first")
  ]

  result$max_prop <- apply(
    proportions, 1, max
  )

  result$second_prop <- apply(
    proportions,
    1,
    function(x) sort(x, decreasing = TRUE)[2]
  )

  result$specificity <- (
    result$max_prop - result$second_prop
  )

  result$x <- result$Right + 0.5 * result$Top

  result$y <- sqrt(3) / 2 * result$Top

  result
}




### plot
plot_ternary <- function(
    comparison_name,
    show_labels = TRUE,
    n_labels = 8,
    show_vertex_labels = FALSE,
    use_immune_genes = FALSE
) {

  plot_df <- ternary_df %>%
    dplyr::filter(Comparison == comparison_name)

  genes_source <- if (use_immune_genes) {
    top_vertex_genes_final
  } else {
    top_vertex_genes
  }

  label_df <- genes_source %>%
    dplyr::filter(Comparison == comparison_name) %>%
    dplyr::group_by(cell_type, dominant) %>%
    dplyr::slice_head(n = n_labels) %>%
    dplyr::ungroup()

  groups_use <- comparisons[[comparison_name]]

  left_label <- unname(group_labels[groups_use[1]])
  right_label <- unname(group_labels[groups_use[2]])
  top_label <- unname(group_labels[groups_use[3]])

  label_left <- label_df %>%
    dplyr::filter(dominant == left_label)

  label_right <- label_df %>%
    dplyr::filter(dominant == right_label)

  label_top <- label_df %>%
    dplyr::filter(dominant == top_label) %>%
    dplyr::group_by(cell_type) %>%
    dplyr::arrange(gene, .by_group = TRUE) %>%
    dplyr::mutate(
      x_label = seq(
        from = 0.14,
        to = 0.86,
        length.out = dplyr::n()
      ),
      y_label = 1.05
    ) %>%
    dplyr::ungroup()

  p <- ggplot2::ggplot() +
    ggplot2::geom_path(
      data = triangle,
      ggplot2::aes(x, y),
      color = "grey40",
      linewidth = 0.5
    ) +
    ggplot2::geom_point(
      data = plot_df,
      ggplot2::aes(
        x = x,
        y = y,
        color = dominant
      ),
      size = 0.55,
      alpha = 0.6
    ) +
    ggplot2::scale_color_manual(
      values = group_colors,
      name = "Dominant condition"
    ) +
    ggplot2::facet_wrap(
      ~ cell_type,
      nrow = 1
    ) +
    ggplot2::coord_equal(
      clip = "off",
      xlim = c(-0.28, 1.28),
      ylim = c(-0.22, 1.35)
    ) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(
        face = "bold",
        size = 13
      ),
      legend.position = "bottom",
      plot.margin = ggplot2::margin(
        18, 25, 18, 25
      )
    )

  if (show_vertex_labels) {
    p <- p +
      ggplot2::annotate(
        "text",
        x = 0,
        y = -0.08,
        label = left_label,
        size = 3.5
      ) +
      ggplot2::annotate(
        "text",
        x = 1,
        y = -0.08,
        label = right_label,
        size = 3.5
      ) +
      ggplot2::annotate(
        "text",
        x = 0.5,
        y = 0.95,
        label = top_label,
        size = 3.5
      )
  }

  if (show_labels) {

    if (nrow(label_left) > 0) {
      p <- p +
        ggrepel::geom_text_repel(
          data = label_left,
          ggplot2::aes(
            x = x,
            y = y,
            label = gene
          ),
          color = "black",
          size = 3.5,
          direction = "y",
          nudge_x = -0.12,
          nudge_y = -0.02,
          hjust = 1,
          box.padding = 0.25,
          point.padding = 0.15,
          segment.color = "grey55",
          segment.size = 0.35,
          min.segment.length = 0,
          force = 1,
          max.overlaps = Inf,
          seed = 123,
          show.legend = FALSE
        )
    }

    if (nrow(label_right) > 0) {
      p <- p +
        ggrepel::geom_text_repel(
          data = label_right,
          ggplot2::aes(
            x = x,
            y = y,
            label = gene
          ),
          color = "black",
          size = 3.5,
          direction = "y",
          nudge_x = 0.12,
          nudge_y = -0.02,
          hjust = 0,
          box.padding = 0.25,
          point.padding = 0.15,
          segment.color = "grey55",
          segment.size = 0.35,
          min.segment.length = 0,
          force = 1,
          max.overlaps = Inf,
          seed = 123,
          show.legend = FALSE
        )
    }

    if (nrow(label_top) > 0) {
      p <- p +
        ggplot2::geom_segment(
          data = label_top,
          ggplot2::aes(
            x = x,
            y = y,
            xend = x_label,
            yend = y_label - 0.025
          ),
          color = "grey60",
          linewidth = 0.3,
          inherit.aes = FALSE
        ) +
        ggplot2::geom_text(
          data = label_top,
          ggplot2::aes(
            x = x_label,
            y = y_label,
            label = gene
          ),
          color = "black",
          size = 4,
          angle = 45,
          hjust = 0,
          vjust = 0,
          inherit.aes = FALSE
        )
    }
  }

  return(p)
}





library(dplyr)
library(ggplot2)
library(ggrepel)

immune_genes <- list(

  "Macrophages" = c(
    "Cd1d1", "Cd74", "H2-Aa", "H2-Ab1",
    "Adgre1", "Cd68", "Csf1r", "Mrc1",
    "Lgals3", "Il1b", "Tnf", "Ccl2",
    "Ccl3", "Ccl4", "Cxcl10", "Ccr2",
    "Stat1", "Irf1", "Irf7", "Isg15",
    "Gbp2", "Gbp4", "Nos2", "Arg1",
    "Gpr183", "Rnf128", "Trim47",
    "Abi2", "Vat1"
  ),

  "Neutrophils" = c(
    "S100a8", "S100a9", "Ly6g", "Csf3r",
    "Cxcr2", "Cxcl2", "Cxcl3",
    "Il1b", "Tnf", "Lcn2", "Mmp9",
    "Retnlg", "Camp", "Ngp",
    "Tnfaip6", "Csf2rb2", "Tifa",
    "Tgfbr1", "Ppard", "Capg",
    "Hk3", "Dok3", "Gpr183",
    "Atg5", "Pdpk1"
  ),

  "T cells" = c(
    "Il17a", "Il17f", "Rorc", "Rora",
    "Foxp3", "Il2ra", "Ctla4", "Il10",
    "Ccr5", "Ccr6", "Cish", "Nfil3",
    "Il6st", "Tob1", "Orai1",
    "Cd69", "Cd44", "Icos", "Tigit",
    "Tcf7", "Lef1", "Il7r",
    "Stat3", "Stat5a", "Batf",
    "Maf", "Irf4"
  )
)