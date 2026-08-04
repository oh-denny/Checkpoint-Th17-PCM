

### group levels
group_levels = c(
  "naive",
  "pbs",
  "igg",
  "antipd1",
  "antictla4"
)

cfu_group_levels = c(
  "pbs",
  "igg",
  "antipd1",
  "antictla4"
)

cluster_colors = c(
  "0"  = "#F06C45",
  "1"  = "#296DB0",
  "2"  = "#E93E3F",
  "3"  = "#88D3D3",
  "4"  = "#2C8C61",
  "5"  = "#155C42",
  "6"  = "#41704D",
  "7"  = "#0A6E75",
  "8"  = "#84C868",
  "9"  = "#B15929",
  "10" = "#A99099"
)

pseudocount = 1e-5


get_cfu_benefit <- function(
    contrasts_df,
    treatment,
    reference = "pbs"
) {

  forward_contrast <- paste(
    treatment,
    "-",
    reference
  )

  reverse_contrast <- paste(
    reference,
    "-",
    treatment
  )

  result <- contrasts_df |>
    filter(
      contrast %in% c(
        forward_contrast,
        reverse_contrast
      )
    )

  if (nrow(result) != 1) {
    stop(
      paste0(
        "Contraste de CFU não encontrado ou duplicado: ",
        treatment,
        " versus ",
        reference
      )
    )
  }

  if (result$contrast == forward_contrast) {


    benefit <- -result$estimate

  } else {

    benefit <- result$estimate
  }

  as.numeric(benefit)
}

create_concordance_ranking <- function(
    cluster_wide,
    treatment,
    comparison_groups,
    cfu_benefit
) {

  required_columns <- c(
    treatment,
    comparison_groups
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(cluster_wide)
  )

  if (length(missing_columns) > 0) {
    stop(
      paste(
        "Colunas ausentes:",
        paste(missing_columns, collapse = ", ")
      )
    )
  }

  treatment_log2fc <- cluster_wide[[treatment]]

  comparison_matrix <- cluster_wide |>
    select(
      all_of(comparison_groups)
    ) |>
    as.matrix()

  comparison_mean <- rowMeans(
    comparison_matrix,
    na.rm = TRUE
  )

  specificity <- treatment_log2fc -
    comparison_mean

  cluster_wide |>
    mutate(
      cluster = as.character(
        RNA_snn_res.0.2
      ),

      treatment = treatment,

      treatment_log2fc =
        treatment_log2fc,

      comparison_mean =
        comparison_mean,

      specificity =
        specificity,

      cfu_benefit =
        cfu_benefit,

      abundance_component =
        pmax(treatment_log2fc, 0),

      specificity_component =
        pmax(specificity, 0),

      concordance_score =
        abundance_component *
        specificity_component *
        cfu_benefit
    ) |>
    arrange(
      desc(concordance_score)
    ) |>
    mutate(
      rank = row_number()
    )
}



plot_concordance_ranking <- function(
    ranking_data,
    treatment_label,
    cluster_colors,
    show_legend = TRUE
) {

  max_score <- max(
    ranking_data$concordance_score,
    na.rm = TRUE
  )

  # Evita problema caso todos os scores sejam zero
  if (!is.finite(max_score) || max_score <= 0) {
    max_score <- 1
  }

  ggplot(
    ranking_data,
    aes(
      x = specificity,
      y = treatment_log2fc,
      color = cluster
    )
  ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.4,
      color = "gray50"
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.4,
      color = "gray50"
    ) +
    geom_point(
      aes(
        size = concordance_score
      ),
      alpha = 0.9
    ) +
    scale_color_manual(
      values = cluster_colors,
      drop = FALSE,
      name = "Cluster"
    ) +
    scale_size_continuous(
      range = c(2.5, 10),
      limits = c(0, max_score),
      name = "Concordance\nscore"
    ) +
    labs(
      x = bquote(
        "Specificity for " *
          .(treatment_label) *
          " (" *
          Delta * log[2] * "FC)"
      ),
      y = bquote(
        log[2] * "FC under " *
          .(treatment_label)
      )
    ) +
    theme_classic(
      base_size = 18
    ) +
    theme(
      legend.position = if (
        show_legend
      ) {
        "bottom"
      } else {
        "none"
      },

      axis.title = element_text(
        color = "black"
      ),

      axis.text = element_text(
        color = "black"
      )
    )
}


create_prioritization_ranking <- function(
    cluster_wide,
    treatment,
    comparison_groups
) {

  required_columns <- c(
    "RNA_snn_res.0.2",
    treatment,
    comparison_groups
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(cluster_wide)
  )

  if (length(missing_columns) > 0) {
    stop(
      paste0(
        "Missing columns: ",
        paste(
          missing_columns,
          collapse = ", "
        )
      )
    )
  }

  treatment_log2fc <- cluster_wide[[treatment]]

  comparison_mean <- rowMeans(
    as.matrix(
      cluster_wide[
        comparison_groups
      ]
    ),
    na.rm = TRUE
  )

  specificity <-
    treatment_log2fc -
    comparison_mean

  cluster_wide |>
    dplyr::mutate(
      cluster = as.character(
        RNA_snn_res.0.2
      ),

      treatment =
        treatment,

      treatment_log2fc =
        treatment_log2fc,

      comparison_mean =
        comparison_mean,

      specificity =
        specificity,

      abundance_component =
        pmax(
          treatment_log2fc,
          0
        ),

      specificity_component =
        pmax(
          specificity,
          0
        ),

      prioritization_score =
        abundance_component *
        specificity_component
    ) |>
    dplyr::arrange(
      dplyr::desc(
        prioritization_score
      )
    ) |>
    dplyr::mutate(
      rank = dplyr::row_number()
    )
}


plot_prioritization_ranking <- function(
    ranking_data,
    treatment_label,
    cluster_colors,
    show_legend = TRUE
) {

  ggplot2::ggplot(
    ranking_data,
    ggplot2::aes(
      x = specificity,
      y = treatment_log2fc,
      color = cluster
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.4,
      color = "gray50"
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.4,
      color = "gray50"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        size = prioritization_score
      ),
      alpha = 0.9
    ) +
    ggplot2::scale_color_manual(
      values = cluster_colors,
      drop = FALSE,
      name = "Cluster"
    ) +
    ggplot2::scale_size_continuous(
      range = c(
        2.5,
        9
      ),
      name = "Prioritization\nscore"
    ) +
    ggplot2::labs(
      x = bquote(
        "Specificity for " *
          .(treatment_label) *
          " (" *
          Delta *
          log[2] *
          "FC)"
      ),

      y = bquote(
        log[2] *
          "FC under " *
          .(treatment_label)
      )
    ) +
    ggplot2::theme_classic(
      base_size = 18
    ) +
    ggplot2::theme(
      legend.position = if (
        show_legend
      ) {
        "bottom"
      } else {
        "none"
      },

      axis.text = ggplot2::element_text(
        color = "black"
      ),

      axis.title = ggplot2::element_text(
        color = "black"
      )
    )
}



### trajectory
select_de_genes = function(
    de_table,
    n_up = 50,
    n_down = 50,
    padj_cutoff = 0.05,
    logfc_cutoff = 0.25
) {

  de_df = de_table |>
    as.data.frame() |>
    tibble::rownames_to_column("gene")

  fc_column = if ("avg_log2FC" %in% colnames(de_df)) {
    "avg_log2FC"
  } else if ("avg_logFC" %in% colnames(de_df)) {
    "avg_logFC"
  } else {
    stop("Não foi encontrada a coluna avg_log2FC ou avg_logFC.")
  }

  de_filtered = de_df |>
    filter(
      !is.na(p_val),
      p_val < padj_cutoff,
      abs(.data[[fc_column]]) > logfc_cutoff
    ) |>
    filter(
      !grepl(
        "^mt-|^Rpl|^Rps|^Hba|^Hbb|^Igh|^Igk|^Igl",
        gene,
        ignore.case = TRUE
      )
    )

  genes_up = de_filtered |>
    filter(.data[[fc_column]] > 0) |>
    arrange(desc(.data[[fc_column]])) |>
    slice_head(n = n_up) |>
    pull(gene)

  genes_down = de_filtered |>
    filter(.data[[fc_column]] < 0) |>
    arrange(.data[[fc_column]]) |>
    slice_head(n = n_down) |>
    pull(gene)

  unique(c(genes_up, genes_down))
}