# Shared helper: REffect wrapper; caller supplies expression table and comparisons.
# Inputs: caller arguments/session objects; no automatic data loading.
# Outputs: functions and documented constants; plots only when called.
# Dependencies/prerequisites: see docs/SCRIPT_INVENTORY.md and SOFTWARE_REQUIREMENTS.md.


library(Seurat)
library(REffect)
library(dplyr)
library(ggplot2)
library(patchwork)

run_relative_effect = function(
    data,
    reference,
    treatment,
    genes,
    colors,
    n_boot = 1000
) {

  df_use <- data %>%
    filter(groups %in% c(reference, treatment)) %>%
    mutate(
      condition = ifelse(groups == reference, 0, 1)
    )

  result <- REffect(
    data = df_use,
    dep_vars = genes,
    indep_var = "condition",
    n_boot = n_boot,
    plot = TRUE,
    plot_colors = colors,
    plot_labels = list(
      y = "Relative effect",
      x = NULL
    ),
    plot_theme = theme_classic(),
    plot_legend_position = "none",
    plot_legend_labels = c(reference, treatment)
  )

  return(result)
}