# Figure/panel: 5C–D
# Purpose: Ligand-target chord diagrams.
# Inputs (under data/ unless session objects): session plot_df.
# Packages: circlize, dplyr.
# Outputs: ./Figure_pbs_naive.tiff; ./Figure_antictla_pbs.tiff.
# Prerequisites: Run figure_5A; functions/circos.r.
# Execution: repository root; see README and docs/SCRIPT_INVENTORY.md.

if (!exists("plot_df", inherits = TRUE)) stop("Missing plot_df. Run Figure 5A in this R session first; plot_df must contain ligand-target predictions.")



source("./functions/circos.r")

library(dplyr)
library(circlize)

mac_col = "#12544F"
neu_col = "#E36A6A"
tar_col = "#ABD2FA"

tiff("./Figure_pbs_naive.tiff", width = 9, height = 9, units = "in", res = 600, compression = "lzw", bg = "white")
draw_circos_plot(plot_df, comparison_name = "PBS vs Naive", main_text = "Figure 5B. PBS vs Naive", top_ligands_per_sender = 20, top_targets = 30, min_weight = 0.01)
dev.off()

tiff("./Figure_antictla_pbs.tiff", width = 9, height = 9, units = "in", res = 600, compression = "lzw", bg = "white")
draw_circos_plot(plot_df, comparison_name = "anti-CTLA-4 vs PBS", main_text = "Figure 5C. anti-CTLA-4 vs PBS", top_ligands_per_sender = 20, top_targets = 30, min_weight = 0.01)
dev.off()
