


source("./functions/circos.r")

library(dplyr)
library(circlize)

mac_col = "#4E9C81"
neu_col = "#E36A6A"
tar_col = "#9DD9DF"

tiff("./Figure_pbs_naive.tiff", width = 9, height = 9, units = "in", res = 600, compression = "lzw", bg = "white")
draw_circos_plot(plot_df, comparison_name = "PBS vs Naive", main_text = "Figure 5B. PBS vs Naive", top_ligands_per_sender = 20, top_targets = 30, min_weight = 0.01)
dev.off()

tiff("./Figure_antictla_pbs.tiff", width = 9, height = 9, units = "in", res = 600, compression = "lzw", bg = "white")
draw_circos_plot(plot_df, comparison_name = "anti-CTLA-4 vs PBS", main_text = "Figure 5C. anti-CTLA-4 vs PBS", top_ligands_per_sender = 20, top_targets = 30, min_weight = 0.01)
dev.off()
