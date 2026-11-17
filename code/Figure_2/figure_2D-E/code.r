# Figure/panel: 2D–E
# Purpose: Sample-averaged relative-expression ternary plots.
# Inputs (under data/ unless session objects): combined_data_1.rds.
# Packages: grDevices, httpgd, purrr.
# Outputs: Ternary_Naive_PBS_PD1.tiff; Ternary_Naive_PBS_CTLA4.tiff.
# Prerequisites: functions/relative_expression.r.
# Execution: repository root; see README and docs/SCRIPT_INVENTORY.md.



source("functions/relative_expression.r")

rds_path = file.path("data", "combined_data_1.rds")
combined_data_integrated = readRDS(rds_path)

obj = combined_data_integrated
DefaultAssay(obj) = "RNA"

obj$cell_class = case_when(
  obj$subtype_cell == "Macrophage" ~ "Macrophages",
  obj$subtype_cell == "Neutrophils" ~ "Neutrophils",
  grepl("^T cells", obj$subtype_cell) ~ "T cells",
  TRUE ~ NA_character_
)

expr = GetAssayData(obj,assay = "RNA",layer = "data")
meta = obj@meta.data
meta$cell_id = rownames(meta)

comparisons = list(
  "Naive / PBS / anti-PD-1" = c(
    "naive", "pbs", "antipd1"
  ),
  "Naive / PBS / anti-CTLA-4" = c(
    "naive", "pbs", "antictla4"
  )
)

celltypes = c(
  "Macrophages",
  "Neutrophils",
  "T cells"
)

group_labels = c(
  naive = "Naive",
  pbs = "PBS",
  antipd1 = "anti-PD-1",
  antictla4 = "anti-CTLA-4"
)

group_colors = c(
  "Naive" = "#F15B83",
  "PBS" = "#788A64",
  "anti-PD-1" = "#26361E",
  "anti-CTLA-4" = "#99C4EF"
)

ternary_list = list()
for (comparison_name in names(comparisons)) {
  groups_use = comparisons[[comparison_name]]
  for (celltype in celltypes) {
    message(
      "Processing: ",
      comparison_name,
      " | ",
      celltype
    )
    result = get_ternary_data(
      celltype = celltype,
      groups_use = groups_use,
      comparison_name = comparison_name
    )
    key = paste(
      comparison_name,
      celltype,
      sep = "_"
    )
    ternary_list[[key]] = result
  }
}

ternary_df = bind_rows(ternary_list)
ternary_df %>%
  count(Comparison, cell_type)




top_vertex_genes_final <- ternary_df %>%
  mutate(
    immune_candidate = purrr::map2_lgl(
      cell_type, gene,
      ~ .y %in% immune_gene_dict[[.x]]
    ),
    forced_gene = ifelse(
      cell_type == "T cells" & gene %in% force_keep_tcells,
      TRUE, FALSE
    ),
    bad_gene = grepl(blacklist_pattern, gene)
  ) %>%
  filter(!bad_gene) %>%
  mutate(
    priority_score =
      specificity +
      ifelse(immune_candidate, 1, 0) +
      ifelse(forced_gene, 2, 0)
  ) %>%
  group_by(Comparison, cell_type, dominant) %>%
  arrange(desc(priority_score), desc(specificity), desc(max_prop)) %>%
  slice_head(n = 10) %>%
  ungroup()



triangle = data.frame(
  x = c(0, 1, 0.5, 0),
  y = c(0, 0, sqrt(3) / 2, 0)
)

p_pd1 = plot_ternary(
  comparison_name = "Naive / PBS / anti-PD-1",
  show_labels = TRUE,
  n_labels = 7,
  show_vertex_labels = FALSE,
  use_immune_genes = TRUE
)

p_ctla4 = plot_ternary(
  comparison_name = "Naive / PBS / anti-CTLA-4",
  show_labels = TRUE,
  n_labels = 7,
  show_vertex_labels = FALSE,
  use_immune_genes = TRUE
)


print(p_pd1)
print(p_ctla4)

if (interactive()) httpgd::hgd()


graphics.off()

grDevices::tiff(
  filename = "Ternary_Naive_PBS_PD1.tiff",
  width = 10,
  height = 6,
  units = "in",
  res = 600,
  compression = "lzw",
  type = "cairo",
  bg = "white"
)

print(p_pd1)

dev.off()

graphics.off()

grDevices::tiff(
  filename = "Ternary_Naive_PBS_CTLA4.tiff",
  width = 10,
  height = 6,
  units = "in",
  res = 600,
  compression = "lzw",
  type = "cairo",
  bg = "white"
)

print(p_ctla4)

dev.off()







