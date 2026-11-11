

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

top_vertex_genes = ternary_df %>%
  group_by(
    Comparison,
    cell_type,
    dominant
  ) %>%
  slice_max(
    order_by = specificity,
    n = 8,
    with_ties = FALSE
  ) %>%
  ungroup()

# Aux table
top_vertex_table = top_vertex_genes %>%
  select(
    Comparison,
    cell_type,
    dominant,
    gene,
    Left,
    Right,
    Top,
    max_prop,
    specificity
  ) %>%
  arrange(
    Comparison,
    cell_type,
    dominant,
    desc(specificity)
  )

print(top_vertex_table, n = Inf)

exclude_pattern = paste(
  c(
    "^Gm[0-9]",
    "^Hba-",
    "^Hbb-",
    "^Rpl[0-9]",
    "^Rps[0-9]",
    "^Mrpl",
    "^Mrps",
    "^mt-"
  ),
  collapse = "|"
)
immune_table = bind_rows(
  lapply(names(immune_genes), function(ct) {

    tibble(
      cell_type = ct,
      gene = immune_genes[[ct]]
    )
  })
)

ternary_immune = ternary_df %>%
  filter(
    !grepl(exclude_pattern, gene)
  ) %>%
  inner_join(
    immune_table,
    by = c("cell_type", "gene")
  ) %>%
  mutate(
    max_prop = pmax(Left, Right, Top),
    second_prop = apply(
      cbind(Left, Right, Top),
      1,
      function(z) sort(z, decreasing = TRUE)[2]
    ),
    specificity = max_prop - second_prop
  )


top_vertex_genes_immune = ternary_immune %>%
  group_by(
    Comparison,
    cell_type,
    dominant
  ) %>%
  slice_max(
    order_by = specificity,
    n = 4,
    with_ties = FALSE
  ) %>%
  ungroup()

print(
  top_vertex_genes_immune %>%
    select(
      Comparison, cell_type, dominant,
      gene, max_prop, specificity
    ),
  n = Inf
)


triangle = data.frame(
  x = c(0, 1, 0.5, 0),
  y = c(0, 0, sqrt(3) / 2, 0)
)

p_pd1 = plot_ternary(
  comparison_name = "Naive / PBS / anti-PD-1",
  show_labels = TRUE,
  n_labels = 5,
  show_vertex_labels = FALSE,
  use_immune_genes = FALSE
)

p_ctla4 = plot_ternary(
  comparison_name = "Naive / PBS / anti-CTLA-4",
  show_labels = TRUE,
  n_labels = 5,
  show_vertex_labels = FALSE,
  use_immune_genes = FALSE
)


print(p_pd1)
print(p_ctla4)

httpgd::hgd()


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







