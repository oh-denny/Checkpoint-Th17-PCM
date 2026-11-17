# Figure/panel: 5A
# Purpose: Pseudobulk edgeR and NicheNet ligand-target inference.
# Inputs (under data/ unless session objects): combined_data_1.rds; t_cells_CD4.rds; nichenet/*.rds; session pb_th17; live OmniPath.
# Packages: dplyr, grid, ragg, readr.
# Outputs: ./data/OmniNetworks_NNformat/lr_Network_Omnipath.rds; Figure5a_p_nichenet_final.tiff.
# Prerequisites: BLOCKED: pb_th17 construction absent; functions/nichenet.r.
# Execution: repository root; see README and docs/SCRIPT_INVENTORY.md.

if (!exists("pb_th17", inherits = TRUE)) stop("Missing pb_th17. Restore the author-approved pb_th17 pseudobulk list (counts and metadata); its construction is absent. See data/README.md.")


source("./functions/nichenet.r")
set.seed(123)

rds_path = file.path("data", "combined_data_1.rds")
combined_data_integrated = readRDS(rds_path)

t_cells_CD4 = readRDS("./data/t_cells_CD4.rds")

t_cells_TH17 = subset(
  t_cells_CD4,
  subset = seurat_clusters %in% c("1")
)

table(combined_data_integrated$subtype_cell)

mac_neu_seuObj = subset(
  combined_data_integrated,
  subset = subtype_cell %in% c("Macrophage", "Neutrophils")
)

t_cells_TH17$cell_nichenet = "Th17/Treg"

mac_neu_seuObj$cell_nichenet = as.character(
  mac_neu_seuObj$subtype_cell
)

DefaultAssay(t_cells_TH17) = "RNA"
DefaultAssay(mac_neu_seuObj) = "RNA"

obj_nichenet = merge(
  x = mac_neu_seuObj,
  y = t_cells_TH17,
  merge.data = FALSE,
  merge.dr = FALSE,
  project = "NicheNet_Paracoccidioides"
)

DefaultAssay(obj_nichenet) = "RNA"

if (inherits(obj_nichenet[["RNA"]], "Assay5")) {

  obj_nichenet = JoinLayers(
    obj_nichenet,
    assay = "RNA"
  )
}

obj_nichenet = NormalizeData(
  obj_nichenet,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

Idents(obj_nichenet) = "cell_nichenet"

obj_nichenet

lr_network_raw = readRDS("./data/nichenet/weighted_networks_nsga2r_final_mouse.rds")
lr_Network_Nichenet = readRDS("./data/nichenet/lr_network_mouse_21122021.rds")
ligand_target_matrix = readRDS("./data/nichenet/ligand_target_matrix_nsga2r_final_mouse.rds")

lr_Network_Nichenet_Unique = lr_Network_Nichenet %>%
    dplyr::distinct(from, to)
nrow(lr_Network_Nichenet_Unique)

options(timeout = 1800)

url_omnipath = paste0(
  "https://omnipathdb.org/interactions?",
  "datasets=ligrecextra",
  "&organisms=10090",
  "&genesymbols=1",
  "&fields=sources,references"
)

lr_Network_Omnipath = readr::read_tsv(
  url_omnipath,
  show_col_types = FALSE
)

head(lr_Network_Omnipath)
dim(lr_Network_Omnipath)
colnames(lr_Network_Omnipath)

lr_Network_Omnipath_Unique = lr_Network_Omnipath %>%
  transmute(
    from = source_genesymbol,
    to = target_genesymbol
  ) %>%
  filter(
    !is.na(from),
    !is.na(to),
    from != "",
    to != ""
  ) %>%
  distinct()

head(lr_Network_Omnipath_Unique)
nrow(lr_Network_Omnipath_Unique)


lr_Network_Omnipath_Unique = lr_Network_Omnipath %>%
    dplyr::distinct(source_genesymbol,target_genesymbol) %>%
    dplyr::rename(from=source_genesymbol, to=target_genesymbol) %>%
    dplyr::filter(from != to)
nrow(lr_Network_Omnipath_Unique)

Venn_plot = draw.pairwise.venn(nrow(lr_Network_Nichenet_Unique), 
    nrow(lr_Network_Omnipath_Unique), 
    nrow(dplyr::intersect(lr_Network_Nichenet_Unique,
    lr_Network_Omnipath_Unique)), 
    category = c("NicheNet: Ligand-Receptor", "Omnipath: LigRecExtra"), 
    lty = rep("blank", 2), fill = c("light blue", "red"), alpha = rep(0.4, 2), 
    cat.pos = c(0, 0), cat.dist = rep(0.025, 2), cex= 1.5, cat.cex=1.5, 
    verbose = FALSE)
grid.draw(Venn_plot)

options(timeout = 1800)

url_omnipath = paste0(
  "https://omnipathdb.org/interactions?",
  "datasets=ligrecextra",
  "&organisms=10090",
  "&genesymbols=1",
  "&fields=sources,references"
)

lr_Network_Omnipath_raw = readr::read_tsv(
  url_omnipath,
  show_col_types = FALSE
)

head(lr_Network_Omnipath_raw)
dim(lr_Network_Omnipath_raw)

lr_Network_Omnipath = lr_Network_Omnipath_raw %>%
  transmute(
    from = source_genesymbol,
    to = target_genesymbol
  ) %>%
  filter(
    !is.na(from),
    !is.na(to),
    from != "",
    to != ""
  ) %>%
  distinct()

head(lr_Network_Omnipath)
nrow(lr_Network_Omnipath)

dir.create(
  "./data/OmniNetworks_NNformat",
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  lr_Network_Omnipath,
  "./data/OmniNetworks_NNformat/lr_Network_Omnipath.rds"
)

obj = obj_nichenet
DefaultAssay(obj) = "RNA"
obj$cell_nichenet = as.character(obj$cell_nichenet)
obj$cell_nichenet[obj$cell_nichenet == "Macrophage"] <- "Macrophages"
Idents(obj) = "cell_nichenet"
expr = GetAssayData(obj,assay = "RNA",layer = "data")
meta = obj[[]]
expression_pct = 0.10

ligands = unique(lr_Network_Nichenet$from)
receptors = unique(lr_Network_Nichenet$to)

## ligands expression in Macrophages and Neutrophils

MAC_CTLA4_lr = get_expressed_lr(
  "Macrophages",
  "antictla4"
)

MAC_PD1_lr = get_expressed_lr(
  "Macrophages",
  "antipd1"
)

MAC_PBS_lr = get_expressed_lr(
  "Macrophages",
  "pbs"
)

NEU_CTLA4_lr = get_expressed_lr(
  "Neutrophils",
  "antictla4"
)

NEU_PD1_lr = get_expressed_lr(
  "Neutrophils",
  "antipd1"
)

NEU_PBS_lr = get_expressed_lr(
  "Neutrophils",
  "pbs"
)

MAC_CTLA4_ligands = unique(MAC_CTLA4_lr$from)
MAC_PD1_ligands = unique(MAC_PD1_lr$from)
MAC_PBS_ligands = unique(MAC_PBS_lr$from)
NEU_CTLA4_ligands = unique(NEU_CTLA4_lr$from)
NEU_PD1_ligands = unique(NEU_PD1_lr$from)
NEU_PBS_ligands = unique(NEU_PBS_lr$from)
sender_ligands = list(
  MAC_CTLA4 = MAC_CTLA4_ligands,
  MAC_PD1 = MAC_PD1_ligands,
  MAC_PBS = MAC_PBS_ligands,
  NEU_CTLA4 = NEU_CTLA4_ligands,
  NEU_PD1 = NEU_PD1_ligands,
  NEU_PBS = NEU_PBS_ligands
)

lengths(sender_ligands)

## Target genes in Th17/Treg cells - pseudobulk

DDS_TH17_CTLA4 = run_pseudobulk_edger(pb_th17,condition = "antictla4",reference = "pbs")
DDS_TH17_PD1 = run_pseudobulk_edger(pb_th17,condition = "antipd1",reference = "pbs")
DDS_TH17_PBS = run_pseudobulk_edger(pb_th17,condition = "pbs",reference = "naive")

geneset_TH17_CTLA4 = get_target_genes(DDS_TH17_CTLA4)
geneset_TH17_PD1 = get_target_genes(DDS_TH17_PD1)
geneset_TH17_PBS = get_target_genes(DDS_TH17_PBS)

receiver_genesets = list(
  CTLA4 = geneset_TH17_CTLA4,
  PD1 = geneset_TH17_PD1,
  PBS = geneset_TH17_PBS
)

lengths(receiver_genesets)
background_TH17_CTLA4 = DDS_TH17_CTLA4$Gene
background_TH17_CTLA4 = intersect(background_TH17_CTLA4,rownames(ligand_target_matrix))
geneset_TH17_CTLA4 = intersect(geneset_TH17_CTLA4,background_TH17_CTLA4)

ligand_activities_MAC_CTLA4 = predict_ligand_activities(
  geneset = geneset_TH17_CTLA4,
  background_expressed_genes = background_TH17_CTLA4,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = MAC_CTLA4_ligands
) %>%
  arrange(desc(aupr_corrected)) %>%
  mutate(rank = row_number())

head(ligand_activities_MAC_CTLA4, 20)

best_ligands_MAC_CTLA4 = ligand_activities_MAC_CTLA4 %>%
  filter(aupr_corrected > 0) %>%
  slice_head(n = 20) %>%
  pull(test_ligand)

ligand_target_MAC_CTLA4 = map_dfr(
  best_ligands_MAC_CTLA4,
  function(ligand) {
    get_weighted_ligand_target_links(
      ligand = ligand,
      geneset = geneset_TH17_CTLA4,
      ligand_target_matrix = ligand_target_matrix,
      n = 250
    )
  }

) %>%
  drop_na() %>%
  filter(weight > 0)

ggplot(ligand_target_MAC_CTLA4,aes(x = target,y = ligand,size = weight)) +
  geom_point(shape = 21,fill = "#519E89",colour = "grey25",stroke = 0.25, aes(alpha = weight)) +
  scale_size_area(
    max_size = 7,
    name = "Regulatory potential"
  ) +
  labs(
    title = "Macrophages → Th17/Treg",
    subtitle = "anti-CTLA-4 vs PBS",
    x = "Predicted target genes",
    y = "Prioritized ligands"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      face = "italic"
    ),
    axis.text.y = element_text(
      face = "italic"
    )
  )


### Run all 
design_nichenet = tibble(
  sender = rep(
    c("Macrophages", "Neutrophils"),
    each = 3
  ),
  comparison = rep(
    c(
      "anti-CTLA-4 vs PBS",
      "anti-PD-1 vs PBS",
      "PBS vs Naive"
    ),
    times = 2
  ),
  potential_ligands = list(
    MAC_CTLA4_ligands,
    MAC_PD1_ligands,
    MAC_PBS_ligands,
    NEU_CTLA4_ligands,
    NEU_PD1_ligands,
    NEU_PBS_ligands
  ),
  geneset = rep(
    list(
      geneset_TH17_CTLA4,
      geneset_TH17_PD1,
      geneset_TH17_PBS
    ),
    times = 2
  ),
  background = rep(
    list(
      DDS_TH17_CTLA4$Gene,
      DDS_TH17_PD1$Gene,
      DDS_TH17_PBS$Gene
    ),
    times = 2
  )
)

design_nichenet
### run nichnet
nichenet_results = pmap(design_nichenet,run_nichenet)
names(nichenet_results) = paste(design_nichenet$sender,design_nichenet$comparison,sep = "_")

map_int(
  nichenet_results,
  function(x) {
    if (is.null(x)) {
      return(0L)
    }
    nrow(x$ligand_target)
  }
)
plot_df = map_dfr(compact(nichenet_results),"ligand_target")
head(plot_df)
dim(plot_df)

ligand_activities_all = imap_dfr(
  compact(nichenet_results),
  function(x, name) {
    x$ligand_activity %>%
      mutate(analysis = name)
  }
)
comparison_order = c(
  "anti-CTLA-4 vs PBS",
  "anti-PD-1 vs PBS",
  "PBS vs Naive"
)
plot_df = plot_df %>%
  mutate(
    sender = factor(
      sender,
      levels = c(
        "Macrophages",
        "Neutrophils"
      )
    ),
    comparison = factor(
      comparison,
      levels = comparison_order
    ),
    target_id = paste(
      comparison,
      target,
      sep = "__"
    )
  )

target_order = plot_df %>%
  group_by(comparison, target_id) %>%
  summarise(
    max_weight = max(weight),
    .groups = "drop"
  ) %>%
  arrange(
    comparison,
    desc(max_weight)
  ) %>%
  pull(target_id)

plot_df$target_id = factor(plot_df$target_id,levels = target_order)
ligand_order = plot_df %>%
  group_by(ligand) %>%
  summarise(
    activity = max(aupr_corrected),
    .groups = "drop"
  ) %>%
  arrange(activity) %>%
  pull(ligand)
plot_df$ligand = factor(plot_df$ligand,levels = ligand_order)

colors_comparison = c(
  "anti-CTLA-4 vs PBS" = "#519E89",
  "anti-PD-1 vs PBS"   = "#E9A6AE",
  "PBS vs Naive"      = "#5885B5"
)

p_nichenet_final = ggplot(plot_df,aes(x = target_id,y = ligand,size = weight,fill = comparison)) +
  geom_point(shape = 21,colour = "grey25",stroke = 0.25, alpha = .95) +
  facet_grid(sender ~ ., scales = "free_y",space = "free_y") +
  scale_x_discrete(
    labels = function(x) {
      sub("^.*__", "", x)
    }
  ) +
  scale_fill_manual(values = colors_comparison) +
  scale_size_area(max_size = 7,breaks = c(0.05, 0.10, 0.20),name = "Regulatory potential") +
  labs(
    x = "Predicted target genes in Th17/Treg cells",
    y = "Prioritized ligands",
    fill = "Comparisons"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 8,
      face = "italic",
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 9,
      face = "italic",
      colour = "black"
    ),
    strip.background = element_blank(),
    strip.text.y.right = element_text(angle = -90,size = 15,face = "bold"),
    panel.spacing.y = grid::unit(0.5,"cm"),
    axis.title.x = element_text(margin = margin(t = 15)),
    legend.position = "bottom",
    legend.box = "horizontal"
  ) +
  guides(fill = guide_legend(override.aes = list(size = 4)))



p_nichenet_final
ggsave(
  "Figure5a_p_nichenet_final.tiff",
  plot = p_nichenet_final,
  device = ragg::agg_tiff,
  width = 9,
  height = 8,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)
