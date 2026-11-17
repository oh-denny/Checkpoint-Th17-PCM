
library(ggalluvial)
library(ggrepel)

# run figure 4H
go_alluvial = go_neut_focus |>
  dplyr::select(
    Treatment,
    Description,
    geneID,
    p.adjust
  ) |>
  separate_rows(
    geneID,
    sep = "/"
  )

colnames(go_alluvial)[colnames(go_alluvial) == "Description"] = "Pathway"
colnames(go_alluvial)[colnames(go_alluvial) == "geneID"] = "Gene"

gene_frequency = go_alluvial |>
  count(Gene, sort = TRUE)

head(gene_frequency, 30)
top_genes = gene_frequency |>
  slice_head(n = 20) |>
  pull(Gene)

go_alluvial_plot = go_alluvial |>
  filter(Gene %in% top_genes)

go_alluvial_plot = go_alluvial_plot |>
  distinct(
    Treatment,
    Pathway,
    Gene
  ) |>
  mutate(
    Weight = 1
  )

group_colors = c(
  "anti-PD-1"   = "#28331A", 
  "anti-CTLA-4" = "#9BC5F3"  
) 

p_alluvial <- ggplot(
  go_alluvial_plot,
  aes(
    axis1 = Pathway,
    axis2 = Gene,
    y = Weight
  )
) +
  geom_alluvium(
    aes(fill = Treatment),
    width = 0.08,
    alpha = 0.65,
    knot.pos = 0.45
  ) +
  geom_stratum(
    width = 0.10,
    fill = "white",
    color = "black",
    linewidth = 0.6
  ) +
  theme_minimal(base_size = 16) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  scale_x_discrete(limits = c("Description", "Gene"), expand = c(.1, .1)) +
  labs(x = "", y = "") +
  theme(
    legend.position = 'none'
  ) +
  scale_fill_manual(
    values = group_colors
  )

pb = ggplot_build(p_alluvial)
strata_data = pb$data[[2]] 

tiff("enrichment_alluvial.tiff",units = 'in', width = 9, height = 5, res = 300)
p_alluvial + 
  geom_text_repel(
    data = strata_data |> filter(x == 1),
    aes(x = x, y = y, label = stratum),
    inherit.aes = FALSE,
    force = 1,
    hjust = 1,
    nudge_x = -1.4,
    min.segment.length = 1,
    segment.curvature = 2,
    box.padding = unit(.02, "lines"),
    direction = "y",
    size = 4,
    segment.size = .3,
    color = 'black'
  ) +
  geom_text_repel(
    data = strata_data |> filter(x == 2),
    aes(x = x, y = y, label = stratum),
    inherit.aes = FALSE,
    force = 1,
    nudge_x = .4,
    color = 'black',
    #nudge_y = 2,
    min.segment.length = 1,
    segment.curvature = 2,
    direction = "y",
    box.padding = unit(.02, "lines"),
    size = 4,
    segment.size = .3,
    max.overlaps = getOption("ggrepel.max.overlaps", default = 1000)
  ) 
dev.off()
