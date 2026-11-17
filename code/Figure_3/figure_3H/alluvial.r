

### run heatmap first
library(ggrepel)
library(ggalluvial)
program_levels = c(
  "Less differentiated",
  "T-cell activation",
  "Migration and adhesion",
  "Th17-associated",
  "Treg-associated"
)
contrast_levels = c(
  "PBS vs Naive",
  "IgG vs PBS",
  "Anti-PD-1 vs PBS",
  "Anti-CTLA-4 vs PBS"
)

alluvial_data = program_summary_complete |>
  dplyr::mutate(
    Program = factor(
      as.character(Program),
      levels = program_levels
    ),
    Contrast = factor(
      as.character(Contrast),
      levels = contrast_levels
    ),
    Enriched = (
      Best_term != "Not enriched" &
      best_padj < 0.05 &
      minus_log10_padj > 0
    ),
    raw_weight = dplyr::if_else(
      Enriched,
      minus_log10_padj,
      0
    )
  ) |>
  dplyr::filter(
    Enriched
  )

alluvial_data = alluvial_data |>
  dplyr::group_by(Program) |>
  dplyr::mutate(
    program_total = sum(
      raw_weight,
      na.rm = TRUE
    ),
    Flow_weight = dplyr::if_else(
      program_total > 0,
      raw_weight / program_total,
      0
    )
  ) |>
  dplyr::ungroup()

alluvial_data |>
  dplyr::group_by(Program) |>
  dplyr::summarise(
    Total_node_width = sum(Flow_weight),
    .groups = "drop"
  )

program_colors = c(
  "Less differentiated"   = "#D99B7F",
  "T-cell activation"     = "#A56F63",
  "Migration and adhesion" = "#464858",
  "Th17-associated"       = "#618764",
  "Treg-associated"       = "#9CB080"
)

contrast_colors = c(
  "PBS vs Naive"        = "#6F7F5A",
  "IgG vs PBS"          = "#4C2A78",
  "Anti-PD-1 vs PBS"    = "#28331A",
  "Anti-CTLA-4 vs PBS"  = "#9BC5F3"
)
gap_factor = 0.80
alluvial_data = alluvial_data |>
  dplyr::mutate(
    flow_weight_plot = Flow_weight * gap_factor
  )

alluvial_data$Program = factor(alluvial_data$Program, levels = c("Th17-associated","Migration and adhesion", "T-cell activation" ,"Less differentiated"))
alluvial_data$Program = factor(alluvial_data$Program, levels = c("Less differentiated","T-cell activation","Migration and adhesion", "Th17-associated"))

p_program_alluvial <- ggplot(
  alluvial_data,
  aes(axis1 = Best_term,axis2 = Program,y = Flow_weight)) +
  geom_alluvium(
    aes(fill = Program),
    width = 0.10,
    alpha = 0.75,
    knot.pos = 0.45,
    color = NA
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
  scale_fill_manual(values = program_colors)

pb = ggplot_build(p_program_alluvial)
strata_data = pb$data[[2]] 

p_program_alluvial = p_program_alluvial + 
  geom_text_repel(
    data = strata_data |> filter(x == 1),
    aes(x = x, y = y, label = stratum),
    inherit.aes = FALSE,
    force = 1,
    hjust = 1,
    nudge_x = -1,
    min.segment.length = 1,
    segment.curvature = 2,
    box.padding = unit(.02, "lines"),
    direction = "y",
    size = 3,
    segment.size = .3,
    color = 'black'
  ) +
  geom_text_repel(
    data = strata_data |> filter(x == 2),
    aes(x = x, y = y, label = stratum),
    inherit.aes = FALSE,
    force = 1,
    nudge_x = .6,
    color = 'white',
    #nudge_y = 2,
    min.segment.length = 0,
    segment.curvature = 0,
    direction = "y",
    box.padding = unit(.02, "lines"),
    size = 0,
    segment.size = 0,
    max.overlaps = getOption("ggrepel.max.overlaps", default = 1000)
  ) 


tiff("alluvial.tiff",units = 'in', width = 7, height = 4, res = 300)
p_program_alluvial
dev.off()






