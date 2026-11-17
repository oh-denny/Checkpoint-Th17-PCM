

library(tidyr)
library(dplyr)
library(ggplot2)

# Run figure 4-A

cellchat_summary_clean = data.frame(

  Interaction = c(
    "Macrophage → Treg/Th17",
    "Treg/Th17 → Macrophage",
    "Neutrophils → Treg/Th17",
    "Treg/Th17 → Neutrophils"
  ),

  PBS = c(
    get_comm_weight(
      cc_pbs_clean,
      "Macrophage",
      "Treg and Th17"
    ),

    get_comm_weight(
      cc_pbs_clean,
      "Treg and Th17",
      "Macrophage"
    ),

    get_comm_weight(
      cc_pbs_clean,
      "Neutrophils",
      "Treg and Th17"
    ),

    get_comm_weight(
      cc_pbs_clean,
      "Treg and Th17",
      "Neutrophils"
    )
  ),

  anti_PD1 = c(
    get_comm_weight(
      cc_pd1_clean,
      "Macrophage",
      "Treg and Th17"
    ),

    get_comm_weight(
      cc_pd1_clean,
      "Treg and Th17",
      "Macrophage"
    ),

    get_comm_weight(
      cc_pd1_clean,
      "Neutrophils",
      "Treg and Th17"
    ),

    get_comm_weight(
      cc_pd1_clean,
      "Treg and Th17",
      "Neutrophils"
    )
  ),

  anti_CTLA4 = c(
    get_comm_weight(
      cc_ctl_clean,
      "Macrophage",
      "Treg and Th17"
    ),

    get_comm_weight(
      cc_ctl_clean,
      "Treg and Th17",
      "Macrophage"
    ),

    get_comm_weight(
      cc_ctl_clean,
      "Neutrophils",
      "Treg and Th17"
    ),

    get_comm_weight(
      cc_ctl_clean,
      "Treg and Th17",
      "Neutrophils"
    )
  )
)

cellchat_long = cellchat_summary_clean |>
  pivot_longer(
    cols = c(
      PBS,
      anti_PD1,
      anti_CTLA4
    ),
    names_to = "group",
    values_to = "weight"
  ) |>
  mutate(
    group = factor(
      group,
      levels = c(
        "PBS",
        "anti_PD1",
        "anti_CTLA4"
      ),
      labels = c(
        "PBS",
        "anti-PD-1",
        "anti-CTLA-4"
      )
    ),

    Myeloid = case_when(
      grepl("Macrophage", Interaction) ~ "Macrophage",
      grepl("Neutrophils", Interaction) ~ "Neutrophils"
    )
  )

interaction_colors <- c(
  "Macrophage → Treg/Th17"   = "#396544",
  "Treg/Th17 → Macrophage"   = "#73A681",
  "Neutrophils → Treg/Th17"  = "#E53637",
  "Treg/Th17 → Neutrophils"  = "#E98989"
)

p_line = ggplot(
  cellchat_long,
  aes(
    x = group,
    y = weight,
    group = Interaction,
    color = Interaction
  )
) +

  geom_line(
    linewidth = 1
  ) +

  geom_point(
    size = 3.5
  ) +

  facet_wrap(
    ~ Myeloid,
    scales = "free_y"
  ) +

  scale_color_manual(
    values = interaction_colors
  ) +

  labs(
    x = NULL,
    y = "CellChat interaction strength",
    color = NULL
  ) +

  theme_classic(base_size = 14) +

  theme(
    legend.position = "top",
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    )
  )


tiff("communication_panel.tiff",units = 'in', width = 6, height = 4.5, res = 300)
p_line
dev.off()