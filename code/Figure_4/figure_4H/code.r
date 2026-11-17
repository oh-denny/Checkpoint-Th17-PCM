# Figure/panel: 4H
# Purpose: Neutrophil differential expression and GO enrichment.
# Inputs (under data/ unless session objects): session obj_cellchat_clean.
# Packages: Seurat, clusterProfiler, dplyr, ggplot2, org.Mm.eg.db, tidyr.
# Outputs: enrichment.tiff.
# Prerequisites: Run figure_4A.
# Execution: repository root; see README and docs/SCRIPT_INVENTORY.md.

library(Seurat)
if (!exists("obj_cellchat_clean", inherits = TRUE)) stop("Missing obj_cellchat_clean. Run Figure 4A in this R session first.")


library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(tidyr)
library(ggplot2)

Idents(obj_cellchat_clean) = "comm_celltype"

neut = subset(obj_cellchat_clean,idents = "Neutrophils")
Idents(neut) = "groups"
de_neut_pd1 = FindMarkers(
  neut,
  ident.1 = "antipd1",
  ident.2 = "pbs",
  assay = "RNA",
  min.pct = 0.10,
  logfc.threshold = 0.25
)
de_neut_ctl = FindMarkers(
  neut,
  ident.1 = "antictla4",
  ident.2 = "pbs",
  assay = "RNA",
  min.pct = 0.10,
  logfc.threshold = 0.25
)

genes_pd1_neut = rownames(
  de_neut_pd1[
    de_neut_pd1$p_val < 0.05 &
    de_neut_pd1$avg_log2FC > 0,
  ]
)
genes_ctl_neut = rownames(
  de_neut_ctl[
    de_neut_ctl$p_val < 0.05 &
    de_neut_ctl$avg_log2FC > 0,
  ]
)
go_pd1_neut = enrichGO(
  gene = genes_pd1_neut,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP"
)
go_ctl_neut = enrichGO(
  gene = genes_ctl_neut,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP"
)

grep(
  "chemotaxis|migration|chemokine|neutrophil|granulocyte",
  go_pd1_neut@result$Description,
  ignore.case = TRUE,
  value = TRUE
)

grep(
  "chemotaxis|migration|chemokine|neutrophil|granulocyte",
  go_ctl_neut@result$Description,
  ignore.case = TRUE,
  value = TRUE
)

pathways_focus = c(
  "neutrophil chemotaxis",
  "neutrophil migration",
  "granulocyte chemotaxis",
  "granulocyte migration",
  "leukocyte chemotaxis",
  "leukocyte migration involved in inflammatory response",
  "response to chemokine",
  "chemokine-mediated signaling pathway",
  "regulation of neutrophil chemotaxis",
  "regulation of neutrophil migration"
)

go_pd1_df = go_pd1_neut@result |>
  as.data.frame() |>
  filter(Description %in% pathways_focus) |>
  mutate(Treatment = "anti-PD-1")

go_ctl_df = go_ctl_neut@result |>
  as.data.frame() |>
  filter(Description %in% pathways_focus) |>
  mutate(Treatment = "anti-CTLA-4")

go_neut_focus = bind_rows(
  go_pd1_df,
  go_ctl_df
)

go_neut_focus = go_neut_focus |>
  mutate(
    GeneRatio_num = sapply(
      GeneRatio,
      function(x) {
        z <- strsplit(x, "/")[[1]]
        as.numeric(z[1]) / as.numeric(z[2])
      }
    ),

    minus_log10_padj = -log10(p.adjust)
  )

p_go_neut = ggplot(
  go_neut_focus,
  aes(
    x = Treatment,
    y = Description
  )
) +
  geom_point(
    aes(
      size = GeneRatio_num,
      color = minus_log10_padj
    )
  ) +
  scale_color_gradient(
    low = "#F4B6B6",
    high = "#B2182B"
  ) +
  scale_size(
    range = c(3, 10)
  ) +
  labs(
    x = NULL,
    y = NULL,
    size = "Gene ratio",
    color = expression(-log[10]("adjusted p-value"))
  ) +
  theme_classic(base_size = 18) +
  theme(
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 18
    )
  )

tiff("enrichment.tiff",units = 'in', width = 12, height = 6, res = 300)
p_go_neut
dev.off()