# Manuscript figure mapping

| Figure panel | Script | Analysis | Required inputs | Main outputs |
|---|---|---|---|---|
| 1B | [code/Figure_1/figure_1B/code.r](../code/Figure_1/figure_1B/code.r) | CFU group linear model and estimated marginal means | combined_data_1.rds; CFU/cfu_exp1.csv; CFU/mapping_exp1.csv | forestplot.tiff |
| 1C | [code/Figure_1/figure_1C/code.r](../code/Figure_1/figure_1C/code.r) | Annotated lung-cell UMAP | combined_data_1.rds | Figure1B_UMAP_600dpi.tiff |
| 1D | [code/Figure_1/figure_1D/code.r](../code/Figure_1/figure_1D/code.r) | Animal-level composition and Dirichlet regression | combined_data_1.rds | heatmap_dirichlet.tiff |
| 1E–F | [code/Figure_1/Figure_1E_and_1F/code.r](../code/Figure_1/Figure_1E_and_1F/code.r) | Treatment-specific cellular prioritization from abundance | combined_data_1.rds; CFU/cfu_exp1.csv (loaded) | cellular_prioritization_antictla4.tiff; cellular_prioritization_antipd1.tiff (under results/figure1_audit) |
| 1G | [code/Figure_1/Figure_1G/code.r](../code/Figure_1/Figure_1G/code.r) | Cluster count fold changes and CFU association | combined_data_1.rds; CFU/mapping_exp1.csv; CFU/cfu_exp1.csv | barplot_logFC.tiff; pro_cell.tiff (under results/figure1_audit) |
| 2B | [code/Figure_2/figure_2b/code.r](../code/Figure_2/figure_2b/code.r) | Cell-type marker dot plot | combined_data_1.rds | Figure2_DotPlot_BonsaiColors.tiff |
| 2C input | [code/Figure_2/figure_2C/export_mtx/code_sketch.r](../code/Figure_2/figure_2C/export_mtx/code_sketch.r) | 5,000-cell leverage-score sketch | combined_data_1.rds | sketched_data_1.rds |
| 2C input | [code/Figure_2/figure_2C/export_mtx/export.r](../code/Figure_2/figure_2C/export_mtx/export.r) | Export sparse counts, identifiers and metadata | sketched_data_1.rds | data/input/{matrix.mtx, genes.tsv, cells.tsv, cell_metadata.csv} |
| 2D–E | [code/Figure_2/figure_2D-E/code.r](../code/Figure_2/figure_2D-E/code.r) | Sample-averaged relative-expression ternary plots | combined_data_1.rds | Ternary_Naive_PBS_PD1.tiff; Ternary_Naive_PBS_CTLA4.tiff |
| 2F | [code/Figure_2/figure_2F/code.r](../code/Figure_2/figure_2F/code.r) | Animal-averaged T-cell module scores; Kruskal–Wallis and Dunn tests | combined_data_1.rds | Figure2_Tcell_ModuleScores_Dunn.tiff |
| Supplementary 2A | [code/Figure_2/supp_figure_2A/code.r](../code/Figure_2/supp_figure_2A/code.r) | Marker heatmap, ternary-gene dot plots and cell-class scores | combined_data_1.rds; session top_vertex_genes | Figure2A_UMAP_600dpi.tiff |
| 3A | [code/Figure_3/figure_3A/code.r](../code/Figure_3/figure_3A/code.r) | T-cell selection, CD4 reclustering and UMAP | combined_data_1.rds | DIMT_supp_per group.tiff; ./data/t_cells_CD4.rds; DIMT.tiff |
| 3C | [code/Figure_3/figure_3C/code.r](../code/Figure_3/figure_3C/code.r) | Annotation-selected 1,908-cell population; contrast-specific DE ordering and DDRTree | data/combined_data_1.rds | results/trajectory/: four CDS objects, four ordering lists, plots and gene workbook |
| 3D–G | [code/Figure_3/figure_3D-G/pseudotime_plot.r](../code/Figure_3/figure_3D-G/pseudotime_plot.r) | Fixed signature scores along trajectories | Four CDS objects and ordering lists; original annotation-selected cells from combined_data_1.rds | results/trajectory_scores/: figure and four scored CDS objects |
| 3D–G companion | [code/Figure_3/figure_3D-G/relative_effect.r](../code/Figure_3/figure_3D-G/relative_effect.r) | Bootstrap relative effects for selected CD4 genes | t_cells_CD4.rds | CD4_relative_effects.tiff |
| 3H | [code/Figure_3/figure_3H/heatmap.r](../code/Figure_3/figure_3H/heatmap.r) | GO enrichment and functional-program classification | four cds_condition_*.rds | GO_programs_gene_supported_heatmap.tiff |
| 3H companion | [code/Figure_3/figure_3H/alluvial.r](../code/Figure_3/figure_3H/alluvial.r) | Enrichment-to-program alluvial plot | session program_summary_complete | alluvial.tiff |
| 3I | [code/Figure_3/figure_3I/code.r](../code/Figure_3/figure_3I/code.r) | CD4 marker feature plots | t_cells_CD4.rds | feature_plot.pdf |
| 3J | [code/Figure_3/figure_3J/code.r](../code/Figure_3/figure_3J/code.r) | Th17 module score split by experimental group | t_cells_CD4.rds | feature_plot.pdf |
| 4A | [code/Figure_4/figure_4A/code.r](../code/Figure_4/figure_4A/code.r) | Myeloid–T-cell communication population preparation | combined_data_1.rds; t_cells_CD4.rds | DIMT_sub.tiff |
| 4B–D | [code/Figure_4/figure_4B-D/code.r](../code/Figure_4/figure_4B-D/code.r) | Group-specific CellChat networks | session obj_cellchat_clean | communication_pbs.tiff; communication_pd1.tiff; communication_ctla4.tiff |
| 4E | [code/Figure_4/figure_4E/code.r](../code/Figure_4/figure_4E/code.r) | CellChat interaction-strength heatmaps | session cc_pbs_clean, cc_pd1_clean, cc_ctl_clean | communication_heatmap.tiff |
| 4F–G | [code/Figure_4/figure_4F-G/code.r](../code/Figure_4/figure_4F-G/code.r) | Myeloid–T-cell interaction-strength comparison | session cc_*_clean and get_comm_weight | communication_panel.tiff |
| 4H | [code/Figure_4/figure_4H/code.r](../code/Figure_4/figure_4H/code.r) | Neutrophil differential expression and GO enrichment | session obj_cellchat_clean | enrichment.tiff |
| 4I | [code/Figure_4/figure_4I/code.r](../code/Figure_4/figure_4I/code.r) | Neutrophil pathway-to-gene alluvial plot | session go_neut_focus | enrichment_alluvial.tiff |
| 5A | [code/Figure_5/figure_5A/code.r](../code/Figure_5/figure_5A/code.r) | Pseudobulk edgeR and NicheNet ligand-target inference | combined_data_1.rds; t_cells_CD4.rds; nichenet/*.rds; session pb_th17; live OmniPath | ./data/OmniNetworks_NNformat/lr_Network_Omnipath.rds; Figure5a_p_nichenet_final.tiff |
| 5B | [code/Figure_5/figure_5B/code.r](../code/Figure_5/figure_5B/code.r) | Sender ligand and receiver receptor expression | session obj_nichenet | results/nichenet/Figure5B.pdf; results/nichenet/Figure5B.png |
| 5C–D | [code/Figure_5/figure_5C-D/code.r](../code/Figure_5/figure_5C-D/code.r) | Ligand-target chord diagrams | session plot_df | ./Figure_pbs_naive.tiff; ./Figure_antictla_pbs.tiff |
| 2C | [code/Figure_2/figure_2C/sanity_bonsai/01_run_sanity.sh](../code/Figure_2/figure_2C/sanity_bonsai/01_run_sanity.sh) | Sanity inference | See [2C guide](../code/Figure_2/figure_2C/README.md) | results/sanity, results/bonsai, config or viewer (stage-specific) |
| 2C | [code/Figure_2/figure_2C/sanity_bonsai/02_run_bonsai.sh](../code/Figure_2/figure_2C/sanity_bonsai/02_run_bonsai.sh) | MPI Bonsai reconstruction | See [2C guide](../code/Figure_2/figure_2C/README.md) | results/sanity, results/bonsai, config or viewer (stage-specific) |
| 2C | [code/Figure_2/figure_2C/sanity_bonsai/03_prepare_scout.sh](../code/Figure_2/figure_2C/sanity_bonsai/03_prepare_scout.sh) | Scout annotation preparation | See [2C guide](../code/Figure_2/figure_2C/README.md) | results/sanity, results/bonsai, config or viewer (stage-specific) |
| 2C | [code/Figure_2/figure_2C/sanity_bonsai/05_launch_scout.sh](../code/Figure_2/figure_2C/sanity_bonsai/05_launch_scout.sh) | Launch local Scout viewer | See [2C guide](../code/Figure_2/figure_2C/README.md) | results/sanity, results/bonsai, config or viewer (stage-specific) |
| 2C | [code/Figure_2/figure_2C/sanity_bonsai/common.sh](../code/Figure_2/figure_2C/sanity_bonsai/common.sh) | Shared path/configuration checks | See [2C guide](../code/Figure_2/figure_2C/README.md) | results/sanity, results/bonsai, config or viewer (stage-specific) |
| 2C visualization | [code/Figure_2/figure_2C/sanity_bonsai/04_set_palette.py](../code/Figure_2/figure_2C/sanity_bonsai/04_set_palette.py) | Explicit palette utility | Bonsai JSON settings | Updated settings and .bak backup |


Panel labels are inferred from existing directories. Scripts for 1A, 2A and 3B are absent from the current tree; no scripts or results were invented. Figure 5C–D plot titles and supplementary 2A output names require author reconciliation with the manuscript.
