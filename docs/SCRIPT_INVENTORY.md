# Script inventory and execution dependencies

Run R analyses from the repository root. Each script header states purpose, inputs, packages, outputs and prerequisites. The complete [figure mapping](FIGURE_MAPPING.md) lists every available panel and output. The [main README](../README.md) describes reproduction order and scientific limitations.

## Figure analysis scripts

| Script | Purpose | Required upstream state |
|---|---|---|
| [code/Figure_1/Figure_1E_and_1F/code.r](../code/Figure_1/Figure_1E_and_1F/code.r) | Treatment-specific cellular prioritization from abundance. | functions/functions.r. |
| [code/Figure_1/Figure_1G/code.r](../code/Figure_1/Figure_1G/code.r) | Cluster count fold changes and CFU association. | None after explicit CFU loading. |
| [code/Figure_1/figure_1B/code.r](../code/Figure_1/figure_1B/code.r) | CFU group linear model and estimated marginal means. | None; processed inputs required. |
| [code/Figure_1/figure_1C/code.r](../code/Figure_1/figure_1C/code.r) | Annotated lung-cell UMAP. | None; annotations and UMAP must already exist. |
| [code/Figure_1/figure_1D/code.r](../code/Figure_1/figure_1D/code.r) | Animal-level composition and Dirichlet regression. | None; processed inputs required. |
| [code/Figure_2/figure_2C/export_mtx/code_sketch.r](../code/Figure_2/figure_2C/export_mtx/code_sketch.r) | 5,000-cell leverage-score sketch. | None; overwrites sketched_data_1.rds if run. |
| [code/Figure_2/figure_2C/export_mtx/export.r](../code/Figure_2/figure_2C/export_mtx/export.r) | Export sparse counts, identifiers and metadata. | code_sketch.r or verified historical sketch. |
| [code/Figure_2/figure_2D-E/code.r](../code/Figure_2/figure_2D-E/code.r) | Sample-averaged relative-expression ternary plots. | functions/relative_expression.r. |
| [code/Figure_2/figure_2F/code.r](../code/Figure_2/figure_2F/code.r) | Animal-averaged T-cell module scores; Kruskal–Wallis and Dunn tests. | None; processed inputs required. |
| [code/Figure_2/figure_2b/code.r](../code/Figure_2/figure_2b/code.r) | Cell-type marker dot plot. | None; processed inputs required. |
| [code/Figure_2/supp_figure_2A/code.r](../code/Figure_2/supp_figure_2A/code.r) | Marker heatmap, ternary-gene dot plots and cell-class scores. | REVIEW: top_vertex_genes has no verified producer. |
| [code/Figure_3/figure_3A/code.r](../code/Figure_3/figure_3A/code.r) | T-cell selection, CD4 reclustering and UMAP. | None; overwrites historical t_cells_CD4.rds if run. |
| [code/Figure_3/figure_3C/code.r](../code/Figure_3/figure_3C/code.r) | DE ordering genes and DDRTree trajectories. | functions/functions.r, trajectory_input.r and patch/patch.r; fresh output directory. |
| [code/Figure_3/figure_3D-G/pseudotime_plot.r](../code/Figure_3/figure_3D-G/pseudotime_plot.r) | Fixed mean-expression module scores along pseudotime. | Historical trajectories; functions/trajectory.r; provenance review. |
| [code/Figure_3/figure_3D-G/relative_effect.r](../code/Figure_3/figure_3D-G/relative_effect.r) | Bootstrap relative effects for selected CD4 genes. | functions/rel.r; REffect installed in advance. |
| [code/Figure_3/figure_3H/alluvial.r](../code/Figure_3/figure_3H/alluvial.r) | Enrichment-to-program alluvial plot. | Run figure_3H/heatmap.r first in the same session. |
| [code/Figure_3/figure_3H/heatmap.r](../code/Figure_3/figure_3H/heatmap.r) | GO enrichment and functional-program classification. | Historical trajectories; functions/classify_enrich.r. |
| [code/Figure_3/figure_3I/code.r](../code/Figure_3/figure_3I/code.r) | CD4 marker feature plots. | None; historical CD4 object required. |
| [code/Figure_3/figure_3J/code.r](../code/Figure_3/figure_3J/code.r) | Th17 module score split by experimental group. | None; historical CD4 object required. |
| [code/Figure_4/figure_4A/code.r](../code/Figure_4/figure_4A/code.r) | Myeloid–T-cell communication population preparation. | None; review hard-coded cluster selections before regeneration. |
| [code/Figure_4/figure_4B-D/code.r](../code/Figure_4/figure_4B-D/code.r) | Group-specific CellChat networks. | Run figure_4A; functions/cellchat.r. |
| [code/Figure_4/figure_4E/code.r](../code/Figure_4/figure_4E/code.r) | CellChat interaction-strength heatmaps. | Run figure_4A then figure_4B-D. |
| [code/Figure_4/figure_4F-G/code.r](../code/Figure_4/figure_4F-G/code.r) | Myeloid–T-cell interaction-strength comparison. | Run figure_4A then figure_4B-D. |
| [code/Figure_4/figure_4H/code.r](../code/Figure_4/figure_4H/code.r) | Neutrophil differential expression and GO enrichment. | Run figure_4A. |
| [code/Figure_4/figure_4I/code.r](../code/Figure_4/figure_4I/code.r) | Neutrophil pathway-to-gene alluvial plot. | Run figure_4A then figure_4H. |
| [code/Figure_5/figure_5A/code.r](../code/Figure_5/figure_5A/code.r) | Pseudobulk edgeR and NicheNet ligand-target inference. | BLOCKED: pb_th17 construction absent; functions/nichenet.r. |
| [code/Figure_5/figure_5B/code.r](../code/Figure_5/figure_5B/code.r) | Sender ligand and receiver receptor expression. | Run figure_5A or provenance-verified historical obj_nichenet.rds. |
| [code/Figure_5/figure_5C-D/code.r](../code/Figure_5/figure_5C-D/code.r) | Ligand-target chord diagrams. | Run figure_5A; functions/circos.r. |

## Original trajectory population

`functions/trajectory_input.r::load_trajectory_population()` selects `T cells CD3+` and `T cells CD4+/CD3+` from the annotated integrated input `data/combined_data_1.rds`. This yields 1,908 cells, whose identifiers match all four stored historical CDS objects. The distinct 6,219-cell `t_cells_CD4.rds` is used for CD4 characterization rather than as the trajectory source.

Figure 3C retains the current `FindMarkers` settings (`min.pct=0.10`, `logfc.threshold=0.10`, `only.pos=FALSE`) and existing `select_de_genes` defaults and up/down counts. A separate CDS is fit for each ordering list using seed 1234, `lowerDetectionLimit=1`, VGAM negative-binomial expression family, DDRTree with two components and `ncenter=30`, and `orderCells()` without an explicit root state. No threshold or root was changed.

Regenerated outputs go to a new `results/trajectory/` directory. Use `options(paracocco.trajectory_data_dir="results/trajectory")` for downstream scoring and enrichment. The default remains historical inputs in `data/`. Score transfer requires the exact original population and writes separate outputs under `results/trajectory_scores/`. Historical ordering-gene discrepancies remain unresolved until numerical historical agreement is established.

## Shared helpers

| Helper | Responsibility |
|---|---|
| [functions/cellchat.r](../functions/cellchat.r) | CellChat inference and matrix lookup; caller supplies Seurat object. |
| [functions/circos.r](../functions/circos.r) | Chord plot helpers; caller supplies plot_df and mac_col, neu_col, tar_col. |
| [functions/classify_enrich.r](../functions/classify_enrich.r) | GO program signatures and annotation helpers; annotate_th17_terms is explicit. |
| [functions/functions.r](../functions/functions.r) | Composition prioritization and DE ordering-gene selection. |
| [functions/nichenet.r](../functions/nichenet.r) | NicheNet and edgeR helpers; caller supplies expr, meta, expression_pct, network resources and pseudobulk list. |
| [functions/patch/patch.r](../functions/patch/patch.r) | Monocle 2 compatibility patches for dplyr and igraph; apply patch_monocle2 explicitly. |
| [functions/rel.r](../functions/rel.r) | REffect wrapper; caller supplies expression table and comparisons. |
| [functions/relative_expression.r](../functions/relative_expression.r) | Ternary expression helpers; caller supplies expr, meta, comparisons, palette and plot tables. |
| [functions/trajectory.r](../functions/trajectory.r) | Fixed gene signatures, mean-expression scoring and pseudotime plotting. |
| [functions/trajectory_input.r](../functions/trajectory_input.r) | Original trajectory population, selected from the annotated integrated lung object. |

Monocle compatibility aliases and implementations are preserved. Some helpers intentionally depend on caller-created expression matrices, metadata, network models, palettes or tables; load the documented parent analysis first. `functions/classify_enrich.r` exposes GO annotation functions and signatures without requiring a table at source time.

## Sanity/Bonsai stages

All stages are under [code/Figure_2/figure_2C/sanity_bonsai](../code/Figure_2/figure_2C/sanity_bonsai/); follow the [2C guide](../code/Figure_2/figure_2C/README.md).

| Script | Prerequisite | Result |
|---|---|---|
| common.sh | Local software configuration | Shared paths and directory preparation |
| 01_run_sanity.sh | Exported Matrix Market inputs; Sanity installation | Sanity delta/d_delta and identifiers |
| 02_run_bonsai.sh | Sanity results; Bonsai/MPI environment | Configuration and hierarchy |
| 03_prepare_scout.sh | Bonsai results and metadata | Scout-ready annotations |
| 04_set_palette.py | Scout settings JSON | Explicit sequential palette update with backup |
| 05_launch_scout.sh | Prepared results/settings | Local interactive viewer |

## Retained validation tools

- [validate_repository.R](../scripts/validate_repository.R): parse R files and record dependency/I/O expressions without executing analyses; writes to a caller-selected or temporary directory.
- [validate_repository.py](../scripts/validate_repository.py): check relative Markdown links, Python/Bash syntax, literal R source paths, Git path spelling and large-file exclusions.
- [smoke_test_helpers.R](../scripts/smoke_test_helpers.R): load helpers and apply Monocle compatibility patches; does not reconstruct scientific results.
- [smoke_test_panel.R](../scripts/smoke_test_panel.R): execute reviewed output-only panels in isolated temporary directories.
- [smoke_test_trajectory.R](../scripts/smoke_test_trajectory.R): verify the original trajectory population against historical CDS identities and confirm complete fixed-signature score transfer without fitting trajectories or writing scientific objects.

## Remaining prerequisites and limitations

`pb_th17` (sample-level counts and metadata) and `top_vertex_genes` have no supplied construction. Do not replace them with guessed aggregation or a differently filtered ranking. Figure 4 session objects and Figure 5 predictions are not automatically persisted. Historical `cellchat_list.rds` and `obj_nichenet.rds` must be checked before reuse. Some output filenames conflict with directory panel labels; final panel definitions require manuscript confirmation. Raw preprocessing and public deposition identifiers remain unavailable.
