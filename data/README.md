# Scientific data inputs

Files below are relative to `data/`. Large scientific objects are preserved locally and excluded individually from regular Git tracking. No verified study sequencing accession or processed-data deposit is supplied. Obtain a provenance-verified bundle from the authors; local availability is not public availability.

## Input populations and execution

- `combined_data_1.rds` is the annotated integrated lung Seurat object (26,200 cells). Original `subtype_cell` annotations, RNA layers, `cartig`, `groups`, cluster metadata and existing embeddings are required. Raw preprocessing and annotation provenance are not supplied.
- The **trajectory population** is selected directly from that integrated object: `T cells CD3+` (123) and `T cells CD4+/CD3+` (1,785), totaling **1,908 cells**. All four historical CDS objects contain exactly these cell identifiers. `functions/trajectory_input.r` implements the selection.
- `t_cells_CD4.rds` contains a **distinct 6,219-cell CD4 characterization population**. It is used for characterization and communication selections, not as the trajectory source. Figure 3A is its candidate regeneration workflow; historical cluster correspondence remains to be confirmed.
- `cds_condition_{PBS_naive,IGG_pbs,PD1_pbs,CTLA4_pbs}.rds` and `genes_{infection,igg,pd1,ctla4}.rds` are historical trajectory intermediates. Stored CDS ordering genes differ from the separate ordering-list files. Figure 3C now fits each contrast separately from the original population and writes a new bundle under `results/trajectory/`; it does not replace these files. Use one consistent bundle for downstream scoring/enrichment.
- `sketched_data_1.rds` is the 5,000-cell historical sketch. The Figure 2C exporter produces `data/input/matrix.mtx`, `genes.tsv`, `cells.tsv` and `cell_metadata.csv`. Local export reproduced 15,936 expressed genes and 1,609,362 nonzero entries, preserving fractional values and identifier order. Downstream Sanity/Bonsai reconstruction was not run.
- `CFU/cfu_exp1.csv` and `CFU/mapping_exp1.csv` hold experiment-1 CFU observations and animal mapping for Figures 1B/G; CFU is also loaded by 1E–F. The CSV template is not measured data. Original workbooks, Prism files and slides remain local under `CFU/EXP */`, excluded pending distribution review.
- NicheNet mouse resources under `nichenet/` support Figure 5A. Historical resource references include Zenodo record 7074291; model identity/hashes and redistribution terms must be verified. That record is not the study sequencing accession.
- `OmniNetworks_NNformat/lr_Network_Omnipath.rds` is a generated network snapshot. Figure 5A currently downloads and rewrites it. Preserve the historical snapshot and retrieval provenance before exact reproduction.

## Expected files

| File | Format / size | Role | Consumers and generation |
|---|---|---|---|
| `CFU/cfu_exp1.csv` | csv (473 bytes) | Processed data | 1B; 1E–F; 1G |
| `CFU/cfu_scrna_input_template.csv` | csv (529 bytes) | Input template; not observations | No current verified consumer; author provenance/distribution review |
| `CFU/mapping_exp1.csv` | csv (323 bytes) | Processed data | 1B; 1G |
| `OmniNetworks_NNformat/lr_Network_Omnipath.rds` | rds (18,494 bytes) | Generated external network snapshot | 5A writes snapshot; no current read |
| `cds_condition_CTLA4_pbs.rds` | rds (6,169,546 bytes) | Intermediate scientific result | 3D–G; 3H; historical trajectory; 3C regeneration writes a separate results/trajectory bundle |
| `cds_condition_IGG_pbs.rds` | rds (6,170,792 bytes) | Intermediate scientific result | 3D–G; 3H; historical trajectory; 3C regeneration writes a separate results/trajectory bundle |
| `cds_condition_PBS_naive.rds` | rds (6,174,200 bytes) | Intermediate scientific result | 3D–G; 3H; historical trajectory; 3C regeneration writes a separate results/trajectory bundle |
| `cds_condition_PD1_pbs.rds` | rds (6,169,061 bytes) | Intermediate scientific result | 3D–G; 3H; historical trajectory; 3C regeneration writes a separate results/trajectory bundle |
| `cellchat_list.rds` | rds (31,193,022 bytes) | Intermediate scientific result | No current verified consumer; author provenance/distribution review |
| `combined_data2.rds` | rds (55,433,436 bytes) | Processed data; producer unavailable | No current verified consumer; author provenance/distribution review |
| `combined_data_1.rds` | rds (124,477,358 bytes) | Processed data; producer unavailable | 1B; 1C; 1D; 1E–F; 1G; 2B; 2C input; 2D–E; 2F; Supplementary 2A; 3A; 3C; trajectory scores; 4A; 5A |
| `combined_data_integrated.rds` | rds (124,462,450 bytes) | Processed data; producer unavailable | No current verified consumer; author provenance/distribution review |
| `genes_ctla4.rds` | rds (492 bytes) | Intermediate scientific result | 3D–G; 3C regeneration writes a separate results/trajectory bundle |
| `genes_igg.rds` | rds (493 bytes) | Intermediate scientific result | 3D–G; 3C regeneration writes a separate results/trajectory bundle |
| `genes_infection.rds` | rds (546 bytes) | Intermediate scientific result | 3D–G; 3C regeneration writes a separate results/trajectory bundle |
| `genes_pd1.rds` | rds (441 bytes) | Intermediate scientific result | 3D–G; 3C regeneration writes a separate results/trajectory bundle |
| `nichenet/ligand_target_matrix_nsga2r_final_mouse.rds` | rds (191,048,057 bytes) | External model resource | 5A external model; original snapshot provenance required |
| `nichenet/lr_network_mouse_21122021.rds` | rds (32,993 bytes) | External model resource | 5A external model; original snapshot provenance required |
| `nichenet/weighted_networks_nsga2r_final_mouse.rds` | rds (77,390,723 bytes) | External model resource | 5A external model; original snapshot provenance required |
| `obj_nichenet.rds` | rds (12,267,107 bytes) | Intermediate scientific result | 5B |
| `sketched_data_1.rds` | rds (24,247,482 bytes) | Intermediate scientific result | 2C input; 2C input |
| `t_cells_CD4.rds` | rds (121,178,512 bytes) | Processed data | 3A; relative effects; 3I; 3J; 4A; 5A |

## Missing inputs and remaining limitations

`pb_th17` has no supplied producer or saved input. `functions/nichenet.r::run_pseudobulk_edger` expects `$counts` (genes × biological samples) and `$metadata` (rownames matching count columns and a `groups` column, with at least two samples per contrast group). The original aggregation, count provenance, filtering and sample definitions must be obtained; no replacement aggregation is inferred.

Supplementary 2A requires `top_vertex_genes`; the supplied `top_vertex_genes_final` ranking is not an automatic substitute. Communication objects (`obj_cellchat_clean`, `cc_*_clean`), enrichment summaries and NicheNet predictions are shared-session objects; follow the main README's execution order. Local `cellchat_list.rds` and `obj_nichenet.rds` are not automatically interchangeable with current cleaned populations or missing predictions.

RNA count layers contain fractional values. Do not round or redefine counts without confirming their provenance and methodological suitability. Exact historical ordering genes, trajectory root provenance, and numerical/visual agreement with manuscript outputs remain unresolved.

## Large-file distribution

Four local RDS files exceed 100 MiB and are excluded from ordinary Git tracking:

| File | Bytes |
|---|---:|
| combined_data_1.rds | 124,477,358 |
| combined_data_integrated.rds | 124,462,450 |
| t_cells_CD4.rds | 121,178,512 |
| nichenet/ligand_target_matrix_nsga2r_final_mouse.rds | 191,048,057 |

Other bulk processed/intermediate RDS files are also individually ignored. Small ordering lists, CFU CSVs and resource tables are not blanket-excluded. Preserve raw sequencing in a suitable sequencing archive and large processed objects/model snapshots in a versioned scientific archive, with persistent identifiers, file checksums and provenance. Author approval of data/resource distribution and a repository license are still needed. No files have been uploaded or permanently deleted.
