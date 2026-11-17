# Figure 2C: Seurat sketch → Sanity → Bonsai → Scout

Run from the repository root. This branch reconstructs exploratory transcriptomic similarity/hierarchy, not a causal treatment model. Counts in local Seurat objects include fractional values; their provenance and suitability for external methods require author review. No rounding or scientific parameters were changed during the audit.

## Inputs and software

The historical `data/combined_data_1.rds` is the input to [code_sketch.r](export_mtx/code_sketch.r), which saves `data/sketched_data_1.rds`. Alternatively use the verified historical sketch. [export.r](export_mtx/export.r) reads that saved object and writes:

| File under data/input/ | Meaning |
|---|---|
| matrix.mtx | Genes × cells, RNA counts; fractional values preserved |
| genes.tsv | Gene identifiers in matrix row order |
| cells.tsv | Cell identifiers in matrix column order |
| cell_metadata.csv | cell_id, cartig, groups, seurat_clusters, subtype_cell |

The local audit identified a 5,000-cell historical sketch. A temporary export during this audit reproduced the previously described 15,936 expressed genes × 5,000 cells and 1,609,362 stored entries. Reloaded Matrix Market values and gene/cell/metadata ordering matched the historical sketch. This validates export, not upstream sketch regeneration or downstream Sanity/Bonsai computation.

Install external [Sanity](https://github.com/jmborghans/Sanity) and [Bonsai-data-representation](https://github.com/dhdegroot/Bonsai-data-representation) checkouts with their upstream dependencies. The shell wrappers require Bash, Python, MPI (`mpiexec`) and `mpi4py`. The historical guide specified Python >=3.7; a complete working external environment and source revisions remain to be captured. The repository does not vendor those projects.

## Configuration

```bash
cp config/pipeline.env.example config/pipeline.env
# Edit BONSAI_REPO and SANITY_HOME to your installed checkouts.
```

`PROJECT_DIR` defaults to the repository root, derived from the wrapper location, so every stage uses the same root-level data/config/results directories. Use an absolute `PROJECT_DIR` for a separate prepared run directory. The example preserves existing defaults: Sanity 8 threads; Bonsai 4 MPI processes, z-score cutoff 1.0, KNN 10; Scout port 8765. Local config is ignored by Git.

## Execution order

Protect existing objects before regenerating a sketch; the R producer writes the same filename. Export also overwrites existing exported inputs. Use a fresh working copy/output area.

```r
source("code/Figure_2/figure_2C/export_mtx/code_sketch.r", print.eval = TRUE)
source("code/Figure_2/figure_2C/export_mtx/export.r", print.eval = TRUE)
```

Then, from the repository root with the external environment activated:

```bash
bash code/Figure_2/figure_2C/sanity_bonsai/01_run_sanity.sh
bash code/Figure_2/figure_2C/sanity_bonsai/02_run_bonsai.sh
bash code/Figure_2/figure_2C/sanity_bonsai/03_prepare_scout.sh
# Optional explicit visualization operation; preserves original categorical palettes.
python code/Figure_2/figure_2C/sanity_bonsai/04_set_palette.py results/bonsai/bonsai_vis_settings.json --to-map RdPu
bash code/Figure_2/figure_2C/sanity_bonsai/05_launch_scout.sh
```

Sanity produces `results/sanity/delta.txt`, `d_delta.txt`, geneID.txt and cellID.txt. Bonsai creates a YAML config under `config/` and reconstructed results under `results/bonsai/`. Scout preprocessing copies metadata to `annotations/`. Stage logs go to `logs/`. The viewer uses localhost port 8765 by default. No reconstruction or viewer was run during this audit.

Stages 01 and 02 guard against certain existing outputs; other intermediate/config files may still be overwritten. Never use these wrappers to validate against the only copy of historical results. The palette utility changes only matching sequential-map strings, writes a backup once, and must be invoked explicitly. Its execution on a synthetic JSON fixture was tested separately from scientific reproduction.

## Provenance and interpretation

Record source checkout commits, external environment versions, exported count provenance, input/output checksums and the final panel export process. Retain `cartig` sample identities and confirm biological replication. The code currently ends at an interactive viewer; it does not automate a final Figure 2C image export. Tree structure and group mixing alone do not establish cell-state transitions or causal treatment effects.
