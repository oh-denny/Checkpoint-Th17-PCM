# Paracoccidioidomycosis scRNA-seq — Sanity / Bonsai pipeline

Reproducible scripts for preprocessing a Seurat sketch with **Sanity**, reconstructing a cell-state hierarchy with **Bonsai**, and exploring it with **Bonsai-scout**. The input is a **5,000-cell leverage-score sketch** of a pulmonary *Paracoccidioides brasiliensis* scRNA-seq dataset.

> **Scope:** exploratory transcriptomic hierarchy, **not** proof of cell-state transitions or causal treatment effects. The source RNA count matrix contained fractional values whose provenance was not confirmed; confirm statistical suitability before using derived results as definitive manuscript evidence.

## Repository structure

```text
.
├── config/
│   └── pipeline.env.example    # environment-specific paths (copy to pipeline.env)
├── data/input/                 # local inputs; ignored by Git
├── scripts/
│   ├── common.sh
│   ├── 01_run_sanity.sh
│   ├── 02_run_bonsai.sh
│   ├── 03_prepare_scout.sh
│   ├── 04_set_palette.py
│   └── 05_launch_scout.sh
├── results/                    # ignored by Git
├── annotations/                # ignored by Git
├── logs/                       # ignored by Git
└── .gitignore
```

## Requirements

- Linux/bash; Python >=3.7 (prefer the preconfigured Conda `bonsai` environment)
- Installed [Sanity](https://github.com/jmborghans/Sanity) and [Bonsai-data-representation](https://github.com/dhdegroot/Bonsai-data-representation) repositories
- MPI (`mpiexec`) and `mpi4py` for Bonsai; dependencies from the upstream Bonsai requirements
- Seurat-generated sparse Matrix Market input plus matching gene/cell identifiers

The upstream projects remain external dependencies and are **not vendored** into this repository. Pin the commit hashes of both upstream projects and Conda package versions for full reproducibility.

## Configure

```bash
cp config/pipeline.env.example config/pipeline.env
# Edit BONSAI_REPO and SANITY_HOME (absolute local install paths)
# Activate your Conda environment: conda activate bonsai
```

Put the following files under `data/input/`:

- `matrix.mtx` — genes × cells; **raw available RNA counts**, including fractional values, exported as Matrix Market `real general` without rounding
- `genes.tsv` — gene names in matrix row order
- `cells.tsv` — cell IDs in matrix column order
- `cell_metadata.csv` — includes `cell_id`, `cartig` (biological replicate), `groups`, and `subtype_cell`

In the demonstrated execution, the final exported matrix had **15,936 expressed genes × 5,000 cells** (1,609,362 stored nonzero entries). The original Seurat sketch sampled 5,000 of 26,200 cells using `SketchData(method="LeverageScore", seed=123)`; these scripts begin *after* exporting that sketch.

## Run

```bash
# From repository root; activate the environment first
bash scripts/01_run_sanity.sh
bash scripts/02_run_bonsai.sh
bash scripts/03_prepare_scout.sh
python scripts/04_set_palette.py results/bonsai/bonsai_vis_settings.json --to-map RdPu
bash scripts/05_launch_scout.sh
```

For long-running reconstruction, use tmux, e.g. `tmux new -s bonsai-paracocco`, then run step 02 inside it; detach with Ctrl+B, D. Avoid rerunning step 02 when a final tree already exists.

If scout is running on a remote host, open an SSH tunnel **on your local computer**:

```bash
ssh -N -L 8765:localhost:8765 USER@SERVER
```

Then open `http://localhost:8765` in a local browser. Choose an unused port if needed.

## Analysis settings used

| Stage | Parameters |
|---|---|
| Seurat SketchData | 5,000 cells, `LeverageScore`, seed 123 |
| Input matrix | 15,936 genes × 5,000 cells; fractional values preserved |
| Sanity | `-n 8 -e 1 -v_m MAP` |
| Bonsai | 4 MPI processes; `zscore_cutoff=1.0`; `use_knn=10`; `pickup_intermediate=False` |
| Bonsai-scout | `take_all_genes=True` |
| Visualization | Sequential colormap `RdPu` (original `YlOrRd`) |

## Reproducibility and interpretation

- Keep the original Seurat object and sample sheet outside the repository; do not upload sample-level data to public GitHub without authorization.
- Retain replicate metadata (`cartig`) and use **animals**, not cells, as independent experimental units for downstream treatment statistics.
- Bonsai visualizes transcriptomic similarity and hierarchy; group mixing is not evidence of no treatment effect.
- Aesthetic palettes do not change inferred expression values. Verify feature scaling before comparing colors across genes.
- Regeneration overwrites configuration and may consume significant compute time; use a fresh output directory for independent runs.

## Acknowledgments

Sanity and Bonsai are third-party research software. Cite their original publications and check repository licenses before redistributing code.
