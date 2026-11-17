# Software requirements

Versions below were observed in the local environment. They are not a validated reconstruction of the historical manuscript environment.

- R: **4.5.0** (2025-04-11), aarch64 macOS. Some installed packages warn that they were built under R 4.5.1–4.5.3; this mismatch was not repaired by upgrading the environment.
- Bioconductor: installed `BiocVersion` **3.21.1** (release series 3.21); no claim that all packages form a coherent validated Bioconductor installation.
- Python: **3.14.4** for local syntax and palette-utility smoke tests. This does not validate external Sanity/Bonsai with that Python version. The existing 2C guide states Python >=3.7; full external compatibility remains untested.
- Bash and Git are used for stage wrappers and repository checks. No `mpiexec` or `Sanity` executable was found on PATH. Conda exists locally, but no Bonsai environment or upstream checkout was established.

## Observed direct R requirements

The table is extracted from `library()`, namespace calls and relevant implicit algorithm backends. Base/recommended packages are included where explicitly used. Versions are observations, not newly selected pins. UCell is not called by current scripts.

| Package | Observed installed version |
|---|---|
| Biobase | 2.70.0 |
| CellChat | 1.6.1 |
| circlize | 0.4.18 |
| clusterProfiler | 4.16.0 |
| ComplexHeatmap | 2.24.1 |
| DDRTree | 0.1.6 |
| DirichletReg | 0.7-2 |
| dplyr | 1.2.1 |
| edgeR | 4.6.3 |
| emmeans | 2.0.3 |
| enrichplot | 1.28.4 |
| forcats | 1.0.1 |
| ggalluvial | 0.12.6 |
| ggplot2 | 4.0.3 |
| ggpubr | 0.6.3 |
| ggrepel | 0.9.8 |
| grDevices | 4.5.0 |
| grid | 4.5.0 |
| httpgd | 2.1.4 |
| igraph | 2.3.2 |
| Matrix | 1.7-5 |
| mgcv | 1.9-4 |
| monocle | 2.36.0 |
| nichenetr | 2.2.1 |
| OmnipathR | 3.16.2 |
| openxlsx | 4.2.8.1 |
| org.Mm.eg.db | 3.21.0 |
| patchwork | 1.3.2 |
| purrr | 1.2.2 |
| readr | 2.2.0 |
| REffect | 0.1.0 |
| rlang | 1.2.0 |
| rstatix | 0.7.3 |
| scales | 1.4.0 |
| Seurat | 5.5.0 |
| SeuratObject | 5.4.0 |
| stringr | 1.6.0 |
| tibble | 3.3.1 |
| tidyr | 1.3.2 |
| tidyverse | 2.0.0 |
| VennDiagram | 1.8.2 |
| VGAM | 1.1-14 |

Gene-table export uses `openxlsx::write.xlsx`; writexl and Augur are not requirements of current analyses.

## Package sources and compatibility

Installed DESCRIPTION metadata records these GitHub sources and exact revisions:

| Package | Repository | Observed revision |
|---|---|---|
| nichenetr 2.2.1 | saeyslab/nichenetr | `038f80e3eb7acf3264ba336bb9a8d2097f4e0c2a` |
| REffect 0.1.0 | oh-denny/relative_effect | `e5f01d1d9bd35e7c6680dd773a5172544b2bdcb1` |
| CellChat 1.6.1 | sqjin/CellChat | `e4f68625b074247d619c2e488d33970cc531e17c` |

`monocle` is **Monocle 2**, not monocle3. Retain [functions/patch/patch.r](../functions/patch/patch.r) and call `patch_monocle2()` after sourcing it when regenerating or plotting patched trajectories. The patch alters namespace functions for igraph DFS/neighbors/MST and dplyr selection compatibility. Repeated application was smoke-tested. The corrected four-contrast trajectory script also completed using the locally observed versions; equivalence to historical trajectories has not been established. Original aliases were preserved.

Seurat scripts use v5 `layer=` APIs, `JoinLayers`, `SketchData`, and `SeuratObject::LayerData`. Older Seurat versions and unjoined Assay5 layers may fail. CellChat 1.6.1 must be tested with these objects before Figure 4 is accepted. GO enrichment depends on the installed annotation snapshot `org.Mm.eg.db`; using a new snapshot can change terms without changing code. NicheNet resource hashes and the OmniPath retrieval date must be preserved separately from package versions.

## Installation pathway

Prefer the original author environment or an isolated copy of it. Do not run installation commands inside analytical scripts. A clean bootstrap, **not a validated reconstruction**, can use CRAN packages, Bioconductor packages and the observed GitHub revisions below. Network installation and restoration have not been tested.

```r
# Bootstrap tooling in an isolated library.
install.packages(c("BiocManager", "remotes"))
# For each CRAN package in the observed table, retain its observed version:
# remotes::install_version("PACKAGE", version = "OBSERVED_VERSION", upgrade = "never")
# Bioconductor release aligned with locally observed metadata:
BiocManager::install(version = "3.21", ask = FALSE, update = FALSE)
BiocManager::install(c("monocle", "DDRTree", "clusterProfiler", "org.Mm.eg.db",
                      "ComplexHeatmap", "edgeR", "enrichplot", "Biobase"),
                    ask = FALSE, update = FALSE)
remotes::install_github("saeyslab/nichenetr@038f80e3eb7acf3264ba336bb9a8d2097f4e0c2a", upgrade = "never")
remotes::install_github("oh-denny/relative_effect@e5f01d1d9bd35e7c6680dd773a5172544b2bdcb1", upgrade = "never")
remotes::install_github("sqjin/CellChat@e4f68625b074247d619c2e488d33970cc531e17c", upgrade = "never")
```

Bioconductor bootstrap can resolve versions different from the table; inspect and compare before scientific execution. The table and GitHub revisions alone do not pin transitive dependencies. Source builds may require C/C++/Fortran toolchains, BLAS/LAPACK and package-specific graphics/text libraries (notably ragg), plus MPI for Bonsai. Exact OS library versions and full external compilation requirements remain uncaptured.

## renv decision

`renv` is not installed locally. The shared R library is large, there is no historical lockfile, source restoration has not been tested, and a full pipeline cannot pass because upstream scientific inputs are missing. Therefore no fabricated or unvalidated `renv.lock` was generated. A dedicated renv project is appropriate **after** restoring the original author environment and resolving scientific blockers. Snapshot the actual working versions and GitHub commits, restore into a second clean library, and validate representative and full analyses before labeling a lockfile validated. A successful local execution alone does not establish historical reproducibility.

## Python and external tools

The only local Python analysis utility, `04_set_palette.py`, imports `argparse`, `json`, `pathlib` and `shutil` from the standard library; it needs no pip requirements file. The wrapper invokes external `mpi4py` through `python -m mpi4py`, Sanity, Bonsai and Bonsai-scout scripts. Their source trees, commits and full dependency versions are not supplied. Install their own declared environments and record exact revisions; an incomplete project `environment.yml` would misrepresent reproducibility. See the [2C guide](../code/Figure_2/figure_2C/README.md).

The workflow diagram was rendered from editable SVG using locally installed `rsvg` 2.7.0. This is a documentation-generation dependency, not an analytical requirement.
