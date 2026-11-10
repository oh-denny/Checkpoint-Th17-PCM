#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/common.sh"
require_config
for f in matrix.mtx genes.tsv cells.tsv; do require_file "$INPUT_DIR/$f"; done
require_file "$SANITY_HOME/scripts/sort_mtx_by_row.py"
require_file "$SANITY_HOME/bin/Sanity"
# The Python inside the active Conda environment must be >=3.7.
if [[ ! -s "$INPUT_DIR/matrix_sorted.mtx.gz" ]]; then
  python "$SANITY_HOME/scripts/sort_mtx_by_row.py" \
    -i "$INPUT_DIR/matrix.mtx" -o "$INPUT_DIR/matrix_sorted.mtx.gz"
fi
# Do not silently reuse or overwrite an incomplete Sanity run.
if [[ -e "$SANITY_OUT/delta.txt" || -e "$SANITY_OUT/d_delta.txt" ]]; then
  echo "Sanity output already exists at $SANITY_OUT; inspect/remove it before re-running." >&2
  exit 1
fi
"$SANITY_HOME/bin/Sanity" \
  -n "$SANITY_THREADS" \
  -f "$INPUT_DIR/matrix_sorted.mtx.gz" \
  -d "$SANITY_OUT" -e 1 -v_m MAP \
  -mtx_genes "$INPUT_DIR/genes.tsv" \
  -mtx_cells "$INPUT_DIR/cells.tsv" \
  2>&1 | tee "$LOG_DIR/sanity.log"
require_file "$SANITY_OUT/delta.txt"
require_file "$SANITY_OUT/d_delta.txt"
printf 'Genes: '; wc -l < "$SANITY_OUT/geneID.txt"
printf 'Cells: '; wc -l < "$SANITY_OUT/cellID.txt"
