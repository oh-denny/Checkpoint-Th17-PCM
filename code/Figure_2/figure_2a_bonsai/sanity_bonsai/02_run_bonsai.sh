#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/common.sh"
require_config
require_file "$SANITY_OUT/delta.txt"
require_file "$SANITY_OUT/d_delta.txt"
require_file "$BONSAI_REPO/bonsai/create_config_file.py"
require_file "$BONSAI_REPO/bonsai/bonsai_main.py"
if [[ -d "$BONSAI_OUT/final_bonsai_zscore${BONSAI_ZSCORE}" ]]; then
  echo 'Final Bonsai output already exists; refusing to reconstruct/overwrite.' >&2
  exit 1
fi
cd "$BONSAI_REPO"
python bonsai/create_config_file.py \
  --new_yaml_path "$CONFIG_OUT" \
  --dataset paracocco_sketch5000 \
  --data_folder "$SANITY_OUT" \
  --results_folder "$BONSAI_OUT" \
  --input_is_sanity_output True \
  --zscore_cutoff "$BONSAI_ZSCORE" \
  --pickup_intermediate False \
  --use_knn "$BONSAI_KNN"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
mpiexec -n "$BONSAI_MPI_PROCESSES" python -m mpi4py bonsai/bonsai_main.py \
  --config_filepath "$CONFIG_OUT" --step all \
  2>&1 | tee "$LOG_DIR/bonsai.log"
