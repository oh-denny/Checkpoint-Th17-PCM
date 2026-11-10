#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/common.sh"
require_config
require_file "$INPUT_DIR/cell_metadata.csv"
require_file "$BONSAI_REPO/bonsai_scout/bonsai_scout_preprocess.py"
cp "$INPUT_DIR/cell_metadata.csv" "$ANNOT_DIR/cell_metadata.csv"
cd "$BONSAI_REPO"
python bonsai_scout/bonsai_scout_preprocess.py \
  --results_folder "$BONSAI_OUT" \
  --annotation_path "$ANNOT_DIR" \
  --take_all_genes True \
  --config_filepath ''
