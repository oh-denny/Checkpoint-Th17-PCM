#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/common.sh"
require_config
require_file "$BONSAI_REPO/bonsai_scout/run_bonsai_scout_app.py"
require_file "$BONSAI_OUT/bonsai_vis_settings.json"
cd "$BONSAI_REPO"
python bonsai_scout/run_bonsai_scout_app.py \
  --results_folder "$BONSAI_OUT" \
  --settings_filename bonsai_vis_settings.json \
  --port "$SCOUT_PORT"
