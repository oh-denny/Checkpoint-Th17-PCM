#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd -- "$SCRIPT_DIR/.." && pwd)"
# Optional machine-specific configuration; environment variables may override it.
if [[ -f "$REPO_DIR/config/pipeline.env" ]]; then
  # shellcheck disable=SC1091
  source "$REPO_DIR/config/pipeline.env"
fi
PROJECT_DIR="${PROJECT_DIR:-$REPO_DIR}"
BONSAI_REPO="${BONSAI_REPO:-}"
SANITY_HOME="${SANITY_HOME:-}"
SANITY_THREADS="${SANITY_THREADS:-8}"
BONSAI_MPI_PROCESSES="${BONSAI_MPI_PROCESSES:-4}"
BONSAI_ZSCORE="${BONSAI_ZSCORE:-1.0}"
BONSAI_KNN="${BONSAI_KNN:-10}"
SCOUT_PORT="${SCOUT_PORT:-8765}"
INPUT_DIR="$PROJECT_DIR/data/input"
SANITY_OUT="$PROJECT_DIR/results/sanity"
BONSAI_OUT="$PROJECT_DIR/results/bonsai"
CONFIG_OUT="$PROJECT_DIR/config/bonsai_sketch5000.yaml"
ANNOT_DIR="$PROJECT_DIR/annotations"
LOG_DIR="$PROJECT_DIR/logs"
mkdir -p "$INPUT_DIR" "$SANITY_OUT" "$BONSAI_OUT" "$ANNOT_DIR" "$LOG_DIR"
require_file() { [[ -f "$1" ]] || { echo "ERROR: missing file: $1" >&2; exit 1; }; }
require_config() {
  [[ -n "$BONSAI_REPO" && -d "$BONSAI_REPO" ]] || { echo 'ERROR: set BONSAI_REPO in config/pipeline.env' >&2; exit 1; }
  [[ -n "$SANITY_HOME" && -d "$SANITY_HOME" ]] || { echo 'ERROR: set SANITY_HOME in config/pipeline.env' >&2; exit 1; }
}
