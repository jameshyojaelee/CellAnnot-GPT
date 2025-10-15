#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VENV_BIN="${ROOT_DIR}/.venv/bin/python"

if [[ ! -x "${VENV_BIN}" ]]; then
  echo "[gpt-cell-annotator] Missing virtual environment. Run scripts/setup_env.sh first." >&2
  exit 1
fi

export PYTHONPATH=""
exec "${VENV_BIN}" -m pytest "$@"
