#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
UV_BIN="${UV_BIN:-${HOME}/.local/bin/uv}"

if ! command -v "${UV_BIN}" >/dev/null 2>&1; then
  echo "[gpt-cell-annotator] Installing uv into ~/.local/bin (override with UV_BIN)" >&2
  python3 -m pip install --user --upgrade uv
fi

echo "[gpt-cell-annotator] Creating Python 3.11 virtual environment under ${ROOT_DIR}/.venv"
"${UV_BIN}" venv --python 3.11 "${ROOT_DIR}/.venv"

echo "[gpt-cell-annotator] Installing project requirements from requirements.txt"
"${UV_BIN}" pip install --python "${ROOT_DIR}/.venv/bin/python" -r "${ROOT_DIR}/requirements.txt"

cat <<'EOF'
[gpt-cell-annotator] Environment ready.

Activate it with:
  source .venv/bin/activate

Run commands with a clean interpreter path (recommended on HPC modules):
  PYTHONPATH= .venv/bin/python -m pytest
EOF
