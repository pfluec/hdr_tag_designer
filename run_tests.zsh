#!/usr/bin/env zsh
set -euo pipefail

ENV_NAME="hdr-tag-designer"
SCRIPT_DIR="${0:A:h}"
cd "$SCRIPT_DIR"

if ! conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo "Missing conda environment: $ENV_NAME"
  echo "Create it with: conda env create -f environment.yml"
  exit 1
fi

eval "$(conda shell.zsh hook)"
conda activate "$ENV_NAME"

PYTHONPATH=. python -m unittest discover -s tests -v
PYTHONPATH=. python scripts/generate_tubb5_report.py
