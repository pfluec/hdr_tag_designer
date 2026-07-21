#!/usr/bin/env zsh
set -euo pipefail

ENV_NAME="hdr-tag-designer"
SCRIPT_DIR="${0:A:h}"

cd "$SCRIPT_DIR"

if ! conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo "Creating conda environment: $ENV_NAME"
  conda create -n "$ENV_NAME" python=3.11 -y
fi

eval "$(conda shell.zsh hook)"
conda activate "$ENV_NAME"

python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -m streamlit run app.py
