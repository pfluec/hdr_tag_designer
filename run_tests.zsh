#!/usr/bin/env zsh
set -eu

SCRIPT_DIR="${0:A:h}"
cd "$SCRIPT_DIR"
PYTHONPATH=. python -m unittest discover -s tests -v
PYTHONPATH=. python scripts/generate_tubb5_report.py
