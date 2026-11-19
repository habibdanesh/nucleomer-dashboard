#!/usr/bin/env bash
# One-time setup: Python environment in .venv, packages, and the data bundle.  Usage: ./setup.sh [bundle.tar.gz]
set -euo pipefail
cd "$(dirname "$0")"

PY=""
for c in python3.12 python3.13 python3.11 python3.10 python3 python; do
  if command -v "$c" >/dev/null 2>&1 && "$c" -c 'import sys; sys.exit(sys.version_info < (3, 10))' 2>/dev/null; then
    PY="$c"; break
  fi
done
if [ -z "$PY" ]; then
  echo "Python 3.10 or newer was not found. Install one (e.g. 'brew install python@3.12' or https://www.python.org/downloads/) and run ./setup.sh again."
  exit 1
fi
echo "Using $($PY --version) ($(command -v $PY))"

[ -x .venv/bin/python ] || "$PY" -m venv .venv
.venv/bin/python -m pip install --quiet --upgrade pip
echo "Installing packages (first time takes a minute)"
.venv/bin/python -m pip install --quiet -r requirements.txt
.venv/bin/python get_data.py "$@"
echo
echo "Done. Start the dashboard with:  ./run.sh"
