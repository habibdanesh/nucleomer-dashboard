#!/usr/bin/env bash
# Start the dashboard and open it in the browser.  Usage: ./run.sh [--port 8050]
set -euo pipefail
cd "$(dirname "$0")"
[ -x .venv/bin/python ] || { echo "Run ./setup.sh first."; exit 1; }
PORT=8050
for ((i = 1; i <= $#; i++)); do [ "${!i}" = "--port" ] && { j=$((i + 1)); PORT="${!j}"; }; done
URL="http://127.0.0.1:$PORT"
( for _ in $(seq 60); do curl -s -o /dev/null "$URL" && break; sleep 1; done
  if command -v open >/dev/null; then open "$URL"; elif command -v xdg-open >/dev/null; then xdg-open "$URL"; fi ) >/dev/null 2>&1 &
echo "Nucleomer dashboard at $URL  (Ctrl+C to stop)"
exec .venv/bin/python app.py "$@"
