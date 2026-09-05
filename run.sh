#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

BIN=build/teleop-hands-latency-analysis

if [ $# -ne 2 ]; then
  echo "usage: run.sh INPUT.csv OUTPUT.jpg" >&2
  exit 1
fi
[ -f "$1" ] || { echo "no such file: $1" >&2; exit 1; }

make
exec "$BIN" "$1" "$2"
