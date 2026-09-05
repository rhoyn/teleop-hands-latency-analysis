#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

make build/teleop-hands-recorder
exec build/teleop-hands-recorder "$@"
