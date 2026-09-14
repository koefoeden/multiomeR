#!/usr/bin/env bash
set -euo pipefail
export PATH="$HOME/.pixi/bin:$HOME/.local/bin:/usr/local/bin:/usr/bin:/bin:$PATH"
exec python3 "$(dirname "$0")/watch.py" "$@"
