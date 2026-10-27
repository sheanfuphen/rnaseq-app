#!/usr/bin/env bash
# RNAseq Analysis App — Windows Subsystem for Linux (WSL).
#
# Usage:
#   ./start-wsl.sh
#   ./start-wsl.sh --setup-only
#
# Opens the app in your Windows default browser from WSL.
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=scripts/common.sh
source "$SCRIPT_DIR/scripts/common.sh"

cd "$SCRIPT_DIR"

VENV_DIR="$SCRIPT_DIR/.venv"
REQUIREMENTS="$SCRIPT_DIR/requirements.txt"
APP_FILE="$SCRIPT_DIR/app.py"
PORT="${STREAMLIT_PORT:-8501}"

SETUP_ONLY=false
if [[ "${1:-}" == "--setup-only" ]]; then
  SETUP_ONLY=true
elif [[ -n "${1:-}" ]]; then
  echo "Usage: $0 [--setup-only]" >&2
  exit 1
fi

if ! grep -qi microsoft /proc/version 2>/dev/null; then
  echo "Note: This does not look like WSL. You can use ./start.sh on macOS/Linux instead."
fi

# Recommend WSL utilities for opening the browser
if ! command -v wslview >/dev/null 2>&1 && ! command -v powershell.exe >/dev/null 2>&1; then
  echo "Tip: install wslu for reliable browser launch from WSL:"
  echo "  sudo apt install -y wslu"
  echo ""
fi

PYTHON="$(pick_python)" || {
  echo "Error: Python 3.9+ is required but was not found in WSL." >&2
  print_python_install_help
  exit 1
}

echo "==> Using $($PYTHON --version) at $(command -v "$PYTHON" || echo "$PYTHON")"

setup_venv_unix "$PYTHON" "$VENV_DIR" "$REQUIREMENTS"

if $SETUP_ONLY; then
  exit 0
fi

start_streamlit_unix "$APP_FILE" "$PORT" "$VENV_DIR"
