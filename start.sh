#!/usr/bin/env bash
# RNAseq Analysis App — macOS / Linux startup (Homebrew Python not required).
#
# Usage:
#   ./start.sh                 Install/update deps, start app, open browser
#   ./start.sh --setup-only    Create or refresh .venv only
#   ./start.sh --windows       Print Windows / WSL launcher instructions
#   ./start.sh -windows        Same as --windows
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
SHOW_WINDOWS_HELP=false

for arg in "$@"; do
  case "$arg" in
    --setup-only) SETUP_ONLY=true ;;
    --windows|-windows) SHOW_WINDOWS_HELP=true ;;
    -h|--help)
      sed -n '2,10p' "$0" | sed 's/^# \{0,1\}//'
      exit 0
      ;;
    *)
      echo "Usage: $0 [--setup-only | --windows]" >&2
      exit 1
      ;;
  esac
done

if $SHOW_WINDOWS_HELP; then
  echo "RNAseq app — Windows launchers"
  echo ""
  echo "  Native Windows (CMD):"
  echo "    start-windows.bat"
  echo ""
  echo "  Native Windows (PowerShell):"
  echo "    .\\start-windows.ps1"
  echo ""
  echo "  Windows Subsystem for Linux (WSL):"
  echo "    ./start-wsl.sh"
  echo ""
  echo "  From Git Bash on Windows, you can also run:"
  echo "    cmd.exe /c start-windows.bat"
  exit 0
fi

# If Git Bash on Windows, offer to delegate to the batch launcher
if [[ "$(uname -s)" == MINGW* ]] || [[ "$(uname -s)" == MSYS* ]]; then
  if [[ -f "$SCRIPT_DIR/start-windows.bat" ]]; then
    echo "==> Detected Git Bash on Windows — launching start-windows.bat"
    exec cmd.exe /c "\"$SCRIPT_DIR\\start-windows.bat\""
  fi
fi

PYTHON="$(pick_python)" || {
  echo "Error: Python 3.9+ is required but was not found." >&2
  print_python_install_help
  exit 1
}

echo "==> Using $($PYTHON --version) at $(command -v "$PYTHON" || echo "$PYTHON")"

setup_venv_unix "$PYTHON" "$VENV_DIR" "$REQUIREMENTS"

if $SETUP_ONLY; then
  exit 0
fi

start_streamlit_unix "$APP_FILE" "$PORT" "$VENV_DIR"
