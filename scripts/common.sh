# Shared helpers for start.sh and start-wsl.sh (source, do not execute directly).

common_script_dir() {
  cd "$(dirname "${BASH_SOURCE[1]}")/.." && pwd
}

# Collect Python 3.9+ candidates without requiring Homebrew.
# Order: explicit versioned binaries, framework/pyenv, then generic python3.
_pick_python_candidates() {
  local v cmd path

  # Homebrew (optional — not required)
  for v in 12 11 10 9; do
    for path in \
      "/opt/homebrew/bin/python${v}" \
      "/opt/homebrew/bin/python3.${v}" \
      "/usr/local/bin/python${v}" \
      "/usr/local/bin/python3.${v}"; do
      [[ -x "$path" ]] && echo "$path"
    done
  done

  # python.org macOS installer
  for v in 12 11 10 9; do
    path="/Library/Frameworks/Python.framework/Versions/${v}/bin/python3"
    [[ -x "$path" ]] && echo "$path"
  done

  # pyenv / asdf shims
  if command -v pyenv >/dev/null 2>&1; then
    path="$(pyenv which python 2>/dev/null || true)"
    [[ -n "$path" && -x "$path" ]] && echo "$path"
  fi

  # Standard PATH names
  for cmd in python3.12 python3.11 python3.10 python3.9 python3; do
    command -v "$cmd" 2>/dev/null || true
  done

  # macOS system Python (often 3.9+ on recent macOS)
  [[ -x /usr/bin/python3 ]] && echo /usr/bin/python3
}

pick_python() {
  local cmd
  while IFS= read -r cmd; do
    [[ -z "$cmd" ]] && continue
    if "$cmd" -c 'import sys; raise SystemExit(0 if sys.version_info >= (3, 9) else 1)' 2>/dev/null; then
      echo "$cmd"
      return 0
    fi
  done < <(_pick_python_candidates | awk '!seen[$0]++')
  return 1
}

print_python_install_help() {
  local os
  os="$(uname -s 2>/dev/null || echo unknown)"
  echo ""
  echo "Python 3.9+ was not found. Install Python using one of these options:"
  echo ""
  if [[ "$os" == "Darwin" ]]; then
    echo "  macOS (no Homebrew required):"
    echo "    • Download the macOS installer from https://www.python.org/downloads/"
    echo "    • Run the installer, then re-run this script."
    echo ""
    echo "  macOS (optional, if you use Homebrew):"
    echo "    brew install python@3.12"
  elif grep -qi microsoft /proc/version 2>/dev/null; then
    echo "  WSL (Ubuntu/Debian):"
    echo "    sudo apt update"
    echo "    sudo apt install -y python3 python3-venv python3-pip"
    echo ""
    echo "  Or install Python on Windows, then use WSL:"
    echo "    https://www.python.org/downloads/windows/"
  else
    echo "  Linux:"
    echo "    sudo apt install -y python3 python3-venv python3-pip   # Debian/Ubuntu"
    echo "    sudo dnf install -y python3 python3-pip                   # Fedora"
  fi
  echo ""
  echo "  Windows (native, outside WSL):"
  echo "    Run start-windows.bat or:  .\\start-windows.ps1"
  echo ""
}

check_venv_module() {
  local py="$1"
  if "$py" -m venv --help >/dev/null 2>&1; then
    return 0
  fi
  echo "Error: '$py' cannot create virtual environments (venv module missing)." >&2
  if grep -qi microsoft /proc/version 2>/dev/null; then
    echo "  On WSL/Ubuntu, run: sudo apt install -y python3-venv" >&2
  fi
  return 1
}

open_browser_url() {
  local url="$1"
  if [[ "$(uname -s)" == "Darwin" ]]; then
    open "$url" >/dev/null 2>&1 || true
  elif grep -qi microsoft /proc/version 2>/dev/null; then
    if command -v wslview >/dev/null 2>&1; then
      wslview "$url" >/dev/null 2>&1 || true
    elif command -v powershell.exe >/dev/null 2>&1; then
      powershell.exe -NoProfile -Command "Start-Process '$url'" >/dev/null 2>&1 || true
    elif command -v cmd.exe >/dev/null 2>&1; then
      cmd.exe /c start "" "$url" >/dev/null 2>&1 || true
    fi
  elif command -v xdg-open >/dev/null 2>&1; then
    xdg-open "$url" >/dev/null 2>&1 || true
  fi
}

verify_imports() {
  python -c "
import pandas, numpy, streamlit, plotly, sklearn, scipy, statsmodels, gseapy
print('  OK: core packages')
try:
    import pydeseq2
    print('  OK: pydeseq2 (DGE)')
except ImportError as e:
    print('  WARN: pydeseq2 not available —', e)
try:
    import kaleido
    print('  OK: kaleido (PNG/PDF export)')
except ImportError:
    print('  WARN: kaleido missing — PNG/PDF export may fail')
"
}

is_unix_venv() {
  local venv_dir="$1"
  [[ -f "$venv_dir/bin/activate" ]]
}

is_windows_venv() {
  local venv_dir="$1"
  [[ -f "$venv_dir/Scripts/activate.bat" ]] \
    || [[ -f "$venv_dir/Scripts/Activate.ps1" ]] \
    || [[ -f "$venv_dir/Scripts/activate" ]]
}

on_windows_mount() {
  local path="$1"
  [[ "$path" == /mnt/* ]]
}

warn_wsl_windows_drive() {
  local project_dir="$1"
  if grep -qi microsoft /proc/version 2>/dev/null && on_windows_mount "$project_dir"; then
    echo "Note: Project is on a Windows drive ($project_dir)."
    echo "  A Linux .venv will be created with --copies (slower but reliable on /mnt/c)."
    echo "  For best performance, copy the repo to Linux home, e.g.:"
    echo "    cp -r \"$project_dir\" ~/RNAseq_app_anu && cd ~/RNAseq_app_anu && ./start-wsl.sh"
    echo ""
  fi
}

ensure_unix_venv() {
  local py="$1"
  local venv_dir="$2"
  local project_dir
  project_dir="$(dirname "$venv_dir")"

  check_venv_module "$py" || return 1

  if [[ -d "$venv_dir" ]] && ! is_unix_venv "$venv_dir"; then
    if is_windows_venv "$venv_dir"; then
      echo "==> Removing Windows .venv (Scripts/) — WSL needs bin/activate"
      echo "    On native Windows, use start-windows.bat instead."
    else
      echo "==> Removing incomplete .venv"
    fi
    rm -rf "$venv_dir"
  fi

  if [[ ! -d "$venv_dir" ]]; then
    echo "==> Creating virtual environment in .venv"
    if on_windows_mount "$project_dir"; then
      echo "    Using --copies (required for venvs on /mnt/c)"
      "$py" -m venv --copies "$venv_dir"
    else
      "$py" -m venv "$venv_dir"
    fi
  fi

  if ! is_unix_venv "$venv_dir"; then
    echo "Error: Could not create a WSL/Linux venv at $venv_dir" >&2
    echo "  Move the project off /mnt/c and try again, e.g.:" >&2
    echo "    cp -r \"$project_dir\" ~/RNAseq_app_anu && cd ~/RNAseq_app_anu" >&2
    echo "    ./start-wsl.sh" >&2
    return 1
  fi
}

activate_venv_unix() {
  local venv_dir="$1"
  if ! is_unix_venv "$venv_dir"; then
    echo "Error: No Linux venv at $venv_dir/bin/activate" >&2
    echo "  Run: ./start-wsl.sh --setup-only   (or delete .venv and re-run)" >&2
    return 1
  fi
  # shellcheck source=/dev/null
  source "$venv_dir/bin/activate"
}

setup_venv_unix() {
  local py="$1"
  local venv_dir="$2"
  local requirements="$3"

  warn_wsl_windows_drive "$(dirname "$venv_dir")"
  ensure_unix_venv "$py" "$venv_dir" || return 1
  activate_venv_unix "$venv_dir" || return 1

  echo "==> Upgrading pip"
  python -m pip install --upgrade pip wheel -q

  echo "==> Installing dependencies from requirements.txt"
  python -m pip install -r "$requirements" -q

  echo "==> Verifying imports"
  verify_imports

  echo "==> Setup complete."
}

start_streamlit_unix() {
  local app_file="$1"
  local port="$2"
  local venv_dir="$3"

  if [[ ! -f "$app_file" ]]; then
    echo "Error: app.py not found at $app_file" >&2
    return 1
  fi

  if command -v lsof >/dev/null 2>&1 && lsof -ti:"$port" >/dev/null 2>&1; then
    echo "Warning: port $port is already in use."
    echo "  Stop the other process or run: STREAMLIT_PORT=8502 $0"
  fi

  activate_venv_unix "$venv_dir" || return 1

  echo ""
  echo "==> Starting Streamlit on http://localhost:${port}"
  echo "    Press Ctrl+C to stop."
  echo ""

  (
    sleep 2
    open_browser_url "http://localhost:${port}"
  ) &

  exec streamlit run "$app_file" \
    --server.port="$port" \
    --server.headless=false \
    --browser.gatherUsageStats=false
}
