@echo off
REM RNAseq Analysis App — native Windows (CMD). Uses Python from python.org / Microsoft Store (not Homebrew).
REM
REM Usage:
REM   start-windows.bat
REM   start-windows.bat setup-only
REM
setlocal EnableExtensions
cd /d "%~dp0"

set "VENV_DIR=%CD%\.venv"
set "REQUIREMENTS=%CD%\requirements.txt"
set "APP_FILE=%CD%\app.py"
if not defined STREAMLIT_PORT set "STREAMLIT_PORT=8501"

set "SETUP_ONLY=0"
if /I "%~1"=="setup-only" set "SETUP_ONLY=1"
if /I "%~1"=="--setup-only" set "SETUP_ONLY=1"

REM --- Find Python 3.9+ (py launcher preferred on Windows) ---
set "PYTHON_CMD="
where py >nul 2>&1
if not errorlevel 1 (
  py -3.12 -c "import sys; assert sys.version_info[:2] >= (3, 9)" 2>nul && set "PYTHON_CMD=py -3.12"
  if not defined PYTHON_CMD py -3.11 -c "import sys; assert sys.version_info[:2] >= (3, 9)" 2>nul && set "PYTHON_CMD=py -3.11"
  if not defined PYTHON_CMD py -3.10 -c "import sys; assert sys.version_info[:2] >= (3, 9)" 2>nul && set "PYTHON_CMD=py -3.10"
  if not defined PYTHON_CMD py -3.9  -c "import sys; assert sys.version_info[:2] >= (3, 9)" 2>nul && set "PYTHON_CMD=py -3.9"
  if not defined PYTHON_CMD py -3    -c "import sys; assert sys.version_info[:2] >= (3, 9)" 2>nul && set "PYTHON_CMD=py -3"
)

if not defined PYTHON_CMD (
  where python >nul 2>&1
  if not errorlevel 1 (
    python -c "import sys; assert sys.version_info[:2] >= (3, 9)" 2>nul && set "PYTHON_CMD=python"
  )
)

if not defined PYTHON_CMD (
  echo.
  echo ERROR: Python 3.9+ was not found.
  echo.
  echo Install Python for Windows:
  echo   1. https://www.python.org/downloads/windows/
  echo   2. Check "Add python.exe to PATH" during setup
  echo   3. Re-run start-windows.bat
  echo.
  echo Or from Microsoft Store: search for "Python 3.12"
  echo.
  echo For WSL instead, run:  bash start-wsl.sh
  echo.
  pause
  exit /b 1
)

echo ==^> Using:
%PYTHON_CMD% --version

REM --- Virtual environment ---
if not exist "%VENV_DIR%\Scripts\activate.bat" (
  echo ==^> Creating virtual environment in .venv
  %PYTHON_CMD% -m venv "%VENV_DIR%"
  if errorlevel 1 (
    echo ERROR: Failed to create venv. Try: %PYTHON_CMD% -m pip install virtualenv
    pause
    exit /b 1
  )
)

call "%VENV_DIR%\Scripts\activate.bat"

echo ==^> Upgrading pip
python -m pip install --upgrade pip wheel -q

echo ==^> Installing dependencies
python -m pip install -r "%REQUIREMENTS%" -q

echo ==^> Verifying imports
python -c "import pandas, numpy, streamlit, plotly, sklearn, scipy, statsmodels, gseapy; print('  OK: core packages')"
python -c "import pydeseq2; print('  OK: pydeseq2')" 2>nul || echo   WARN: pydeseq2 unavailable
python -c "import kaleido; print('  OK: kaleido')" 2>nul || echo   WARN: kaleido unavailable

echo ==^> Setup complete.

if "%SETUP_ONLY%"=="1" exit /b 0

if not exist "%APP_FILE%" (
  echo ERROR: app.py not found
  pause
  exit /b 1
)

echo.
echo ==^> Starting Streamlit on http://localhost:%STREAMLIT_PORT%
echo     Press Ctrl+C to stop.
echo.

REM Open default browser after a short delay
start "" cmd /c "timeout /t 2 /nobreak >nul && start http://localhost:%STREAMLIT_PORT%"

streamlit run "%APP_FILE%" --server.port=%STREAMLIT_PORT% --server.headless=false --browser.gatherUsageStats=false

endlocal
