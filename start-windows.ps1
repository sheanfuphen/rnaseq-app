# RNAseq Analysis App — native Windows (PowerShell).
#
# Usage:
#   .\start-windows.ps1
#   .\start-windows.ps1 -SetupOnly
#
param(
    [switch]$SetupOnly
)

$ErrorActionPreference = "Stop"
Set-Location $PSScriptRoot

$VenvDir = Join-Path $PSScriptRoot ".venv"
$Requirements = Join-Path $PSScriptRoot "requirements.txt"
$AppFile = Join-Path $PSScriptRoot "app.py"
$Port = if ($env:STREAMLIT_PORT) { $env:STREAMLIT_PORT } else { "8501" }

function Find-Python {
    $candidates = @(
        @("py", @("-3.12")),
        @("py", @("-3.11")),
        @("py", @("-3.10")),
        @("py", @("-3.9")),
        @("py", @("-3")),
        @("python", @())
    )
    foreach ($c in $candidates) {
        $exe = $c[0]
        $args = $c[1]
        if (-not (Get-Command $exe -ErrorAction SilentlyContinue)) { continue }
        try {
            & $exe @args -c "import sys; assert sys.version_info[:2] >= (3, 9)" 2>$null
            if ($LASTEXITCODE -eq 0) {
                return @{ Exe = $exe; Args = $args }
            }
        } catch { }
    }
    return $null
}

$py = Find-Python
if (-not $py) {
    Write-Host ""
    Write-Host "ERROR: Python 3.9+ was not found." -ForegroundColor Red
    Write-Host ""
    Write-Host "Install from https://www.python.org/downloads/windows/"
    Write-Host 'Enable "Add python.exe to PATH", then re-run this script.'
    Write-Host ""
    Write-Host "For WSL, run:  ./start-wsl.sh"
    exit 1
}

$pyCmd = $py.Exe
$pyArgs = $py.Args
Write-Host "==> Using:" -NoNewline
& $pyCmd @pyArgs --version

if (-not (Test-Path "$VenvDir\Scripts\Activate.ps1")) {
    Write-Host "==> Creating virtual environment in .venv"
    & $pyCmd @pyArgs -m venv $VenvDir
}

& "$VenvDir\Scripts\Activate.ps1"

Write-Host "==> Upgrading pip"
python -m pip install --upgrade pip wheel -q

Write-Host "==> Installing dependencies"
python -m pip install -r $Requirements -q

Write-Host "==> Verifying imports"
python -c "import pandas, numpy, streamlit, plotly, sklearn, scipy, statsmodels, gseapy; print('  OK: core packages')"

if ($SetupOnly) {
    Write-Host "==> Setup complete."
    exit 0
}

if (-not (Test-Path $AppFile)) {
    throw "app.py not found at $AppFile"
}

Write-Host ""
Write-Host "==> Starting Streamlit on http://localhost:$Port"
Write-Host "    Press Ctrl+C to stop."
Write-Host ""

Start-Job -ScriptBlock {
    param($url)
    Start-Sleep -Seconds 2
    Start-Process $url
} -ArgumentList "http://localhost:$Port" | Out-Null

streamlit run $AppFile --server.port=$Port --server.headless=false --browser.gatherUsageStats=false
