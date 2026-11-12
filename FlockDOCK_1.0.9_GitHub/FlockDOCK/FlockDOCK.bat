@echo off
setlocal
cd /d "%~dp0"
set "ROOT=%~dp0"
set "PY=%ROOT%.venv\Scripts\python.exe"
echo ================================================================
echo                         FlockDOCK 1.0.9
echo ================================================================
if not exist "%PY%" (echo Environment missing. Run setup_FlockDOCK.bat first.& pause& exit /b 1)
"%PY%" -c "import pandas,PySide6,requests,Bio,rdkit,scipy,meeko,gemmi,rich" >nul 2>&1
if errorlevel 1 (echo Python dependencies are incomplete. Run setup_FlockDOCK.bat.& pause& exit /b 1)
"%PY%" -m flockdock_gui.main
set "RC=%ERRORLEVEL%"
if not "%RC%"=="0" (echo FlockDOCK exited with code %RC%.& pause)
exit /b %RC%
