@echo off
setlocal EnableExtensions
cd /d "%~dp0"
set "ROOT=%~dp0"
set "VENV=%ROOT%.venv"
set "PY=%VENV%\Scripts\python.exe"
set "BASEPY="

echo ================================================================
echo                    FlockDOCK 1.0.9 SETUP
echo ================================================================
echo.

REM Prefer supported Python versions with reliable Windows wheels.
where py >nul 2>&1
if not errorlevel 1 (
    py -3.12 --version >nul 2>&1 && set "BASEPY=py -3.12"
    if not defined BASEPY py -3.11 --version >nul 2>&1 && set "BASEPY=py -3.11"
    if not defined BASEPY py -3.13 --version >nul 2>&1 && set "BASEPY=py -3.13"
)
if not defined BASEPY (
    where python >nul 2>&1 || (echo ERROR: Install Python 3.11, 3.12, or 3.13 from python.org. Enable the Python Launcher and/or Add Python to PATH.& goto :error)
    python -c "import sys; raise SystemExit(0 if (3,11) <= sys.version_info[:2] <= (3,13) else 1)" || (echo ERROR: Python 3.11, 3.12, or 3.13 is required. Detected version:& python --version& goto :error)
    set "BASEPY=python"
)

echo Using base interpreter: %BASEPY%
%BASEPY% --version || goto :error

REM An existing venv without pip is usually incomplete/corrupt. Rebuild it cleanly.
if exist "%VENV%" (
    if not exist "%PY%" goto :rebuild
    "%PY%" -m pip --version >nul 2>&1
    if errorlevel 1 goto :rebuild
    goto :install
)
:rebuild
echo.
echo Creating a clean virtual environment...
if exist "%VENV%" rmdir /s /q "%VENV%"
%BASEPY% -m venv "%VENV%"
if errorlevel 1 (
    echo ERROR: Could not create the virtual environment.
    goto :error
)

:install
if not exist "%PY%" (echo ERROR: Virtual-environment Python was not created.& goto :error)

REM Some Python distributions omit or fail to seed pip in new virtual environments.
echo Bootstrapping pip...
"%PY%" -m ensurepip --upgrade
if errorlevel 1 (
    echo ERROR: ensurepip could not install pip. Your base Python installation may be incomplete.
    echo Repair/reinstall Python 3.11-3.13 from python.org, including pip and the Python Launcher.
    goto :error
)
"%PY%" -m pip --version || goto :error

echo.
echo Installing FlockDOCK dependencies. This may take several minutes...
"%PY%" -m pip install --upgrade --retries 10 --timeout 60 pip setuptools wheel
if errorlevel 1 goto :network_error
"%PY%" -m pip install --prefer-binary --retries 10 --timeout 60 -r "%ROOT%requirements.txt"
if errorlevel 1 goto :network_error

echo.
echo Verifying imports...
"%PY%" -c "import pandas,PySide6,requests,Bio,rdkit,scipy,meeko,gemmi,rich" || goto :error

echo.
echo ================================================================
echo FlockDOCK Windows environment installed successfully.
echo ================================================================
echo GPU docking additionally requires WSL2 Ubuntu and an NVIDIA driver with WSL CUDA support.
echo P2Rank and Java 17 are checked when automatic pocket prediction is run.
echo Run FlockDOCK.bat to launch the application.
pause
exit /b 0

:network_error
echo.
echo ERROR: Python packages could not be downloaded from PyPI.
echo Check that this computer has a working internet connection and that
echo a firewall, proxy, or VPN is not blocking pypi.org or files.pythonhosted.org.
echo Then rerun setup_FlockDOCK.bat. The installer already retried the download.
pause
exit /b 1

:error
echo.
echo SETUP FAILED.
echo If pip bootstrap failed, repair/reinstall Python 3.11, 3.12, or 3.13 from python.org.
echo Ensure the installer includes pip, venv, and the Python Launcher, then rerun this script.
echo Full setup log is shown above.
pause
exit /b 1
