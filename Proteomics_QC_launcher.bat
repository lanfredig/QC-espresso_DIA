@echo off
cd /d "%~dp0"
setlocal EnableDelayedExpansion
title Proteomics QC espresso 0.2.7

:: ---------------------------------------------------------------------------
:: Banner / Header
:: ---------------------------------------------------------------------------
echo ===========================================================
echo   STARTING RENDER (PORTABLE MODE) Proteomics QC espresso
echo ===========================================================
echo   author: GL                                     v.0.2.7
echo ===========================================================

:: ---------------------------------------------------------------------------
:: PATHS (FIXED BASE_DIR ONLY)
:: ---------------------------------------------------------------------------
set "BASE_DIR=%~dp0"
for %%I in ("%BASE_DIR%") do set "BASE_DIR=%%~fI"

set "R_SCRIPTS_DIR=%BASE_DIR%app\Rscripts"
set "OUTPUT_DIR=%BASE_DIR%output\Results"
set "QMD_FILE=%R_SCRIPTS_DIR%\DIA_QC.qmd"
set "DIANN_REPORT=%BASE_DIR%output\DIANN-out\report.parquet"
set "FASTA_DIR=%BASE_DIR%data\FASTA"

if not exist "%OUTPUT_DIR%" mkdir "%OUTPUT_DIR%"

:: ---------------------------------------------------------------------------
:: 0. INPUT VALIDATION
:: ---------------------------------------------------------------------------
echo [INFO] Validating required input files...

if not exist "%QMD_FILE%" (
    echo [ERROR] QC script not found: %QMD_FILE%
    pause & exit /b 1
)

if not exist "%DIANN_REPORT%" (
    echo [ERROR] DIA-NN report not found: %DIANN_REPORT%
    pause & exit /b 1
)

if not exist "%FASTA_DIR%" (
    echo [ERROR] FASTA directory not found: %FASTA_DIR%
    pause & exit /b 1
)

echo [INFO] All required inputs found.

:: ===========================================================================
:: 1. QUARTO RESOLUTION
:: ===========================================================================
echo [INFO] Resolving Quarto...

set "QUARTO_CMD="

if exist "%BASE_DIR%app\runtime\quarto\bin\quarto.exe" (
    set "QUARTO_CMD=%BASE_DIR%app\runtime\quarto\bin\quarto.exe"
)

if "%QUARTO_CMD%"=="" (
    for /f "delims=" %%i in ('where quarto.exe 2^>nul') do set "QUARTO_CMD=%%i"
)

if "%QUARTO_CMD%"=="" (
    for /f "delims=" %%i in ('where quarto.cmd 2^>nul') do set "QUARTO_CMD=%%i"
)

if "%QUARTO_CMD%"=="" (
    echo [ERROR] Quarto not found.
    pause & exit /b 1
)

:: ===========================================================================
:: 2. R RESOLUTION
:: ===========================================================================
echo [INFO] Resolving R...

set "RSCRIPT_EXE="

for /d %%D in ("%BASE_DIR%app\runtime\R\R-*") do (
    if exist "%%D\bin\Rscript.exe" set "RSCRIPT_EXE=%%D\bin\Rscript.exe"
)

if "%RSCRIPT_EXE%"=="" (
    for /f "delims=" %%i in ('where Rscript.exe 2^>nul') do set "RSCRIPT_EXE=%%i"
)

if "%RSCRIPT_EXE%"=="" (
    echo [ERROR] Rscript.exe not found.
    pause & exit /b 1
)

echo [INFO] Using Quarto:  "%QUARTO_CMD%"
echo [INFO] Using Rscript: "%RSCRIPT_EXE%"

:: ===========================================================================
:: 3. ENVIRONMENT VARIABLES
:: ===========================================================================
set "QUARTO_EXEC=%QUARTO_CMD%"
set "QUARTO_R=%RSCRIPT_EXE%"

:: ===========================================================================
:: 4. RENDERING
:: ===========================================================================
echo [INFO] Rendering report...

pushd "%R_SCRIPTS_DIR%" >nul

"%QUARTO_CMD%" render "DIA_QC.qmd" ^
  --to html ^
  --output-dir "%OUTPUT_DIR%" ^
  -P report_path="%DIANN_REPORT%" ^
  -P fasta_path="%FASTA_DIR%" ^
  -P output_dir="%OUTPUT_DIR%"

set "RENDER_RC=%ERRORLEVEL%"
popd >nul

if not "%RENDER_RC%"=="0" (
    echo [ERROR] Rendering failed.
    pause & exit /b 1
)

:: ===========================================================================
:: 5. OUTPUT HANDLING
:: ===========================================================================
for /f %%I in ('powershell -NoProfile -Command "Get-Date -Format yyyyMMdd_HHmmss"') do set "TS=%%I"

set "ORIGINAL_HTML=%OUTPUT_DIR%\DIA_QC.html"
set "FINAL_HTML=%OUTPUT_DIR%\DIA_QC_%TS%.html"

if exist "%ORIGINAL_HTML%" (
    move /Y "%ORIGINAL_HTML%" "%FINAL_HTML%" >nul
    echo [SUCCESS] Report saved: "%FINAL_HTML%"
    start "" "%FINAL_HTML%"
) else (
    echo [WARNING] Rendered HTML file not found.
)

:: ---------------------------------------------------------------------------
:: Closing banner
:: ---------------------------------------------------------------------------
echo ===========================================================
echo   PROCESS FINISHED
echo ===========================================================
pause