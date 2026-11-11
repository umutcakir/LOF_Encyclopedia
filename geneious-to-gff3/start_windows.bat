@echo off
rem Starts Geneious to GFF3 version 1.1. Double-click this file.
rem Nothing needs to be installed except Python: Flask is included in the "vendor" folder.
setlocal
cd /d "%~dp0"

rem This version keeps its jobs separate from version 1.0, so the two can live side by side.
if not defined GFF3TOOL_DATA set "GFF3TOOL_DATA=%LOCALAPPDATA%\Geneious_to_GFF3_v1_1\jobs"
if not defined GFF3TOOL_PORT set "GFF3TOOL_PORT=8775"

set "PY="
where py >nul 2>nul && set "PY=py -3"
if not defined PY where python >nul 2>nul && set "PY=python"
if not defined PY goto nopython
%PY% -c "import sys; sys.exit(sys.version_info < (3, 9))" >nul 2>nul || goto nopython

%PY% app.py
pause
exit /b 0

:nopython
echo Python 3.9 or newer was not found on this computer.
echo Install it from https://www.python.org/downloads/ and tick "Add python.exe to PATH"
echo during the installation. Then double-click this file again.
pause
exit /b 1
