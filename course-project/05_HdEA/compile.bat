@echo off
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul 2>&1
if %ERRORLEVEL% NEQ 0 (echo vcvars failed & pause & exit /b 1)

set INC=C:\Program Files (x86)\Microsoft SDKs\MPI\Include
set LIBP=C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x64

echo === Building HdEA Ring Individual ===
cl /EHsc /I"%INC%" /O2 "%~dp0TSP_HdEA_ring.cpp" /Fe"%~dp0TSP_HdEA_ring.exe" /link /LIBPATH:"%LIBP%" msmpi.lib
if %ERRORLEVEL% EQU 0 (echo BUILD SUCCESS) else (echo BUILD FAILED & pause)
