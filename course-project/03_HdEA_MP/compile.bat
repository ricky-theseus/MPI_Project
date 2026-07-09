@echo off
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul 2>&1
set "INC=C:\Program Files (x86)\Microsoft SDKs\MPI\Include\"
set "LIBP=C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x64\"
cl /EHsc /I"%INC%" "%~dp0TSP_HdEA_MP.cpp" /link /LIBPATH:"%LIBP%" msmpi.lib /OUT:"%~dp0TSP_HdEA_MP.exe"
if %ERRORLEVEL% EQU 0 (echo BUILD SUCCESS) else (echo BUILD FAILED)
