@echo off
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul 2>&1
if %ERRORLEVEL% NEQ 0 (echo vcvars failed & pause & exit /b 1)

set INC=C:\Program Files (x86)\Microsoft SDKs\MPI\Include
set LIBP=C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x64
set SRC=D:\MPI_Project\course-project

echo === Building dEA ===
cl /EHsc /I"%INC%" "%SRC%\01_分布式进化算法\TSP_dEA.cpp" /Fe"%SRC%\TSP_dEA.exe" /link /LIBPATH:"%LIBP%" msmpi.lib
if %ERRORLEVEL% NEQ 0 (echo dEA FAILED) else (echo dEA OK)

echo === Building HdEA ===
cl /EHsc /I"%INC%" "%SRC%\02_分层分布式进化算法\TSP_HdEA.cpp" /Fe"%SRC%\TSP_HdEA.exe" /link /LIBPATH:"%LIBP%" msmpi.lib
if %ERRORLEVEL% NEQ 0 (echo HdEA FAILED) else (echo HdEA OK)

echo === Building HdEA-MP ===
cl /EHsc /I"%INC%" "%SRC%\03_HdEA_MP\TSP_HdEA_MP.cpp" /Fe"%SRC%\TSP_HdEA_MP.exe" /link /LIBPATH:"%LIBP%" msmpi.lib
if %ERRORLEVEL% NEQ 0 (echo HdEA-MP FAILED) else (echo HdEA-MP OK)

echo === All builds done ===
