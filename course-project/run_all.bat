@echo off
set RUNS=10
set GENS=200000
set NP=4

echo ========================================
echo TSP Parallel EA Experiment Runner
echo Runs=%RUNS% MaxGen=%GENS% Procs=%NP%
echo Start: %DATE% %TIME%
echo ========================================

echo.
echo --- Running Algorithm 1: dEA ---
mpiexec -n %NP% "01_???_??_??\TSP_dEA.exe" %RUNS% %GENS%

echo.
echo --- Running Algorithm 2: HdEA ---
mpiexec -n %NP% "02_???_???_??_??\TSP_HdEA.exe" %RUNS% %GENS%

echo.
echo --- Running Algorithm 3: HdEA-MP ---
mpiexec -n %NP% "03_????_???_???_??_??_??\TSP_HdEA_MP.exe" %RUNS% %GENS%

echo.
echo ========================================
echo All experiments complete!
echo End: %DATE% %TIME%
echo ========================================
pause