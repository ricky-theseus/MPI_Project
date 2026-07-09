@echo off
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul 2>&1
cl /EHsc /I"%MSMPI_INC%" "D:\MPI_Project\course-project\03_HdEA_MP\TSP_HdEA_MP.cpp" /link /LIBPATH:"%MSMPI_LIB64%" msmpi.lib /OUT:"D:\MPI_Project\course-project\03_HdEA_MP\TSP_HdEA_MP.exe"
if %ERRORLEVEL% EQU 0 (echo BUILD SUCCESS) else (echo BUILD FAILED)
pause
