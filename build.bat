@echo off
:: -------------------------------------------------------
:: build.bat  -  Windows build script for KSROP
:: Requires Intel Fortran (ifort) on PATH
::
:: Usage:
::   build.bat             build propagator
::   build.bat tests       build unit test binary
::   build.bat test        build + run all tests
::   build.bat clean       delete build artefacts
:: -------------------------------------------------------

set FC=ifort
set FFLAGS=-O2 -warn:all
set SRC=driver_KS.F Subrouts.F Legendre.F
set TARGET=driver_KS.exe
set TEST_SRC=test_subrouts.F Subrouts.F Legendre.F
set TEST_BIN=test_subrouts.exe

if "%1"=="clean" goto :clean
if "%1"=="tests" goto :tests
if "%1"=="test"  goto :runtest

:: Default: build propagator
:build
echo Building %TARGET% ...
%FC% %FFLAGS% %SRC% -o %TARGET%
if errorlevel 1 (
    echo BUILD FAILED
    exit /b 1
)
echo Done: %TARGET%
goto :eof

:tests
echo Building %TEST_BIN% ...
%FC% %FFLAGS% %TEST_SRC% -o %TEST_BIN%
if errorlevel 1 (
    echo BUILD FAILED
    exit /b 1
)
echo Done: %TEST_BIN%
goto :eof

:runtest
call build.bat
if errorlevel 1 goto :eof
call build.bat tests
if errorlevel 1 goto :eof
echo.
echo === Unit tests ===
%TEST_BIN%
echo.
echo === Integration test ===
python test_driver.py %TARGET%
goto :eof

:clean
echo Cleaning build artefacts ...
del /q *.obj *.mod *.smod *.ilk *.pdb *.optrpt 2>nul
del /q %TARGET% %TEST_BIN% 2>nul
echo Done.
goto :eof
