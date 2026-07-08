@echo off
REM run_tc.bat TC<n>  [nrev]
REM   Copies benchmark\TC<n>.opm to input\input.opm, runs driver_KS and
REM   ksbench, then calls compare_oem.py on the two output OEM files.
REM   Optional second argument overrides the revolution count in input.dat
REM   (default: uses whatever is already in input\input.dat).
REM
REM Example: run_tc.bat TC1
REM          run_tc.bat TC3 50

if "%1"=="" (
    echo Usage: run_tc.bat TC^<n^>
    exit /b 1
)

set TC=%1
set OPM=benchmark\%TC%.opm

if not exist "%OPM%" (
    echo ERROR: %OPM% not found.
    exit /b 1
)

echo === Benchmark run: %TC% ===
copy /Y "%OPM%" input\input.opm > nul
echo Copied %OPM% to input\input.opm

echo Running driver_KS ...
driver_KS.exe
if errorlevel 1 ( echo driver_KS FAILED & exit /b 1 )

echo Running ksbench ...
ksbench.exe
if errorlevel 1 ( echo ksbench FAILED & exit /b 1 )

REM Find the most recently written KSROP OEM
for /f "delims=" %%F in ('dir /b /o-d output\KSROP_*.oem 2^>nul') do (
    set KSROP_OEM=output\%%F
    goto :found
)
echo ERROR: no KSROP_*.oem found in output\
exit /b 1
:found

set BENCH_OEM=output\ksbench.oem
echo.
echo Comparing: %KSROP_OEM%  vs  %BENCH_OEM%
python benchmark\compare_oem.py %KSROP_OEM% %BENCH_OEM% --csv output\%TC%_error.csv

echo.
echo Done.  Error CSV: output\%TC%_error.csv
