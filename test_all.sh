#!/bin/bash
# test_all.sh — Run lint checks + all test suites for KSROP
#
# Usage: bash test_all.sh
#
# Runs:
#   1. Fortran lint check (lint_check.sh)
#   2. All compiled test executables
#   3. Reports overall pass/fail
#
# Prerequisites: all test executables must be compiled

cd "$(dirname "$0")"

TOTAL_PASS=0
TOTAL_FAIL=0
TOTAL_TESTS=0

echo "######################################################"
echo "#          KSROP Full Test & Lint Suite               #"
echo "######################################################"
echo ""

# --- Phase 1: Lint ---
echo "====== Phase 1: Code Lint ======"
if bash lint_check.sh; then
    echo ""
    echo "  >> Lint: PASSED"
else
    echo ""
    echo "  >> Lint: FAILED (errors found)"
    TOTAL_FAIL=$((TOTAL_FAIL + 1))
fi
echo ""

# --- Phase 2: Fortran tests ---
echo "====== Phase 2: Unit & Regression Tests ======"

TESTS="test_subrouts test_bugs test_tle test_tle2sv test_tle2opm"

for t in $TESTS; do
    if [ -f "./${t}.exe" ]; then
        EXE="./${t}.exe"
    elif [ -f "./${t}" ]; then
        EXE="./${t}"
    else
        echo "  [SKIP]  $t — executable not found"
        continue
    fi
    echo ""
    echo "  --- $t ---"
    OUTPUT=$($EXE 2>&1)
    SUMMARY=$(echo "$OUTPUT" | grep "Total:")
    PASSED=$(echo "$SUMMARY" | grep -o 'Passed: *[0-9]*' | grep -o '[0-9]*')
    FAILED=$(echo "$SUMMARY" | grep -o 'Failed: *[0-9]*' | grep -o '[0-9]*')
    TOTAL=$(echo "$SUMMARY" | grep -o 'Total: *[0-9]*' | grep -o '[0-9]*')

    if [ -n "$PASSED" ] && [ -n "$FAILED" ]; then
        TOTAL_PASS=$((TOTAL_PASS + PASSED))
        TOTAL_FAIL=$((TOTAL_FAIL + FAILED))
        TOTAL_TESTS=$((TOTAL_TESTS + TOTAL))
        if [ "$FAILED" -eq 0 ]; then
            echo "  >> $t: $TOTAL tests, ALL PASSED"
        else
            echo "  >> $t: $TOTAL tests, $FAILED FAILED"
            # Show failures
            echo "$OUTPUT" | grep '\[FAIL\]'
        fi
    else
        echo "  >> $t: could not parse output"
        TOTAL_FAIL=$((TOTAL_FAIL + 1))
    fi
done

# --- Phase 3: Python integration tests ---
echo ""
echo "====== Phase 3: Integration Tests ======"

if [ -f "driver_KS.exe" ]; then
    DRIVER="driver_KS.exe"
elif [ -f "driver_KS" ]; then
    DRIVER="driver_KS"
else
    DRIVER=""
fi

if [ -n "$DRIVER" ]; then
    if command -v python &>/dev/null || command -v python3 &>/dev/null; then
        PY=$(command -v python3 || command -v python)
        for ptest in test_driver.py test_initial_conditions.py; do
            if [ -f "$ptest" ]; then
                echo ""
                echo "  --- $ptest ---"
                OUTPUT=$($PY "$ptest" "$DRIVER" 2>&1) || true
                P=$(echo "$OUTPUT" | grep -o 'Passed: *[0-9]*' | grep -o '[0-9]*' | tail -1)
                F=$(echo "$OUTPUT" | grep -o 'Failed: *[0-9]*' | grep -o '[0-9]*' | tail -1)
                T=$(echo "$OUTPUT" | grep -o 'Total[a-z ]*: *[0-9]*' | grep -o '[0-9]*' | tail -1)
                if [ -n "$P" ] && [ -n "$F" ] && [ -n "$T" ]; then
                    TOTAL_PASS=$((TOTAL_PASS + P))
                    TOTAL_FAIL=$((TOTAL_FAIL + F))
                    TOTAL_TESTS=$((TOTAL_TESTS + T))
                    if [ "$F" -eq 0 ]; then
                        echo "  >> $ptest: $T checks, ALL PASSED"
                    else
                        echo "  >> $ptest: $T checks, $F FAILED"
                    fi
                else
                    echo "  >> $ptest: FAILED (could not parse output)"
                    echo "$OUTPUT" | tail -20
                    TOTAL_FAIL=$((TOTAL_FAIL + 1))
                fi
            fi
        done
    else
        echo "  [SKIP]  Python not available"
    fi
else
    echo "  [SKIP]  driver_KS.exe not compiled"
fi

# --- Summary ---
echo ""
echo "######################################################"
echo "#  GRAND TOTAL: $TOTAL_TESTS tests, $TOTAL_PASS passed, $TOTAL_FAIL failed"
echo "######################################################"

if [ $TOTAL_FAIL -gt 0 ]; then
    exit 1
else
    echo "#  ALL TESTS PASSED"
    exit 0
fi
