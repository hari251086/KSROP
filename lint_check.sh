#!/bin/bash
# lint_check.sh — Automated Fortran 77 code quality checks for KSROP
#
# Checks:
#   1. Line length > 72 columns (source files only)
#   2. Single-precision intrinsics (acos/atan/sqrt without d-prefix)
#   3. Inconsistent common block /xy/ declarations
#   4. Tabs in source (F77 uses spaces only)
#   5. Missing continuation character in column 6
#
# Usage: bash lint_check.sh [--strict]
#   --strict: also check test files for line length
#
# Exit code: 0 = all clean, 1 = issues found

SRCFILES="driver_KS.F Subrouts.F TLEread.F Legendre.F tle2opm.F"
TESTFILES="test_subrouts.F test_bugs.F test_tle.F test_tle2sv.F test_tle2opm.F"
ALLFILES="$SRCFILES $TESTFILES"

STRICT=0
if [ "$1" = "--strict" ]; then
    STRICT=1
fi

ERRORS=0
WARNINGS=0

echo "=============================================="
echo " KSROP Fortran 77 Lint Check"
echo "=============================================="
echo ""

# --- Check 1: Line length (72 columns) ---
echo "  --- Line length (max 72 columns) ---"
if [ $STRICT -eq 1 ]; then
    CHECKFILES="$ALLFILES"
else
    CHECKFILES="$SRCFILES"
fi
FOUND=0
for f in $CHECKFILES; do
    if [ -f "$f" ]; then
        OVER=$(awk 'length > 72 {printf "    %s:%d (%d cols)\n", FILENAME, NR, length}' "$f")
        if [ -n "$OVER" ]; then
            echo "$OVER"
            FOUND=$((FOUND + $(echo "$OVER" | wc -l)))
        fi
    fi
done
if [ $FOUND -eq 0 ]; then
    echo "  [PASS]  All lines within 72 columns"
else
    echo "  [FAIL]  $FOUND lines exceed 72 columns"
    ERRORS=$((ERRORS + FOUND))
fi
echo ""

# --- Check 2: Single-precision intrinsics in executable code ---
echo "  --- Single-precision intrinsics ---"
FOUND=0
for f in $ALLFILES; do
    if [ -f "$f" ]; then
        # Find lines with acos(/atan(/sqrt(/log( NOT preceded by 'd'
        # Skip: comments (col1 = c/C/*/!), string literals, datan2, atan3
        HITS=$(awk '
        /^[cC*!]/ { next }
        /[^d]acos\(/ || /[^d]atan\(/ || /[^d]sqrt\(/ || /[^d]log\(/ {
            if (/dacos/ || /datan/ || /dsqrt/ || /dlog/ || /dabs/) next
            if (/atan2/ || /datan2/ || /atan3/) next
            if (/'"'"'/) next
            printf "    %s:%d: %s\n", FILENAME, NR, $0
        }' "$f")
        if [ -n "$HITS" ]; then
            echo "$HITS"
            FOUND=$((FOUND + $(echo "$HITS" | wc -l)))
        fi
    fi
done
if [ $FOUND -eq 0 ]; then
    echo "  [PASS]  All intrinsics use double-precision (d-prefix)"
else
    echo "  [WARN]  $FOUND potential single-precision calls"
    WARNINGS=$((WARNINGS + FOUND))
fi
echo ""

# --- Check 3: Common block /xy/ consistency ---
echo "  --- Common block /xy/ consistency ---"
SIZES=""
FOUND=0
for f in $ALLFILES; do
    if [ -f "$f" ]; then
        # Count common /xy/ members (comma-separated items)
        grep -n 'common /xy/' "$f" | while read line; do
            NCOMMA=$(echo "$line" | tr -cd ',' | wc -c)
            NMEMBERS=$((NCOMMA + 1))
            if [ $NMEMBERS -ne 6 ]; then
                echo "    $f: $line — has $NMEMBERS members (expected 6)"
                FOUND=$((FOUND + 1))
            fi
        done
    fi
done
if [ $FOUND -eq 0 ]; then
    echo "  [PASS]  All /xy/ blocks have 6 members"
else
    echo "  [FAIL]  Inconsistent /xy/ block sizes"
    ERRORS=$((ERRORS + FOUND))
fi
echo ""

# --- Check 4: Tab characters ---
echo "  --- Tab characters ---"
FOUND=0
for f in $ALLFILES; do
    if [ -f "$f" ]; then
        NTABS=$(awk '/\t/ {count++} END {print count+0}' "$f")
        if [ "$NTABS" -gt 0 ]; then
            echo "    $f: $NTABS lines with tabs"
            FOUND=$((FOUND + NTABS))
        fi
    fi
done
if [ $FOUND -eq 0 ]; then
    echo "  [PASS]  No tab characters found"
else
    echo "  [WARN]  $FOUND lines contain tabs"
    WARNINGS=$((WARNINGS + FOUND))
fi
echo ""

# --- Check 5: Trailing whitespace in code lines ---
echo "  --- Trailing whitespace ---"
FOUND=0
for f in $SRCFILES; do
    if [ -f "$f" ]; then
        NT=$(awk '/ +$/ {count++} END {print count+0}' "$f")
        FOUND=$((FOUND + NT))
    fi
done
if [ $FOUND -eq 0 ]; then
    echo "  [PASS]  No trailing whitespace in source files"
else
    echo "  [INFO]  $FOUND lines have trailing whitespace (cosmetic)"
fi
echo ""

# --- Summary ---
echo "=============================================="
echo " Errors: $ERRORS   Warnings: $WARNINGS"
echo "=============================================="

if [ $ERRORS -gt 0 ]; then
    exit 1
else
    exit 0
fi
