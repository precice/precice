#!/bin/bash
# Smoke test for the mocked preCICE library: runs a curated set of short
# tutorials and the course solutions against libpreciceMocked and compares the
# normalized output with the committed baselines in baselines/.
#
# The baselines are generated (and validated against real preCICE) with
# smoke_update.sh; rerun that whenever preCICE, the tutorials, or the course
# material changes intentionally.
#
# Target runtime on the development machine: < 5 minutes.
#
# usage: smoke.sh [--tutorials-only|--courses-only]
set -u
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=smoke_lib.sh
source "$SCRIPT_DIR/smoke_lib.sh"

RUN_TUTORIALS=1
RUN_COURSES=1
CASE_TIMEOUT=90
COURSE_TIMEOUT=120

case "${1:-}" in
    --tutorials-only) RUN_COURSES=0 ;;
    --courses-only) RUN_TUTORIALS=0 ;;
    "") ;;
    *) echo "usage: $0 [--tutorials-only|--courses-only]"; exit 2 ;;
esac

if [ ! -f "$MOCK_LIB" ]; then
    echo "Error: mock library not found at $MOCK_LIB" >&2
    exit 1
fi
if [ ! -d "$BASELINE_DIR" ]; then
    echo "Error: no baselines found in $BASELINE_DIR — run smoke_update.sh first." >&2
    exit 1
fi

smoke_setup_env

PASS=()
DIFF=()
suite_start=$SECONDS

check_against_baseline() { # <label> <baseline file> <actual file>
    local label="$1" baseline="$2" actual="$3"
    if [ ! -f "$baseline" ]; then
        echo "  ✗ $label: missing baseline $baseline"
        DIFF+=("$label (missing baseline)")
        return
    fi
    if diff -u "$baseline" "$actual" > /tmp/smoke_diff.$$; then
        echo "  ✓ $label"
        PASS+=("$label")
    else
        echo "  ✗ $label: output differs from baseline:"
        head -40 /tmp/smoke_diff.$$ | sed 's/^/    /'
        DIFF+=("$label")
    fi
    rm -f /tmp/smoke_diff.$$
}

# The four courses are independent directory trees, so they run as parallel
# background jobs while the tutorial cases run sequentially in the foreground;
# the wall time is roughly max(courses) instead of the sum.
COURSE_PIDS=()
COURSE_OUTS=()
COURSE_NAMES=()
if [ "$RUN_COURSES" -eq 1 ]; then
    for entry in "${SMOKE_COURSES[@]}"; do
        num="${entry%%:*}"; cname="${entry##*:}"
        out=$(mktemp)
        smoke_run_mock_course "$num" "$COURSE_TIMEOUT" > "$out" &
        COURSE_PIDS+=($!)
        COURSE_OUTS+=("$out")
        COURSE_NAMES+=("$cname")
    done
fi

if [ "$RUN_TUTORIALS" -eq 1 ]; then
    echo "== Tutorials (mock) =="
    for group in "${SMOKE_GROUPS[@]}"; do
        IFS='|' read -r slug tut cases <<< "$group"
        for c in $cases; do
            smoke_clean_case "$tut" "$c"
            actual=$(mktemp)
            smoke_run_mock_case "$tut" "$c" "$CASE_TIMEOUT" > "$actual"
            check_against_baseline "$slug/$c" "$BASELINE_DIR/tutorials/$slug--$c.out" "$actual"
            rm -f "$actual"
        done
    done
fi

if [ "$RUN_COURSES" -eq 1 ]; then
    echo "== Courses (mock) =="
    for i in "${!COURSE_PIDS[@]}"; do
        wait "${COURSE_PIDS[$i]}"
        check_against_baseline "course-${COURSE_NAMES[$i]}" \
            "$BASELINE_DIR/courses/${COURSE_NAMES[$i]}.out" "${COURSE_OUTS[$i]}"
        rm -f "${COURSE_OUTS[$i]}"
    done
fi

echo ""
echo "=================================================="
echo "SMOKE SUMMARY  (took $((SECONDS - suite_start))s)"
echo "=================================================="
echo "Matched baselines: ${#PASS[@]}"
echo "Deviations:        ${#DIFF[@]}"
if [ -f "$BASELINE_DIR/real-status.txt" ]; then
    echo ""
    echo "Baselines last validated against real preCICE:"
    sed 's/^/  /' "$BASELINE_DIR/real-status.txt" | head -3
fi
if [ ${#DIFF[@]} -gt 0 ]; then
    echo ""
    echo "DEVIATING CASES:"
    for d in "${DIFF[@]}"; do echo "  - $d"; done
    exit 1
fi
exit 0
