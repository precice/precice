#!/bin/bash
# Regenerates the smoke-test baselines. Run this whenever preCICE, the
# tutorials, or the course material changes intentionally.
#
# Two phases:
#   1. real validation — each tutorial group is run as an actual coupled pair
#      against the real libprecice, and each course with concurrent
#      participants; the outcomes are recorded in baselines/real-status.txt.
#      This proves the cases (still) work outside the mock. Skip with
#      --skip-real when only the mock changed.
#   2. mock baselines — the same cases run against libpreciceMocked; their
#      normalized outputs become the new baselines that smoke.sh compares to.
#
# usage: smoke_update.sh [--skip-real]
set -u
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=smoke_lib.sh
source "$SCRIPT_DIR/smoke_lib.sh"

SKIP_REAL=0
REAL_TIMEOUT=600
CASE_TIMEOUT=90
COURSE_TIMEOUT=120

case "${1:-}" in
    --skip-real) SKIP_REAL=1 ;;
    "") ;;
    *) echo "usage: $0 [--skip-real]"; exit 2 ;;
esac

if [ ! -f "$MOCK_LIB" ]; then
    echo "Error: mock library not found at $MOCK_LIB" >&2
    exit 1
fi

smoke_setup_env
mkdir -p "$BASELINE_DIR/tutorials" "$BASELINE_DIR/courses"

# ---------------------------------------------------------------------------
if [ "$SKIP_REAL" -eq 0 ]; then
    echo "== Phase 1: validating cases against real preCICE =="
    {
        echo "validated: $(date -u +%Y-%m-%dT%H:%MZ) on $(hostname)"
        echo "precice: $("$REPO_ROOT/build/precice-version" 2>/dev/null || echo unknown)"
    } > "$BASELINE_DIR/real-status.txt"

    for group in "${SMOKE_GROUPS[@]}"; do
        IFS='|' read -r slug tut cases <<< "$group"
        echo "-- $slug ($tut: $cases)"
        for c in $cases; do smoke_clean_case "$tut" "$c"; done
        pids=()
        rcs=()
        declare -A pid_case=()
        for c in $cases; do
            timeout "$REAL_TIMEOUT" bash -c "cd '$TUTORIALS_DIR/$tut/$c' && ./run.sh" \
                > "/tmp/smoke_real_$c.$$.log" 2>&1 &
            pids+=($!)
            pid_case[$!]="$c"
        done
        status="pass"
        for pid in "${pids[@]}"; do
            wait "$pid"; rc=$?
            rcs+=("$rc")
            c=${pid_case[$pid]}
            if [ "$rc" -eq 124 ]; then
                status="TIMEOUT(>${REAL_TIMEOUT}s:$c)"
            elif [ "$rc" -ne 0 ]; then
                status="FAIL(rc=$rc:$c)"
                echo "   participant $c exited with $rc; last output:"
                tail -8 "/tmp/smoke_real_$c.$$.log" | sed 's/^/     /'
            fi
            rm -f "/tmp/smoke_real_$c.$$.log"
        done
        echo "tutorial $slug: $status" >> "$BASELINE_DIR/real-status.txt"
        echo "   -> $status"
    done

    for entry in "${SMOKE_COURSES[@]}"; do
        num="${entry%%:*}"; cname="${entry##*:}"
        echo "-- course $cname (real, concurrent participants)"
        raw=$(mktemp)
        bash "$SCRIPT_DIR/run_course_solutions.sh" --course-dir "$COURSE_DIR" \
            --course-number "$num" --timeout "$REAL_TIMEOUT" > "$raw" 2>&1
        summary=$(grep -E "^(Passed|Failed|Timeout|Skipped):" "$raw" | tr '\n' ' ')
        echo "course $cname: $summary" >> "$BASELINE_DIR/real-status.txt"
        echo "   -> $summary"
        rm -f "$raw"
    done
    echo ""
else
    echo "== Phase 1 skipped (--skip-real): real-status.txt left unchanged =="
    echo ""
fi

# ---------------------------------------------------------------------------
echo "== Phase 2: regenerating mock baselines =="
# Courses run as parallel background jobs (independent trees) while the
# tutorial cases run sequentially, mirroring smoke.sh.
course_pids=()
for entry in "${SMOKE_COURSES[@]}"; do
    num="${entry%%:*}"; cname="${entry##*:}"
    smoke_run_mock_course "$num" "$COURSE_TIMEOUT" > "$BASELINE_DIR/courses/$cname.out" &
    course_pids+=($!)
done

for group in "${SMOKE_GROUPS[@]}"; do
    IFS='|' read -r slug tut cases <<< "$group"
    for c in $cases; do
        smoke_clean_case "$tut" "$c"
        out="$BASELINE_DIR/tutorials/$slug--$c.out"
        smoke_run_mock_case "$tut" "$c" "$CASE_TIMEOUT" > "$out"
        echo "  wrote $(basename "$out") ($(head -1 "$out"))"
    done
done

for pid in "${course_pids[@]}"; do wait "$pid"; done
for entry in "${SMOKE_COURSES[@]}"; do
    cname="${entry##*:}"
    echo "  wrote courses/$cname.out ($(head -1 "$BASELINE_DIR/courses/$cname.out"))"
done

echo ""
echo "Baselines updated. Review the diff (git diff $BASELINE_DIR) before committing."
