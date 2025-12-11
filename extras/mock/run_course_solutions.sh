#!/bin/bash


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# repo root is two levels up from extras/mock
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

MOCK_LIB="$REPO_ROOT/build/libpreciceMocked.so.3.3.0"
COURSE_DIR="$REPO_ROOT/examples/course"
TIMEOUT=300
USE_MOCK=0

FAILED=()
PASSED=()
TIMED_OUT=()

print_usage() {
    cat <<EOF
Usage: $0 [--mock] [--course-dir PATH] [--timeout SECONDS] [--help]

Options:
  --mock            Run solutions with mocked precice (LD_PRELOAD)
  --mock-lib PATH   Path to mocked precice library (default: relative to repo)
  --course-number N  Select course to run: 1=B,2=T,3=CHT,4=MAP (default: all)
  --course-dir PATH Path to the course directory (default: $COURSE_DIR)
  --timeout N       Timeout per solution in seconds (default: $TIMEOUT)
  --help            Show this help

This script only runs solution scripts. It searches for files matching:
  - run_solution.sh
  - */solutions/run.sh
  - */solution*/run.sh
If your course uses a different pattern, pass a different course dir or modify the script.
EOF
}

# Parse args
while [[ $# -gt 0 ]]; do
    case "$1" in
        --mock)
            USE_MOCK=1
            shift
            ;;
        --mock-lib)
            MOCK_LIB="$2"
            shift 2
            ;;
        --course-number|-n)
            COURSE_NUMBER_ARG="$2"
            shift 2
            ;;
        --course-dir)
            COURSE_DIR="$2"
            shift 2
            ;;
        --timeout)
            TIMEOUT="$2"
            shift 2
            ;;
        --help|-h)
            print_usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            print_usage
            exit 2
            ;;
    esac
done

if [ "$USE_MOCK" -eq 1 ]; then
    if [ ! -f "$MOCK_LIB" ]; then
        echo "Error: mock library not found at $MOCK_LIB"
        exit 1
    fi
    echo "Running solutions with mock library: $MOCK_LIB"
else
    echo "Running solutions with real precice (no LD_PRELOAD)"
fi

echo "Course directory: $COURSE_DIR"
echo "Timeout per solution: ${TIMEOUT}s"
echo "=================================================="


if [ ! -d "$COURSE_DIR" ]; then
    echo "Error: course directory not found: $COURSE_DIR"
    exit 1
fi

# Allow selecting a specific course by number
# 1 -> B, 2 -> T, 3 -> CHT, 4 -> MAP
COURSE_NUMBER=0
TASK_DIRS=()

# If user passed --course-number earlier, parse it from env (we added parsing below)
if [ -n "${COURSE_NUMBER_ARG:-}" ]; then
    COURSE_NUMBER="$COURSE_NUMBER_ARG"
fi

if [ "$COURSE_NUMBER" -ne 0 ]; then
    case "$COURSE_NUMBER" in
        1) COURSE_NAME="B" ;;
        2) COURSE_NAME="T" ;;
        3) COURSE_NAME="CHT" ;;
        4) COURSE_NAME="MAP" ;;
        *) echo "Unknown course number: $COURSE_NUMBER"; exit 2 ;;
    esac
    # special-case CHT which has solution/all
    if [ "$COURSE_NAME" = "CHT" ]; then
        SOL_ROOT="$COURSE_DIR/$COURSE_NAME/solution/all"
    else
        SOL_ROOT="$COURSE_DIR/$COURSE_NAME/solution"
    fi
    if [ ! -d "$SOL_ROOT" ]; then
        echo "Course solution directory not found: $SOL_ROOT"; exit 1
    fi
    # collect task directories sorted (T1, T2, ... or any subdir)
    mapfile -t TASK_DIRS < <(find "$SOL_ROOT" -mindepth 1 -maxdepth 1 -type d -printf '%p\n' | sort -V)
else
    # No specific course — run all solutions found under examples/course in deterministic order
    mapfile -t COURSE_NAMES < <(find "$COURSE_DIR" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort -V)
    for cname in "${COURSE_NAMES[@]}"; do
        if [ "$cname" = "CHT" ]; then
            SOL_ROOT="$COURSE_DIR/$cname/solution/all"
        else
            SOL_ROOT="$COURSE_DIR/$cname/solution"
        fi
        if [ -d "$SOL_ROOT" ]; then
            mapfile -t tmp < <(find "$SOL_ROOT" -mindepth 1 -maxdepth 1 -type d -printf '%p\n' | sort -V)
            TASK_DIRS+=("${tmp[@]}")
        fi
    done
fi

if [ ${#TASK_DIRS[@]} -eq 0 ]; then
    echo "No task directories found to run under $COURSE_DIR"
    exit 0
fi

echo "Found ${#TASK_DIRS[@]} task directories to run"
echo ""

for dir in "${TASK_DIRS[@]}"; do
    name=$(realpath --relative-to="$COURSE_DIR" "$dir" 2>/dev/null || echo "$dir")
    echo "Running task: $name"

    # run clean.sh before in task, generator, and propagator if present
    run_clean_if_exists() {
        d="$1"
        if [ -x "$d/clean.sh" ]; then
            (cd "$d" && ./clean.sh) || true
        elif [ -f "$d/clean.sh" ]; then
            (cd "$d" && bash clean.sh) || true
        fi
    }

    run_clean_if_exists "$dir"

    # determine generator(s) and propagator(s) locations for this task
    GEN_DIRS=()
    PROP_DIRS=()

    # find generator-like dirs or files: generator, generator_*, root generator.py
    if [ -d "$dir/generator" ]; then
        GEN_DIRS+=("$dir/generator")
    fi
    for g in "$dir"/generator_*; do
        [ -e "$g" ] || continue
        if [ -d "$g" ]; then GEN_DIRS+=("$g"); fi
    done
    if [ -f "$dir/generator.py" ]; then GEN_DIRS+=("$dir"); fi

    # find propagator dirs
    if [ -d "$dir/propagator" ]; then
        PROP_DIRS+=("$dir/propagator")
    fi
    for p in "$dir"/propagator_*; do
        [ -e "$p" ] || continue
        if [ -d "$p" ]; then PROP_DIRS+=("$p"); fi
    done
    if [ -f "$dir/propagator.py" ]; then PROP_DIRS+=("$dir"); fi

    # Temporary outputs arrays
    GEN_OUTS=()
    PROP_OUTS=()
    PIDS=()
    RCODES=()

    if [ ${#GEN_DIRS[@]} -gt 0 ] || [ ${#PROP_DIRS[@]} -gt 0 ]; then
        echo "  Detected ${#GEN_DIRS[@]} generator(s) and ${#PROP_DIRS[@]} propagator(s)"

        for gd in "${GEN_DIRS[@]}"; do run_clean_if_exists "$gd"; done
        for pd in "${PROP_DIRS[@]}"; do run_clean_if_exists "$pd"; done

        if [ "$USE_MOCK" -eq 1 ]; then
            echo "  Running sequentially (mock mode)"
            for gd in "${GEN_DIRS[@]}"; do
                out=$(mktemp /tmp/course_gen_out.XXXXXX)
                GEN_OUTS+=("$out")
                # pick a python file to run
                if [ -f "$gd/generator.py" ]; then
                    timeout "$TIMEOUT" bash -c "cd '$gd' && LD_PRELOAD=$MOCK_LIB python3 generator.py" > "$out" 2>&1 || true
                    rc=$?
                else
                    pyfile=$(ls -1 "$gd"/*.py 2>/dev/null | head -n1 || true)
                    if [ -n "$pyfile" ]; then
                        fname=$(basename "$pyfile")
                        timeout "$TIMEOUT" bash -c "cd '$gd' && LD_PRELOAD=$MOCK_LIB python3 $fname" > "$out" 2>&1 || true
                        rc=$?
                    else rc=0; fi
                fi
                RCODES+=("$rc")
            done
            for pd in "${PROP_DIRS[@]}"; do
                out=$(mktemp /tmp/course_prop_out.XXXXXX)
                PROP_OUTS+=("$out")
                if [ -f "$pd/propagator.py" ]; then
                    timeout "$TIMEOUT" bash -c "cd '$pd' && LD_PRELOAD=$MOCK_LIB python3 propagator.py" > "$out" 2>&1 || true
                    rc=$?
                else
                    pyfile=$(ls -1 "$pd"/*.py 2>/dev/null | head -n1 || true)
                    if [ -n "$pyfile" ]; then
                        fname=$(basename "$pyfile")
                        timeout "$TIMEOUT" bash -c "cd '$pd' && LD_PRELOAD=$MOCK_LIB python3 $fname" > "$out" 2>&1 || true
                        rc=$?
                    else rc=0; fi
                fi
                RCODES+=("$rc")
            done
        else
            echo "  Starting generator(s) and propagator(s) concurrently (separate outputs)"
            for gd in "${GEN_DIRS[@]}"; do
                out=$(mktemp /tmp/course_gen_out.XXXXXX)
                GEN_OUTS+=("$out")
                if [ -f "$gd/generator.py" ]; then
                    timeout "$TIMEOUT" bash -c "cd '$gd' && python3 generator.py" > "$out" 2>&1 & pid=$!
                else
                    pyfile=$(ls -1 "$gd"/*.py 2>/dev/null | head -n1 || true)
                    if [ -n "$pyfile" ]; then
                        fname=$(basename "$pyfile")
                        timeout "$TIMEOUT" bash -c "cd '$gd' && python3 $fname" > "$out" 2>&1 & pid=$!
                    else
                        pid=0
                    fi
                fi
                PIDS+=("$pid")
            done
            for pd in "${PROP_DIRS[@]}"; do
                out=$(mktemp /tmp/course_prop_out.XXXXXX)
                PROP_OUTS+=("$out")
                if [ -f "$pd/propagator.py" ]; then
                    timeout "$TIMEOUT" bash -c "cd '$pd' && python3 propagator.py" > "$out" 2>&1 & pid=$!
                else
                    pyfile=$(ls -1 "$pd"/*.py 2>/dev/null | head -n1 || true)
                    if [ -n "$pyfile" ]; then
                        fname=$(basename "$pyfile")
                        timeout "$TIMEOUT" bash -c "cd '$pd' && python3 $fname" > "$out" 2>&1 & pid=$!
                    else
                        pid=0
                    fi
                fi
                PIDS+=("$pid")
            done

            for pid in "${PIDS[@]}"; do
                if [ "$pid" != "0" ]; then
                    if wait $pid; then
                        rc=0
                    else
                        rc=$?
                    fi
                    RCODES+=("$rc")
                fi
            done
        fi

        echo "  ----- Generator outputs -----"
        for out in "${GEN_OUTS[@]}"; do sed 's/^/    /' "$out" || true; done
        echo "  ----- Propagator outputs -----"
        for out in "${PROP_OUTS[@]}"; do sed 's/^/    /' "$out" || true; done

        for gd in "${GEN_DIRS[@]}"; do run_clean_if_exists "$gd"; done
        for pd in "${PROP_DIRS[@]}"; do run_clean_if_exists "$pd"; done

        any_fail=0; any_to=0
        for r in "${RCODES[@]}"; do
            if [ "$r" -eq 124 ] 2>/dev/null; then any_to=1; fi
            if [ "$r" -ne 0 ] 2>/dev/null; then any_fail=1; fi
        done
        if [ "$any_fail" -eq 0 ]; then
            echo "  ✓ PASSED"
            PASSED+=("$name")
        else
            if [ "$any_to" -eq 1 ]; then
                echo "  ⏱ TIMEOUT (>${TIMEOUT}s)"
                TIMED_OUT+=("$name")
            else
                echo "  ✗ FAILED"
                FAILED+=("$name")
            fi
        fi

        for out in "${GEN_OUTS[@]}" "${PROP_OUTS[@]}"; do rm -f "$out"; done
        echo ""
    else
        echo "  No runnable candidate found in $dir — skipping"
        FAILED+=("$name")
    fi
done
echo "=================================================="
echo "SUMMARY"
echo "=================================================="
echo "Passed: ${#PASSED[@]}"
echo "Failed: ${#FAILED[@]}"
echo "Timeout: ${#TIMED_OUT[@]}"
echo "Total: $((${#PASSED[@]} + ${#FAILED[@]} + ${#TIMED_OUT[@]}))"

if [ ${#FAILED[@]} -gt 0 ]; then
    echo "FAILED SOLUTIONS:";
    for t in "${FAILED[@]}"; do echo "  - $t"; done
    echo ""
fi

if [ ${#TIMED_OUT[@]} -gt 0 ]; then
    echo "TIMED OUT SOLUTIONS (>${TIMEOUT}s):";
    for t in "${TIMED_OUT[@]}"; do echo "  - $t"; done
    echo ""
fi

echo "PASSED SOLUTIONS:";
for t in "${PASSED[@]}"; do echo "  - $t"; done

if [ ${#FAILED[@]} -gt 0 ]; then
    exit 1
else
    exit 0
fi
