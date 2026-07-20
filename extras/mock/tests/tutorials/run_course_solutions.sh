#!/bin/bash

# Runs all course solution tasks, either with real precice or with the mocked
# library (LD_PRELOAD).
# run with
# ./run_course_solutions.sh --mock &> mock_course.log

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# repo root is four levels up from extras/mock/tests/tutorials
REPO_ROOT="$(cd "$SCRIPT_DIR/../../../.." && pwd)"

MOCK_LIB="$REPO_ROOT/build/libpreciceMocked.so.3.3.0"
COURSE_DIR="$REPO_ROOT/examples/course"
# Python venv with a working numpy/pyprecice combo (same reasoning as in
# run_all_tutorials.sh): plain `python3` resolves to ~/.local where pyprecice's
# binary is ABI-incompatible with numpy 2.x; putting this venv on PATH makes
# `python3` resolve to it instead.
VENV_DIR="$REPO_ROOT/.venv"
TIMEOUT=300
USE_MOCK=0

FAILED=()
PASSED=()
TIMED_OUT=()
SKIPPED=()

print_usage() {
    cat <<EOF
Usage: $0 [--mock] [--course-number N] [--course-dir PATH] [--timeout SECONDS] [--help]

Options:
  --mock             Run solutions with mocked precice (LD_PRELOAD). Participants
                     run sequentially, since the mock needs no coupling partner.
  --mock-lib PATH    Path to mocked precice library (default: relative to repo)
  --course-number N  Select course to run: 1=B,2=T,3=CHT,4=MAP (default: all)
  --course-dir PATH  Path to the course directory (default: $COURSE_DIR)
  --timeout N        Timeout per solution in seconds (default: $TIMEOUT)
  --help             Show this help

A task is any directory under <course>/solution/ that contains a
precice-config.xml. Participants of a task are:
  - generator*/ and propagator*/ subdirs, or generator.py/propagator.py in the
    task dir itself (run with python3)
  - any other subdir containing a run.sh (e.g. CHT's fluid-openfoam/solid-nutils)
Tasks without any runnable participant (e.g. MAP T5, which is driven manually
with ASTE) are reported as skipped.
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

# Activate the working Python venv so participants that call `python3` use a
# numpy/pyprecice that actually imports (see comment at VENV_DIR).
if [ -x "$VENV_DIR/bin/python3" ]; then
    export PATH="$VENV_DIR/bin:$PATH"
    export VIRTUAL_ENV="$VENV_DIR"
    echo "Python venv: $VENV_DIR"
fi

# Constrain pip (including its isolated build envs) for participants that build
# their own venv (e.g. CHT's solid-nutils): nutils needs numpy<2 and pyprecice's
# sdist must build against the same numpy it runs against.
if [ -f "$SCRIPT_DIR/pip-constraints.txt" ]; then
    export PIP_CONSTRAINT="$SCRIPT_DIR/pip-constraints.txt"
fi

# Headless matplotlib: several scripts end in plt.show(), which blocks until
# the plot window is closed; Agg makes that a no-op.
export MPLBACKEND=Agg

echo "Course directory: $COURSE_DIR"
echo "Timeout per solution: ${TIMEOUT}s"
echo "=================================================="


if [ ! -d "$COURSE_DIR" ]; then
    echo "Error: course directory not found: $COURSE_DIR"
    exit 1
fi

# A task dir is any dir under solution/ containing a precice-config.xml. This
# finds B/T/MAP solution/T* as well as CHT's solution/all (whose participants
# are one level deeper) without special-casing.
collect_tasks() {
    local sol_root="$1" cfg
    [ -d "$sol_root" ] || return 0
    while IFS= read -r cfg; do
        TASK_DIRS+=("$(dirname "$cfg")")
    done < <(find "$sol_root" -maxdepth 2 -name precice-config.xml | sort -V)
}

# Allow selecting a specific course by number
# 1 -> B, 2 -> T, 3 -> CHT, 4 -> MAP
COURSE_NUMBER=0
TASK_DIRS=()

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
    if [ ! -d "$COURSE_DIR/$COURSE_NAME/solution" ]; then
        echo "Course solution directory not found: $COURSE_DIR/$COURSE_NAME/solution"; exit 1
    fi
    collect_tasks "$COURSE_DIR/$COURSE_NAME/solution"
else
    # No specific course — run all solutions found under the course dir in
    # deterministic order
    mapfile -t COURSE_NAMES < <(find "$COURSE_DIR" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort -V)
    for cname in "${COURSE_NAMES[@]}"; do
        collect_tasks "$COURSE_DIR/$cname/solution"
    done
fi

if [ ${#TASK_DIRS[@]} -eq 0 ]; then
    echo "No task directories found to run under $COURSE_DIR"
    exit 0
fi

echo "Found ${#TASK_DIRS[@]} task directories to run"
echo ""

run_clean_if_exists() {
    d="$1"
    if [ -x "$d/clean.sh" ]; then
        (cd "$d" && ./clean.sh) || true
    elif [ -f "$d/clean.sh" ]; then
        (cd "$d" && bash clean.sh) || true
    fi
}

# Echo the command to run inside participant dir $1 (empty if none found):
# run.sh if present, else a python file named after the dir (generator_left/
# holds generator_left.py), else the first *.py in the dir.
resolve_cmd() {
    local d="$1" base pyfile
    if [ -f "$d/run.sh" ]; then
        echo "bash ./run.sh"
        return
    fi
    base=$(basename "$d")
    if [ -f "$d/$base.py" ]; then
        echo "python3 $base.py"
        return
    fi
    pyfile=$(ls -1 "$d"/*.py 2>/dev/null | head -n1 || true)
    if [ -n "$pyfile" ]; then
        echo "python3 $(basename "$pyfile")"
    fi
}

for dir in "${TASK_DIRS[@]}"; do
    name=$(realpath --relative-to="$COURSE_DIR" "$dir" 2>/dev/null || echo "$dir")
    echo "Running task: $name"

    run_clean_if_exists "$dir"

    # Collect participants: dir + command pairs (parallel arrays)
    PART_DIRS=()
    PART_CMDS=()

    for g in "$dir/generator" "$dir"/generator_*; do
        [ -d "$g" ] || continue
        cmd=$(resolve_cmd "$g")
        [ -n "$cmd" ] || continue
        PART_DIRS+=("$g"); PART_CMDS+=("$cmd")
    done
    if [ -f "$dir/generator.py" ]; then
        PART_DIRS+=("$dir"); PART_CMDS+=("python3 generator.py")
    fi
    for p in "$dir/propagator" "$dir"/propagator_*; do
        [ -d "$p" ] || continue
        cmd=$(resolve_cmd "$p")
        [ -n "$cmd" ] || continue
        PART_DIRS+=("$p"); PART_CMDS+=("$cmd")
    done
    if [ -f "$dir/propagator.py" ]; then
        PART_DIRS+=("$dir"); PART_CMDS+=("python3 propagator.py")
    fi
    # run.sh-based participants (e.g. CHT): any other immediate subdir with run.sh
    for s in "$dir"/*/; do
        s="${s%/}"
        [ -d "$s" ] && [ -f "$s/run.sh" ] || continue
        case "$(basename "$s")" in generator*|propagator*) continue ;; esac
        PART_DIRS+=("$s"); PART_CMDS+=("bash ./run.sh")
    done

    if [ ${#PART_DIRS[@]} -eq 0 ]; then
        echo "  ⤼ SKIPPED (no runnable participants found)"
        SKIPPED+=("$name")
        echo ""
        continue
    fi

    echo "  Detected ${#PART_DIRS[@]} participant(s):"
    for i in "${!PART_DIRS[@]}"; do
        prel=$(realpath --relative-to="$dir" "${PART_DIRS[$i]}" 2>/dev/null || echo "${PART_DIRS[$i]}")
        echo "    - $prel: ${PART_CMDS[$i]}"
    done

    for pdir in "${PART_DIRS[@]}"; do
        [ "$pdir" != "$dir" ] && run_clean_if_exists "$pdir"
    done

    OUTS=()
    RCODES=()

    if [ "$USE_MOCK" -eq 1 ]; then
        echo "  Running sequentially (mock mode)"
        for i in "${!PART_DIRS[@]}"; do
            out=$(mktemp /tmp/course_part_out.XXXXXX)
            OUTS+=("$out")
            timeout "$TIMEOUT" bash -c "cd '${PART_DIRS[$i]}' && LD_PRELOAD='$MOCK_LIB' ${PART_CMDS[$i]}" > "$out" 2>&1
            RCODES+=($?)
        done
    else
        echo "  Starting participants concurrently (separate outputs)"
        PIDS=()
        for i in "${!PART_DIRS[@]}"; do
            out=$(mktemp /tmp/course_part_out.XXXXXX)
            OUTS+=("$out")
            timeout "$TIMEOUT" bash -c "cd '${PART_DIRS[$i]}' && ${PART_CMDS[$i]}" > "$out" 2>&1 &
            PIDS+=($!)
        done
        for pid in "${PIDS[@]}"; do
            wait "$pid"
            RCODES+=($?)
        done
    fi

    for i in "${!PART_DIRS[@]}"; do
        prel=$(realpath --relative-to="$dir" "${PART_DIRS[$i]}" 2>/dev/null || echo "${PART_DIRS[$i]}")
        echo "  ----- output: $prel (rc=${RCODES[$i]}) -----"
        sed 's/^/    /' "${OUTS[$i]}" || true
    done

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

    for out in "${OUTS[@]}"; do rm -f "$out"; done
    echo ""
done
echo "=================================================="
echo "SUMMARY"
echo "=================================================="
echo "Passed: ${#PASSED[@]}"
echo "Failed: ${#FAILED[@]}"
echo "Timeout: ${#TIMED_OUT[@]}"
echo "Skipped: ${#SKIPPED[@]}"
echo "Total: $((${#PASSED[@]} + ${#FAILED[@]} + ${#TIMED_OUT[@]} + ${#SKIPPED[@]}))"

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

if [ ${#SKIPPED[@]} -gt 0 ]; then
    echo "SKIPPED SOLUTIONS (no runnable participants):";
    for t in "${SKIPPED[@]}"; do echo "  - $t"; done
    echo ""
fi

echo "PASSED SOLUTIONS:";
for t in "${PASSED[@]}"; do echo "  - $t"; done

if [ ${#FAILED[@]} -gt 0 ]; then
    exit 1
else
    exit 0
fi
