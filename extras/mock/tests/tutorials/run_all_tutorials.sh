#!/bin/bash

# SET CORRECT PATHS BEFORE RUNNING
# Usage: ./run_all_tutorials.sh [--mock-lib PATH] [--tutorials-dir PATH] [--timeout SECONDS]
# recommended to run from the extras/mock/ directory with the command:
# ./run_all_tutorials.sh &> mock_all_tutorials.log

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../../.." && pwd)"

MOCK_LIB="$REPO_ROOT/build/libpreciceMocked.so.3.3.0"
TUTORIALS_DIR="$REPO_ROOT/examples/tutorials"
TIMEOUT=60  # 1 minute timeout per tutorial
USE_MOCK=1
MOCK_CONFIG=""
CLEAN=1        # run each tutorial's clean.sh before running it
SKIP_MISSING=0 # skip tutorials whose solver/tooling is not installed
SKIP_LONG=0    # skip tutorials known to exceed any reasonable timeout

# Cases that run correctly against the mock but simulate so many time windows
# (or have such slow solver startup) that they exceed any reasonable timeout.
# flow-around-controlled-moving-cylinder: 12000 windows of 0.00025 s
# turek-hron-fsi3: 15000 windows of 0.001 s (~2h)
LONG_TUTORIALS=(
    "channel-transport/fluid-nutils"
    "channel-transport/transport-nutils"
    "flow-around-controlled-moving-cylinder/fluid-openfoam"
    "heat-exchanger/fluid-outer-openfoam"
    "perpendicular-flap/fluid-nutils"
    "perpendicular-flap/solid-nutils"
    "turek-hron-fsi3/fluid-openfoam"
    "two-scale-heat-conduction/micro-nutils"
)
# Python venv with a working numpy/pyprecice combo. Plain `python3` in the
# tutorials picks up ~/.local where pyprecice's binary is ABI-incompatible with
# numpy 2.x; putting this venv on PATH makes `python3` resolve to it instead.
VENV_DIR="$REPO_ROOT/.venv"

print_usage() {
    cat <<EOF
Usage: $0 [--mock-lib PATH] [--tutorials-dir PATH] [--timeout SECONDS] [--mock-config PATH] [--help]

Options:
  --mock-lib PATH      Path to mocked precice library (default: $MOCK_LIB)
  --tutorials-dir      Path to tutorials dir (default: $TUTORIALS_DIR)
  --timeout N          Timeout per tutorial in seconds (default: $TIMEOUT)
  --mock-config PATH   Path to mock config file to copy to each tutorial
  --skip-missing       Skip tutorials whose solver/tooling is not installed
  --skip-long          Skip tutorials known to run far longer than any timeout
  --no-clean           Do not run each tutorial's clean.sh before running it
  --help, -h           Show this help
EOF
}

# Parse args
while [[ $# -gt 0 ]]; do
    case "$1" in
        --mock-lib)
            MOCK_LIB="$2"; shift 2 ;;
        --tutorials-dir)
            TUTORIALS_DIR="$2"; shift 2 ;;
        --timeout)
            TIMEOUT="$2"; shift 2 ;;
        --mock-config)
            MOCK_CONFIG="$2"; shift 2 ;;
        --skip-missing)
            SKIP_MISSING=1; shift ;;
        --skip-long)
            SKIP_LONG=1; shift ;;
        --no-clean)
            CLEAN=0; shift ;;
        --help|-h)
            print_usage; exit 0 ;;
        *) echo "Unknown option: $1"; print_usage; exit 2 ;;
    esac
done

FAILED_TUTORIALS=()
PASSED_TUTORIALS=()
TIMEOUT_TUTORIALS=()
SKIPPED_TUTORIALS=()

# Print a reason and return 0 if the tutorial's solver/tooling is not available
# on this machine (used by --skip-missing).
missing_requirement() {
    local dir="$1" name app
    name=$(basename "$dir")

    # OpenFOAM cases: the required solver is the application in controlDict.
    # This also covers custom solvers (heatTransfer) and solids4Foam.
    if [ -f "$dir/system/controlDict" ]; then
        app=$(sed -n 's/^[[:space:]]*application[[:space:]]*\([A-Za-z0-9_]*\);.*/\1/p' "$dir/system/controlDict" | head -1)
        if [ -n "$app" ] && ! command -v "$app" > /dev/null 2>&1; then
            echo "OpenFOAM solver '$app' not found"
            return 0
        fi
    fi

    case "$name" in
        *calculix*)
            command -v ccx_preCICE > /dev/null 2>&1 || { echo "ccx_preCICE not found"; return 0; } ;;
        *fenics*)
            python3 -c 'import fenicsprecice' > /dev/null 2>&1 || { echo "python module fenicsprecice not importable"; return 0; } ;;
        *su2*)
            command -v SU2_preCICE_FSI.py > /dev/null 2>&1 || { echo "SU2_preCICE_FSI.py not found"; return 0; } ;;
        *codeaster*)
            command -v as_run > /dev/null 2>&1 || { echo "code_aster (as_run) not found"; return 0; } ;;
        *dune*)
            find "$dir" -maxdepth 1 -name 'dune-*' -type f -executable 2> /dev/null | grep -q . \
                || { echo "DUNE adapter binary not built"; return 0; } ;;
        *dumux*)
            if ! find "$dir" -maxdepth 1 -name '*dumux*' -type f -executable 2> /dev/null | grep -q . \
               && ! python3 -c 'import dumux' > /dev/null 2>&1; then
                echo "DuMuX not available"
                return 0
            fi ;;
        aste-turbine)
            command -v precice-aste-evaluate > /dev/null 2>&1 || { echo "ASTE (precice-aste-evaluate) not found"; return 0; } ;;
    esac
    return 1
}

# Verify mock library exists
if [ ! -f "$MOCK_LIB" ]; then
    echo "Error: Mock library not found at $MOCK_LIB"
    exit 1
fi

# Verify mock config exists if provided and convert to absolute path
if [ -n "$MOCK_CONFIG" ]; then
    if [ ! -f "$MOCK_CONFIG" ]; then
        echo "Error: Mock config file not found at $MOCK_CONFIG"
        exit 1
    fi
    # Convert to absolute path
    MOCK_CONFIG="$(cd "$(dirname "$MOCK_CONFIG")" && pwd)/$(basename "$MOCK_CONFIG")"
fi

# Activate the working Python venv so tutorials that call `python3` directly use
# a numpy/pyprecice that actually imports. Tutorials that create their own venv
# (e.g. nutils) are unaffected. FEniCS tutorials still need system dolfin and
# may fail independently of this.
if [ -x "$VENV_DIR/bin/python3" ]; then
    export PATH="$VENV_DIR/bin:$PATH"
    export VIRTUAL_ENV="$VENV_DIR"
fi

# Make pkg-config find the locally-built libprecice (3.3.0). The Rust `precice`
# crate constrains `libprecice >= 3.3, < 3.4`; a system-wide 3.2.x .pc does not
# satisfy it, so cargo builds fail before reaching the mock. (Note: the Rust
# tutorial sources may still need updating for the 3.3 Result-based binding API.)
PRECICE_PC_DIR="$REPO_ROOT/build/lib/pkgconfig"
if [ -f "$PRECICE_PC_DIR/libprecice.pc" ]; then
    export PKG_CONFIG_PATH="$PRECICE_PC_DIR${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
fi

# Constrain pip (including its isolated build envs) for tutorials that build
# their own venv. nutils still relies on numpy<2 APIs, and pyprecice's sdist must
# build against the same numpy it runs against. A requirements.txt pin alone does
# not reach the build env; PIP_CONSTRAINT does.
if [ -f "$SCRIPT_DIR/pip-constraints.txt" ]; then
    export PIP_CONSTRAINT="$SCRIPT_DIR/pip-constraints.txt"
fi

# Headless matplotlib: several scripts end in plt.show(), which blocks until
# the plot window is closed; Agg makes that a no-op.
export MPLBACKEND=Agg

echo "Running all tutorials with mock library..."
echo "Mock library: $MOCK_LIB"
if [ -x "$VENV_DIR/bin/python3" ]; then
    echo "Python venv: $VENV_DIR"
fi
if [ -n "$MOCK_CONFIG" ]; then
    echo "Mock config: $MOCK_CONFIG"
fi
echo "=================================================="
echo ""

# Find all run.sh files and run them
while IFS= read -r run_script; do
    tutorial_dir=$(dirname "$run_script")
    tutorial_name=$(echo "$tutorial_dir" | sed "s|$TUTORIALS_DIR/||")

    echo "Running: $tutorial_name"

    if [ "$SKIP_MISSING" -eq 1 ] && reason=$(missing_requirement "$tutorial_dir"); then
        echo "  ⤼ SKIPPED ($reason)"
        echo ""
        SKIPPED_TUTORIALS+=("$tutorial_name")
        continue
    fi

    if [ "$SKIP_LONG" -eq 1 ]; then
        skip_this=0
        for long_name in "${LONG_TUTORIALS[@]}"; do
            case "$tutorial_dir" in
            */"$long_name")
                skip_this=1
                break
                ;;
            esac
        done
        if [ "$skip_this" -eq 1 ]; then
            echo "  ⤼ SKIPPED (known long-running)"
            echo ""
            SKIPPED_TUTORIALS+=("$tutorial_name")
            continue
        fi
    fi

    # Change to tutorial directory and run with mock
    cd "$tutorial_dir"

    # Clean leftover state from previous (possibly killed) runs: stale result
    # dirs otherwise poison later runs, e.g. via OpenFOAM's startFrom latestTime.
    if [ "$CLEAN" -eq 1 ] && [ -x ./clean.sh ]; then
        ./clean.sh > /dev/null 2>&1 || echo "  (clean.sh failed, continuing)"
    fi

    # Find precice-config.xml location and copy mock config there if provided
    precice_config_dir=""
    if [ -n "$MOCK_CONFIG" ]; then
        # Look for precice-config.xml in current and parent directories
        if [ -f "./precice-config.xml" ]; then
            precice_config_dir="."
        elif [ -f "../precice-config.xml" ]; then
            precice_config_dir=".."
        fi

        if [ -n "$precice_config_dir" ]; then
            cp "$MOCK_CONFIG" "$precice_config_dir/precice-mock-config.xml"
        fi
    fi

    # Run with timeout and capture both stdout and stderr
    start_time=$SECONDS
    timeout $TIMEOUT bash -c "LD_PRELOAD=$MOCK_LIB ./run.sh" > /tmp/tutorial_output.log 2>&1
    exit_code=$?
    elapsed=$((SECONDS - start_time))

    # Remove mock config if it was copied
    if [ -n "$precice_config_dir" ] && [ -f "$precice_config_dir/precice-mock-config.xml" ]; then
        rm "$precice_config_dir/precice-mock-config.xml"
    fi

    if [ $exit_code -eq 0 ]; then
        echo "  ✓ PASSED (${elapsed}s)"
        PASSED_TUTORIALS+=("$tutorial_name (${elapsed}s)")
    elif [ $exit_code -eq 124 ]; then
        echo "  ⏱ TIMEOUT (>$TIMEOUT seconds)"
        TIMEOUT_TUTORIALS+=("$tutorial_name")
        # Print last part of output for timeout cases
        echo "  Last output:"
        tail -20 /tmp/tutorial_output.log | sed 's/^/    /'
    else
        echo "  ✗ FAILED (exit code: $exit_code)"
        FAILED_TUTORIALS+=("$tutorial_name")
        # Print errors
        echo "  Error output:"
        tail -30 /tmp/tutorial_output.log | sed 's/^/    /'
    fi
    echo ""

done < <(find "$TUTORIALS_DIR" -name "run.sh" -type f | sort)

# Print summary
echo "=================================================="
echo "SUMMARY"
echo "=================================================="
echo "Passed: ${#PASSED_TUTORIALS[@]}"
echo "Failed: ${#FAILED_TUTORIALS[@]}"
echo "Timeout: ${#TIMEOUT_TUTORIALS[@]}"
echo "Skipped: ${#SKIPPED_TUTORIALS[@]}"
echo "Total: $((${#PASSED_TUTORIALS[@]} + ${#FAILED_TUTORIALS[@]} + ${#TIMEOUT_TUTORIALS[@]} + ${#SKIPPED_TUTORIALS[@]}))"
echo ""

if [ ${#SKIPPED_TUTORIALS[@]} -gt 0 ]; then
    echo "SKIPPED TUTORIALS (missing solver/tooling):"
    for tutorial in "${SKIPPED_TUTORIALS[@]}"; do
        echo "  - $tutorial"
    done
    echo ""
fi

if [ ${#FAILED_TUTORIALS[@]} -gt 0 ]; then
    echo "FAILED TUTORIALS:"
    for tutorial in "${FAILED_TUTORIALS[@]}"; do
        echo "  - $tutorial"
    done
    echo ""
fi

if [ ${#TIMEOUT_TUTORIALS[@]} -gt 0 ]; then
    echo "TIMEOUT TUTORIALS (>$TIMEOUT seconds):"
    for tutorial in "${TIMEOUT_TUTORIALS[@]}"; do
        echo "  - $tutorial"
    done
    echo ""
fi

echo "PASSED TUTORIALS:"
for tutorial in "${PASSED_TUTORIALS[@]}"; do
    echo "  - $tutorial"
done

# Return exit code based on failures
if [ ${#FAILED_TUTORIALS[@]} -gt 0 ]; then
    exit 1
else
    exit 0
fi
