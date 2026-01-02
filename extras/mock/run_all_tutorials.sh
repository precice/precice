#!/bin/bash

# SET CORRECT PATHS BEFORE RUNNING
# Usage: ./run_all_tutorials.sh [--mock-lib PATH] [--tutorials-dir PATH] [--timeout SECONDS]
# recommended to run from the extras/mock/ directory with the command:
# ./run_all_tutorials.sh &> mock_all_tutorials.log

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

MOCK_LIB="$REPO_ROOT/build/libpreciceMocked.so.3.3.0"
TUTORIALS_DIR="$REPO_ROOT/examples/tutorials"
TIMEOUT=300  # 5 minute timeout per tutorial
USE_MOCK=1
MOCK_CONFIG=""

print_usage() {
    cat <<EOF
Usage: $0 [--mock-lib PATH] [--tutorials-dir PATH] [--timeout SECONDS] [--mock-config PATH] [--help]

Options:
  --mock-lib PATH      Path to mocked precice library (default: $MOCK_LIB)
  --tutorials-dir      Path to tutorials dir (default: $TUTORIALS_DIR)
  --timeout N          Timeout per tutorial in seconds (default: $TIMEOUT)
  --mock-config PATH   Path to mock config file to copy to each tutorial
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
        --help|-h)
            print_usage; exit 0 ;;
        *) echo "Unknown option: $1"; print_usage; exit 2 ;;
    esac
done

FAILED_TUTORIALS=()
PASSED_TUTORIALS=()
TIMEOUT_TUTORIALS=()

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

echo "Running all tutorials with mock library..."
echo "Mock library: $MOCK_LIB"
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

    # Change to tutorial directory and run with mock
    cd "$tutorial_dir"

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
    timeout $TIMEOUT bash -c "LD_PRELOAD=$MOCK_LIB ./run.sh" > /tmp/tutorial_output.log 2>&1
    exit_code=$?

    # Remove mock config if it was copied
    if [ -n "$precice_config_dir" ] && [ -f "$precice_config_dir/precice-mock-config.xml" ]; then
        rm "$precice_config_dir/precice-mock-config.xml"
    fi

    if [ $exit_code -eq 0 ]; then
        echo "  ✓ PASSED"
        PASSED_TUTORIALS+=("$tutorial_name")
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
echo "Total: $((${#PASSED_TUTORIALS[@]} + ${#FAILED_TUTORIALS[@]} + ${#TIMEOUT_TUTORIALS[@]}))"
echo ""

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
