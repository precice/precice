#!/bin/bash

# SET CORRECT PATHS BEFORE RUNNING
# recommended to run from the extras/mock/ directory with the command:
# ./run_all_tutorials.sh &> mock_all_tutorials.log
MOCK_LIB="/home/julian/Documents/precice/build/libpreciceMocked.so.3.3.0"
TUTORIALS_DIR="/home/julian/Documents/precice/examples/tutorials"
TIMEOUT=300  # 5 minutes timeout per tutorial
FAILED_TUTORIALS=()
PASSED_TUTORIALS=()
TIMEOUT_TUTORIALS=()

# Verify mock library exists
if [ ! -f "$MOCK_LIB" ]; then
    echo "Error: Mock library not found at $MOCK_LIB"
    exit 1
fi

echo "Running all tutorials with mock library..."
echo "Mock library: $MOCK_LIB"
echo "=================================================="
echo ""

# Find all run.sh files and run them
while IFS= read -r run_script; do
    tutorial_dir=$(dirname "$run_script")
    tutorial_name=$(echo "$tutorial_dir" | sed "s|$TUTORIALS_DIR/||")

    echo "Running: $tutorial_name"

    # Change to tutorial directory and run with mock
    cd "$tutorial_dir"

    # Run with timeout and capture both stdout and stderr
    timeout $TIMEOUT bash -c "LD_PRELOAD=$MOCK_LIB ./run.sh" > /tmp/tutorial_output.log 2>&1
    exit_code=$?

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
