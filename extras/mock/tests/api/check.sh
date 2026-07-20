#!/bin/bash
# Golden-output test for the mocked preCICE library.
#
# Runs one driver scenario against libpreciceMocked (the driver binary is linked
# directly against it) and compares the observable "T> " lines with the
# checked-in expected output. The mock is deterministic, so any diff is a
# behavior change.
#
# usage:  check.sh <driver> <tests-dir> <name> <participant> <cfg-dir> <scenario...>
# Set UPDATE_EXPECTED=1 to (re)generate the expected file instead of comparing.
set -u

DRIVER="$1"
TESTS_DIR="$2"
NAME="$3"
PART="$4"
CFGDIR="$5"
shift 5

CONFIG="$TESTS_DIR/$CFGDIR/precice-config.xml"
EXPECTED="$TESTS_DIR/expected/$NAME.txt"

tmp=$(mktemp) || exit 1
trap 'rm -f "$tmp" "$tmp.filtered"' EXIT

"$DRIVER" "$PART" "$CONFIG" "$@" > "$tmp" 2>&1
rc=$?
if [ $rc -ne 0 ]; then
    echo "driver exited with code $rc; full output:" >&2
    sed 's/^/  /' "$tmp" >&2
    exit 1
fi

# Only the T> lines are the test contract; normalize the config path, which
# appears in some error messages.
grep "^T> " "$tmp" | sed "s|$CONFIG|<CONFIG>|g" > "$tmp.filtered"

if [ "${UPDATE_EXPECTED:-0}" = "1" ]; then
    mkdir -p "$(dirname "$EXPECTED")"
    cp "$tmp.filtered" "$EXPECTED"
    echo "updated $EXPECTED"
    exit 0
fi

if [ ! -f "$EXPECTED" ]; then
    echo "Missing expected output: $EXPECTED" >&2
    echo "Generate it with: UPDATE_EXPECTED=1 $0 $DRIVER $TESTS_DIR $NAME $PART $CFGDIR $*" >&2
    exit 1
fi

if ! diff -u "$EXPECTED" "$tmp.filtered"; then
    echo "Mock behavior differs from expected output ($NAME)" >&2
    exit 1
fi
