#!/bin/bash
# Differential test: real preCICE (coupled pair) vs the mock (LD_PRELOAD, standalone).
#
# Usage:  ./run.sh [all|lifecycle|errors|emptymesh|diff]
#
# Builds the driver against libprecice, runs each scenario against the real library
# (two coupled processes) and against libpreciceMocked.so, and stores the outputs in
# out/. "diff" compares the "T> " control-flow lines of all recorded runs; lines
# starting with "T> V " carry data values and are expected to differ (the mock does
# not couple or interpolate data).
#
# Environment: PRECICE_BUILD_DIR (default: <repo>/build)
set -u
ROOT="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$ROOT/../../../.." && pwd)"
BUILD=${PRECICE_BUILD_DIR:-$REPO/build}
MOCK=$BUILD/libpreciceMocked.so
TMO=90

if [ ! -f "$MOCK" ] || [ ! -f "$BUILD/libprecice.so" ]; then
  echo "libprecice.so / libpreciceMocked.so not found in $BUILD" >&2
  echo "Build them first or set PRECICE_BUILD_DIR." >&2
  exit 1
fi

build_driver() {
  if [ ! -x "$ROOT/driver" ] || [ "$ROOT/driver.cpp" -nt "$ROOT/driver" ]; then
    echo "building driver..."
    g++ -std=c++17 -O1 "$ROOT/driver.cpp" -I"$REPO/src" -I"$BUILD/src" \
        -L"$BUILD" -lprecice -Wl,-rpath,"$BUILD" -o "$ROOT/driver" || exit 1
  fi
}

run_real_pair() { # <label> <cfgdir> <scenario...>
  local label=$1 cfg=$2; shift 2
  local dir=$ROOT/out/real-$label
  rm -rf "$dir"; mkdir -p "$dir"; cd "$dir" || exit 1
  timeout $TMO "$ROOT/driver" SolverOne "$ROOT/$cfg/precice-config.xml" "$@" >one.out 2>&1 &
  local p1=$!
  sleep 0.2
  timeout $TMO "$ROOT/driver" SolverTwo "$ROOT/$cfg/precice-config.xml" "$@" >two.out 2>&1 &
  local p2=$!
  wait $p1; local r1=$?
  wait $p2; local r2=$?
  echo "[real-$label] SolverOne=$r1 SolverTwo=$r2"
}

run_real_solo() { # <label> <cfgdir> <participant> <scenario...>
  local label=$1 cfg=$2 part=$3; shift 3
  local dir=$ROOT/out/real-$label
  rm -rf "$dir"; mkdir -p "$dir"; cd "$dir" || exit 1
  local out=one.out; [ "$part" = SolverTwo ] && out=two.out
  timeout $TMO "$ROOT/driver" "$part" "$ROOT/$cfg/precice-config.xml" "$@" >"$out" 2>&1
  echo "[real-$label] $part=$?"
}

run_real_victim() { # <label> <cfgdir> <victim-scenario...>  (SolverTwo runs lifecycle)
  local label=$1 cfg=$2; shift 2
  local dir=$ROOT/out/real-$label
  rm -rf "$dir"; mkdir -p "$dir"; cd "$dir" || exit 1
  timeout 30 "$ROOT/driver" SolverOne "$ROOT/$cfg/precice-config.xml" "$@" >one.out 2>&1 &
  local p1=$!
  timeout 30 "$ROOT/driver" SolverTwo "$ROOT/$cfg/precice-config.xml" lifecycle >two.out 2>&1
  wait $p1
  echo "[real-$label] SolverOne=$?"
}

run_mock() { # <label> <cfgdir> <participant> <scenario...>
  local label=$1 cfg=$2 part=$3; shift 3
  local dir=$ROOT/out/mock-$label
  mkdir -p "$dir"; cd "$dir" || exit 1
  local out=one.out; [ "$part" = SolverTwo ] && out=two.out
  LD_PRELOAD=$MOCK timeout $TMO "$ROOT/driver" "$part" "$ROOT/$cfg/precice-config.xml" "$@" >"$out" 2>&1
  echo "[mock-$label] $part=$?"
}

diff_all() {
  local rc=0
  for rdir in "$ROOT"/out/real-*; do
    [ -d "$rdir" ] || continue
    local label=${rdir##*/real-}
    local mdir=$ROOT/out/mock-$label
    [ -d "$mdir" ] || continue
    for f in one two; do
      [ -f "$rdir/$f.out" ] && [ -f "$mdir/$f.out" ] || continue
      if diff <(grep "^T> " "$rdir/$f.out" | grep -v "^T> V ") \
              <(grep "^T> " "$mdir/$f.out" | grep -v "^T> V ") >/dev/null; then
        echo "OK   $label/$f"
      else
        echo "DIFF $label/$f"
        diff <(grep "^T> " "$rdir/$f.out" | grep -v "^T> V ") \
             <(grep "^T> " "$mdir/$f.out" | grep -v "^T> V ") | head -20
        rc=1
      fi
    done
  done
  return $rc
}

build_driver
mkdir -p "$ROOT/out"

case "${1:-all}" in
lifecycle|all)
  run_real_pair explicit explicit lifecycle
  run_mock explicit explicit SolverOne lifecycle
  run_mock explicit explicit SolverTwo lifecycle
  run_real_pair explicit-sub explicit lifecycle-sub
  run_mock explicit-sub explicit SolverOne lifecycle-sub
  run_mock explicit-sub explicit SolverTwo lifecycle-sub
  run_real_pair implicit implicit lifecycle
  run_mock implicit implicit SolverOne lifecycle
  run_mock implicit implicit SolverTwo lifecycle
  run_real_pair implicit-sub implicit lifecycle-sub
  run_mock implicit-sub implicit SolverOne lifecycle-sub
  run_mock implicit-sub implicit SolverTwo lifecycle-sub
  run_real_pair maxtime maxtime lifecycle
  run_mock maxtime maxtime SolverOne lifecycle
  run_mock maxtime maxtime SolverTwo lifecycle
  run_real_pair firstpart firstpart lifecycle-fixed
  run_mock firstpart firstpart SolverOne lifecycle-fixed
  run_mock firstpart firstpart SolverTwo lifecycle-fixed
  ;;&
errors|all)
  run_real_solo err-ctor explicit SolverOne errors ctor-badrank
  run_mock err-ctor explicit SolverOne errors ctor-badrank
  run_real_solo err-preinit-one explicit SolverOne errors preinit
  run_real_solo err-preinit-two explicit SolverTwo errors preinit
  run_mock err-preinit-one explicit SolverOne errors preinit
  run_mock err-preinit-two explicit SolverTwo errors preinit
  run_real_pair err-postinit explicit errors postinit
  run_mock err-postinit explicit SolverOne errors postinit
  run_mock err-postinit explicit SolverTwo errors postinit
  for probe in advzero advneg advinf advbig; do
    run_real_victim err-$probe explicit errors $probe
    run_mock err-$probe explicit SolverOne errors $probe
  done
  run_real_victim err-nocheckpoint implicit errors nocheckpoint
  run_mock err-nocheckpoint implicit SolverOne errors nocheckpoint
  run_real_pair err-afterloop explicit errors afterloop
  run_mock err-afterloop explicit SolverOne errors afterloop
  run_mock err-afterloop explicit SolverTwo errors afterloop
  ;;&
emptymesh|all)
  run_real_victim emptymesh explicit errors emptymesh
  run_mock emptymesh explicit SolverOne errors emptymesh
  ;;&
diff|all)
  echo "==== control-flow diff (real vs mock) ===="
  diff_all
  ;;
esac
echo DONE
