#! /bin/sh
#
# Invoke with
# runsolverdummies.sh path/to/solverdummy1 path/to/solverdummy2 path/to/config

$1 $3 SolverOne MeshOne  &
PID=$!

$2 $3 SolverTwo MeshTwo
RET=$?

set -e
wait $PID
exit $RET
