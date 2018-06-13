#!/bin/bash
set -e

# This script uses SCons to build preCICE in an already set-up environment.

# Clean any previous log files
rm -f scons_clean.log scons.log

# Clean any previous build
(cd ${PRECICE_ROOT}; scons --clean ${SCONS_OPTIONS}) &> scons_clean.log 2>&1 || exit

# Build preCICE using SCons

(cd ${PRECICE_ROOT}; scons --config=force ${SCONS_OPTIONS})  &> scons.log 2>&1 || exit

# Create a spinner indicating progress.
# See https://stackoverflow.com/a/12498305/5158031
pid=$! # Process Id of the previous running command
spin='-\|/'
i=0
while kill -0 ${pid} 2>/dev/null
do
  i=$(( (i+1) %4 ))
  printf "\rBuilding preCICE (see scons.log): ${spin:$i:1}"
  sleep .1
done
printf "\nBuilding done."

# Run the tests
printf "\nRunning the tests:"
(cd ${PRECICE_ROOT}; ./tools/compileAndTest.py -t) || exit

# Refresh activate
# TODO: Why?
source ${CONDA_ROOT}/bin/activate precice
