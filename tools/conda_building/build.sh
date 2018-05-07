#!/bin/bash

source activate.sh

# TODO: make parallel jobs variable
SCONS_OPTIONS="petsc=no compiler=mpic++"

# clean
(cd $PRECICE_ROOT; scons --clean $SCONS_OPTIONS) &> scons_clean.log || exit &
# create a spinner indicating progress. See https://stackoverflow.com/a/12498305/5158031
pid=$! # Process Id of the previous running command

spin='-\|/'

i=0
while kill -0 $pid 2>/dev/null
do
  i=$(( (i+1) %4 ))
  printf "\rcleaning last build: ${spin:$i:1}"
  sleep .1
done

printf "\ncleaning done."
echo ""

# building
(cd $PRECICE_ROOT; scons --config=force $SCONS_OPTIONS -j8)  &> scons.log || exit &
# create a spinner indicating progress. See https://stackoverflow.com/a/12498305/5158031
pid=$! # Process Id of the previous running command

spin='-\|/'

i=0
while kill -0 $pid 2>/dev/null
do
  i=$(( (i+1) %4 ))
  printf "\rbuilding precice: ${spin:$i:1}"
  sleep .1
done

printf "\nbuilding done."

# testing
(cd $PRECICE_ROOT; ./tools/compileAndTest.py -b) || exit

# refresh activate
source $CONDA_ROOT/bin/activate precice
