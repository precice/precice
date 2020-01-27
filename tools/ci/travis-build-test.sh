#!/bin/bash

set -x
set -e

# CMake build
export PRECICE_BUILD_DIR=$TRAVIS_BUILD_DIR/build
mkdir $PRECICE_BUILD_DIR && cd $PRECICE_BUILD_DIR
cmake -DPRECICE_PETScMapping=$PETSC -DPRECICE_MPICommunication=$MPI -DPRECICE_PythonActions=on -DBUILD_SHARED_LIBS=ON -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_COMPILER_LAUNCHER=ccache $TRAVIS_BUILD_DIR
cmake --build . -- -j $(nproc)
# Run tests with moderate console output
ctest -LE mpiports --output-on-failure -O $TRAVIS_BUILD_DIR/tests/boost-test-output
# Run tests which should be muted
ctest -L mpiports
