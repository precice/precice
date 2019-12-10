#!/bin/bash

set -x
set -e

if [ "$BUILD_TYPE" = "SCONS" ]; then
    # SCons build
    export PRECICE_BUILD_DIR=$TRAVIS_BUILD_DIR/build/last
    cd $TRAVIS_BUILD_DIR
    scons -j $(nproc) petsc=$PETSC mpi=$MPI python=on compiler=$CXX
    cd $TRAVIS_BUILD_DIR/tests
    if [ "$MPI" = "on" ]; then
        mpirun.openmpi -n 4 --output-filename boost-test-output  $PRECICE_BUILD_DIR/testprecice -r detailed
    else
        $PRECICE_BUILD_DIR/testprecice -x -r detailed > boost-test-output
    fi
else
    # CMake build
    export PRECICE_BUILD_DIR=$TRAVIS_BUILD_DIR/build
    mkdir $PRECICE_BUILD_DIR && cd $PRECICE_BUILD_DIR
    cmake -DPETSC=$PETSC -DMPI=$MPI -DPYTHON=on -DBUILD_SHARED_LIBS=ON -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_COMPILER_LAUNCHER=ccache $TRAVIS_BUILD_DIR
    cmake --build . -- -j $(nproc)
    # Run tests with moderate console output
    ctest -LE mpiports --output-on-failure -O $TRAVIS_BUILD_DIR/tests/boost-test-output
    # Run tests which should be muted
    ctest -L mpiports
fi
