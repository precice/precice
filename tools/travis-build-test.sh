#!/bin/bash

if [ "$BUILD_TYPE" = "SCONS" ]; then
    export PRECICE_BUILD_DIR=$TRAVIS_BUILD_DIR/build/last
    cd $TRAVIS_BUILD_DIR
    scons -j 2 petsc=$PETSC mpi=$MPI python=on compiler=$CXX staticlib bin solib tests symlink
    cd $TRAVIS_BUILD_DIR/tests
    if [ "$MPI" = "on" ]; then
        mpirun.openmpi -n 4 --output-filename boost-test-output  $PRECICE_BUILD_DIR/testprecice -r detailed
    else
        $PRECICE_BUILD_DIR/testprecice -x -r detailed > boost-test-output
    fi
else
    mkdir build && cd build
    cmake -DPETSC=$PETSC -DMPI=$MPI -DPYTHON=on -DBUILD_SHARED_LIBS=ON -DCMAKE_BUILD_TYPE=Debug ..
    cmake --build . -- -j 2
    ctest -V > $TRAVIS_BUILD_DIR/tests/boost-test-output
fi
