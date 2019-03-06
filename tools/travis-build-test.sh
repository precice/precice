#!/bin/bash

set -x
set -e


echo "LOCAL_INSTALL: $LOCAL_INSTALL"
ls -l $LOCAL_INSTALL
ls -l $LOCAL_INSTALL/include
ls -l $LOCAL_INSTALL/lib
ls -l $LOCAL_INSTALL/eigen3
echo "EIGEN3_ROOT_DIR: $EIGEN3_ROOT_DIR"

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
    mkdir $TRAVIS_BUILD_DIR/build && cd $TRAVIS_BUILD_DIR/build
    cmake -DPETSC=$PETSC -DMPI=$MPI -DPYTHON=on -DBUILD_SHARED_LIBS=ON -DCMAKE_BUILD_TYPE=Debug $TRAVIS_BUILD_DIR
    cmake --build . -- -j 2
    ctest --output-on-failure -O $TRAVIS_BUILD_DIR/tests/boost-test-output
fi
