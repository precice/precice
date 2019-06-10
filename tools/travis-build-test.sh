#!/bin/bash

set -x
set -e

if [ "$BUILD_TYPE" = "SCONS" ]; then
    # SCons build
    export PRECICE_BUILD_DIR=$TRAVIS_BUILD_DIR/build/last
    cd $TRAVIS_BUILD_DIR
    scons -j $(nproc) petsc=$PETSC mpi=$MPI python=on compiler=$CXX staticlib bin solib tests symlink
    cd $TRAVIS_BUILD_DIR/tests
    if [ "$MPI" = "on" ]; then
        mpirun.openmpi -n 4 --output-filename boost-test-output  $PRECICE_BUILD_DIR/testprecice -r detailed
    else
        $PRECICE_BUILD_DIR/testprecice -x -r detailed > boost-test-output
    fi
else
    # CMake build
    mkdir $TRAVIS_BUILD_DIR/build && cd $TRAVIS_BUILD_DIR/build
    cmake -DPETSC=$PETSC -DMPI=$MPI -DPYTHON=on -DBUILD_SHARED_LIBS=ON -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_COMPILER_LAUNCHER=ccache $TRAVIS_BUILD_DIR
    cmake --build . -- -j $(nproc)
    ctest --output-on-failure -O $TRAVIS_BUILD_DIR/tests/boost-test-output
fi

export FUTURE_PYTHON_BINDINGS_DIR=$TRAVIS_BUILD_DIR/src/precice/bindings/python_future
cd $FUTURE_PYTHON_BINDINGS_DIR

# test bindings
pip3 install --user setuptools cython numpy
python3 setup.py test

# install bindings
python3 setup.py config -I${TRAVIS_BUILD_DIR}/src -L${TRAVIS_BUILD_DIR}/build -lprecice
python3 setup.py install --user

export PYTHON_BINDINGS_DIR=$TRAVIS_BUILD_DIR/src/precice/bindings/python
cd $PYTHON_BINDINGS_DIR

# install bindings
pip3 install --user .
