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
    ctest --output-on-failure -O $TRAVIS_BUILD_DIR/tests/boost-test-output
fi

export FUTURE_PYTHON_BINDINGS_DIR=$TRAVIS_BUILD_DIR/src/precice/bindings/python_future
cd $FUTURE_PYTHON_BINDINGS_DIR

# test bindings
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:$TRAVIS_BUILD_DIR/src
python3 setup.py test

# install bindings
python3 setup.py build_ext --include-dirs=$TRAVIS_BUILD_DIR/src --library-dirs=$PRECICE_BUILD_DIR
python3 setup.py install --user

export PYTHON_BINDINGS_DIR=$TRAVIS_BUILD_DIR/src/precice/bindings/python
cd $PYTHON_BINDINGS_DIR

# install bindings
pip3 install --user .
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PRECICE_BUILD_DIR
python3 -c "import precice_future"
