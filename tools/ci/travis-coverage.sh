#!/bin/bash

set -x
set -e

# CMake build
mkdir $TRAVIS_BUILD_DIR/build && cd $TRAVIS_BUILD_DIR/build
CXX=g++
CXXFLAGS="--coverage"
LDFLAGS="--coverage"
cmake -DPRECICE_MPICommunication=on -DPRECICE_PETScMapping=on -DPRECICE_PythonActions=on -DBUILD_SHARED_LIBS=ON -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_CXX_COMPILER_LAUNCHER=ccache $TRAVIS_BUILD_DIR
cmake --build . -- -j $(nproc)
ctest
