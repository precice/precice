#!/bin/bash

set -x
set -e

export PYTHON_BINDINGS_DIR=$TRAVIS_BUILD_DIR/src/precice/bindings/python_future
cd $PYTHON_BINDINGS_DIR

cythonize -i -E TEST=True precice_future.pyx test/test_bindings_module.pyx
python3 -m unittest --verbose

