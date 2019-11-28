#!/bin/bash

set -x
set -e

export PYTHON_BINDINGS_DIR=$TRAVIS_BUILD_DIR/src/precice/bindings/python
cd $PYTHON_BINDINGS_DIR

export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:$TRAVIS_BUILD_DIR/src

python3 setup.py test

