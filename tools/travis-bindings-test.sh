#!/bin/bash

set -x
set -e

export PYTHON_BINDINGS_DIR=$TRAVIS_BUILD_DIR/src/precice/bindings/python_future
cd $PYTHON_BINDINGS_DIR

export PYTHONPATH=/usr/local/lib/python3/

pip3 install --user setuptools cython
python3 setup.py test

