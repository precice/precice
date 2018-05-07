#!/bin/bash
source ./config.sh

############

export PRECICE_GIT_BRANCH
export PRECICE_ROOT
export CONDA_ROOT
export SCONS_PARALLELJOBS
export PRECICE_MPI_IMPLEMENTATION

CONDA_ENV=precice
CONDA_ENV_ROOT=$CONDA_PREFIX

export PKG_CONFIG_PATH=$CONDA_ENV_ROOT/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$CONDA_ENV_ROOT/lib

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PRECICE_ROOT/build/last

export PATH=$PRECICE_ROOT/bin:$CONDA_ROOT/bin:$PATH

############

if [ -z $PRECICE_ROOT ]; then
    echo "please define PRECICE_ROOT"
fi

if [ -z $CONDA_ROOT ]; then
    echo "please define CONDA_ROOT"
fi

source $CONDA_ROOT/bin/activate precice
