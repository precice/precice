#!/bin/bash
source ./precice_config.sh

############

export PRECICE_GIT_BRANCH
export PRECICE_ROOT
export ANACONDA_ROOT
export SCONS_PARALLELJOBS
export PRECICE_MPI_IMPLEMENTATION

CONDA_ENV=precice
CONDA_ENV_ROOT=$CONDA_PREFIX

export PKG_CONFIG_PATH=$CONDA_ENV_ROOT/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$CONDA_ENV_ROOT/lib

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PRECICE_ROOT/build/last

export PATH=$PRECICE_ROOT/bin:$ANACONDA_ROOT/bin:$PATH

############

if [ -z $PRECICE_ROOT ]; then
    echo "please define PRECICE_ROOT"
fi

if [ -z $ANACONDA_ROOT ]; then
    echo "please define ANACONDA_ROOT"
fi

source $ANACONDA_ROOT/bin/activate precice
