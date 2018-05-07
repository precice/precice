#!/bin/bash
set -e

# This script activates an already created Conda environment for preCICE.
# It exports the Conda and preCICE paths and calls Conda activate.

# Load the configuration
source config.sh

# Check if PRECICE_ROOT and CONDA_ROOT are defined
# TODO: Check for existence?
if [ -z ${PRECICE_ROOT} ]; then
    printf "\nError: Please define PRECICE_ROOT."
    exit
fi
if [ -z ${CONDA_ROOT} ]; then
    printf "\nError: Please define CONDA_ROOT."
    exit
fi

# TODO: Do we actually need to export these?
export PRECICE_ROOT
export CONDA_ROOT
export PRECICE_MPI_IMPLEMENTATION

# Add Conda to the system paths
# TODO: Why?
echo "TODO: CONDA_ENV_ROOT:"
echo "${CONDA_ENV_ROOT}"
export PKG_CONFIG_PATH="${CONDA_ENV_ROOT}/lib/pkgconfig:${PKG_CONFIG_PATH}"
export LD_LIBRARY_PATH="${CONDA_ENV_ROOT}/lib"
export PATH="${CONDA_ROOT}/bin:${PATH}"

# Add preCICE to the system paths
export LD_LIBRARY_PATH="${PRECICE_ROOT}/build/last:${LD_LIBRARY_PATH}"
export PATH="${PRECICE_ROOT}/bin:${PATH}"

# Activate the precice environment
source ${CONDA_ROOT}/bin/activate precice
