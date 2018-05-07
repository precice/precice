#!/bin/bash
set -e

# This script creates a Conda environment named precice.

source config.sh

# Create the Conda environment from the precice.yml file.
# This will overwrite any previous Conda environments named precice.
# TODO: Rename precice.yml to precice_conda.yml to not confuse with simulations
conda env create --force -f ${PRECICE_ROOT}/tools/conda_building/precice.yml
