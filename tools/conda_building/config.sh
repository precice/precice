# Path to the Conda installation (top-level, which includes bin/)
CONDA_ROOT=''

# Path to the preCICE sources (top-level, which includes src/)
PRECICE_ROOT=''

# Name and prefix of the Conda environment
CONDA_ENV="precice"
CONDA_ENV_ROOT=${CONDA_PREFIX}

# Option for building with SCons
SCONS_OPTIONS="petsc=no python=yes mpi=yes compiler=mpic++ -j 2"

# MPI implementation, used in the Python API ("openmpi" or "mpich")
# TODO: Update also in precice.yml (maybe make this configurable?)
PRECICE_MPI_IMPLEMENTATION="openmpi"
