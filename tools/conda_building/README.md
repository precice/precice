# Setup and run preCICE using Conda

You can use the [Conda](https://conda.io/) package manager to Build precice, including all the required dependencies. This will create an environment which you have to activate before you are able to use preCICE. Helper scripts for building and using preCICE with Conda are provided in this folder.

## Build preCICE using Conda

1. Get [Conda](https://conda.io/), e.g. by installing the (minimal) [Miniconda](https://conda.io/miniconda.html) or the (full) [Anaconda](https://www.anaconda.com/) Python distribution.
2. Go to `tools/conda_building` (this folder).
3. Define `CONDA_ROOT=<path/to/conda/installation>` and `PRECICE_ROOT=<path/to/precice/folder>` in `config.sh`.
4. Run `./initialize.sh` to create your Conda environment, named `precice`. This will install the required dependencies and set the necessary environment variables.
5. Run `source activate.sh` to activate the environment `precice`.
6. Run `./build.sh` to (clean-up the potentially existing preCICE build) and build preCICE. Output of the cleaning and building process is written to `scons_clean.log` and `scons.log`, respectively.
