# MATLAB bindings

These bindings allow to use preCICE with MATLAB based on the C++ mex and C data API. They are still in an experimental state, so please utilize them with care. Any feedback is welcome.

## Requirements

MATLAB R2018a or later is required. The bindings were tested on R2018b and R2019a.

## Restrictions

- Unfortunately, the bindings can only be used if preCICE was compiled without MPI. This is a drawback especially if the user wants to couple a MATLAB solver to a parallel solver. The reason is that MATLAB crashes instantly if `MPI_Init()` is called in a MEX file. I am trying to understand and fix this, but I do not have a solution yet.
- The bindings were not tested with any toolboxes. In particular, they are not geared towards the parallel computing toolbox up to now.
- Several API functions are still missing.

## Compilation

The MATLAB script `compile_matlab_bindings_for_precice.m` located in this folder compiles the bindings. Simply running it from MATLAB should do.
The script uses `pkg-config` to determine the necessary flags. If `pkg-config` is unable to find the flags, the script will throw an error. Please refer to the [Linking to preCICE](https://github.com/precice/precice/wiki/Linking-to-preCICE) page in the preCICE wiki for details regarding `pkg-config`.

If using the script fails for some reason, please let us know.

## Add to path

To use the bindings, you have to add this folder to your MATLAB path using [`addpath`](https://de.mathworks.com/help/matlab/ref/addpath.html?searchHighlight=addpath&s_tid=doc_srchtitle). Adding this folder alone is sufficient (you don't have to use `genpath`). You can let MATLAB do this at startup by adding the respective line to your startup, see [here](https://de.mathworks.com/help/matlab/matlab_env/add-folders-to-matlab-search-path-at-startup.html).

## Usage

The API is accessible in MATLAB as `precice.\<API function\>`. E.g., a solver interface is created using ``interface = preCICE.SolverInterface(config_file_name,0,1)``.
The function syntax is mostly identical to the syntax of the C++ API.
_More details follow._
