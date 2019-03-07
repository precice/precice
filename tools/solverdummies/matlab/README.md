# MATLAB solverdummy - experimental

This dummy marks the first successful attempts to get preCICE working with MATLAB. Right now, it is not much more than a copy of the cpp solverdummy refactored to work with the [C data API](https://de.mathworks.com/help/matlab/cc-mx-matrix-library.html).

## Compilation

First, you need a C++ compiler that is compatible with MATLAB mex. Use the `mex -setup c++` command to check the configuration of your MEX. If necessary, change the compiler. Find more details in [MATLABs documentation](https://de.mathworks.com/help/matlab/matlab_external/choose-c-or-c-compilers.html).

Then, you can compile the solverdummy by calling

`mex solverdummy.cpp -lprecice -output solverdummy`

in MATLAB. This will create a file named `solverdummy.mexa64` which can be called like every other MATLAB function.

## Run

To run the dummy, open two MATLAB instances and call

* `solverdummy ../precice-config.xml SolverOne MeshOne`
* `solverdummy ../precice-config.xml SolverTwo MeshTwo`

Since `solverdummy` is a MATLAB function and MATLAB treats blank space separated argumenst as char arrays, you can equivalently call

`solverdummy('../precice-config.xml','SolverOne','MeshOne')`

and analogously for the second call. 
Naturally, you may also couple the MATLAB dummy with another dummy instead.
