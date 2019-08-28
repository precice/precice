# MATLAB bindings

These bindings allow to use preCICE with MATLAB based on the C++ MEX and C data API. They are still in an experimental state, so please utilize them with care. Any feedback is welcome.

## Requirements

MATLAB R2018a or later is required. The bindings were tested on R2018b and R2019a.

## Restrictions

- An issue causes MATLAB to crash upon SolverInterface initialization if precice was compiled with openmpi. This issue can be resolved by installing openmpi from source using the option `-disable-dlopen`. For reference, see e.g. [here](https://stackoverflow.com/questions/26901663/error-when-running-openmpi-based-library). Alternatively, the user can switch to a different MPI implementation, e.g. MPICH (other implemenations were not tested).
- The API function `getMeshHandle` is missing. Also, there are no wrapper classes for MeshHandle etc.
- Currently, only one instance of the `SolverInterface` class can exist at the same time in a single MATLAB instance. If the user wishes to couple multiple participants based on MATLAB, he is supposed to start them in different MATLAB instances. If, for some reason, the user needs multiple instances of `SolverInterface`, he should use the OOP variant (Multiple instances of `SolverInterfaceOOP` can exist at the same time).

## Compilation

The MATLAB script `compile_matlab_bindings_for_precice.m` located in this folder compiles the bindings. Simply running it from MATLAB should do.

The script uses `pkg-config` to determine the necessary flags. If `pkg-config` is unable to find the flags, the script will throw an error. Please refer to the [Linking to preCICE](https://github.com/precice/precice/wiki/Linking-to-preCICE) page in the preCICE wiki for details regarding `pkg-config`.

If using the script fails for some reason, please let us know.

## Add to path

To use the bindings, you have to add this folder (`$PRECICE_ROOT/src/bindings/matlab`) to your MATLAB path using [`addpath`](https://de.mathworks.com/help/matlab/ref/addpath.html?searchHighlight=addpath&s_tid=doc_srchtitle). Adding this folder alone is sufficient (you don't have to use `genpath`). You can let MATLAB do this at startup by adding the respective line to your `startup.m`, see [here](https://de.mathworks.com/help/matlab/matlab_env/add-folders-to-matlab-search-path-at-startup.html).

## Usage

The API introduces MATLAB wrapper classes for the `SolverInterface` class and the preCICE constants. They are accessible in MATLAB as `precice.SolverInterface` and `precice.Constants`.

The function syntax is mostly identical to the syntax of the C++ API. The following things should be noted:
- C++ `int`s correspond to MATLAB `int32`s.
- Wherever the C++ API expects pointers, the MATLAB API expects a matrix/vector instead. If the user wants to pass vector data (e.g. vertex coordinates) for multiple vertices, the shape of the corresponding matrix must be `[dim numVertices]`, where `dim` is the problem dimension. Thus, each **column** must correspond to a vertex, and each line must correspond to a coordinate - **not** vice versa. Users should try to respect this in their MATLAB code from the start, because transposing can be costly for huge matrices.
- Output arguments which are pointers passed input arguments to the C++ preCICE API are replaced by output matrices. E.g., the C++ API function
```
readBlockScalarData(int dataID, int size, const int *valueIndices, double *values)
```
is found in the MATLAB bindings as
```
values = readBlockScalarData(dataID, size, valueIndices)
```

## Out of process variant

The C++ MEX API supports [out of process execution](https://de.mathworks.com/help/matlab/matlab_external/out-of-process-execution-of-c-mex-functions.html) of MEX functions. This feature is implemented in the class `precice.SolverInterfaceOOP`. This class works exactly like `precice.SolverInterface`. Internally, however, the gateway function that calls the preCICE routines is run on a `mexHost` object.
This has the following advantages:
- Multiple instances of `SolverInterfaceOOP` can exist at the same time.
- If the gateway function crashes, then MATLAB will not crash. Only the mexHost object will crash.
However, using the OOP variant is **significantly** slower than the normal in process 
