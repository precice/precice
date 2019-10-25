# MATLAB solverdummy

This dummy illustrates the use of precice with MATLAB using the MATLAB bindings based on the C data API.

**Please note:** Currently, the MATLAB solverdummy can't be coupled with other solvers unless an exchange directory is specified with its full path in the config file. So, if you wish to couple with a different solver, please make sure to modify the config file accordingly.

## Compilation

No compiling is necessary. You only have to compile the MATLAB bindings.

## Run

To run the dummy, open two MATLAB instances and call

* `solverdummy ../precice-config.xml SolverOne MeshOne`
* `solverdummy ../precice-config.xml SolverTwo MeshTwo`

Since `solverdummy` is a MATLAB function and MATLAB treats blank space separated argumenst as char arrays, you can equivalently call

`solverdummy('../precice-config.xml','SolverOne','MeshOne')`

and analogously for the second call. 
Naturally, you may also couple the MATLAB dummy with another dummy instead.
