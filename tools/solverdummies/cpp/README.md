# Compilation

## With CMake

* preCICE has to be built using CMake! Otherwise it does not work.
* Postpone for now **TODO**

## Without CMake

Add command `export LIBRARY_PATH=$PRECICE_ROOT/build/last:$LIBARY_PATH` and
`export LIBRARY_PATH=$PRECICE_ROOT/build/last:$LIBARY_PATH` to `~/.bashrc` and run `source ~/.bashrc`
before running compilation command  

# Run

**TODO arguments?...**

You can test the dummy solver by coupling two instances with each other. Open two terminals and run
 * `./solverdummy precice-config.xml SolverOne MeshOne`
 * `./solverdummy precice-config.xml SolverTwo MeshTwo`

# Next Steps

If you want to couple any other solver against the dummy solver be sure to adjust the preCICE configuration (participant names, mesh names, data names etc.) to the needs of your solver, compare our [step-by-step guide for new adapters](https://github.com/precice/precice/wiki/Adapter-Example).
