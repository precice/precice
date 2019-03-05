# Compilation

## With CMake

**preCICE has to be installed using the provided binaries or built using CMake! Otherwise this approach does not work. For more information refer to [our wiki](https://github.com/precice/precice/wiki/Get-preCICE)**

You can use the provided `CMakeLists.txt` to build with CMake.

1. run `cmake .` in this folder
2. run `make`

## Without CMake (deprecated)

**If you used scons to build preCICE, continue here**

For building the C++ solverdummies you first have to modify your `LIBRARY_PATH` and `LD_LIBRARY_PATH` such that `libprecice.so` can be found. You can, for example add the following lines to your `~/.bashrc`:

```
export LIBRARY_PATH=$PRECICE_ROOT/build/last:$LIBRARY_PATH
export LD_LIBRARY_PATH=$PRECICE_ROOT/build/last:$LD_LIBRARY_PATH
export CPLUS_INCLUDE_PATH=$PRECICE_ROOT/src:$CPLUS_INCLUDE_PATH
```

After that, don't forget to `source ~/.bashrc`.

Now run `g++ solverdummy.cpp -lprecice -o solverdummy`.

# Run

You can test the dummy solver by coupling two instances with each other. Open two terminals and run

* `./solverdummy precice-config.xml SolverOne MeshOne`
* `./solverdummy precice-config.xml SolverTwo MeshTwo`

# Next Steps

If you want to couple any other solver against the dummy solver be sure to adjust the preCICE configuration (participant names, mesh names, data names etc.) to the needs of your solver, compare our [step-by-step guide for new adapters](https://github.com/precice/precice/wiki/Adapter-Example).
