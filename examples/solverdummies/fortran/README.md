# Compilation

## With CMake

**preCICE has to be installed using the provided binaries or built using CMake! Otherwise this approach does not work. For more information refer to [our wiki](https://github.com/precice/precice/wiki/Get-preCICE)**

You can use the provided `CMakeLists.txt` to build with CMake.

1. run `cmake .` in this folder
2. run `make`

## Make

Simply type `make`.

Assumption: `pkg-config` can find preCICE (try `pkg-config --modversion libprecice`). If you installed preCICE in a user directory, you may need to set your `PKG_CONFIG_PATH`. See the preCICE wiki on "Linking to preCICE" for more.

# Run

You can test the dummy solver by coupling two instances with each other. Open two terminals and run
 * `./solverdummy ../precice-config.xml SolverOne MeshOne`
 * `./solverdummy ../precice-config.xml SolverTwo MeshTwo`

# Next Steps

If you want to couple any other solver against the dummy solver be sure to adjust the preCICE configuration (participant names, mesh names, data names etc.) to the needs of your solver, compare our [step-by-step guide for new adapters](https://github.com/precice/precice/wiki/Adapter-Example).
