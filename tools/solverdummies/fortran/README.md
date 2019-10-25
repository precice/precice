# Compilation

Simply type `make`.

Assumption: `pkg-config` can find preCICE (try `pkg-config --modversion libprecice`). If you installed preCICE in a user directory, you may need to set your `PKG_CONFIG_PATH`. See the preCICE wiki on "Linking to preCICE" for more.

# Compilation with SCons (deprecated)

Simply type `scons`. Assumption: preCICE is available as a shared library in your `LD_LIBRARY_PATH`.

# Run

You can test the dummy solver by coupling two instances with each other. Open two terminals and run
 * `./solverdummy ../precice-config.xml SolverOne MeshOne`
 * `./solverdummy ../precice-config.xml SolverTwo MeshTwo`

Alternatively, to test the case where two solvers exchange data, run: 
 * `./solverdummy_A`
 * `./solverdummy_B`

# Next Steps

If you want to couple any other solver against the dummy be sure to adjust the preCICE configuration (participant names, mesh names, data names etc.) to the needs of your solver, compare our [step-by-step guide for new adapters](https://github.com/precice/precice/wiki/Adapter-Example).
 
