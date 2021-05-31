# Compilation

**preCICE has to be installed using the provided binaries or built using CMake! Otherwise this approach does not work. For more information refer to [our documentation](https://www.precice.org/docs.html)**

You can use the provided `CMakeLists.txt` to build with CMake.

1. run `cmake .` in this folder
2. run `make`

# Run

You can test the dummy solver by coupling two instances with each other. Open two terminals and run
 * `./solverdummy ../precice-config.xml SolverOne MeshOne`
 * `./solverdummy ../precice-config.xml SolverTwo MeshTwo`

# Next Steps

If you want to couple any other solver against the dummy solver be sure to adjust the preCICE configuration (participant names, mesh names, data names etc.) to the needs of your solver, compare our [step-by-step guide for new adapters](https://www.precice.org/couple-your-code-overview.html).
