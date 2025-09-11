# Compilation

**preCICE has to be installed using the provided binaries or built using CMake! Otherwise this approach does not work. For more information refer to [our documentation](https://www.precice.org/docs.html)**

You can use the provided `CMakeLists.txt` to build with CMake.

1. run `cmake .` in this folder
2. run `make`

**Alternatively, you can directly include the precice module into your source code, for more information refer to the README.md in https://github.com/precice/fortran-module**

# Run

You can test the dummy solver by coupling two instances with each other. Open two terminals and run
 * `./solverdummy ../precice-config.xml SolverOne`
 * `./solverdummy ../precice-config.xml SolverTwo`

# Next Steps

If you want to couple any other solver against the dummy solver be sure to adjust the preCICE configuration (participant names, mesh names, data names etc.) to the needs of your solver, compare our [step-by-step guide for new adapters](https://www.precice.org/couple-your-code-overview.html).
