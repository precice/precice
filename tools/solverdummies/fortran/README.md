# Compilation

Simply type `scons`. Assumption: preCICE is available as a shared library in your `LD_LIBRARY_PATH`.

# Run

You can test the dummy solver by coupling two instances with each other. Open two terminals and run
 * `./SolverDummy precice-config.xml SolverOne MeshOne`
 * `./SolverDummy precice-config.xml SolverTwo MeshTwo`

# Next Steps

If you want to couple any other solver against the dummy be sure to adjust the preCICE configuration (participant names, mesh names, data names etc.) to the needs of your solver, compare our [step-by-step guide for new adapters](https://github.com/precice/precice/wiki/Adapter-Example).
 
