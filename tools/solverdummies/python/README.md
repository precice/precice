# Compilation

Don't forget to build the Python bindings! Go [here](https://github.com/precice/precice/blob/develop/src/precice/bindings/python/README.md) for instructions.

# Run

You can test the dummy solver by coupling two instances with each other. Open two terminals and run
 * `python solverdummy.py precice-config.xml SolverOne MeshOne`
 * `python solverdummy.py precice-config.xml SolverTwo MeshTwo`

# Next Steps

If you want to couple any other solver against the dummy solver be sure to adjust the preCICE configuration (participant names, mesh names, data names etc.) to the needs of your solver, compare our [step-by-step guide for new adapters](https://github.com/precice/precice/wiki/Adapter-Example).
 
