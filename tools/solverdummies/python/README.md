# Compilation

Don't forget to build the python bindings! See `precice/src/precice/bindings/python` for instructions

# Run

You can test the proxy by coupling two instances with each other. Open two terminals and run
 * `python solverproxy.py precice-config.xml SolverOne MeshOne`
 * `python solverproxy.py precice-config.xml SolverTwo MeshTwo`

# Next Steps

If you want to couple any other solver against the proxy be sure to adjust the preCICE configuration (participant names, mesh names, data names etc.) to the needs of your solver, compare our [step-by-step guide for new adapters](https://github.com/precice/precice/wiki/Adapter-Example).
 
