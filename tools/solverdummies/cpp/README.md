# Compilation

**TODO: How do we build it with cmake? What do we have to set?...**

Do some cmake things...

# Run

**TODO arguments?...**

You can test the dummy solver by coupling two instances with each other. Open two terminals and run
 * `./solverdummy precice-config.xml SolverOne MeshOne`
 * `./solverdummy precice-config.xml SolverTwo MeshTwo`

# Next Steps

If you want to couple any other solver against the dummy solver be sure to adjust the preCICE configuration (participant names, mesh names, data names etc.) to the needs of your solver, compare our [step-by-step guide for new adapters](https://github.com/precice/precice/wiki/Adapter-Example).
