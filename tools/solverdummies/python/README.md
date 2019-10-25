# Install Dependencies

Run in this directory:
```
pip3 install --user -r requirements.txt
```

Don't forget to install the precice [python bindings](../../../src/precice/bindings/python/README.md).

# Run

You can test the dummy solver by coupling two instances with each other. Open two terminals and run
 * `python3 solverdummy.py ../precice-config.xml SolverOne MeshOne`
 * `python3 solverdummy.py ../precice-config.xml SolverTwo MeshTwo`

Alternatively, to test the case where two solvers exchange data, run: 
 * `python3 solverdummy_A.py`
 * `python3 solverdummy_B.py`

# Next Steps

If you want to couple any other solver against the dummy solver be sure to adjust the preCICE configuration (participant names, mesh names, data names etc.) to the needs of your solver, compare our [step-by-step guide for new adapters](https://github.com/precice/precice/wiki/Adapter-Example).
 
