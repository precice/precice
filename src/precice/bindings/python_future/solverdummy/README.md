# Install Dependencies

Run in this directory:
```
pip3 install --user -r requirements.txt
```

Don't forget to install the precice [python bindings](../../../src/precice/bindings/python/README.md).

# Run

You can test the dummy solver by coupling two instances with each other. Open two terminals and run
 * `python3 solverdummy.py precice-config.xml SolverOne MeshOne`
 * `python3 solverdummy.py precice-config.xml SolverTwo MeshTwo`

