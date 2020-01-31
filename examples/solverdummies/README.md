# Solverdummies

The `solverdummies` are minimal working examples for using preCICE with different languages. Please refer to the corresponding subfolders for information on compilation and usage.

# Combining different solverdummies

preCICE allows to couple codes written in different programming languages. The solverdummies can be used to demonstrate this feature. First, you have to build the solverdummies you want to use. Then run the following commands from this folder; each command in one terminal:

```
./c/solverdummy precice-config.xml SolverOne MeshOne
./cpp/solverdummy precice-config.xml SolverTwo MeshTwo
```

This combines the c solverdummy with the cpp solverdummy. Note that both solverdummies use the same `precice-config.xml`. Feel free to experiment with other combinations of solverdummies or couple your own solver with one of the solverdummies. Note that there are also solverdummies for [Fortran 2003 and above](https://github.com/precice/fortran-module/tree/master/examples/solverdummy), [python](https://github.com/precice/python-bindings/tree/develop/solverdummy), and [MATLAB](https://github.com/precice/matlab-bindings/tree/develop/solverdummy).
