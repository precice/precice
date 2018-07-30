Python Interface for preCICE
----------------------------

# Dependencies

Download and install Cython from http://cython.org/#download or install it using your package manager.

# Building

1. Open terminal in this folder.
2. Define the environment variable `PRECICE_ROOT` as the path to preCICE on your system.
3. Define the environment variable `PRECICE_MPI_IMPLEMENTATION` to either `openmpi` or `mpich`, depending on your MPI implementation. You can determine your implementation by typing `mpic++ -v` or `mpic++ --showme::version`.
4. Execute the following command:

```
$ python setup.py build_ext --inplace
```
(works with python 2 and 3)


This generates `PySolverInterface.so` and `PySolverInterface.cpp`. If you use preCICE as a static library, you have to manually add all third-party libraries (boost, python, petsc) that you use to setup.py. It is recommended to use preCICE as a shared library here.

# Using

1. Add the path `$PRECICE_ROOT/src/precice/bindings/python` to python's import sources. This can be done by including the following lines in your python code:

```
precice_root = os.getenv('PRECICE_ROOT')
precice_python_adapter_root = precice_root + "/src/precice/bindings/python"
sys.path.insert(0, precice_python_adapter_root)
```

2. Import `PySolverInterface` into your code:

```
import PySolverInterface
from PySolverInterface import *
```

3. If you use preCICE with MPI, you also have to add
   
```   
from mpi4py import MPI
```

To install mpi4py, you can e.g. (you need MPI already installed, e.g. libopenmpi-dev) 

```
sudo apt-get install python-pip
sudo pip install mpi4py
```


NOTE: 
- For an example of how the `PySolverInterface` can be used, refer to the 1D tutorial. (https://github.com/precice/elastictube1d/tree/master/PythonTube)
- In case the compilation fails with `shared_ptr.pxd not found` messages, check if you use the latest version of Cython.
