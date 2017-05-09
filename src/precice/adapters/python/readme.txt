Python Interface for preCICE
----------------------------

Dependencies:

1. Download and install the version 0.23.4 of Cython (latest as of Feb 2016) from http://cython.org/#download

Steps:

1. Open terminal in this folder.
2. Define the environment variable PRECICE_ROOT as the path to preCICE on your system.
3. Execute the following command:

   $ python setup.py build_ext --inplace

This generates the PySolverInterface.so library. If you use preCICE as a static library, you have to manually add all thrid-party libraries (boost, python, petsc) that you use to setup.py. It is recommended to use preCICE as a shared library here. 

4. Add the path PRECICE_ROOT/src/precice/adapters/python to python's import sources. This can be done by including the following lines in your python code:

   precice_root = os.getenv('PRECICE_ROOT')
   precice_python_adapter_root = precice_root+"/src/precice/adapters/python"
   sys.path.insert(0, precice_python_adapter_root)

5. Import PySolverInterface into your code:

   import PySolverInterface
   from PySolverInterface import *

6. If you use preCICE with MPI, you also have to add
   
   from mpi4py import MPI

To install mpi4py, you can e.g. (you need MPI already installed, e.g. libopenmpi-dev) 

   sudo apt-get install python-pip
   sudo pip install mpi4py

7. Use the PySolverInterface!

NOTE: 
- For an example of how the PySolverInterface can be used, refer to the 1D tutorial. (https://github.com/precice/elastictube1d/tree/master/PythonDummy)
- In case the compilation fails with `'shared_ptr.pxd' not found` messages, check if `.local` contains the latest version of Cython.
