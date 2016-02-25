Python Interface for preCICE
----------------------------

Dependencies:

1. Download and install the version 0.23.4 of Cython (latest as of Feb 2016) from http://cython.org/#download

Steps:

1. Open terminal in this folder.
2. Define the environment variable PRECICE_ROOT as the path to preCICE on your system.
3. Execute the following command:

   $ python setup.py build_ext --inplace

3. This generates the PySolverInterface.so library. Add the path PRECICE_ROOT/src/precice/adapters/python to python's import sources. This can be done by including the following lines in your python code:

   precice_root = os.getenv('PRECICE_ROOT')
   precice_python_adapter_root = precice_root+"/src/precice/adapters/python"
   sys.path.insert(0, precice_python_adapter_root)

4. Import PySolverInterface into your code:

   import PySolverInterface
   from PySolverInterface import *

5. Use the PySolverInterface!

NOTE: for an example of how the PySolverInterface can be used, refer to the 1D tutorial.
