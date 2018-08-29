Python Interface for preCICE
----------------------------

# Dependencies

Download and install Cython from http://cython.org/#download or install it using your package manager.

# Building

1. Open terminal in this folder.
2. Execute the following command:

```
$ python setup.py build
```
This creates a folder `build` with the binaries.
3. Run 
```
$ python setup.py install
```
to install the module on your system. You might need `sudo`, depending on the how you have installed python. You can use the option `--prefix=your/default/path` to install the module at an arbitrary path of your choice (for example, if you cannot or don't want to use `sudo`).
4. Clean
```
$ python setup.py clean --all
```
again: depending on you system, you might need `sudo`.

It is recommended to use preCICE as a shared library here. `mpic++` is used as default compiler, if you want to use a different compiler, this can be done with the option `--mpicompiler=<yourcompiler>`. Example:
```
$ python setup.py build --mpicompiler=mpicc
```

# Using

1. Import `PySolverInterface` into your code:

```
import PySolverInterface
from PySolverInterface import *
```

2. If you use preCICE with MPI, you also have to add
   
```   
from mpi4py import MPI
```

To install mpi4py, you can e.g. (you need MPI already installed, e.g. libopenmpi-dev) 

```
sudo apt-get install python-pip
sudo pip install mpi4py
```


NOTE: 
- For an example of how the `PySolverInterface` can be used, refer to the [1D elastic tube example](https://github.com/precice/precice/wiki/Using-the-Python-API).
- In case the compilation fails with `shared_ptr.pxd not found` messages, check if you use the latest version of Cython.
