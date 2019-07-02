**DEPRECATED!!!**

`PySolverInterface` are the old python language bindings for preCICE and will be removed in version 2.0.0. Please use the new python language bindings [`precice`](../python/README.md).

Python language bindings for preCICE
----------------------------

# Dependencies

Install the [precice python bindings](../python/README.md).

# Installation

In this directory, run:
```
$ pip3 install --user .
```

# Using

1. Import `PySolverInterface` in your code:

```python
import PySolverInterface
from PySolverInterface import *
```

2. If you compiled preCICE with MPI, you also have to install mpi4py:

```   
$ pip3 install --user mpi4py
```   

Then add the following import to your python file:
```python
from mpi4py import MPI # Initialize MPI
```
