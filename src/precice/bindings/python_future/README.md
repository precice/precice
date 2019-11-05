Python language bindings for preCICE
------------------------------------

These are the python bindings for preCICE.

# Installing the package

We recommend [using pip3](https://github.com/precice/precice/blob/develop/src/precice/bindings/python/README.md#using-pip3) for the sake of simplicity.

## Using pip3

### preCICE system installs

For system installs of preCICE, this works out of the box.

In this directory, execute:
```
$ pip3 install --user .
```
*note the dot at the end of the line*

This will fetch cython, compile the bindings and finally install the precice package.

### preCICE at custom location (setting PATHS)

If preCICE was installed in a custom prefix, or not installed at all, you have to extend the following environment variables:
- `LIBRARY_PATH`, `LD_LIBRARY_PATH` to the library location, or `$prefix/lib`
- `CPATH` either to the `src` directory or the `$prefix/include`

## Using setup.py

### preCICE system installs

In this directory, execute:
```
$ python3 setup.py install --user
```

### preCICE at custom location (setting PATHS)

see above. Then run
```
$ python3 setup.py install --user
```

### preCICE at custom location (explicit include path, library path, or mpicompiler)

1. Install cython via pip3
```
$ pip3 install --user cython
```
2. Open terminal in this folder.
3. Build the bindings

```
$ python3 setup.py build_ext --mpicompiler=mpicc --include-dirs=$PRECICE_ROOT/src --library-dirs=$PRECICE_ROOT/build/last
```

**Options:**
- `--include-dirs=`, default: `''` 
  Path to the headers of preCICE, point to the sources `$PRECICE_ROOT/src`, or the your custom install prefix `$prefix/include`.
- `--library-dirs=`, default: `''` 
  Path to the libary of preCICE, point to the build directory (scons: `$PRECICE_ROOT/build/last`, cmake: wherever you configured the build), or to the custom install prefix `$prefix/lib`.
- `--mpicompiler=`, default: `mpic++` 
  MPI compiler wrapper of choice.

**NOTES:**

- If you build preCICE using CMake, you can pass the path to the CMake binary directory using `--library-dirs`.
- It is recommended to use preCICE as a shared library here.
- If you used scons for building precice and `PRECICE_ROOT` is defined, you can also use the script `build_and_install.sh`.

4. Install the bindings
```
$ python3 setup.py install --user
```

5. Clean-up _optional_
```
$ python3 setup.py clean --all
```

# Test the installation

Update `LD_LIBRARY_PATH` such that python can find `precice.so`
```
$ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PRECICE_ROOT/build/last
```

Run the following to test the installation:
```
$ python3 -c "import precice_future as precice"
```

## Unit tests

1. Clean-up __mandatory__ (because we must not link against the real `precice.so`, but we use a mocked version)
```
$ python3 setup.py clean --all
```

2. Set `CPLUS_INCLUDE_PATH` (we cannot use `build_ext` and the `--include-dirs` option here)
```
$ export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:$PRECICE_ROOT/src
```

3. Run tests with
```
python3 setup.py test
```

# Using with MPI

If precice was compiled with MPI, you have to initialize MPI prior to configuring a SolverInterface.
To do so, install the package `mpi4py` and import `MPI` in your python module.

```
$ pip3 install --user mpi4py
```

```python
from mpi4py import MPI # Initialize MPI 
```

**NOTE:**
- For an example of how `precice` can be used, refer to the [1D elastic tube example](https://github.com/precice/precice/wiki/1D-elastic-tube-using-the-Python-API).
- In case the compilation fails with `shared_ptr.pxd not found` messages, check if you use the latest version of Cython.
- If you want to use the old interface (precice version < 1.4.0), please also install the corresponding wrapper [`PySolverInterface`](https://github.com/precice/precice/tree/changingNameOfPySolverInterface/src/precice/bindings/PySolverInterface).
