# preCICE Change Log

All notable changes to this project will be documented in this file. For future plans, see our [Roadmap](https://github.com/precice/precice/wiki/Roadmap).

## develop
- Build system:
  - Make `python=off` default.

## 1.1.1
- Fix SConstruct symlink build target failing when using lowercase build (debug, release) names.

## 1.1.0
- Build system:
  - Remove the `staticlib` and `bin` from the default targets to reduce the building time and storage requirements.
  - Change build types to mixed case, i.e. ```Debug``` and ```Release```. Old versions are retained for backward compatibility.
  - Make `mpicxx` default setting for compiler.
  - Experiemental support for building with Conda, see `tools/anaconda_building`
  - Use NumPy to figure out Python include paths.
  - Search for PETSc in more paths
  - Add experimental CMake build control files.

- Add a job command file for the SuperMUC HPC system in `tools`.
- `compileAndTest.py` Change the `-b` option to `-t`, do not crash when ./tests do no exist, make `mpirun` command configurable
- Use `libxml2` for XML parsing, this makes `libxml2-dev` a dependency.
- Update EventTimings framework.
- Add python script to plot Events on a timeline.
- PETSc RBF mapping now supports conservative mapping with a separated polynomial
- Converted all tests to the new, boost test based, unit testing framework.
- Removed the `tarch` legacy library.
- Use `boost::signal2` for implement observer pattern for the Mesh class.
- Add contributer guidelines.


## 1.0.3
- Fix compilation for boost 1.66, see issue #93.

## 1.0.2
- Fix bug in the mesh repartitioning for plane-like coupling interfaces and small gaps between both sides.

## 1.0.1
- Fix compilation issue with the python interface.
