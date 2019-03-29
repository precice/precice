# preCICE Change Log

All notable changes to this project will be documented in this file. For future plans, see our [Roadmap](https://github.com/precice/precice/wiki/Roadmap).

## develop
- The SolverInterface is now hardened against invalid IDs and misconfiguration using a consitent mechanism to express requirements.
- The SolverInterface now keeps track of the Mesh states, which results more informative error messages for mesh related functions.
- Sending data between participants is now fully asynchronous. This is relevant in one-way coupling scenarios, where the sending participant doesn't need to wait for the receiving one.

## 1.4.0
- The python modules are now proper packages tracking dependencies etc.
- Fix CMake now importable from binary directory.
- The Python module for the preCICE bindings `PySolverInterface` is renamed to `precice`. This change does not break old code. Please refer to [`src/precice/bindings/python/README.md`](src/precice/bindings/python/README.md) for more information.
- Add a pkg-config for preCICE (`libprecice.pc`).
- Use the Boost stacktrace library for cross-platform stacktrace printing. **This requires Boost 1.65.1**.
- Added explicit linking to `libdl` for `boost::stacktrace`.
- Reimplemented the internals of the nearest-projection mapping to significantly reduce its initialization time.
- The EventTimings now do a time normalization among all ranks, i.e., the first event is considered to happen at t=0, all other events are adapted thereto.
- The old CSV format of the EventTimings log files, split among two files was replaced by a single file, in structured JSON format.
- Fixed memory leaks in the `xml::XMLAttributes` and `xml::Validator*`.
- Removed the `xml::Validator*` classes and replaced them in `xml::XMLAttribute` with a set of "options".
- Made `xml::XMLAttribute` and `xml::XMLTag` chainable.
- Added manpages for binprecice and testprecice.
- Fixed memory leaks in `mesh::Mesh`.
- Fixed mapping classes not flushing the underlying caches on `clear()`.
- Fixed format of version logging (first preCICE log message).
- Added tolerance of iterations in quasi-Newton tests.
- CMake overhaul:
  - Converted to target-based system: precice, testprecice, binprecice
  - New options:
    - `PRECICE_Packages` (default ON) to configure CPack,
    - `PRECICE_InstallTest` (default ON) to configure installation of tests.  
      This includes the binary `testprecice` and necessary files.
      Use `PREFIX/share/precice` as `PRECICE_ROOT`.
  - Moved CMake files from `tools/cmake-modules` to `cmake/` (general scripts) and `cmake/modules` (find modules).
  - Migrated from file-globing to explicit source/interface/test-file lists.  
    Use `tools/updateSourceFiles.py` from project-root to update all necessary files.
  - `install` target installs:
     - the library `PREFIX/lib`.
     - the binaries `PREFIX/bin` and their manfiles into `PREFIX/share/man/man1`.
     - the CMake configuration files into `PREFIX/lib/cmake/precice`.
     - the pkg-config configuration files into `PREFIX/lib/pkgconfig`
     - the necessary files to run testprecice into `PREFIX/share/precice`. Use this as `PRECICE_ROOT` on installed system.
  - CTest definition of tests run in isolated working directories:
    - `precice.Base` for the base test suite
    - `precice.MPI2` run on 2 MPI ranks
    - `precice.MPI4` run on 4 MPI ranks
  - CPack configuration of target `package` to generate binary debian, tar and zip packages.
  - Added `CMakeLists.txt` to `tools/solverdummy/cpp`. It is an example of how to link to preCICE with CMake.
  - Extended the displayed information when configuring.
- Extended `updateSourceFiles.py` to verify the sources using `git ls-files --full-name` if available.
- Fixed the `io::VTKXMLExporter` not to write VertexNormals.
- Improved the user-friendliness of the tests.
  - `make test` will run all tests.
  - `make test_base` only a unproblematic base-set.
  - A timeout will kill hanging tests.
  - All tests sets run in isolated working directories.
- Added an (experimental) Fortran 2003 solver dummy.

## 1.3.0
- Update of build procedure for python bindings (see [`precice/src/bindings/python/README.md`](https://github.com/precice/precice/blob/develop/src/precice/bindings/python/README.md) for instructions). Note: you do not have to add `PySolverInterface.so` to `PYTHONPATH` manually anymore, if you want to use it in your adapter. Python should be able to find it automatically.   
- Make naming of log files consistent, following the pattern `precice-SOLVERNAME-logtype.log`, example: `precice-FLUID-eventTimings.log`
- Enable boost.geometry based preallocation. Speeds up initialization of PetRBF based mapping.
- Actions can now specify a `MeshRequirement`, such as the `ScaleByAreaAction`.
- Many events have been reworked and are now uniformly named. 
- There is a `syncMode` for events (for detailed performance measurements), configurable and off by default. 

## 1.2.0
- Make `polynomial=separate` the default setting for PetRBF.
- Removed ExportVRML functionality
- Build system:
  - Make `python=off` default.
- Building with Conda:
  - The helper scripts are now placed in the directory `tools/conda_building`. All the terms refering to `Anaconda` have been changed to `Conda`.
- Sending data exchange is now fully asynchronous, so that the sending participant never waits for the receiving one.
- Rename `src/precice/adapters` to `src/precice/bindings`
- adding `libprefix` option in scons build process to allow for non-standard library paths

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
