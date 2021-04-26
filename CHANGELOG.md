# preCICE Change Log

All notable changes to this project will be documented in this file. For future plans, see our [Roadmap](https://www.precice.org/fundamentals-roadmap.html).

## 2.2.0

- Added a file sink to the test runner, which additionally writes the test log output to the file `test.log`.
- Added a verbose log file `test.debug.log` to tests. Submitting this file alongside bug reports will significantly simplify debugging the tests.
- Added check for user-defined python actions.
- Added check to ensure `advance()` isn't called with invalid timestep size.
- Added error message for incorrectly configured acceleration data in a serial implicit coupling scheme.
- Added missing documentation to XML tags and attributes.
- Added support for `BUILD_TESTING`, which may be used to disable the compilation and execution of the tests.
- Added support for escaped characters in XML.
- Added the build configuration (Release/Debug) to the preCICE startup statement and inform the user when `Debug` an `Trace` logs are not available.
- Added watch-integral functionality to, for instance, compute total force, total stress, or flow rate at coupling mesh.
- Changed an error to a warning when all sub-vectors of the residual-sum preconditioner are numerically zero.
- Changed the CMake FindPETSc module to a robust wrapper based on pkg-config. The new method is more robust and simplifies the compilation of preCICE in e.g. SLURM jobs.
- Changed the log format for tests. They now contain the participant name and are `|`-separated.
- Clarified the wording of errors messages in data access functions. They now refer to vertex IDs instead of indices.
- Deprecated API method `getMeshIDs`. Call `getMeshID` for specific mesh names instead. This method will be removed in a future release.
- Extended configuration and repartitioning to allow the user to define multiple mappings from and/or to the same mesh.
- Fixed MacOS compilation and test errors.
- Fixed a bug in ReceivedPartition which led to problems when coupling multiple participants.
- Fixed a bug in the configuration of Aitken underrelexation.
- Fixed an issue when running the tests with Intel MPI.
- Fixed an occasional issue solving the system matrix in PETSc-based RBF mappings.
- Fixed boost log_level issues on MacOS.
- Fixed compilation error emitted by intel compilers.
- Fixed cryptic assertion for forgetting the max-iterations tag. Now max-iterations is enforced in the configuration. 
- Fixed indexing bug in solverdummies.
- Fixed input checks for data access functions.
- Fixed interleaved assertion output.
- Fixed parameter type of operator"" in TestContext, which was a system-dependant error.
- Fixed parsing error on systems without locales installed. This fixes issues when running preCICE in minimal docker containers.
- Fixed syntax of the Fortran function `precicef_get_mesh_vertex_size_`, which lead to incorrect name de-mangling.
- Fixed the data type and precision in exported VTK files.
- Fixed two wrong assertions in QR factorization, which did not allow meshes with only a single partition. 
- Improved checks of configuration related to data access.
- Improved compiler compatibility of assertion.
- Improved the error message for not exchanging data over the same mesh used for convergence measures.
- Increased the minimum required C++ version from 11 to 14. This was triggered by `Boost.Geometry` increasing their minimum version to C++14 in Boost `1.75`.
- Removed obsolete `trigger-solver-plot` config option.

## 2.1.1

- Fixed a compilation error emitted by intel compilers when compiling tests.
- Fixed an error message for unique acceleration subtags.
- Fixed parsing error on systems without locales installed. This fixes issues when running preCICE in minimal docker containers.
- Fixed preCICE erroneously expects cyclic communicator for IMVJ with restart.
- Fixed system-dependent compilation error when compiling tests.
- Updated the PKGBUILD for the Arch Linux User Repository to be compatible with preCICE v2.

## 2.1.0

- Added MPI version detection and print it during the CMake configuration.
- Added a new action to sum up data values, `action:summation`.
- Added check for only allowing watchpoints on a provided meshes.
- Added check for the coordinates of configured watchpoints. They have to match the dimensionality of the `<solver-interface>`.
- Added check to prevent `<use-mesh>` from the same participant.
- Added control flow checks to C bingings.
- Added many tests for the communication abstraction.
- Added option to make a convergence measure strict. It has to converge then and leads to a premature simulation stop if not. 
- Added parallel support for Eigen RBF mapping.
- Added platform-specific defaults of the loopback interface name to the `network` attribute of socket connections.
- Added reset of written data 0 in `advance()` to simplify detection of missing write data calls.
- Added several checks to prevent that false API usage leads to complicated quasi-Newton assertions.
- Added support for planar quads based on triangulation, needed, for instance for nearest-projection mapping. Quads will be split along the shortest diagonal to minimize the mapping error. Note that this triangulation will be visible in configured mesh exports.
- Added tests which build solverdummies and tests which run them in various configurations.
- Added tooling to prevent changelog conflicts.
- Added warning when configuring m2n as mpi when preCICE is compiled with OpenMPI. This is known to cause connection issues.
- Allowed multiple convergence measures of same coupling data.
- Changed `com::MPIDirectcommunication` to work only for Master-Slave connections.
- Changed test setup to a simpler and consistent version using the new `testing::TestContext`. Tests now require to run on 4 MPI ranks. They will still compile but not run when `MPICommunication=OFF`.
- Changed the CMake FindNumPy module to only consider information based on the selected python interpreter.
- Changed the minimum required PETSc version to 3.12, which delivers consistent results across platforms.
- Disabled tests based on MPIPorts and MPISinglePorts when using Open MPI.
- Enabled RBF-based tests in partiton unit-tests and serial integration tests.
- Extended iteration logging by total and dropped quasi-Newton columns. 
- Extended title headers of convergence and iteration files by measure abbreviation.
- Fixed MPIPorts and MPISinglePorts not always closing ports.
- Fixed SocketCommunication setting up a port and writing connection info even if there are no requesters.
- Fixed compatibility with Boost 1.73.0.
- Fixed memory leaks and hanging communication when not calling `precice::SolverInterface::finalize()`.
- Fixed memory leaks in C bindings.
- Fixed memory leaks in Fortran bindings in case of missing call to`X_finalize()`.
- Fixed memory leaks when using petrbf mappings in some conditions.
- Fixed occasional errors when using PETRBF with a preCICE-managed MPI Communicator.
- Fixed over-subscription error when executing the tests with recent versions of Open MPI on less than 4 physical cores.
- Fixed silent truncation of numeric values in the config with C locales using decimal comma. The parser now always uses the locale `en_US.UTF-8`.
- Fixed target `test_install` requiring CMake version 1.13.
- Fixed value semantics of `precice::SolverInterface`.
- Improved **all** error messages.
- Improved readability of relative convergence measure INFO logs.
- Refactored `com::Communication` handling of rank adjustments.
- Refactored `cplscheme::BaseCouplingScheme` and derived classes. Introduce `cplscheme::BiCouplingScheme`.
- Refactored `mesh::BoundingBox` into seperate class.
- Removed `m2n:mpi-single`, which never worked outside tests.
- Removed convergence file logging for min-iterations convergence measure.
- Removed deprecated and untested Manifold Mapping. API functions `hasToEvaluateSurrogateModel` and `hasToEvaluateFineModel` remain as nop stubs.
- Removed logging of average convergence rates.
- Removed support for `CPP` and `CPPFLAGS` environment variables.
- Renamed a few functions of the Fortran bindings for consistency with the C++ API. This does not break backwards compatibility, but the previous functions are now considered deprecated and will be removed in v3.
- Removes the summary of event timing when preCICE finalizes and writes it to the file `precice-PARTICIPANT-events-summary.log` instead.
- Restricted the configuration of WatchPoints to provided meshes.
- Split multi-setup integration tests into multiple single-setup tests.

## 2.0.2

- Fixed a critical bug in the testing framework.
- Fixed a critical bug in the partitioning for geometric filter set to `on-master` in `<use-mesh>` tags. The default configuration is `on-slave`.

## 2.0.1

- Fixed broken pkg-config file in some cases due to CMake `GNUInstallDirs`.
- Fixed system-dependent error when displaying mapping distance information for empty partitions.
- Improved RBF error messages by clarifying them and giving a hint to the user on what to do next.

## 2.0.0

- Added CMake build type fallback to `Debug` in case it wasn't provided.
- Added CMake check for C++11 library conformance. This is especially helpful when using Intel Compilers.
- Added CMake options to enable native bindings `PRECICE_ENABLE_C`, `PRECICE_ENABLE_FORTRAN` (on by default).
- Added `examples/` to installation and package
- Added a generator for markdown references `binprecice md`.
- Added caching to the CMake library validation
- Added directional directory level for file-based connection exchange. For each connection, there is now a directory in `precice-run` of the form `Accepter-Requestor`.
- Added distance statistics of nearest-neighbour and nearest-projection mappings between mesh pairs as debug output.
- Added grouped tests by module to CTest.
- Added information to the log of the first written Data values.
- Added log statements to the connection information file writers and listeners including full paths.
- Added step to remove the connection directories after connecting the slaves. `precice-run` will be empty after a successful run.
- Added support for python 3 in python actions.
- Added the mesh name to the information used to generate connection information files, which is required for the two-level initialization.
- Changed CMake to always validate dependencies. Set `PRECICE_ALWAYS_VALIDATE_LIBS=NO` to disable this behaviour.
- Changed the internal handling of meshes by removing sub-meshes, the type hierarchy based on `mesh::PropertyContainer`, and the obsolete `mesh::Group` and `mesh::Merge`. This improves memory consumption, dramatically reduces allocations and improves locality when traversing primitives.
- Changed unit tests to run only on MPI Rank 0.
- Completely removed server mode. Now, the only supported parallelization concept is the peer-to-peer master-slave mode.
- Disabled the installation of the test binary and files by default.
- Dropped official python2 support for python bindings ([remove tests](https://github.com/precice/systemtests/commit/dba924447996574967b2295cf652fb32bec58020)).
- Removed all experimental python bindings [`precice`](https://github.com/precice/precice/tree/v1.6.1/src/precice/bindings/python), [`precice-future`](https://github.com/precice/precice/tree/v1.6.1/src/precice/bindings/python-future), [`PySolverInterface`](https://github.com/precice/precice/tree/v1.6.1/src/precice/bindings/PySolverInterface).
- Fixed a bug in the XML parser which did not correctly checked tag occurrence.
- Fixed a bug in the XML parser which lead to ignored error messages from `libxml2`.
- Fixed the Debian package generation by using `GNUInstallDirs`, providing a correct `changelog` and `SOVERSION`, as well as generating a package name including the `SOVERSION`.
- Improved efficiency of nearest projection mapping of matching meshes using lazy generation of index trees.
- Introduced preCICE-MATLAB bindings (https://github.com/precice/precice/pull/494, https://github.com/precice/precice/pull/580) and provided them in [`precice/matlab-bindings`](https://github.com/precice/matlab-bindings).
- Merged the `SolverInterface::configure()` into the `SolverInterface` constructors. They now have a second parameter for the configuration file.
- Moved Fortan 2003 bindings (`src/precice/bindings/f2003`) and solverdummy (`tools/solverdummy/f2003`) to a separate repository.
- Refactored and made two-level initialization configurable through `"use-two-level-init"` in `m2n`.
- Refactored the XML documentation generation out of the `xml::XMLAttribute` and `xml::XMLTag` classes into `xml/Printer.[c/h]pp`.
- Released finalized version of python bindings in independent repository: [`precice/python-bindings`](https://github.com/precice/python-bindings). Package is named [`pyprecice`](https://github.com/precice/python-bindings/blob/3b9aec6c529814e6904a6a4697cf92388d4c4bf0/setup.py#L18) and supports the preCICE version >= 2.0.0.
- Removed `MeshHandle` from API and replace use in integration tests by `SolverInterfaceImpl::mesh()`.
- Removed an unnecessary assertion in `getMeshVertexIDsFromPositions()`.
- Removed deprecated SCons.
- Removed deprecated `HierarchicalAitkenAcceleration`.
- Removed deprecated `ModifyCoordinatesAction`.
- Removed deprecated voxel queries in `src/query/`.
- Removed packaging files specific to Ubuntu 18.04 as it is covered by CPack.
- Renamed CMake variables `MPI`, `PETSC`, `PYTHON` to `PRECICE_MPICommunication`, `PRECICE_PETScMapping`, `PRECICE_PythonActions`
- Replaced geometric filter option "filter-first" and "broadcast-filter" by "on-master" and "on-slaves", respectively, to generalize to two-level initialization.
- Restructured tools and bindings:
  - Moved developer tools to `tools/`.
  - Moved user tools to `extras/`.
  - Moved native bindings to `extras/bindings/`.
- Simplify parallel configuration:
  - Automatically add `master:mpi-single` for parallel participant if necessary.
  - No longer require `gather-scatter` distribution type for a `m2n` with at least one serial participant.
  - Automatically choose suitable RBF implementation based on whether preCICE was built with PETSc and whether the participant is serial or parallel.
- Sorted out the different meaning of timestep and time window:
  - Renamed API function `isTimestepComplete` to `isTimeWindowComplete`.
  - Renamed C bindings function `precicec_isCouplingTimestepComplete` to `precicec_isTimeWindowComplete`.
  - Renamed `cplscheme` configuration option `timestep-length` to `time-window-size`.
  - Renamed `cplscheme` configuration option `max-timesteps` to `max-time-windows`.
  - Renamed `acceleration` configuration option `timesteps-reused` to `time-windows-reused`.
  - Renamed `acceleration` configuration option `reused-timesteps-at-restart` to `reused-time-windows-at-restart`.
  - Renamed `export` configuration option `timestep-interval` to `every-n-time-windows`.

## 1.6.1

- Fixed a platform-dependent issue with boost.

## 1.6.0

- Added CMake target to uninstall the project.
- Added `INFO` log output of the connection build-up during the instantiation. This dramatically simplifies the identification of connection-related problems.
- Added `precice::getVersionInformation()` to the interface. This returns a semicolon separated list of information containing the version, configuration, compiler and flags used to compile the library.
- Added additional checks for wrong user input in SolverInterface.
- Added basic index validation of user input to the SolverInterface
- Added events for mapping internals such as the computation of the mapping and the index-generation.
- Added more "pythonic" python bindings `precice_future`. More: [here](https://github.com/precice/precice/pull/379) and [here](https://github.com/precice/precice/issues/458)
- Added overrides for dependencies, allowing users to force CMake to use a given version. [Read more](https://github.com/precice/precice/wiki/Building:-Using-CMake#xsdk-compliance)
- Added support for the `CPP` and `CPPFLAGS` environment variables. [Read more](https://github.com/precice/precice/wiki/Building:-Using-CMake#xsdk-compliance)
- Added support for the xsdk default mode. [Read more](https://github.com/precice/precice/wiki/Building:-Using-CMake#xsdk-compliance)
- Added the CMake variable `PRECICE_CTEST_MPI_FLAGS` which can be used to pass additional flags to `mpirun`. A common use-case is to enable oversubscribing on machines with only 2 cores by passing `--oversubscribe`.
- Added the possibility to pass a custom `MPI_Comm` to preCICE via a new SolverInterface constructor. The passed communicator is then used as the global internal communicator. The user has to ensure that the mpi implementations of preCICE and the caller code are consistent and compatible.
- Added the result of `git describe --tags --dirty` to the library, which is now displayed during the configuration of preCICE. This allows quickly check what commit you are actually using and whether there were local changes.
- Added tests for python bindings using a mocked preCICE C++ Interface. [Read more](https://github.com/precice/precice/pull/420)
- Added validation of dependencies in the CMake script.
- Changed the connection publishing to a hash-based approach. This method is faster, more robust and NFS-friendly. The files are rooted in the folder `precice-run`, deleting this folder resolves most connection problems.
- Changed the errors in the configuration stage to throw `std::runtime_error`.
- Changed the id management from `std::set` to `boost::container::flat_set`. This reduces the peak memory consumption of meshes by about 9%.
- Changed the log output for filtered vertices. Log level `INFO` prints the total number of filtered vertices. The detailed information was moved to `DEBUG`.
- Changed the scope of the generated `version.[ch]pp` files to `precice/impl`. This prevents collisions with other random version headers.
- Changed and optimized the index generation for an individual speed up of 10-20x.
- Changed and optimized some complexity issues in the generation of communication maps for an individual speedup of around 1000x for bigger meshes.
- Changed and optimized some complexity issues in `NearestNeighbour::tagFirstRound` leading to an individual speedup of around 100x.
- Deprecated the python bindings `precice`, which will be removed in preCICE Version 2.0.0. If you still want to use them, please install `precice` and `precice_future`. Our recommendation, if you want to use the new bindings: Use `import precice_future as precice`.
- Fixed a bug in the python bindings which ignored the memory layout of numpy array. We now use `numpy.ascontiguousarray` to guarantee a C-compatible layout of data structures. [Read more](https://github.com/precice/precice/pull/448)
- Fixed a major memory issue due to excessive logger instantiations. This reduced peak memory consumption and allocation count of meshes by 50%.
- Fixed compatibility with Eigen versions `>3.3.7`.
- Fixed macro namespace by using the `PRECICE_` prefix. This prevents collisions with foreign macros.
- Fixed the index-based version of the nearest projection mapping and reintegrated it.
- Fixed wrong user input in Triangle creation to trigger assertions instead of a comprehensive user error.
- Fixed a bug in MPI initialization when using PETRBF mappings.
- Fixed a warning message in `ReceivedPartition` displaying a wrong vertex count.

## 1.5.2

- Fixed faulty communication in `master:sockets` leading to crashes.

## 1.5.1

- Fixed the exposure of `boost::asio` implementation details leading to version incompatibilities
- Fixed the CMake linkage of Boost to become compatible with the generated CMake config of Boost 1.70.0

## 1.5.0

- Added CMake alias `precice::precice` which mimics the namespaced library name after calling `find_package(precice)`. This allows seamless use as a subproject.
- Added state-awareness to the SolverInterface to keep track of the Mesh states. This results more informative error messages for mesh related functions.
- Added faster `PetRadialBasisFctMapping` tagging using RTree queries.
- Added potentially missing signal `SIGXCPU`.
- Added sanitization of IDs to the SolverInterface and harden it against misconfiguration using a consistent mechanism to express requirements.
- Added support for using preCICE directly from the binary tree using `precice_DIR` in other projects.
- Added target `test_install`, which can be used to test the install. It configures builds and runs the cpp solverdummy using the installed preCICE library.
- Added the ability to send a mesh to multiple participants.
- Added the attribute `Participant` to the logger, which is now available in format and filter expressions.
- Added the automatic creation of target directories, when exporting meshes with `<export:vtk ...>`.
- Added the data name to the log output of convergence measures.
- Added the generation of a pkg-config file for preCICE, which is usable after installing.
- Added the generation of code coverage reports.
- Changed the tested build-system in CI from SCons to CMake.
- Fixed a **critical bug** in the nearest projection mapping rooting in the variant index tree. We decided to roll-back to the old implementation until there is a fix.
- Fixed PETSc complaining about unused arguments when running tests.
- Fixed `bindings/fortran` to use `precice::SolverInterface` instead of `precice::impl::SolverInterfaceImpl`
- Fixed a potential integer overflow when parsing huge numbers in the configuration code.
- Fixed a potential overflow in the computation of the side length of bounding boxes.
- Fixed bug in preconditioner when we have convergence in the first iteration.
- Fixed integration tests resetting the log-level by calling `configure`.
- Fixed missing includes for Eigen
- Fixed overuse of `std::endl` in code, file generation should be noticeably faster now.
- Fixed the `const`ness of passed pointer arguments as well as some member functions of the `SolverInterface`.
- Fixed the default safety factor for nearest projection mapping by increasing it from 0.1 to 0.5.
- Fixed Python bindings now offer all preCICE API functions (except `get_mesh_handle`)
- Fixed CMake now properly detects PETSc on systems with a directory-scoped MPI installation (`include/openmpi/mpi.h`) when the compiler was not set to the MPI compiler wrapper.
- Refactored `Mesh::computeState` into two logical units, improving the trace information.
- Refactored `computePartitions` using the stl.
- Refactored logger which dramatically reduces compile times.
- Refactored the M2N handling code, which improves the available trace information for issues related to the M2N initialization.
- Removed [`#395`](https://github.com/precice/precice/pull/395) API functions `nameConfiguration(), dataDisplacements(), dataForces(), dataVelocities(), actionPlotOutput(), exportVTK(), exportAll()` from CPP API and all language bindings.
- Removed ancient `PRECIE_NO_SOCKETS` definition.
- Removed dead configuration code.
- Sending data between participants is now fully asynchronous. This is relevant in one-way coupling scenarios, where the sending participant doesn't need to wait for the receiving one.

## 1.4.1

- Bug in re-partitioning fixed, occured for OpenFOAM and empty ranks in parallel.

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
