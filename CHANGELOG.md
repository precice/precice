# preCICE Change Log

All notable changes to this project will be documented in this file. For future plans, see our [Roadmap](https://www.precice.org/fundamentals-roadmap.html).

## 3.3.1

- Added CMake option `Boost_USE_STATIC_LIBS` defaulting to OFF which forces Boost to find shared library by default. (https://github.com/precice/precice/pull/2415)
- Fixed compilation error on non-POSIX systems due to unused inclusion of `<string.h>` (https://github.com/precice/precice/pull/2424)
- Fixed compilation with C++20 required by Kokkos 5.0.0 (https://github.com/precice/precice/pull/2410)
- Fixed incorrect macro for Windows detection. (https://github.com/precice/precice/pull/2424)
- Fixed incorrect rank event output when nesting user-defined profiling sections. (https://github.com/precice/precice/pull/2393)
- Fixed logging of filtered vertices for direct mesh access. (https://github.com/precice/precice/pull/2427)
- Fixed preprocessor if inside function-like macro invocation leading to errors in MSVC. (https://github.com/precice/precice/pull/2424)

## 3.3.0

- Added CMake option `PRECICE_FEATURE_PROFILING_COMPRESSION` to optionally enable LZMA compression of the textual records of the profiling files. (https://github.com/precice/precice/pull/2328)
- Added a clear error message when calling `advance(getMaxTimeStepSize())` on the participant that is supposed to provide the time-window size due to `<time-window-size method="first-participant"/>`. (https://github.com/precice/precice/pull/2347)
- Added all coupling partners of `coupling-scheme:multi` to info statement in compositional coupling schemes. (https://github.com/precice/precice/pull/2353)
- Added detailed CMake, compiler and linker information to the CMake configuration log. (https://github.com/precice/precice/pull/2282)
- Added executables `precice-version`, `precice-config-validate`, and `precice-config-doc`. (https://github.com/precice/precice/pull/2312)
- Added function to create a participant from Fortran with a custom communicator. (https://github.com/precice/precice/pull/2363)
- Added implicit socket-based intra-comms for parallel participants in case preCICE is compiled without MPI. (https://github.com/precice/precice/pull/2336)
- Added missing libraries for compilation on Windows. (https://github.com/precice/precice/pull/2280)
- Added one retry in case opening of connection info files fails during communication buildup. (https://github.com/precice/precice/pull/2276)
- Added option `<export:X update-series="true"/>` to enable keeping the series files of exporters updated after each time window instead of generating them at the end of the simulation. (https://github.com/precice/precice/pull/2272)
- Added support for using `setMeshQuad` and `setMeshQuads` with 2D meshes. (https://github.com/precice/precice/pull/2364)
- Added support for version `5.0.0` of Eigen. (https://github.com/precice/precice/pull/2366)
- Added warning when an MPI communicator is passed to the Participant constructor, but preCICE was build without MPI. (https://github.com/precice/precice/pull/2365)
- Changed Eigen detection to only using the CMake configuration files provided by the Eigen installation. (https://github.com/precice/precice/pull/2366)
- Changed default profiling `mode="fundamental"` to only include steering methods. Use the new `mode="api"` to profile all API methods. (https://github.com/precice/precice/pull/2320)
- Changed output of coupling scheme state in advance to be more compact. (https://github.com/precice/precice/pull/2325)
- Changed profiling output format to a text based format, which is more space efficient than the JSON format. (https://github.com/precice/precice/pull/2328)
- Deprecated the `precice-profiling` script in favor of the `precice-profiling` python package or the `precice-cli`. (https://github.com/precice/precice/pull/2317)
- Deprecated the executable `precice-tools` in favor of `precice-version`, `precice-config-validate`, and `precice-config-doc`. (https://github.com/precice/precice/pull/2312)
- Fixed Kokkos deallocation when solver initializes Kokkos instead of preCICE. (https://github.com/precice/precice/pull/2263)
- Fixed QR3 filter of acceleration not falling back to QR2 for non-constant pre-scaling coefficients. (https://github.com/precice/precice/pull/2370)
- Fixed `precice-profiling merge` mistakenly detecting non-precice json files such as solver or adapter configuration files. (https://github.com/precice/precice/pull/2298)
- Fixed `setMeshTriangle` still creating edges instead of relying on the mesh preprocessing to generate required edges. (https://github.com/precice/precice/pull/2302)
- Fixed `~Participant()` calling `std::abort()`, when coupling scheme expects a requirement. (https://github.com/precice/precice/pull/2352)
- Fixed a bug in the Fortran bindings signature of `precicef_get_mesh_vertex_ids_and_coordinates_`. (https://github.com/precice/precice/pull/2316)
- Fixed compatibility with boost.asio `1.89.0`. (https://github.com/precice/precice/pull/2351 and https://github.com/precice/precice/pull/2354)
- Fixed coupling state info statement of compositional coupling schemes using newlines. (https://github.com/precice/precice/pull/2353)
- Fixed error in configuration for a convergence measure of a multi coupling scheme using data that was exchanged, but not with the controller. (https://github.com/precice/precice/pull/2287)
- Fixed incorrect declaration for Fortran `precicef_set_mesh_triangles_`. (https://github.com/precice/precice/pull/2293)
- Fixed rounding errors computing time-window size to send to the second participant using `<time-window-size method="first-participant" />`. (https://github.com/precice/precice/pull/2344)
- Fixed the use of the deprecated libxml2 header `libxml/SAX.h`, which leads to compile errors for libxml2 `2.14.0` and newer. (https://github.com/precice/precice/pull/2306)
- Fixed user-defined `<intra-comm... />` leading to error when running participant in serial. (https://github.com/precice/precice/pull/2337)
- Improved checks when using incorrectly configured JIT mappings. (https://github.com/precice/precice/pull/2309)
- Improved error messages when requirement functions haven't been called when expected, namely `requiresInitialData()`, `requiresWritingCheckpoint()`, and `requiresReadingCheckpoint()`. (https://github.com/precice/precice/pull/2350)
- Migrated `precice-profiling merge` to the sqlite intermediate output. (https://github.com/precice/precice/pull/2343)
- Moved the generation of the test list to a separate target `precice-test-list`. Use `make testprecice precice-test-list` if you only want to run the tests. (https://github.com/precice/precice/pull/2310)
- Removed all commands from `precice-profiling` except `merge`. (https://github.com/precice/precice/pull/2343)
- Removed support for configuring preCICE as static library, which has been broken since version 3.0.0. (https://github.com/precice/precice/pull/2300)

## 3.2.0

- Added initial data for time=0 in watch integrals and watch points. (https://github.com/precice/precice/pull/2181)
- Added Paraview series files to `vtk`, `vtu`, and `vtp` exporters when exporting time windows. This makes loading exports easier and assigns the correct time stamps to them. (https://github.com/precice/precice/pull/2055)
- Added QR3 filter and use it in default IQN-ILS and IQN-IMVJ configurations (see https://doi.org/10.3390/mca27030040). (https://github.com/precice/precice/pull/1493)
- Added a check to ensure coupled participants use configuration files with the same content. (https://github.com/precice/precice/pull/2217)
- Added a clearer error message when attempting to use a pre-version 3 configuration file. (https://github.com/precice/precice/pull/2007)
- Added a configuration check for missing data samples, which may be caused by missing mappings or `write-data` tags. (https://github.com/precice/precice/pull/2044)
- Added a log info statement showing the amount of selected elements when computing a nearest-projection mapping. (https://github.com/precice/precice/pull/2184)
- Added an error if participants using different preCICE versions attempt coupling. (https://github.com/precice/precice/pull/2202)
- Added an error message for data that is read but never mapped or exchanged. (https://github.com/precice/precice/pull/2064)
- Added current working directory of participant to the startup output of preCICE. (https://github.com/precice/precice/pull/2170)
- Added default initial-relaxation factor of 0.5 to the Aitken acceleration. (https://github.com/precice/precice/pull/2112)
- Added descriptive error when meshes are swapped in read/write mappings. (https://github.com/precice/precice/pull/2124)
- Added error in case of missing samples in origin of read mappings. (https://github.com/precice/precice/pull/2024)
- Added error message when defining multiple acceleration schemes per coupling-scheme. (https://github.com/precice/precice/pull/2111)
- Added experimental support for just-in-time mapping with NN and PU-RBF via the new API functions `mapAndReadData` and `writeAndMapData`. (https://github.com/precice/precice/pull/2099)
- Added experimental support for remeshing at runtime of parallel coupling schemes. Adding `<precice-configuration experimental="true" allow-remeshing="true" />` in the configuration allows calling `resetMesh()` for any provided mesh in the first iteration of a time-window. (https://github.com/precice/precice/pull/2070)
- Added header `precice/Exceptions.hpp` which provides the preCICE exception type `precice::Error`. It is included via `precice/precice.hpp`. (https://github.com/precice/precice/pull/2065)
- Added profiling event for handling exports, which can be a substantial part of `advance`. (https://github.com/precice/precice/pull/2041)
- Added profiling events to most API calls. (https://github.com/precice/precice/pull/2130)
- Added support for Intel MPI version 2021.13. (https://github.com/precice/precice/pull/1826)
- Added support for Intel oneAPI C/DPC++/C++/Fortran Compilers icx, icpx, and ifx version 2024.2. (https://github.com/precice/precice/pull/1826)
- Added support for Kokkos configured with parallel host and device executors. (https://github.com/precice/precice/pull/2169)
- Added support for secondary data in IQN-IMVJ method. (https://github.com/precice/precice/pull/2010)
- Added support for user-defined profiling sections using `startProfilingSection(name)` and `stopLastProfilingSection()`. This is an easy way measure and display profiling of adapter code in the context of the whole participant or simulation. (https://github.com/precice/precice/pull/1657)
- Added support for using unity builds via the CMake option `PRECICE_BUILD_UNITY`, which is enabled by default. (https://github.com/precice/precice/pull/2068)
- Added support for waveform iterations in Aitken acceleration (https://github.com/precice/precice/pull/2163)
- Added support for waveform iterations in Quasi-Newton acceleration (https://github.com/precice/precice/pull/2005)
- Added the synchronization of ranks to the profiling data if `<profiling synchronize="true" />` (https://github.com/precice/precice/pull/2019)
- Added warning when user set larger `max-used-iterations` in acceleration than `max-iterations` in coupling-scheme when no previous time window is reused in acceleration. (https://github.com/precice/precice/pull/2201)
- Added warning when using `<m2n:mpi-multiple-ports />` with Intel MPI due to frequent hangs. (https://github.com/precice/precice/pull/1826)
- Added weight monitoring for residual-sum preconditioner and use it in default IQN-ILS and IQN-IMVJ configurations (see https://doi.org/10.3390/mca27030040). (https://github.com/precice/precice/pull/1493)
- Changed all exchanges of implicit coupling schemes to use `substeps="True"` by default. This will activate exchange of additional data from inside the coupling time window if subcycling is used. Even without subcycling additional data will be communicated (see https://github.com/precice/precice/issues/2226 for details). (https://github.com/precice/precice/pull/2220)
- Changed `m2n:sockets` and `intra-comm:sockets` to disable Nagle's algorithm using `TCP_NODELAY`, which improves data communication delay. This is especially noticeable when using substeps. (https://github.com/precice/precice/pull/2208)
- Changed `setMeshAccessRegion` to be callable once per mesh instead of once per participant. (https://github.com/precice/precice/pull/2072)
- Changed error when using `<time-window-size value="X" method="first-participant"/>` with any value other than -1 to a warning. (https://github.com/precice/precice/pull/2092)
- Changed exporters to skip exporting data of a mesh when the data correctly does not containing any samples. (https://github.com/precice/precice/pull/2132)
- Changed output of `precice-profiling analyze` to show nested events using indentation. (https://github.com/precice/precice/pull/2117)
- Changed scoping of events to be handled in post-processing. This leads to a consistent event naming. (https://github.com/precice/precice/pull/2066)
- Changed the C++ API to throw a `precice::Error` instead of aborting if preCICE detects an error. This allows adapters to clean up on error and makes preCICE easier to use via interactive sessions. (https://github.com/precice/precice/pull/2065)
- Changed the baseline from Ubuntu 20.04 LTS to Ubuntu 22.04 LTS, which lifts minimum required versions of dependencies to: CMake 3.22.1, GCC/stdlibc++ 11.2.0, Eigen 3.4.0, Boost 1.74.0, PETSc 3.15, Python 3.10.6, and Numpy 1.21.5. (https://github.com/precice/precice/pull/1962)
- Changed the behavior of the function `getMeshVertexIDsAndCoordinates`, when combining direct-mesh access with mappings, to return only the vertices within the defined access region (`setMeshAccessRegion`). Setting an access region is mandatory if no other mapping is defined. In parallel, it has been mandatory anyway. (https://github.com/precice/precice/pull/2099)
- Changed the default mode for IMVJ from no-restart to SVD-restart. (https://github.com/precice/precice/pull/2136)
- Fixed a bug in sample retrieval of `readData` which leads to crashes for some combinations of time-window-sizes and time-windows. (https://github.com/precice/precice/pull/2191)
- Fixed a crash when meshes in a mapping neither received nor provided by a participant. (https://github.com/precice/precice/pull/2124)
- Fixed bug resulting in inconsistent watch points and integrals for different coupling schemes. (https://github.com/precice/precice/pull/2188)
- Fixed bugs due to size of gradients if gradient mapping is used and substeps are exchanged. (https://github.com/precice/precice/pull/2223)
- Fixed confusing output of the final advance at the end of the simulation. (https://github.com/precice/precice/pull/2182)
- Fixed crash in finalize after an inter-participant communication fails to connect. (https://github.com/precice/precice/pull/2141)
- Fixed direct mesh access (api access) to correctly read and write data in serial runs, where vertices are filtered through the defined mesh access region. (https://github.com/precice/precice/pull/2224)
- Fixed incompatibility between Cuda and more recent boost.asio versions by using `BOOST_ASIO_USE_TS_EXECUTOR_AS_DEFAULT` if Cuda was enabled. (https://github.com/precice/precice/pull/2218)
- Fixed incorrect escapes in version information for some CMake versions (https://github.com/precice/precice/pull/2030)
- Fixed mismatched MPI communicator sizes to trigger an assertion before their sanitization check. (https://github.com/precice/precice/pull/2065)
- Fixed parallel repartitioning for PU-RBF to be more accurate. (https://github.com/precice/precice/pull/1912)
- Fixed the API functions `setMeshAccessRegion` and `getMeshVertexIDsAndCoordinates` to now require the flag `<receive-mesh ... api-access="true" />` in the configuration file. If not set, the function will throw an error. (https://github.com/precice/precice/pull/2099)
- Fixed unnecessary creation of empty watch point and watch integral files at non-owning participants. (https://github.com/precice/precice/pull/2126)
- Improved convergence log to display mesh and data name in convergence measures. (https://github.com/precice/precice/pull/2032)
- Improved memory footprint when handling non-interpolated data. (https://github.com/precice/precice/pull/2190)
- Improved the runtime of the creation of the Bspline interpolation, by switching the matrix A inside the Bspline class to a sparse matrix. (https://github.com/precice/precice/pull/2021)
- Improved time-interpolation to allow arbitrary high waveform degree. (https://github.com/precice/precice/pull/2039)
- Moved Kokkos initialization into solver class to prevent unnecessary initializations where they are not required. (https://github.com/precice/precice/pull/2080)
- Renamed the `direct-access` flag in `<receive-mesh direct-access="..." />` to `api-access` `<receive-mesh api-access="..." />` and deprecated `direct-access`. (https://github.com/precice/precice/pull/2099)
- Renamed the `develop` CMake preset to `development`. (https://github.com/precice/precice/pull/2241)
- Replaced Ginkgo's triangular solver by cublas and hipblas implementations. Using the Cuda/HIP QR decomposition requires now the Cuda/HIP blas library. (https://github.com/precice/precice/pull/2079)
- Speed-up parallel repartitioning for global RBF mappings. (https://github.com/precice/precice/pull/1912)
- Updated CMake to use `BoostConfig.cmake` instead of CMake's FindBoost module to find boost. (https://github.com/precice/precice/pull/2056)
- Updated Ginkgo integration with Kokkos for accelerated (GPU or OpenMP) RBF data mapping. Kokkos is now a mandatory dependency when using Ginkgo in preCICE. preCICE is now compatible with Ginkgo v1.8.0 or higher with Kokkos version 4.1 or higher. (https://github.com/precice/precice/pull/2051)
- Updated boost asio interface for compatibility with boost 1.86. (https://github.com/precice/precice/pull/2158)
- Updated pre-commit hook versions and applied corresponding changes. (https://github.com/precice/precice/pull/2244)

## 3.1.2

- Fixed incorrect handling of compositional coupling involving an implicit scheme. Explicit schemes now run after the implicit scheme has reached convergence, correctly receive data of the final iteration.

## 3.1.1

- Added missing checks for incorrect Participant names in M2N and coupling-scheme. (https://github.com/precice/precice/pull/1995)
- Fixed skipping initial mapping when data contains only zeros in parallel. (https://github.com/precice/precice/pull/1999)

## 3.1.0

- Added warnings when using invalid options inside log configuration files. (https://github.com/precice/precice/pull/1956)
- Changed the output of `precice-profiling analyze` to be sorted by event name. (https://github.com/precice/precice/pull/1953)
- Fixed a bug where reading from the end of the time window can trigger an assertion. (https://github.com/precice/precice/pull/1982)
- Fixed bug when using log configuration files with invalid options. (https://github.com/precice/precice/pull/1956)
- Fixed oversubscription errors when running tests. (https://github.com/precice/precice/pull/1960)
- Fixed too strict check on allowed mapping types when combining serial and parallel participants. (https://github.com/precice/precice/pull/1964)
- Improved numerical accuracy of time handling for simulations with many time steps or time windows, using separate Kahan sums for time-window start and window progress. (https://github.com/precice/precice/pull/1954)
- Replaced boost.filesystem with `std::filesystem`. (https://github.com/precice/precice/pull/1972)

## 3.0.0

- Added API method `getMaxTimeStepSize()`, replacing return values and simplifying usability of `advance(double dt)` and `initialize()`. (https://github.com/precice/precice/pull/1623)
- Added CMake presets for simpler configuration. (https://github.com/precice/precice/pull/1452)
- Added Wendland's compactly supported C2 and C4 RBF functions. (https://github.com/precice/precice/pull/1378)
- Added `<mapping:rbf>` alias defaulting to the `rbf-pum-direct`. (https://github.com/precice/precice/pull/1616)
- Added `precicef_get_version_information_` to Fortran bindings. (https://github.com/precice/precice/pull/1487)
- Added a `scaled-consistent-volume` mapping constraint in order to scale by volumetric primitives. (https://github.com/precice/precice/pull/1387)
- Added a mapping configuration subtag `executor` to specify where and how RBF data mappings should be computed: `cpu` (default), `cuda`, `hip`, or `openmp`. (https://github.com/precice/precice/pull/1630)
- Added a new mixed convergence measure `absolute-or-relative-convergence-measure`. (https://github.com/precice/precice/pull/1862)
- Added a warning for empty IQN matrix, but keep the simulation going. This allows to start from a zero initial state that only later changes (e.g., an opening valve in a flow simulation) (https://github.com/precice/precice/pull/1895)
- Added an optional preconditioner to the Aitken underrelaxation acceleration. (https://github.com/precice/precice/pull/1379)
- Added check for duplicated `<exchange/>` tags. (https://github.com/precice/precice/pull/1467)
- Added communication of substeps data for waveform coupling. A new attribute `substeps` in the `exchange` tag allows activating and deactivating this additional communication. `substeps` defaults to `false`, since not yet all acceleration schemes support this feature. This default may change once preCICE is fully ported to handling waveforms.  See https://github.com/precice/precice/issues/1825 for details. (https://github.com/precice/precice/pull/1664)
- Added compactly supported Wendland C8 basis function. (https://github.com/precice/precice/pull/1871)
- Added configuration option for participants waiting for each other in finalize: `<precice-configuration wait-in-finalize="1"/>`. (https://github.com/precice/precice/pull/1840)
- Added default values for the IQN-ILS and IQN-IMVJ configuration options `initial-relaxation`, `max-used-iterations`, `time-windows-reused`, `filter`, and `preconditioner`. (https://github.com/precice/precice/pull/1754)
- Added dependency versions to the CMake output. (https://github.com/precice/precice/pull/1430)
- Added descriptive errors when providing incorrect mesh or data names. (https://github.com/precice/precice/pull/1588)
- Added efficient bulk functions `setMeshEdges|Triangles|Quads|Tetrahedra` for setting connectivity based on an array of vertex IDs. (https://github.com/precice/precice/pull/1322)
- Added error message when attempting to initialize data on a directly-accessed mesh, which is not available before initialization. (https://github.com/precice/precice/pull/1592)
- Added error messages for inconsistent argument sizes. (https://github.com/precice/precice/pull/1641)
- Added experimental support for building preCICE with Ginkgo. (https://github.com/precice/precice/pull/1630)
- Added experimental, serial implementations for `axial-geometric-multiscale` and `radial-geometric-multiscale` consistent mapping types for 1D-3D geometric multiscale. (https://github.com/precice/precice/pull/1439)
- Added gradient data to constant and Aitken underrelaxation. (https://github.com/precice/precice/pull/1373)
- Added gradient related API functions to native Fortran and C bindings. (https://github.com/precice/precice/pull/1374)
- Added partition of unity RBF mapping as a more efficient alternative to existing RBF mapping for large cases. The mapping can be configured using `<mapping:rbf-pum-direct ... > <basis-function:... /> </mapping:rbf-pum-direct>`. (https://github.com/precice/precice/pull/1483)
- Added preprocessing of mesh connectivity, which removes duplicates and completes hierarchical connectivity, if required. (https://github.com/precice/precice/pull/1494)
- Added sample generation from write-data of substeps that result from subcycling, allowing higher-order B-spline interpolation. (https://github.com/precice/precice/pull/1671)
- Added support for a single participant to use multiple coupling schemes with different time-window-sizes. This currently doesn't work with actions. (https://github.com/precice/precice/pull/1705)
- Added support for higher-order B-spline interpolation (based on multiple samples in current window). (https://github.com/precice/precice/pull/1422)
- Added support for implicit coupling schemes without convergence measures, which always iterate `max-iterations` times per time window. Not defining convergence measures now requires `max-iterations` to be defined. (https://github.com/precice/precice/pull/1841)
- Added support for multiple data sets in the Aitken underrelaxation acceleration. (https://github.com/precice/precice/pull/1379)
- Added support for using constant underrelaxation with waveforms. (https://github.com/precice/precice/pull/1794)
- Added support for waveform interpolation with multi coupling. (https://github.com/precice/precice/pull/1445)
- Added support for waveform iterations for the serial-explicit coupling. The `second` participant can now use higher-order time interpolation. (https://github.com/precice/precice/pull/1602)
- Added support to include preCICE as CMake subproject. (https://github.com/precice/precice/pull/1615)
- Added tag `min-iterations` to coupling schemes, which allows defining a minimal amount of iterations before convergence measures are taken into account. (https://github.com/precice/precice/pull/1841)
- Added tool `precice-profiling` to post-process the generated profiling files. (https://github.com/precice/precice/pull/1419)
- Added warning if subcycling over multiple windows introduces time drift leading to very small maximum time step sizes. See #1866 for details. (https://github.com/precice/precice/pull/1867)
- Added warnings for `scaled-consistent` mappings unable to uphold their conservation constraint. (https://github.com/precice/precice/pull/1697)
- Added waveform interpolation also for serial coupling. A serial coupling scheme now allows interpolation for both participants. (https://github.com/precice/precice/pull/1352)
- Allowed `initialize=true` for send and write data of both participants in serial coupling schemes (consistently with parallel coupling schemes). (https://github.com/precice/precice/pull/1367)
- Changed data samples containing only 0 to be skipped when mapping samples during `initialize`. This prevents mapping samples when using `<exchange ... initialize="false" />`. (https://github.com/precice/precice/pull/1742)
- Changed API function `getMaxTimeStepSize()` to return `dt = 0.0` if called after `final advance(dt)` (https://github.com/precice/precice/pull/1906)
- Changed API function `readData(...)` such that it may only be called after `initialize()` (https://github.com/precice/precice/pull/1906)
- Changed API function `readData(...)` to only accept `relativeReadTime = 0.0`, if called after final `advance(dt)` (i.e. not `isCouplingOngoing()`) (https://github.com/precice/precice/pull/1906)
- Changed API function `writeData(...)` such that it may only be called before the final `advance(dt)` (i.e. not isCouplingOngoing()) (https://github.com/precice/precice/pull/1906)
- Changed API to guarantee the creation of hierarchical connectivity. For example, `setMeshTriangle` now ensures necessary edges exist internally. (https://github.com/precice/precice/pull/1322)
- Changed API to use names instead of ids for meshes and data. Functions now accept a mesh name instead of a mesh id, and a mesh and data name instead of a data id. (https://github.com/precice/precice/pull/1588)
- Changed C++ API from pointers to `precice::span`, which is forward-compatible with `std::span`. Every provided span contains an associated size, which is required to be correctly sized. (https://github.com/precice/precice/pull/1641)
- Changed CMake to hide non-API symbols by default, significantly reducing binary size. (https://github.com/precice/precice/pull/1498)
- Changed PETSc detection to rely solely on pkg-config. `PETSC_DIR` and `PETSC_ARCH` are no longer used to infer possible locations of the `PETSc.pc` file. (https://github.com/precice/precice/pull/1547)
- Changed `waveform-degree="1"` to be default and time interpolation to non-experimental (no need to use `experimental="true"` anymore). Note on explicit schemes: Depending on the scheme for some participants `waveform-degree="1"` technically leads to a zeroth order interpolation. (https://github.com/precice/precice/pull/1684)
- Changed argument `relativeReadTime` to be mandatory in `readData()`, use `0` to get the current point in time and `getMaxTimeStepSize()` to get the end of the time window. (https://github.com/precice/precice/pull/1622)
- Changed data and mesh names in solverdummy configuration to match naming conventions. (https://github.com/precice/precice/pull/1548)
- Changed default values of `rbf-pum` to 50 vertices per cluster and a relative overlap of 0.15. (https://github.com/precice/precice/pull/1872)
- Changed direct-mesh access to be non-experimental. (https://github.com/precice/precice/pull/1740)
- Changed error to warning when trying to add a zero-value column to `V` matrix in QN-acceleration. (https://github.com/precice/precice/pull/1863)
- Changed execution order of schemes per participant to explicit schemes in configured order, followed by one optional implicit scheme. (https://github.com/precice/precice/pull/1462)
- Changed internal time handling to a Kahan accumulator, preventing issues combining `max-time` with a small `time-windows-size` or when subcycling. (https://github.com/precice/precice/pull/1933 https://github.com/precice/precice/pull/1934)
- Changed read mapping to be conditional during data initialization. This has a minor influence on some events being triggered and avoids performing unnecessary mappings. See #834. (https://github.com/precice/precice/pull/1404)
- Changed storage so that it does not build the B-spline interpolant when reading from existing timestamps. (https://github.com/precice/precice/pull/1837)
- Changed the data API to taking `precice::span`, reducing the API to `readData`, `writeData`, and `writeGradientData`. Sizes are checked at the API boundary and are required to be consistent. (https://github.com/precice/precice/pull/1636)
- Changed the initial guess for global-iterative mappings to reset to zero after each completed time window. (https://github.com/precice/precice/pull/1660)
- Changed the minimum supported C++ version from 14 to 17. The API only requires to be compiled with C++ 11. (https://github.com/precice/precice/pull/1413)
- Changed the partition-of-unity mapping to only create a single cluster, if the desired number of vertices per cluster is smaller than the local number of vertices. (https://github.com/precice/precice/pull/1736)
- Disabled the default PETSc signal handler in case preCICE needs to initialize PETSc. (https://github.com/precice/precice/pull/1434)
- Fixed *init* export written at the end of `initialize`. (https://github.com/precice/precice/pull/1796)
- Fixed CMake to correctly require Eigen version 3.3.7. (https://github.com/precice/precice/pull/1618)
- Fixed `nearest-projection` mapping behavior to map to the closest primitive instead of the highest-dimensional primitive available. (https://github.com/precice/precice/pull/1351)
- Fixed `precice-tools check` not displaying errors if the logger is disabled in the checked config. (https://github.com/precice/precice/pull/1477)
- Fixed a crash when using `precice-tools check` to check configurations that use PETSc-based RBF mappings. (https://github.com/precice/precice/pull/1578)
- Fixed bug for increasing time window size with participant first method and waveform relaxation by adding a corresponding check to prevent wrong usage. See https://github.com/precice/precice/pull/1770 for details. (https://github.com/precice/precice/pull/1789)
- Fixed deadlock in two-level initialization, which occurs if ranks share vertices, but only one of the ranks actually needs the vertices (or any) to compute the mapping. (https://github.com/precice/precice/pull/1847)
- Fixed missing error on mismatched closing XML tags (https://github.com/precice/precice/pull/1574)
- Fixed missing include path for C-bindings when using preCICE via pkg-config directly from the build directory. (https://github.com/precice/precice/pull/1931)
- Fixed preconditioners dividing by zero in case of coupling data always being zero and not changing over time. (https://github.com/precice/precice/pull/1910)
- Fixed remapping of a time windows initial data sample in implicit coupling schemes. (https://github.com/precice/precice/pull/1899)
- Fixed scaled-consistent mapping constraints with an output integral of zero to dividing by zero. (https://github.com/precice/precice/pull/1697)
- Fixed subcycling when combining an implicit coupling scheme with any other coupling scheme. (https://github.com/precice/precice/pull/1705)
- Fixed the Fortran API to use lowercase `id` in functions setting mesh vertices. (https://github.com/precice/precice/pull/1624)
- Fixed the faulty conflict while using `no-restart` mode for restart and `always-build-jacobian` mode at the same time in IQN-IMVJ acceleration. (https://github.com/precice/precice/pull/1856)
- Fixed the formula of the `residual-relative-convergence-measure` in the documentation. (https://github.com/precice/precice/pull/1769)
- Fixed the threshold value for truncation in SVD for IQN-IMVJ using SVD-restart. (https://github.com/precice/precice/pull/1815)
- Improved error messages of unknown tags and attributes in the configuration. (https://github.com/precice/precice/pull/1573)
- Improved parallel single-level initialization significantly by optimizing communication map creation. (https://github.com/precice/precice/pull/1830)
- Improved the computation of the owner rank for shared vertices between ranks in parallel computations, leading to a significant speed-up of the two-level initialization. (https://github.com/precice/precice/pull/1849)
- Improved the runtime of `readData` by implementing a cache for the computations in the BSpline interpolation. (https://github.com/precice/precice/pull/1765)
- Moved `precice::getVersionInformation()` to `precice/Tooling.hpp`. (https://github.com/precice/precice/pull/1473)
- Moved the configuration attribute `dimensions` from the `solver-interface` XML tag to the `mesh` tag. (https://github.com/precice/precice/pull/1742)
- Moved the configuration to enable synchronization from the root tag to profiling: `<profiling synchronize="true" />`. (https://github.com/precice/precice/pull/1743)
- Removed *final* export written in `finalize`. (https://github.com/precice/precice/pull/1796)
- Removed API function `initializeData()`, as data initialization is now performed in `initialize()`. Initial data has to be written before `initialize()`, directly after defining the mesh. (https://github.com/precice/precice/pull/1350)
- Removed API functions `hasMesh` and `hasData`. (https://github.com/precice/precice/pull/1741)
- Removed API functions `isReadDataAvailable()` and is `isWriteDataRequired()` as waveforms eliminate their use-cases. (https://github.com/precice/precice/pull/1362)
- Removed Broyden acceleration. (https://github.com/precice/precice/pull/1735)
- Removed `ComputeCurvatureAction` action. (https://github.com/precice/precice/pull/1614)
- Removed `getMeshID()`, `getDataID()`, and `getDataIDs()` as we directly use names now. (https://github.com/precice/precice/pull/1588)
- Removed `getMeshVertexIDsFromPositions` and `getMeshVertices`. (https://github.com/precice/precice/pull/1556)
- Removed `min-iteration-convergence-measure`, use the new `min-iterations` tag or `max-iterations` instead. (https://github.com/precice/precice/pull/1841)
- Removed `precice::constants` from the C++ API and the C and Fortran bindings. (https://github.com/precice/precice/pull/1487)
- Removed `valid-digits` from `<time-window-size ... />`. preCICE now always uses the internal numerical precision for checking if the end of the time window has been reached. (https://github.com/precice/precice/pull/1882)
- Removed action `ScaleByDtAction`. (https://github.com/precice/precice/pull/1403)
- Removed action timings `read-mapping-prior`, `write-mapping-prior`, and `on-time-window-complete-post`. (https://github.com/precice/precice/pull/1614)
- Removed callback functions `vertexCallback` and `postAction` from `PythonAction` interface. (https://github.com/precice/precice/pull/1614)
- Removed dead-axis option from rbf-pum mappings. (https://github.com/precice/precice/pull/1872)
- Removed deprecated API function `precicef_ongoing_()` in fortran bindings. (https://github.com/precice/precice/pull/1537)
- Removed deprecated API functions `mapWriteDataFrom()` and `mapReadDataTo()` (https://github.com/precice/precice/pull/1222)
- Removed deprecated XML attributes related to normals: `<mesh flip-normals="0">` and `<export:vtk normals="1"/>`. (https://github.com/precice/precice/pull/1475)
- Removed deprecated XML tags using terminology master and slave. (https://github.com/precice/precice/pull/1474)
- Removed deprecated action timings `regular-prior`, `regular-post`, `on-exchange-prior`, and `on-exchange-post`. (https://github.com/precice/precice/pull/1614)
- Removed deprecated aliases `<intra-comm:mpi-single />` and `<m2n:mpi-singleports />`. (https://github.com/precice/precice/pull/1781)
- Removed deprecated mapping timings, as there is no use-case anymore. (https://github.com/precice/precice/pull/1536)
- Removed edge IDs from the API: `setMeshEdge` now returns nothing, all `setMeshXWithEdges` variants were removed from the API, and `setMeshTriangle|Quad` now accepts vertex IDs instead of edge IDs. (https://github.com/precice/precice/pull/1322)
- Removed extrapolation feature from coupling schemes. Initial guess for acceleration is now always constant extrapolation (i.e. use converged results from end of last time window). (https://github.com/precice/precice/pull/1738)
- Removed internally shipped dependency `nlohmann/json`. (https://github.com/precice/precice/pull/1759)
- Removed second order extrapolation in acceleration scheme due to unclear use-case at high maintenance cost. Please contact us if you have a use-case that requires this feature. (https://github.com/precice/precice/pull/1503)
- Removed support for multiple implicit coupling-schemes per participant. Use `<coupling-scheme:multi>` instead. (https://github.com/precice/precice/pull/1462)
- Removed symbolic link `binprecice`. The executable is now called `precice-tools`. (https://github.com/precice/precice/pull/1595)
- Removed the API stubs for `hasToEvaluateSurrogateModel()` and `hasToEvaluateFineModel()`, which were part of the deprecated Manifold mapping. (https://github.com/precice/precice/pull/1472)
- Removed the XML configuration tag `<solver-interface>` and moved all attributes to the parent `<precice-configuration>` tag. (https://github.com/precice/precice/pull/1662)
- Removed the `preallocation` option in the `rbf-global-iterative` mapping configuration. (https://github.com/precice/precice/pull/1758)
- Removed the ability to install tests via `PRECICE_InstallTest`. (https://github.com/precice/precice/pull/1499)
- Removed the deprecated `vtk` exporter fallback for parallel participants. If a parallel participant will try to use the `vtk` exporter, it will cause an error and ask the user to change the exporter to a compatible one (like vtu). (https://github.com/precice/precice/pull/1784)
- Removed trailing space in watchpoint, convergence, and iteration log files. (https://github.com/precice/precice/pull/1620)
- Renamed CMake variables consistently. (https://github.com/precice/precice/pull/1744)
- Renamed `<mapping:rbf... use-qr-decomposition="true" />` to `<mapping:rbf-global-direct ... > <basis-function:... /> </mapping:rbf-global-direct>` (https://github.com/precice/precice/pull/1550)
- Renamed `SolverInterface` to `Participant` in the API. (https://github.com/precice/precice/pull/1643)
- Renamed `isGradientDataRequired()` to `requiresGradientDataFor()` (https://github.com/precice/precice/pull/1487)
- Renamed `isMeshConnectivityRequired()` to `requiresMeshConnectivityFor()` (https://github.com/precice/precice/pull/1487)
- Renamed `waveform-order` to `waveform-degree` and moved this attribute from `read-data` to `<data:scalar ... />`, respectively `<data:vector ... />`. (https://github.com/precice/precice/pull/1714)
- Renamed precice-events to precice-profiling in configuration, tests, and scripts (https://github.com/precice/precice/pull/1787)
- Renamed the `<m2n:... />` attributes `from` and `to` to `acceptor` and `connector`. (https://github.com/precice/precice/pull/1683)
- Renamed the `scaled-consistent` mapping constraint to `scaled-consistent-surface`. (https://github.com/precice/precice/pull/1387)
- Renamed the direct-mesh API function `getMeshVerticesAndIDs` to `getMeshVertexIDsAndCoordinates`. (https://github.com/precice/precice/pull/1740)
- Renamed the primary API header files `SolverInterface.hpp` to `precice.hpp`, `SolverInterfaceC.h` to `preciceC.h`, and `SolverInterfaceFortran.hpp` to `preciceFortran.hpp`. (https://github.com/precice/precice/pull/1654)
- Replaced `getDimensions()` with `getMeshDimensions(meshName)` and `getDataDimensions(meshName, dataName)`. (https://github.com/precice/precice/pull/1631)
- Replaced `isActionRequired()` and `markActionFulfilled()` with explicit calls `requiresInitialData()`, `requiresReadingCheckpoint()`, and `requiresWritingCheckpoint()`. Invoking these functions is a promise to preCICE that the caller fulfills the requirement. (https://github.com/precice/precice/pull/1487)
- Replaced all `<mapping:rbf... />` related tags. RBF mappings are now defined in terms of the applied solver (`<mapping:rbf-global-direct ...` or `<mapping:rbf-global-iterative`) and the applied basis function is a subtag of the solver. Users should use the additionally added auto selection of an appropriate solver as follows: `<mapping:rbf  ...> <basis-function:... /> </mapping:rbf>`. Example: `<mapping:compact-polynomial-c0 direction="read" from= ... support-radius="0.3" />` would become `<mapping:rbf  direction="read" from= ...> <basis-function:compact-polynomial-c0 support-radius="0.3" /> </mapping:rbf>` (https://github.com/precice/precice/pull/1550)
- Replaced tag `<use-mesh />` with dedicated tags for `<receive-mesh />` and `<provide-mesh />`. (https://github.com/precice/precice/pull/1450)
- Replaced the Eigen-based QR decomposition by a Cholesky decomposition for s.p.d. data mapping matrices, resulting in a signifacant speed-up (affecting Gaussian RBF, compact polynomial RBFs and inverse multiquadrics). (https://github.com/precice/precice/pull/1372)
- Replaced the internal profiling tool to improve error-robustness and reduce complexity. (https://github.com/precice/precice/pull/1419)
- Updated fmt to version 10.2.1 (https://github.com/precice/precice/pull/1928)

## 2.5.1

- Changed error to warning in case of an empty IQN matrix, keeping the simulation going. This allows to start from a zero initial state that only later changes (e.g., an opening valve in a flow simulation). (https://github.com/precice/precice/pull/1895)
- Changed error to warning when trying to add a zero-value column to `V` matrix in QN-acceleration (Backport of https://github.com/precice/precice/pull/1863)
- Fixed compatibility with libxml version 2.12.0. (https://github.com/precice/precice/pull/1886)
- Fixed missing include path for C-bindings when using preCICE via pkg-config directly from the build directory. (Backport of https://github.com/precice/precice/pull/1931)

## 2.5.0

- Added 3D support to Linear Cell Interpolation mapping (`<mapping:linear-cell-interpolation >/`) using tetrahedra. (https://github.com/precice/precice/pull/1337)
- Added C and Fortran bindings for setMeshTetrahedron (`precicec_setMeshTetrahedron` and `precicef_set_tetrahedron`). (https://github.com/precice/precice/pull/1382)
- Added a clang-tidy CI and enforced the updated clang-tidy format. (https://github.com/precice/precice/pull/999)
- Added a new interface file `precice/Version.h`, which provides version macros `PRECICE_VERSION_MAJOR`, `PRECICE_VERSION_MINOR`, `PRECICE_VERSION_PATCH` as well as the stringified version `PRECICE_VERSION`. Use the convenience macro `PRECICE_VERSION_GREATER_EQUAL(2,4,0)` to check for compatibility. `SolverInterface.hpp` and `SolverInterfaceC.h` now include `Version.h` to simplify macro usage on older version lacking the header. (https://github.com/precice/precice/pull/1392)
- Added checks during the establishment of inter-participant communication in order to safeguard against common issues. (https://github.com/precice/precice/pull/1290)
- Added checks of MPI ports functions and descriptive error messages. (https://github.com/precice/precice/pull/1292)
- Added default network interface for BSD-flavored operating systems. (https://github.com/precice/precice/pull/1332)
- Added pre-commit configuration, which drastically simplifies contributing. Required formatting tools will now be installed by pre-commit. (https://github.com/precice/precice/pull/1359)
- Added support for "Linear Cell Interpolation" data mapping in 2D. (https://github.com/precice/precice/pull/1297)
- Added success statement to `precice-tools check valid-config.xml`, which used to output nothing. ((https://github.com/precice/precice/pull/1406)
- Added support for polynomial="separate" and polynomial="off" for Eigen based RBF mappings (use-qr-decomposition="true"). (https://github.com/precice/precice/pull/1335)
- Added support for tetrahedral meshes with `setMeshTetrahedron`. (https://github.com/precice/precice/pull/1314)
- Added support for triangles in 2D meshes. (https://github.com/precice/precice/pull/1286)
- Added support in VTK and VTU exporters for exporting tetrahedral elements. (https://github.com/precice/precice/pull/1314)
- Added support to export gradient data in VTK exporting. Scalar gradients are written to a vector `<ScalarDataname>_gradient` (`v0_dx,v0_dy,v0_dz,v1_dx ...`), vector gradients are written to multiple vectors `<VectorDataname>_dx/dy/dz` (`v0x_dx,v0y_dx,v0z_dz, ...`). (https://github.com/precice/precice/pull/1315)
- Added support to export gradient data in XML based exporting, VTU, VTP. (https://github.com/precice/precice/pull/1340)
- Added tetrahedron communication. (https://github.com/precice/precice/pull/1325)
- Added the `PRECICE_BUILD_TOOLS` (ON by default) CMake variable in order to enable/disable the generation of the `precice-tools` executable. (https://github.com/precice/precice/pull/1344)
- Changed behavior of "Watchpoints" to look for cell interpolation instead of boundary interpolation, if possible. (https://github.com/precice/precice/pull/1361)
- Deprecated API functions `isReadDataAvailable()` and `isWriteDataRequired()` to simplify implementation of waveform iteration. (https://github.com/precice/precice/pull/1224)
- Fixed a bug leading to large whitespace sections in logs when using MPI ports. (https://github.com/precice/precice/pull/1292)
- Fixed data type of 'uint' to compatible 'unsigned int' type in RBF test. (https://github.com/precice/precice/pull/1298)
- Fixed erroneous behavior of first-participant timestepping in compositional coupling schemes.  (https://github.com/precice/precice/pull/1307)
- Fixed inconsistency between `writeVectorGradientData` and `writeBlockScalarGradientData` and cleaned-up doxygen documentation of gradient API functions. (https://github.com/precice/precice/pull/1302)
- Fixed the API function `isGradientDataRequired` to behave as documented. (https://github.com/precice/precice/pull/1295)
- Fixed the git revision detection picking up revisions of super-projects. (https://github.com/precice/precice/pull/1398)
- Fixed the triangle-to-point distance calculation leading to erroneous projections in corner cases of the nearest-projection mapping. (https://github.com/precice/precice/pull/1395)
- Improved efficienty of mesh filtering and communication. (https://github.com/precice/precice/pull/1311)
- Improved the basis function implementation of RBF mappings for faster evaluations. (https://github.com/precice/precice/pull/1338)
- Improved the memory footprint of triangles in meshes. (https://github.com/precice/precice/pull/1311)
- Improved the runtime of the Eigen based RBF matrix assembly (use-qr-decomposition="true"). (https://github.com/precice/precice/pull/1320)
- Passing `nullptr` to the `SolverInterface` as `communicator` is now forbidden. Use the `SolverInterface` constructor without the `communicator` argument, if you don't want to provide a custom MPI communicator. (https://github.com/precice/precice/pull/1261)
- Refactored RBF mappings to share a common base class (https://github.com/precice/precice/pull/1279)
- Refactored RBF system assembly and RBF system solving into a dedicated class. (https://github.com/precice/precice/pull/1319)
- Refactored mapping class interface into mapConsistent and mapConservative. (https://github.com/precice/precice/pull/1301)
- Removed argument 'rowsFirst' in `writeVectorGradientData` and `writeBlockVectorGradientData`. The functions have now the signature `writeVectorGradientData( int dataID, int valueIndex, const double *gradientValues)` with the gradient format `(vx_dx, vy_dx, vx_dy, vy_dy)`. `writeBlockVectorGradientData` allows to pass data for multiple vertices point-wise in the same format. (https://github.com/precice/precice/pull/1302)
- Removed restriction to only use 5 iterations in the `RS_LS` restart mode of the IQN-IMVJ acceleration method. (https://github.com/precice/precice/pull/1257)
- Removed the explicit gradient data flag (`<data: ... gradient="on" />`) for nearest-neighbor-gradient mapping. The gradient data requirement is now deduced automatically by preCICE. (https://github.com/precice/precice/pull/1371)
- Renamed preCICE's `master` branch to `main` branch. (https://github.com/precice/precice/pull/1385)
- Upgraded fmt to 9.0.0, json to 3.11.1 and tcbrindle span. (https://github.com/precice/precice/pull/1396)

## 2.4.0

- Added API methods for entering gradient data in the solver interface. Methods available are for scalar gradient data, vector gradient data, block scalar gradient data and block vector gradient data. (https://github.com/precice/precice/pull/1169)
- Added CMake options to enable debug logging, trace logging and assertions in release builds. (https://github.com/precice/precice/pull/1177)
- Added a pkg-config file for using preCICE directly from the build directory. (https://github.com/precice/precice/pull/1238)
- Added a warning if the relative convergence measure is set too low. (https://github.com/precice/precice/pull/1266)
- Added a warning when using `<export:vtk />` in a parallel participant. (https://github.com/precice/precice/pull/1136)
- Added experimental API functions for waveform interpolation. Currently restricted to parallel-implicit coupling. Refer to the [user documentation](https://precice.org/couple-your-code-waveform.html). (https://github.com/precice/precice/pull/1187)
- Added export to CSV using `<export:csv />`. (https://github.com/precice/precice/pull/1144)
- Added export to VTP using `<export:vtp />`. (https://github.com/precice/precice/pull/1137)
- Added export to VTU using `<export:vtu />`. (https://github.com/precice/precice/pull/1136)
- Added the `support-radius` as an additional configuration option to the Gaussian RBF mapping configuration. (https://github.com/precice/precice/pull/1163)
- Added the API header `precice/Tooling.hpp`, which includes utility functions. (https://github.com/precice/precice/pull/1122)
- Added the command `check` to `binprecice`, which checks a given configuration file for correctness. (https://github.com/precice/precice/pull/1132)
- Added the command `version` to `binprecice`, which prints the version string of the used preCICE library. (https://github.com/precice/precice/pull/1122)
- Added tooling to CMake, which simplifies contributing. (https://github.com/precice/precice/pull/1143)
- Changed baseline from Ubuntu 18.04 LTS to Ubuntu 20.04 LTS. preCICE now requires Boost version `1.71.0` and CMake version `3.16.3`. (https://github.com/precice/precice/pull/1259)
- Deprecated API functions `mapWriteDataFrom` and `mapReadDataTo`. Note: compiling the tests will trigger this warning. (https://github.com/precice/precice/pull/859)
- Fixed a bug in ID management, which lead to a crash when reconstructing a SolverInterface. (https://github.com/precice/precice/pull/1190)
- Fixed a bug in the socket communication back-end, which occasionally lead to crashes. (https://github.com/precice/precice/pull/1262)
- Fixed a bug which crashed preCICE when explicitly enforcing gather scatter communication using `<m2n:X enforce-gather-scatter=1 />`. (https://github.com/precice/precice/pull/1268)
- Fixed an incompatibility with Boost 1.79.0 (https://github.com/precice/precice/pull/1250)
- Fixed missing data pieces in the master pvtu file for empty vertex distributions (e.g. a mesh is not exchanged) (https://github.com/precice/precice/pull/1233)
- Fixed missing edge connectivity in VTK exports. (https://github.com/precice/precice/pull/1127)
- Fixed the `format-all-docker` script to change of reformatted files to root. (https://github.com/precice/precice/pull/1210)
- Fixed the names of experimental direct-access API in the Fortran bindings to `precicef_set_mesh_access_region_` and `precicef_get_mesh_vertices_and_IDs_`. (https://github.com/precice/precice/pull/1129)
- Fixed the naming scheme of parallel VTU files (`.pvtu` and `.vtu` pieces). Paraview now correctly detects exports of multiple time-steps as a time series. (https://github.com/precice/precice/pull/1126)
- Fixed the wrong data assignment for multiple read mappings to the same mesh or multiple write mappings from the same mesh. (https://github.com/precice/precice/pull/1267)
- Improved the efficiency of vertex to edge projections. (https://github.com/precice/precice/pull/1226)
- Migrated to the [FindPython3](https://cmake.org/cmake/help/v3.16/module/FindPython3.html) module for locating Python3 and NumPy. (https://github.com/precice/precice/pull/1263)
- Migrated to using the target provided by [FindLibXml2](https://cmake.org/cmake/help/v3.16/module/FindLibXml2.html) (https://github.com/precice/precice/pull/1263)
- Renamed `binprecice` to `precice-tools` and added a symbolic link for backwards-compatibility. The link will be removed in preCICE version 3.0.0. (https://github.com/precice/precice/pull/1175)
- Removed meshName from solverdummies input parameters. (https://github.com/precice/precice/pull/1256)
- Renamed all instances of master to primary and of slave to secondary, along with their plural and joint forms. (https://github.com/precice/precice/pull/1253)
- Renamed class `MasterSlave` to `IntraComm` and changed all corresponding names. Changed the XML tag `master` to `intra-comm` and changed the XML attributes `on-slaves` and `on-master` to `on-secondary-ranks` and `on-primary-rank` respectively. The tag `master` and attributes `on-slaves` and `on-master` are deprecated. (https://github.com/precice/precice/pull/1274)
- Upgraded dependencies fmtlib to 8.1.1, json to 3.10.5 and tcbrindle span. (https://github.com/precice/precice/pull/1241)

## 2.3.0

- Added `isMeshConnectivityRequired(meshID)` to the SolverInterface API. This is useful to generate connectivity information only if required by preCICE.
- Added a configuration option to enable experimental API functions.
- Added bindings to doxygen.
- Added check of user-written data to be finite.
- Added convergence information for linear solvers of PETSc-based RBF mappings. Reports of converged solvers are logged as DEBUG, stopped solvers as WARNING and diverged solvers as ERROR messages.
- Added experimental support for direct mesh access of meshes.
- Added fmtlib for more readable and maintainable message formatting.
- Added support for 3D meshes in `ScaleByAreaAction`.
- Added the `Rank` of partitions in VTK exports.
- Added the associated data name when printing convergence information.
- Added the constraint `scaled-consistent` to mapping methods, which scales a consistent mapping such that the surface integral is equal on both sides of the interface.
- Added the experimental API functions `setMeshAccessRegion` and `getMeshVerticesAndIDs`.
- Added the xml configuration option '<use-mesh ... direct-access=true />' for this feature.
- Adopted the Contributor Covenant code of conduct (see `CODE_OF_CONDUCT.md`).
- Change `m2n:mpi` to use the more efficient single-ports implementation. To use the old implementation, use `m2n:mpi-mulitple-ports`.
- Changed edge and triangle normals to be computed on-demand.
- Changed formatting scripts to use `git ls-files`. This simplifies handling of edge cases.
- Changed linear solvers of PETSc-based RBF mappings now issue a warning when they reach the maximum iterations. This used to result in an ERROR.
- Changed the convergence measures INFO print to a fixed width scientific format. (issue #975)
- Changed the default build type of the library from `static` to `shared`.
- Changed the precice log files to print a fixed width scientific notation using 8 digits after the comma for floating point numbers and a fixed width for integer values. (issue #975)
- Changed the watch-point calculation algorithm to use already built index trees.
- Deprecated the `flip-normals` attribute of meshes. This is not functional anymore.
- Deprecated vertex normals in the `vertexCallback()` of the python actions. preCICE will pass `None` if the normal parameter is defined.
- Fixed FindNumPy using fallback of `find_path`.
- Fixed FindPETSc to also find PETSc if the include directory is in CPATH, which occurred with some versions of `pkg-config`.
- Fixed socket communication to fail without a network connection.
- Fixed wrong error in tightly converging QN coupling. (Issue #976)
- Implement a simplified interface to query index trees of meshes.
- Improved memory consumption and data locality of vertices.
- Improved mesh memory usage of vertices by ~41% and of connectivity by ~55%.
- Improved usabilty and readability of logging macros and assertions.
- Introduced package `time` and moved functionality related to extrapolation there.
- Migrated `PRECICE_XXX` macros to fmtlib.
- Modifying release PR templates according to restructured tutorials.
- Removed limitation of multi-coupling scheme for coupling topologies without a central participant. A non-centric controlling participant still needs to run in serial.
- Removed obsolete `query` package functionality. Interpolation is refactored into `mapping` package.
- Removed the FASTEST fortran bindings.
- Removed vertex normals from the vtk exporters.
- Replaced `prettyprint` with `fmtlib`.

## 2.2.1

- Fixed a bug leading to a freeze when using `sync-mode` with lazy indexing.
- Fixed empty received partitions for filtering on slaves.
- Fixed gather-scatter communication deadlock with empty master ranks.

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
- Enabled RBF-based tests in partition unit-tests and serial integration tests.
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
- Refactored `mesh::BoundingBox` into separate class.
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
- Moved Fortran 2003 bindings (`src/precice/bindings/f2003`) and solverdummy (`tools/solverdummy/f2003`) to a separate repository.
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

- Bug in re-partitioning fixed, occurred for OpenFOAM and empty ranks in parallel.

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
  - The helper scripts are now placed in the directory `tools/conda_building`. All the terms referring to `Anaconda` have been changed to `Conda`.
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
  - Experimental support for building with Conda, see `tools/anaconda_building`
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
- Add contributor guidelines.


## 1.0.3
- Fix compilation for boost 1.66, see issue #93.

## 1.0.2
- Fix bug in the mesh repartitioning for plane-like coupling interfaces and small gaps between both sides.

## 1.0.1
- Fix compilation issue with the python interface.
