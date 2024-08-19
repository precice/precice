cmake_minimum_required (VERSION 3.22.1)
project(preCICE VERSION 3.1.2 LANGUAGES CXX)
set(preCICE_SOVERSION ${preCICE_VERSION_MAJOR})

#
# Overview of this configuration
#
# PREAMBLE
# Setup Options
# Find Mandatory Dependencies
# Find Configurable Dependencies
# Configuration of Target precice
# Configuration of Target precice-tools
# Configuration of Target testprecice
# Install Targets for precice
# CPack
#


#
# PREAMBLE
#

# Make our cmake modules visible to CMake
list (APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
list (APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules")

include(CheckSTL)
include(CopyTargetProperty)
include(XSDKMacros)
include(Validation)

# CMake Policies

# CMP0074: find_package() uses <PackageName>_ROOT variables.
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()
# CMP0075 Include file check macros honor CMAKE_REQUIRED_LIBRARIES
if(POLICY CMP0075)
  cmake_policy(SET CMP0075 NEW)
endif()

#
# Setup Options
#

if(NOT CMAKE_BUILD_TYPE)
  message(WARNING "You did not specify a CMAKE_BUILD_TYPE.

  We will assume you asked for a Debug build.")
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "The type of this build" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

option(PRECICE_FEATURE_MPI_COMMUNICATION "Enables MPI-based communication and running coupling tests." ON)
option(PRECICE_FEATURE_PETSC_MAPPING "Enable use of the PETSc linear algebra library for global-iterative RBF data mappings." ON)
option(PRECICE_FEATURE_PYTHON_ACTIONS "Enable Python support for preCICE actions." ON)
option(PRECICE_CONFIGURE_PACKAGE_GENERATION "Configure package generation." ON)
option(PRECICE_FEATURE_GINKGO_MAPPING "Enable use of Ginkgo for accelerated data mapping." OFF)
option(BUILD_SHARED_LIBS "Build shared instead of static libraries" ON)
option(BUILD_TESTING "Build tests" ON)
option(PRECICE_ALWAYS_VALIDATE_LIBS "Validate libraries even after the validatation succeeded." OFF)
option(PRECICE_BINDINGS_C "Enable the native C bindings" ON)
option(PRECICE_BINDINGS_FORTRAN "Enable the native Fortran bindings" ON)
option(PRECICE_BUILD_TOOLS "Build the \"precice-tools\" executable" ON)
option(PRECICE_BUILD_UNITY "Use unity builds in preCICE" ON)
option(PRECICE_FEATURE_LIBBACKTRACE_STACKTRACES "Enable libbacktrace for stacktrace generation." OFF)

option(CMAKE_INTERPROCEDURAL_OPTIMIZATION "Enable interprocedural optimization for all targets." OFF)

option(PRECICE_RELEASE_WITH_DEBUG_LOG "Enable debug logging in release builds" OFF)
option(PRECICE_RELEASE_WITH_TRACE_LOG "Enable trace logging in release builds" OFF)
option(PRECICE_RELEASE_WITH_ASSERTIONS "Enable assertions in release builds" OFF)

xsdk_tpl_option_override(PRECICE_FEATURE_MPI_COMMUNICATION TPL_ENABLE_MPI)
xsdk_tpl_option_override(PRECICE_FEATURE_PETSC_MAPPING TPL_ENABLE_PETSC)
xsdk_tpl_option_override(PRECICE_FEATURE_PYTHON_ACTIONS TPL_ENABLE_PYTHON)

if(PRECICE_FEATURE_PETSC_MAPPING AND NOT PRECICE_FEATURE_MPI_COMMUNICATION)
  message(FATAL_ERROR "Please enable MPI to use PETSC.")
endif()

set(PRECICE_CTEST_MPI_FLAGS "" CACHE STRING "Add additional flags to mpiexec for running tests via CTest.")

include(XSDKOptions)

# Print information about this configuration
include(PrintHelper)
print_configuration(
  ADDITIONAL
  "PRECICE_CONFIGURE_PACKAGE_GENERATION;Configure package generation"
  "PRECICE_CTEST_MPI_FLAGS;Additional CTest MPI Flags"
  )
print_empty()

include(FeatureSummary)
add_feature_info(PRECICE_FEATURE_MPI_COMMUNICATION PRECICE_FEATURE_MPI_COMMUNICATION
  "Enables the MPI communication back-end.

   This enables the MPI communication back-end which is highly recommended on multi-node systems.
   See the documentation of the CMake module FindMPI to control its functionality.

   This feature can be enabled/disabled by setting the PRECICE_FEATURE_MPI_COMMUNICATION CMake option.
")
add_feature_info(PRECICE_FEATURE_PETSC_MAPPING PRECICE_FEATURE_PETSC_MAPPING
  "Enables the PETSc-powered radial basic function mappings.

   This enables the PETSc-based variant of radial basis function mappings which can run in parallel,
   also across different compute nodes. This is highly recommended for large cases running in parallel.

   This feature can be enabled/disabled by setting the PRECICE_FEATURE_PETSC_MAPPING CMake option.
   Requires PRECICE_FEATURE_MPI_COMMUNICATION.
  ")

  add_feature_info(PRECICE_FEATURE_GINKGO_MAPPING PRECICE_FEATURE_GINKGO_MAPPING
  "Enables the Ginkgo backend which supports the usage of GPUs and OpenMP for (global) radial-basis function mappings.

   Currently, this enables radial-basis function mappings to leverage the compute power of Nvidia/AMD GPUs and OpenMP.
   It uses iterative solvers provided by Ginkgo or direct solvers offered by hip/cuSolver.

   This feature can be enabled/disabled by setting the PRECICE_FEATURE_GINKGO_MAPPING CMake option.
   Requires the Ginkgo library, Kokkos (v4.1 or higher) and Eigen 3.4.
  ")

  add_feature_info(PRECICE_FEATURE_PYTHON_ACTIONS PRECICE_FEATURE_PYTHON_ACTIONS
  "Enables the support for user-defined python actions.

   preCICE allows to manipulate coupling data at runtime using configurable actions.
   This feature enables the support for user-defined actions written in Python based on numpy.

   This feature can be enabled/disabled by setting the PRECICE_FEATURE_PYTHON_ACTIONS CMake option.
  ")
add_feature_info(PRECICE_BINDINGS_C PRECICE_BINDINGS_C
  "Enables the native C bindings.

   preCICE provides native bindings for C, which are compiled into the library.
   This feature enables the compilation and installation of the bindings into the library.
   Note that we strongly recommend to compile with C bindings enabled for compatibility reasons.

   This feature can be enabled/disabled by setting the PRECICE_BINDINGS_C CMake option.
  ")
add_feature_info(PRECICE_BINDINGS_FORTRAN PRECICE_BINDINGS_FORTRAN
  "Enables the native Fortran bindings.

   preCICE provides native bindings for Fortran, which are compiled into the library.
   This feature enables the compilation and installation of the bindings into the library.
   Note that we strongly recommend to compile with Fortran bindings enabled for compatibility reasons.

   This feature can be enabled/disabled by setting the PRECICE_BINDINGS_FORTRAN CMake option.
  ")
add_feature_info(PRECICE_BUILD_TOOLS PRECICE_BUILD_TOOLS
  "Build the \"precice-tools\" executable

   preCICE offers in addition to the core library a variety of tools in order to check configuration
   files or generate xml references. The tools are copmiled into an executable called \"precice-tools\".

   This feature can be enabled/disabled by setting the PRECICE_BUILD_TOOLS CMake option.
  ")
  add_feature_info(PRECICE_FEATURE_LIBBACKTRACE_STACKTRACES PRECICE_FEATURE_LIBBACKTRACE_STACKTRACES
  "Enables libbacktrace for stracktrace generation.

   preCICE assertions are helpful during development, however, the default stacktraces may not be sufficient.
   This feature enables the support for libbacktrace in boost.stacktrace.

   This feature can be enabled/disabled by setting the PRECICE_FEATURE_LIBBACKTRACE_STACKTRACES CMake option.
  ")


feature_summary(WHAT ENABLED_FEATURES  DESCRIPTION "=== ENABLED FEATURES ===" QUIET_ON_EMPTY)
feature_summary(WHAT DISABLED_FEATURES DESCRIPTION "=== DISABLED FEATURES ===" QUIET_ON_EMPTY)



#
# Find Mandatory Dependencies
#
print_section("DEPENDENCIES")

find_package (Threads REQUIRED)

if(TPL_ENABLE_BOOST)
  xsdk_tpl_require(BOOST BOOST_ROOT)
  # Use BOOST_ROOT to set the directory
  set(Boost_NO_SYSTEM_PATHS ON CACHE BOOL "" FORCE)
  unset(ENV{BOOST_ROOT})
endif()

if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()

find_package(Boost 1.74.0 REQUIRED CONFIG
  COMPONENTS log log_setup program_options system thread unit_test_framework
  )
mark_as_advanced(Boost_INCLUDE_DIR Boost_LOG_LIBRARY_RELEASE Boost_LOG_SETUP_LIBRARY_RELEASE Boost_PROGRAM_OPTIONS_LIBRARY_RELEASE Boost_SYSTEM_LIBRARY_RELEASE Boost_THREAD_LIBRARY_RELEASE Boost_UNIT_TEST_FRAMEWORK_LIBRARY_RELEASE)
mark_as_advanced(Boost_INCLUDE_DIR Boost_LOG_LIBRARY_DEBUG Boost_LOG_SETUP_LIBRARY_DEBUG Boost_PROGRAM_OPTIONS_LIBRARY_DEBUG Boost_SYSTEM_LIBRARY_DEBUG Boost_THREAD_LIBRARY_DEBUG Boost_UNIT_TEST_FRAMEWORK_LIBRARY_DEBUG)
message(STATUS "Found Boost ${Boost_VERSION}")

# Eigen
if(TPL_ENABLE_EIGEN3)
  # Use EIGEN3_ROOT to set the directory
  xsdk_tpl_require(EIGEN3 EIGEN3_INCLUDE_DIR)
endif()
find_package(Eigen3 3.4.0 REQUIRED)
message(STATUS "Found Eigen ${EIGEN3_VERSION}")
precice_validate_eigen()

# LibXML2
find_package(LibXml2 REQUIRED)
message(STATUS "Found LibXml2 ${LIBXML2_VERSION_STRING}")
precice_validate_libxml2()

# libfmt
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/fmt)

#
# Find Configurable Dependencies
#

# Option: PRECICE_FEATURE_MPI_COMMUNICATION
if (PRECICE_FEATURE_MPI_COMMUNICATION)
  set(MPI_DETERMINE_LIBRARY_VERSION ON) # Try to detect the library vendor

  find_package(MPI REQUIRED)
  message(STATUS "MPI Version: ${MPI_CXX_LIBRARY_VERSION_STRING}")
endif()

# Option: PETSC
if (PRECICE_FEATURE_PETSC_MAPPING)
  if (TPL_ENABLE_PETSC)
    xsdk_tpl_require(PETSC PETSC_DIR PETSC_ARCH)
    # PETSc detection uses primarily these ENVs
  endif()
  find_package(PETSc 3.15 REQUIRED)
  # No validation required as PETSc does this internally
  message(STATUS "Found PETSc ${PETSc_VERSION}")
else()
  message(STATUS "PETSc support disabled")
endif()

# Option: Ginkgo
if (PRECICE_FEATURE_GINKGO_MAPPING)
  find_package(Ginkgo 1.8.0 REQUIRED)
  message(STATUS "Found Ginkgo ${GINKGO_PROJECT_VERSION} at ${Ginkgo_DIR}")

  find_package(Kokkos 4.1 REQUIRED)
  message(STATUS "Found Kokkos version ${Kokkos_VERSION} at ${Kokkos_DIR}")

  include(CMakeDependentOption)
  cmake_dependent_option(PRECICE_WITH_CUDA "Enable Ginkgo Cuda support in preCICE" ON "PRECICE_FEATURE_GINKGO_MAPPING;Kokkos_ENABLE_CUDA" OFF)
  cmake_dependent_option(PRECICE_WITH_HIP "Enable Ginkgo HIP support in preCICE" ON "PRECICE_FEATURE_GINKGO_MAPPING;Kokkos_ENABLE_HIP" OFF)
  cmake_dependent_option(PRECICE_WITH_OPENMP "Enable Ginkgo OpenMP support in preCICE" ON "PRECICE_FEATURE_GINKGO_MAPPING;Kokkos_ENABLE_OPENMP" OFF)

  # We aim for the same configuration as Ginkgo itself and check immediately for device compiler
  # such that preciceCore inherits all relevant flags directly
  # 1. Check if Ginkgo was built with CUDA and if CUDA is available on the system
  if(PRECICE_WITH_CUDA)
    # Necessary for the unified kernels and the qr decomposition
    include(CheckLanguage)
    check_language(CUDA)
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    set(CMAKE_CUDA_FLAGS "${GINKGO_CUDA_ARCH_FLAGS} --expt-relaxed-constexpr")

    # Similar to Ginkgo, we don't let CMake deduce the architecture, but inherit the Ginkgo flags
    if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.18)
      set(CMAKE_CUDA_ARCHITECTURES OFF)
    endif()
  # Check the Kokkos version
  if (Kokkos_VERSION VERSION_LESS 4.4)
    set(KOKKOS_VERSION_BELOW_4_4 TRUE)
  else()
    set(KOKKOS_VERSION_BELOW_4_4 FALSE)
  endif()
else()
  set(KOKKOS_VERSION_BELOW_4_4 FALSE)
endif()

  # 2. Check if Ginkgo was built with HIP and if HIP is available on the system
  if(PRECICE_WITH_HIP)
    include(CheckLanguage)
    check_language(HIP)
    enable_language(HIP)
    set(CMAKE_HIP_STANDARD 17)
    set(CMAKE_HIP_STANDARD_REQUIRED ON)
    if(NOT DEFINED HIP_PATH)
     if(NOT DEFINED $ENV{HIP_PATH})
       set(HIP_PATH "/opt/rocm/hip" CACHE PATH "Path to which HIP has been installed")
       set(ENV{HIP_PATH} ${HIP_PATH})
     else()
       set(HIP_PATH $ENV{HIP_PATH} CACHE PATH "Path to which HIP has been installed")
     endif()
    endif()

    if(NOT DEFINED ROCM_PATH)
     set(ROCM_PATH "/opt/rocm" CACHE PATH "Path to which ROCM has been installed")
     set(ENV{ROCM_PATH} ${ROCM_PATH})
    endif()

    find_package(hipsolver REQUIRED)
  endif()

  if(PRECICE_WITH_OPENMP)
    find_package(OpenMP 3.0 REQUIRED)
  endif()
  # That's essentially all backends we use from Ginkgo. The reference executor and OMP doesn't need any special treatments and
  # DPCPP has not yet been tested/setup
else()
  message(STATUS "Ginkgo support disabled")
endif()

if (PRECICE_BUILD_TOOLS AND BUILD_TESTING)
  # We need this to run the script precice-profiling
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
  message(STATUS "Found Python ${Python3_VERSION}")
endif()

# Option Python
if (PRECICE_FEATURE_PYTHON_ACTIONS)
  find_package(Python3 REQUIRED COMPONENTS Development NumPy)
  message(STATUS "Found Python ${Python3_VERSION} with NumPy ${Python3_NumPy_VERSION}")
  precice_validate_libpython()
  precice_validate_numpy()
else()
  message(STATUS "Python support disabled")
endif()

# Option ENABLE_LIBBACKTRACE
# Note that the FindBacktrace module does not work reliably to find libbacktrace installations
if (PRECICE_FEATURE_LIBBACKTRACE_STACKTRACES)
  add_subdirectory(thirdparty/libbacktrace)
endif()

#
# Setup miscellaneous features
#

# Enable glibc++ assertions in debug builds
add_compile_definitions("$<$<CONFIG:DEBUG>:_GLIBCXX_ASSERTIONS>")


#
# Configuration of Target precice
#

print_empty()
print_section("TARGETS & PACKAGES")

# Core of preCICE bundling the sources and dependencies.
# The main library and tests link against this.
add_library(preciceCore OBJECT)
set_target_properties(preciceCore PROPERTIES
  CXX_STANDARD 17
  CXX_STANDARD_REQUIRED Yes
  CXX_EXTENSIONS No
  POSITION_INDEPENDENT_CODE YES
  CXX_VISIBILITY_PRESET default
  VISIBILITY_INLINES_HIDDEN NO
  )

if(NOT ${Kokkos_VERSION} VERSION_LESS 4.4)
set_target_properties(preciceCore PROPERTIES
  CXX_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN YES
)
endif()

# Setup release override options
if(PRECICE_RELEASE_WITH_DEBUG_LOG)
  target_compile_definitions(preciceCore PUBLIC PRECICE_RELEASE_WITH_DEBUG_LOG)
endif()

if(PRECICE_RELEASE_WITH_TRACE_LOG)
  target_compile_definitions(preciceCore PUBLIC PRECICE_RELEASE_WITH_TRACE_LOG)
endif()

if(PRECICE_RELEASE_WITH_ASSERTIONS)
  target_compile_definitions(preciceCore PUBLIC PRECICE_RELEASE_WITH_ASSERTIONS)
endif()


# Setup Boost
target_compile_definitions(preciceCore PUBLIC BOOST_ALL_DYN_LINK BOOST_ASIO_ENABLE_OLD_SERVICES BOOST_GEOMETRY_DISABLE_DEPRECATED_03_WARNING)
target_link_libraries(preciceCore PUBLIC
  Boost::boost
  Boost::log
  Boost::log_setup
  Boost::program_options
  Boost::system
  Boost::thread
  Boost::unit_test_framework
  )
if(UNIX OR APPLE OR MINGW)
  target_compile_definitions(preciceCore PUBLIC _GNU_SOURCE)
  target_link_libraries(preciceCore PUBLIC ${CMAKE_DL_LIBS})
endif()
if(PRECICE_FEATURE_LIBBACKTRACE_STACKTRACES)
  target_compile_definitions(preciceCore PRIVATE BOOST_STACKTRACE_USE_BACKTRACE)
  target_link_libraries(preciceCore PRIVATE internal::libbacktrace)
endif()

# Force Intel compiler to honor nan and infinite and use value-preserving floating point mode
target_compile_options(preciceCore PUBLIC $<$<CXX_COMPILER_ID:IntelLLVM>:-fhonor-infinities -fhonor-nans -fp-model=precise>)

# Setup Eigen3
target_link_libraries(preciceCore PUBLIC Eigen3::Eigen)
target_compile_definitions(preciceCore PUBLIC "$<$<CONFIG:DEBUG>:EIGEN_INITIALIZE_MATRICES_BY_NAN>")

# Setup LIBXML2
target_link_libraries(preciceCore PUBLIC LibXml2::LibXml2)

# Setup FMT
target_link_libraries(preciceCore PUBLIC fmtlib-static)

# Setup MPI
if (PRECICE_FEATURE_MPI_COMMUNICATION)
  target_link_libraries(preciceCore PUBLIC MPI::MPI_CXX)
else()
  target_compile_definitions(preciceCore PUBLIC PRECICE_NO_MPI)
endif()

# Setup PETSC
if (PRECICE_FEATURE_PETSC_MAPPING AND PRECICE_FEATURE_MPI_COMMUNICATION)
  target_link_libraries(preciceCore PUBLIC PETSc::PETSc)
else()
  target_compile_definitions(preciceCore PUBLIC PRECICE_NO_PETSC)
endif()

# Setup Ginkgo
if (PRECICE_FEATURE_GINKGO_MAPPING)
    # Next we configure the common kernels for the matrix assembly
    target_sources(preciceCore PRIVATE
    src/mapping/device/GinkgoRBFKernels.cpp
    src/mapping/device/Ginkgo.cpp
    )

  target_link_libraries(preciceCore PRIVATE Ginkgo::ginkgo Kokkos::kokkos)
  target_compile_definitions(preciceCore PUBLIC EIGEN_NO_CUDA)

  # Disable compiler launcher since Kokkos ships its own nvcc wrapper
   set_target_properties(preciceCore PROPERTIES CXX_COMPILER_LAUNCHER "")

  # Now we setup cuSolver and hipSOLVER for preCICE
  if(PRECICE_WITH_CUDA)
    target_compile_definitions(preciceCore PUBLIC PRECICE_WITH_CUDA)
    # Required to make nvcc happy
    target_compile_definitions(preciceCore PRIVATE BOOST_PP_VARIADICS=1)
    target_sources(preciceCore PRIVATE ${PROJECT_SOURCE_DIR}/src/mapping/device/CudaQRSolver.cuh ${PROJECT_SOURCE_DIR}/src/mapping/device/CudaQRSolver.cu)
    target_link_libraries(preciceCore PRIVATE cusolver)
    message(STATUS "Enabled CUDA support for preCICE via the Ginkgo library.")
  endif()

  if(PRECICE_WITH_HIP)
    target_compile_definitions(preciceCore PUBLIC PRECICE_WITH_HIP)
    set(_hiplist ${PROJECT_SOURCE_DIR}/src/mapping/device/HipQRSolver.hip.cpp)
    target_sources(preciceCore PRIVATE ${_hiplist})
    target_link_libraries(preciceCore PRIVATE hip::device Ginkgo::ginkgo roc::hipsolver)

    message(STATUS "Enabled HIP support for preCICE via the Ginkgo library.")
  endif()

  if(PRECICE_WITH_OPENMP)
    target_compile_definitions(preciceCore PUBLIC PRECICE_WITH_OPENMP)
    target_link_libraries(preciceCore PUBLIC OpenMP::OpenMP_CXX)
    # There is nothing special required for OpenMP, looking for OpenMP and setting the flag is handled by the unified kernels
    message(STATUS "Enabled OpenMP support for preCICE via the Ginkgo library.")
  endif()

  # Finally, we link Ginkgo into preCICE
  target_link_libraries(preciceCore PUBLIC Ginkgo::ginkgo)
else()
  target_compile_definitions(preciceCore PUBLIC PRECICE_NO_GINKGO)
endif()

# Option Python
if (PRECICE_FEATURE_PYTHON_ACTIONS)
  target_link_libraries(preciceCore PUBLIC Python3::NumPy Python3::Python)
  target_compile_definitions(preciceCore PUBLIC NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION)
else()
  target_compile_definitions(preciceCore PUBLIC PRECICE_NO_PYTHON)
endif()

# Includes configuration for the core
target_include_directories(preciceCore PUBLIC
  $<BUILD_INTERFACE:${preCICE_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${preCICE_BINARY_DIR}/src>
  $<INSTALL_INTERFACE:include>
  )

# Add precice as an empty target
add_library(precice)
set_target_properties(precice PROPERTIES
  VERSION ${preCICE_VERSION}
  SOVERSION ${preCICE_SOVERSION}
  CXX_VISIBILITY_PRESET default
  VISIBILITY_INLINES_HIDDEN NO
  )

  if(NOT ${Kokkos_VERSION} VERSION_LESS 4.4)
set_target_properties(precice PROPERTIES
  CXX_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN YES
 )
endif()


target_include_directories(precice PUBLIC
  $<BUILD_INTERFACE:${preCICE_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${preCICE_BINARY_DIR}/src>
  $<INSTALL_INTERFACE:include>
  )

target_link_libraries(precice PRIVATE preciceCore)

# Sources Configuration
include(${CMAKE_CURRENT_LIST_DIR}/src/sources.cmake)

# File Configuration
configure_file("${PROJECT_SOURCE_DIR}/src/precice/impl/versions.hpp.in" "${PROJECT_BINARY_DIR}/src/precice/impl/versions.hpp" @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/src/testing/SourceLocation.hpp.in" "${PROJECT_BINARY_DIR}/src/testing/SourceLocation.hpp" @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/src/precice/Version.h.in" "${PROJECT_BINARY_DIR}/src/precice/Version.h" @ONLY)

if(PRECICE_BUILD_UNITY)
  set_target_properties(preciceCore PROPERTIES
    UNITY_BUILD Yes
    UNITY_BUILD_MODE GROUP
  )
  get_target_property(_precice_sources preciceCore SOURCES)
  foreach(SOURCE_FILE IN LISTS _precice_sources)
    if(SOURCE_FILE MATCHES "^src/([a-z0-9]+)/")
      set_source_files_properties(${SOURCE_FILE} PROPERTIES UNITY_GROUP "${CMAKE_MATCH_1}")
    endif()
  endforeach()
  unset(_precice_sources)
endif()

# Configure version file generation
include(GenerateVersionInformation)
# Configure CMake variables only. The result contains a placeholder for the git revision
configure_file("${PROJECT_SOURCE_DIR}/src/precice/impl/versions.cpp.in" "${PROJECT_BINARY_DIR}/src/precice/impl/versions.cpp.in" @ONLY)
# Dynamically replace the git revision information.
include(${CMAKE_CURRENT_LIST_DIR}/cmake/DetectGitRevision.cmake)


# Setup export header
include(GenerateExportHeader)
GENERATE_EXPORT_HEADER(preciceCore
  EXPORT_FILE_NAME "src/precice/export.h"
  BASE_NAME PRECICE
  EXPORT_MACRO_NAME PRECICE_API
  INCLUDE_GUARD_NAME PRECICE_EXPORT
  )

# Add the export header to the public headers of the preCICE lib
set_property(TARGET precice APPEND PROPERTY PUBLIC_HEADER
    ${PROJECT_BINARY_DIR}/src/precice/export.h)

#
# Configuration of Target precice-tools
#
if (PRECICE_BUILD_TOOLS)
  add_executable(precice-tools "src/drivers/main.cpp")
  target_link_libraries(precice-tools
   PRIVATE
   precice
   )
  set_target_properties(precice-tools PROPERTIES
    # precice is a C++17 project
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED Yes
    CXX_EXTENSIONS No
    )
endif()

#
# Configuration of Target testprecice
#
IF (BUILD_TESTING)
  add_executable(testprecice "src/testing/main.cpp")
  target_link_libraries(testprecice
    PRIVATE
    preciceCore
    )
  set_target_properties(testprecice PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED Yes
    CXX_EXTENSIONS No
    )
  target_include_directories(testprecice PRIVATE
    ${preCICE_SOURCE_DIR}/src
    ${preCICE_SOURCE_DIR}/tests
    )

  message(STATUS "Including test sources")
  # Test Sources Configuration
  include(${CMAKE_CURRENT_LIST_DIR}/src/tests.cmake)
  include(${CMAKE_CURRENT_LIST_DIR}/tests/tests.cmake)

  if(PRECICE_BUILD_UNITY)
    set_target_properties(testprecice PROPERTIES
      UNITY_BUILD Yes
      UNITY_BUILD_MODE GROUP
      )
    # Add integration tests from tests/ sources to a single unity group
    get_target_property(_testprecice_sources testprecice SOURCES)
    foreach(SOURCE_FILE IN LISTS _testprecice_sources)
      if(SOURCE_FILE MATCHES "^tests/")
        set_source_files_properties(${SOURCE_FILE} PROPERTIES UNITY_GROUP "integration")
      endif()
    endforeach()
  endif()

else(BUILD_TESTING)
  message(STATUS "Excluding test sources")
endif(BUILD_TESTING)

# Include Native C Bindings
if (PRECICE_BINDINGS_C)
  # include(${CMAKE_CURRENT_LIST_DIR}/extras/bindings/c/CMakeLists.txt)
  add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/extras/bindings/c)
endif()

# Include Native Fortran Bindings
if (PRECICE_BINDINGS_FORTRAN)
  # include((${CMAKE_CURRENT_LIST_DIR}/extras/bindings/fortran/CMakeLists.txt)
  add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/extras/bindings/fortran)
endif()
#
# Install Targets for precice
#

# Setup General Install for:
# precice - the library
# precice-tools - the precice binary
include(GNUInstallDirs)
install(TARGETS precice
  EXPORT preciceTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/precice
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/precice
  )

if(PRECICE_BUILD_TOOLS)
  install(TARGETS precice-tools
    EXPORT preciceTargets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )

  install(PROGRAMS ${preCICE_SOURCE_DIR}/tools/profiling/precice-profiling
    DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

# Install examples
install(DIRECTORY examples
  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/precice
  PATTERN ".gitignore" EXCLUDE
  )

# Export the Targets to install
install(EXPORT preciceTargets
  FILE preciceTargets.cmake
  NAMESPACE precice::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/precice
  )

# Generate a Package Config File for precice
include(CMakePackageConfigHelpers)
write_basic_package_version_file("preciceConfigVersion.cmake"
  VERSION ${preCICE_VERSION}
  COMPATIBILITY SameMajorVersion
  )

# Install the Config and the ConfigVersion files
install(FILES "cmake/preciceConfig.cmake" "${preCICE_BINARY_DIR}/preciceConfigVersion.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/precice
  )

# Setup the config in the build directory
export(EXPORT preciceTargets
  NAMESPACE precice::
  FILE "preciceTargets.cmake")
file(COPY "cmake/preciceConfig.cmake"
  DESTINATION "${preCICE_BINARY_DIR}")

# Add an alias to allow subprojects to use the namespaced name
add_library(precice::precice ALIAS precice)


# Set the directory used to prepare files for packaging
set(PRECICE_PACKAGING_DIR "${CMAKE_CURRENT_BINARY_DIR}/packaging")
mark_as_advanced(PRECICE_PACKAGING_DIR)

# Compress and install the manpages
find_program(GZIP_EXE gzip DOC "The gzip executable")
if(GZIP_EXE)
  # Process manpages for binaries
  file(COPY docs/man/man1/precice-tools.1 DESTINATION packaging/man1)
  file(GLOB PRECICE_MAN_PAGES "${PRECICE_PACKAGING_DIR}/man1/*.1")
  foreach(manpage ${PRECICE_MAN_PAGES})
    message(STATUS "Compressing manpage: ${manpage}")
    execute_process(COMMAND "${GZIP_EXE}" "-9nf" "${manpage}")
  endforeach()

  # Install compressed manpages
  install(DIRECTORY ${PRECICE_PACKAGING_DIR}/man1
    DESTINATION ${CMAKE_INSTALL_MANDIR}
    )
else()
  message(WARNING "Installing uncompressed manpages")
  # Install uncompressed manpages
  install(DIRECTORY docs/man/man1
    DESTINATION ${CMAKE_INSTALL_MANDIR}
    )
endif()

# Configure the pkg-config file for the binary directory
set(_precice_bindings "")
if (PRECICE_BINDINGS_C)
  set(_precice_bindings "-I${PROJECT_SOURCE_DIR}/extras/bindings/c/include")
endif()
configure_file(
  "${PROJECT_SOURCE_DIR}/cmake/binary.pc.in"
  "libprecice.pc"
  @ONLY
  )
unset(_precice_bindings)

# Configure the pkg-config file for the installation
configure_file(
  "${PROJECT_SOURCE_DIR}/cmake/install.pc.in"
  "lib/pkgconfig/libprecice.pc"
  @ONLY
  )
install(DIRECTORY "${preCICE_BINARY_DIR}/lib/pkgconfig"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

#
# Add uninstall
#

include(${CMAKE_CURRENT_LIST_DIR}/cmake/Uninstall.cmake)


#
# CPack
#

if (PRECICE_CONFIGURE_PACKAGE_GENERATION)
  include(${CMAKE_CURRENT_LIST_DIR}/cmake/CPackConfig.cmake)
endif()


#
# CTest
#
if (BUILD_TESTING)
  print_empty()
  print_section("TESTS")

  include(${CMAKE_CURRENT_LIST_DIR}/cmake/CTestConfig.cmake)

  # Add precice-profiling tests
  add_subdirectory(tools/profiling)

  #
  # Add test_install
  #

  include(${CMAKE_CURRENT_LIST_DIR}/cmake/TestInstall.cmake)
endif()


#
# Tooling
#

add_custom_target(
  doxygen
  COMMAND doxygen
  WORKING_DIRECTORY ${preCICE_SOURCE_DIR}
  USES_TERMINAL
  )

add_custom_target(
  format
  COMMAND tools/formatting/format-all
  WORKING_DIRECTORY ${preCICE_SOURCE_DIR}
  USES_TERMINAL
  )

add_custom_target(
  sourcesIndex
  COMMAND tools/building/updateSourceFiles.py
  WORKING_DIRECTORY ${preCICE_SOURCE_DIR}
  USES_TERMINAL
  )

add_custom_target(
  changelog
  COMMAND tools/building/createChangelog
  WORKING_DIRECTORY ${preCICE_SOURCE_DIR}
  USES_TERMINAL
  )

add_custom_target(
  tidy
  COMMAND tools/linting/run_clang_tidy.sh
  WORKING_DIRECTORY ${preCICE_SOURCE_DIR}
  USES_TERMINAL
  )
