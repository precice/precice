cmake_minimum_required (VERSION 3.16.3)
project(preCICE VERSION 2.5.0 LANGUAGES CXX)
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

include(CheckCXX11Library)
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

option(PRECICE_MPICommunication "Enables MPI-based communication and running coupling tests." ON)
option(PRECICE_PETScMapping "Enable use of the PETSc linear algebra library." ON)
option(PRECICE_GinkgoParallelization "Enable use of the Ginkgo linear algebra library for parallelization." OFF)
option(PRECICE_PythonActions "Python support" ON)
option(PRECICE_Packages "Configure package generation." ON)
option(BUILD_SHARED_LIBS "Build shared instead of static libraries" ON)
option(BUILD_TESTING "Build tests" ON)
option(PRECICE_ALWAYS_VALIDATE_LIBS "Validate libraries even after the validatation succeeded." OFF)
option(PRECICE_ENABLE_C "Enable the native C bindings" ON)
option(PRECICE_ENABLE_FORTRAN "Enable the native Fortran bindings" ON)
option(PRECICE_BUILD_TOOLS "Build the \"precice-tools\" executable" ON)
option(PRECICE_ENABLE_LIBBACKTRACE "Enable libbacktrace for stacktrace generation." OFF)

option(PRECICE_RELEASE_WITH_DEBUG_LOG "Enable debug logging in release builds" OFF)
option(PRECICE_RELEASE_WITH_TRACE_LOG "Enable trace logging in release builds" OFF)
option(PRECICE_RELEASE_WITH_ASSERTIONS "Enable assertions in release builds" OFF)

xsdk_tpl_option_override(PRECICE_MPICommunication TPL_ENABLE_MPI)
xsdk_tpl_option_override(PRECICE_PETScMapping TPL_ENABLE_PETSC)
xsdk_tpl_option_override(PRECICE_PythonActions TPL_ENABLE_PYTHON)

if(PRECICE_PETScMapping AND NOT PRECICE_MPICommunication)
  message(FATAL_ERROR "Please enable MPI to use PETSC.")
endif()

set(PRECICE_CTEST_MPI_FLAGS "" CACHE STRING "Add additional flags to mpiexec for running tests via CTest.")

include(XSDKOptions)

# Print information about this configuration
include(PrintHelper)
print_configuration(
  ADDITIONAL
  "PRECICE_Packages;Configure package generation"
  "PRECICE_CTEST_MPI_FLAGS;Additional CTest MPI Flags"
  )
print_empty()

include(FeatureSummary)
add_feature_info(MPICommunication PRECICE_MPICommunication
  "Enables the MPI communication back-end.

   This enables the MPI communication back-end which is highly recommended on multi-node systems.
   See the documentation of the CMake module FindMPI to control its functionality.

   This feature can be enabled/disabled by setting the PRECICE_MPICommunication CMake option.
")
add_feature_info(PETScMapping PRECICE_PETScMapping
  "Enables the PETSc-powered radial basic function mappings.

   This enables the PETSc-based variant of radial basis function mappings which can run in parallel,
   also across different compute nodes. This is highly recommended for large cases running in parallel.

   This feature can be enabled/disabled by setting the PRECICE_PETScMapping CMake option.
   Requires MPICommunication.
  ")
  add_feature_info(GinkgoParallelization PRECICE_GinkgoParallelization
  "Enables the Ginkgo library backend which supports the usage of GPUs and OpenMP in radial basis function mappings.

   Currently, this enables radial basis function mappings to leverage the compute power of Nvidia and AMD GPUs.
   It uses iterative solvers to compute the solution of interpolation problems.

   This feature can be enabled/disabled by setting the PRECICE_GinkgoParallelization CMake option.
   Requires the Ginkgo library.
  ")
add_feature_info(PythonActions PRECICE_PythonActions
  "Enables the support for user-defined python actions.

   preCICE allows to manipulate coupling data at runtime using configurable actions.
   This feature enables the support for user-defined actions written in Python based on numpy.

   This feature can be enabled/disabled by setting the PRECICE_PythonActions CMake option.
  ")
add_feature_info(CBindings PRECICE_ENABLE_C
  "Enables the native C bindings.

   preCICE provides native bindings for C, which are compiled into the library.
   This feature enables the compilation and installation of the bindings into the library.
   Note that we strongly recommend to compile with C bindings enabled for compatibility reasons.

   This feature can be enabled/disabled by setting the PRECICE_ENABLE_C CMake option.
  ")
add_feature_info(FortranBindings PRECICE_ENABLE_FORTRAN
  "Enables the native Fortran bindings.

   preCICE provides native bindings for Fortran, which are compiled into the library.
   This feature enables the compilation and installation of the bindings into the library.
   Note that we strongly recommend to compile with Fortran bindings enabled for compatibility reasons.

   This feature can be enabled/disabled by setting the PRECICE_ENABLE_FORTRAN CMake option.
  ")
add_feature_info(PreciceTools PRECICE_BUILD_TOOLS
  "Build the \"precice-tools\" executable

   preCICE offers in addition to the core library a variety of tools in order to check configuration
   files or generate xml references. The tools are copmiled into an executable called \"precice-tools\".

   This feature can be enabled/disabled by setting the PRECICE_BUILD_TOOLS CMake option.
  ")
  add_feature_info(libbacktrace PRECICE_ENABLE_LIBBACKTRACE
  "Enables libbacktrace for stracktrace generation.

   preCICE assertions are helpful during development, however, the default stacktraces may not be sufficient.
   This feature enables the support for libbacktrace in boost.stacktrace.

   This feature can be enabled/disabled by setting the PRECICE_ENABLE_LIBBACKTRACE CMake option.
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
find_package(Boost 1.71.0 REQUIRED
  COMPONENTS filesystem log log_setup program_options system thread unit_test_framework
  )
mark_as_advanced(Boost_FILESYSTEM_LIBRARY_RELEASE Boost_INCLUDE_DIR Boost_LOG_LIBRARY_RELEASE Boost_LOG_SETUP_LIBRARY_RELEASE Boost_PROGRAM_OPTIONS_LIBRARY_RELEASE Boost_SYSTEM_LIBRARY_RELEASE Boost_THREAD_LIBRARY_RELEASE Boost_UNIT_TEST_FRAMEWORK_LIBRARY_RELEASE)
mark_as_advanced(Boost_FILESYSTEM_LIBRARY_DEBUG Boost_INCLUDE_DIR Boost_LOG_LIBRARY_DEBUG Boost_LOG_SETUP_LIBRARY_DEBUG Boost_PROGRAM_OPTIONS_LIBRARY_DEBUG Boost_SYSTEM_LIBRARY_DEBUG Boost_THREAD_LIBRARY_DEBUG Boost_UNIT_TEST_FRAMEWORK_LIBRARY_DEBUG)
message(STATUS "Found Boost ${Boost_VERSION}")

# Eigen
if(TPL_ENABLE_EIGEN3)
  # Use EIGEN3_ROOT to set the directory
  xsdk_tpl_require(EIGEN3 EIGEN3_INCLUDE_DIR)
endif()
find_package(Eigen3 3.3.7 REQUIRED)
message(STATUS "Found Eigen ${EIGEN3_VERSION}")
precice_validate_eigen()

# LibXML2
find_package(LibXml2 REQUIRED)
message(STATUS "Found LibXml2 ${LIBXML2_VERSION_STRING}")
precice_validate_libxml2()

# nlohmann/JSON
if(TPL_ENABLE_JSON)
  xsdk_tpl_require(JSON JSON_INCLUDE_DIR)
  add_library(JSON INTERFACE IMPORTED)
  set_property(TARGET JSON PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${JSON_INCLUDE_DIR}")
else()
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/json)
endif()
precice_validate_json()

# libfmt
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/fmt)
precice_validate_fmtlib()

#
# Find Configurable Dependencies
#

# Option: PRECICE_MPICommunication
if (PRECICE_MPICommunication)
  set(MPI_DETERMINE_LIBRARY_VERSION ON) # Try to detect the library vendor

  find_package(MPI REQUIRED)
  message(STATUS "MPI Version: ${MPI_CXX_LIBRARY_VERSION_STRING}")

  # We need to disable the MPIPorts tests in case we detect Open MPI
  if(MPI_CXX_LIBRARY_VERSION_STRING MATCHES "Open MPI")
    set(PRECICE_MPI_OPENMPI TRUE CACHE BOOL "Is the MPI implementation OpenMPI?" FORCE)
    mark_as_advanced(PRECICE_MPI_OPENMPI)
  endif()
endif()

# Option: PETSC
if (PRECICE_PETScMapping)
  if (TPL_ENABLE_PETSC)
    xsdk_tpl_require(PETSC PETSC_DIR PETSC_ARCH)
    # PETSc detection uses primarily these ENVs
  endif()
  find_package(PETSc 3.12 REQUIRED)
  # No validation required as PETSc does this internally
  message(STATUS "Found PETSc ${PETSc_VERSION}")
else()
  message(STATUS "PETSc support disabled")
endif()

# Option: Ginkgo
if (PRECICE_GinkgoParallelization)

  # Check if CUDA is available on the system
  include(CheckLanguage)
  check_language(CUDA)

  # Enable CUDA if driver and libraries are available; otherwise it'll be ignored
  if(CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    set(CMAKE_CUDA_FLAGS "-lineinfo --allow-unsupported-compiler --use_fast_math --std=c++17 --expt-extended-lambda --expt-relaxed-constexpr -Xcompiler -fopenmp -Wall")
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
    set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
    message(STATUS "Enabled CUDA support for Ginkgo library.")

    add_compile_definitions("PRECICE_WITH_CUDA")

  endif()

  if(NOT DEFINED HIP_PATH)
      set(HIP_PATH "/opt/rocm/hip" CACHE PATH "Path to which HIP has been installed")
   endif()

  if(NOT DEFINED ROCM_PATH)
     set(ROCM_PATH "/opt/rocm" CACHE PATH "Path to which ROCM has been installed")
  endif()

  if(NOT DEFINED HCC_PATH)
      set(HCC_PATH "${ROCM_PATH}/hcc" CACHE PATH "Path to which HCC has been installed")
  endif()

  if(NOT DEFINED HIP_CLANG_PATH)
      set(HIP_CLANG_PATH "${ROCM_PATH}/llvm/bin" CACHE PATH "Path to which HIP compatible clang binaries have been installed")
  endif()

  find_package(HIP QUIET)

  if(HIP_FOUND)
      set(CMAKE_MODULE_PATH "${HIP_PATH}/cmake" ${CMAKE_MODULE_PATH})
      list(APPEND CMAKE_PREFIX_PATH
     	"${HIP_PATH}/lib/cmake"
      	"${HIP_PATH}/../lib/cmake"
      )

      include(${HIP_PATH}/cmake/FindHIP.cmake)

      set_target_properties(hip::device PROPERTIES INTERFACE_COMPILE_OPTIONS "-fPIC")
      set_target_properties(hip::device PROPERTIES INTERFACE_LINK_LIBRARIES "hip::host")

      set(CMAKE_HIP_FLAGS "-fPIC -Wno-deprecated-register -w -fgpu-rdc")
      set(CMAKE_POSITION_INDEPENDENT_CODE ON)

      set(GINKGO_HIPCC_OPTIONS "-fPIC")

      find_package(hipsolver)

      message(STATUS "Enabled AMD HIP support for Ginkgo library.")

      add_compile_definitions("PRECICE_WITH_HIP")
  endif()

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

  include_directories(${EIGEN3_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/fmt/include)

  # Check if OpenMP is available
  find_package(OpenMP)
  if(OpenMP_FOUND)
    message(STATUS "Enabled OpenMP support for Ginkgo library.")
  endif()

endif()

# Option Python
if (PRECICE_PythonActions)
  find_package(Python3 REQUIRED COMPONENTS Development NumPy)
  message(STATUS "Found Python ${Python3_VERSION} with NumPy ${Python3_NumPy_VERSION}")
  precice_validate_libpython()
  precice_validate_numpy()
else()
  message(STATUS "Python support disabled")
endif()

# Option ENABLE_LIBBACKTRACE
# Note that the FindBacktrace module does not work reliably to find libbacktrace installations
if (PRECICE_ENABLE_LIBBACKTRACE)
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
  CXX_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN YES
  )

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
  Boost::filesystem
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
if(PRECICE_ENABLE_LIBBACKTRACE)
  target_compile_definitions(preciceCore PRIVATE BOOST_STACKTRACE_USE_BACKTRACE)
  target_link_libraries(preciceCore PRIVATE internal::libbacktrace)
endif()

# Setup Eigen3
target_link_libraries(preciceCore PUBLIC Eigen3::Eigen)
target_compile_definitions(preciceCore PUBLIC "$<$<CONFIG:DEBUG>:EIGEN_INITIALIZE_MATRICES_BY_NAN>")

# Setup LIBXML2
target_link_libraries(preciceCore PUBLIC LibXml2::LibXml2)

# Setup JSON
target_link_libraries(preciceCore PUBLIC JSON fmt-header-only)

# Setup MPI
if (PRECICE_MPICommunication)
  target_link_libraries(preciceCore PUBLIC MPI::MPI_CXX)
else()
  target_compile_definitions(preciceCore PUBLIC PRECICE_NO_MPI)
endif()

# Setup PETSC
if (PRECICE_PETScMapping AND PRECICE_MPICommunication)
  target_link_libraries(preciceCore PUBLIC PETSc::PETSc)
else()
  target_compile_definitions(preciceCore PUBLIC PRECICE_NO_PETSC)
endif()

# Setup Ginkgo
if (PRECICE_GinkgoParallelization)
  set(CMAKE_CXX_FLAGS "-fPIC")
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
  set(Ginkgo_Path "" CACHE STRING "Path to the experimental Ginkgo Repository of branch a195f856e.") # TODO: Remove as soon as GPU feature is in main branch
  find_package(Ginkgo REQUIRED HINTS "${Ginkgo_Path}/build" NO_CMAKE_SYSTEM_PATH)
  target_include_directories(preciceCore PUBLIC "${Ginkgo_Path}/include")
  include(${Ginkgo_Path}/cmake/unified_kernels.cmake)
  file(GLOB ginkgo_lib "${Ginkgo_Path}/build/lib/*")

  # This command registers Ginkgo kernels for different execution backends
  ginkgo_add_unified_kernels(preciceCore src/mapping/ginkgo_kernels/GinkgoRBFKernels.cpp)

  # This compiles cuSolver or rocSOLVER into preCICE for GPU QR decomposition
  if(CMAKE_CUDA_COMPILER)
    target_sources(preciceCore PRIVATE src/mapping/cuda_kernels/qr_decomp.cu)
    target_link_libraries(preciceCore PUBLIC ${ginkgo_lib} cusolver cublas)
  elseif(HIP_FOUND)
    hip_add_library(preciceHIPQR src/mapping/cuda_kernels/qr_decomp.hip.cpp HIPCC_OPTIONS "-std=c++14" CLANG_OPTIONS "-fPIC" STATIC)
    target_link_libraries(preciceHIPQR PUBLIC ${ginkgo_lib} hip::host roc::hipsolver)
    target_include_directories(preciceHIPQR PUBLIC "${Ginkgo_Path}/include" "/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi" "/usr/lib/x86_64-linux-gnu/openmpi/include") # TODO: Fix include path
    target_link_libraries(preciceCore PRIVATE preciceHIPQR)
    endif()

  if(OpenMP_FOUND)
    target_compile_options(preciceCore PUBLIC "-fopenmp")
  endif()
else()
  target_compile_definitions(preciceCore PUBLIC PRECICE_NO_GINKGO)
endif()

# Option Python
if (PRECICE_PythonActions)
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
  CXX_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN YES
  )

target_include_directories(precice PUBLIC
  $<BUILD_INTERFACE:${preCICE_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${preCICE_BINARY_DIR}/src>
  $<INSTALL_INTERFACE:include>
  )

target_link_libraries(precice PRIVATE preciceCore)

# Sources Configuration
include(${CMAKE_CURRENT_LIST_DIR}/src/sources.cmake)

# File Configuration
include(GenerateVersionInformation)
include(${CMAKE_CURRENT_LIST_DIR}/cmake/DetectGitRevision.cmake)
configure_file("${PROJECT_SOURCE_DIR}/src/precice/impl/versions.hpp.in" "${PROJECT_BINARY_DIR}/src/precice/impl/versions.hpp" @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/src/testing/SourceLocation.hpp.in" "${PROJECT_BINARY_DIR}/src/testing/SourceLocation.hpp" @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/src/precice/Version.h.in" "${PROJECT_BINARY_DIR}/src/precice/Version.h" @ONLY)

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
else(BUILD_TESTING)
  message(STATUS "Excluding test sources")
endif(BUILD_TESTING)

# Include Native C Bindings
if (PRECICE_ENABLE_C)
  # include(${CMAKE_CURRENT_LIST_DIR}/extras/bindings/c/CMakeLists.txt)
  add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/extras/bindings/c)
endif()

# Include Native Fortran Bindings
if (PRECICE_ENABLE_FORTRAN)
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
configure_file(
  "${PROJECT_SOURCE_DIR}/cmake/binary.pc.in"
  "libprecice.pc"
  @ONLY
  )

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

if (PRECICE_Packages)
  include(${CMAKE_CURRENT_LIST_DIR}/cmake/CPackConfig.cmake)
endif()


#
# CTest
#
if (BUILD_TESTING)
  print_empty()
  print_section("TESTS")

  include(${CMAKE_CURRENT_LIST_DIR}/cmake/CTestConfig.cmake)


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
