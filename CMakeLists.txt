cmake_minimum_required (VERSION 3.10.2)
project(preCICE VERSION 2.2.0 LANGUAGES CXX)
set(preCICE_SOVERSION ${preCICE_VERSION_MAJOR})

#
# Overview of this configuration
#
# PREAMBLE
# Setup Options
# Find Mandatory Dependencies
# Find Configurable Dependencies
# Configuration of Target precice
# Configuration of Target binprecice
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
  set(CMAKE_BUILD_TYPE Debug CACHE STRINGS "The type of this build" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

option(PRECICE_MPICommunication "Enables MPI-based communication and running coupling tests." ON)
option(PRECICE_PETScMapping "Enable use of the PETSc linear algebra library." ON)
option(PRECICE_PythonActions "Python support" ON)
option(PRECICE_Packages "Configure package generation." ON)
option(PRECICE_InstallTest "Add test binary and necessary files to install target." OFF)
option(BUILD_SHARED_LIBS "Build shared libraries by default" OFF)
option(BUILD_TESTING "Build tests" ON)
option(PRECICE_ALWAYS_VALIDATE_LIBS "Validate libraries even after the validatation succeeded." OFF)
option(PRECICE_ENABLE_C "Enable the native C bindings" ON)
option(PRECICE_ENABLE_FORTRAN "Enable the native Fortran bindings" ON)

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
  "PRECICE_InstallTest;Install tests/testfiles"
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
find_package(Boost 1.65.1 REQUIRED
  COMPONENTS filesystem log log_setup program_options system thread unit_test_framework
  )

# Eigen
if(TPL_ENABLE_EIGEN3)
  # Use EIGEN3_ROOT to set the directory
  xsdk_tpl_require(EIGEN3 EIGEN3_INCLUDE_DIR)
endif()
find_package(Eigen3 3.2 REQUIRED)
precice_validate_eigen()

# LibXML2
if(TPL_ENABLE_LIBXML2)
  xsdk_tpl_require(LIBXML2 LIBXML2_LIBRARIES LIBXML2_INCLUDE_DIRS)
endif()
find_package(LibXml2 REQUIRED)
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

# prettyprint
if(TPL_ENABLE_PRETTYPRINT)
  xsdk_tpl_require(PRETTYPRINT PRETTYPRINT_INCLUDE_DIR)
  add_library(prettyprint INTERFACE IMPORTED)
  set_property(TARGET prettyprint PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${PRETTYPRINT_INCLUDE_DIR}")
else()
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/prettyprint)
endif()
precice_validate_prettyprint()


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
else()
  message(STATUS "PETSc support disabled")
endif()

# Option Python
if (PRECICE_PythonActions)
  if (TPL_ENABLE_PYTHON)
    xsdk_tpl_require(PYTHON PYTHON_LIBRARY PYTHON_INCLUDE_DIR NumPy_INCLUDE_DIR)
    find_package(PythonLibs REQUIRED)

    # Override NumPy
    # TODO: Incorporate into the FindNumPy module
    set(NumPy_FOUND True CACHE BOOL "NumPy found?" FORCE)
    if(NOT TARGET NumPy::NumPy)
      add_library(NumPy::NumPy INTERFACE IMPORTED)
      set_property(TARGET NumPy::NumPy PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${NumPy_INCLUDE_DIR}")
    endif()
  else()
    find_package(PythonLibs REQUIRED)
    find_package(NumPy REQUIRED)
  endif()
  if(PYTHONLIBS_VERSION_STRING VERSION_LESS 2.7)
    message(FATAL_ERROR "The selected Python Version is too old: ${PYTHONLIBS_VERSION_STRING} (minimum 2.7)")
  endif()
  precice_validate_libpython()
  precice_validate_numpy()
else()
  message(STATUS "Python support disabled")
endif()


#
# Configuration of Target precice
#

print_empty()
print_section("TARGETS & PACKAGES")

# Add a dummy to silence add_library warning for cmake < 3.11.
set(preCICE_DUMMY "")
if(CMAKE_VERSION VERSION_LESS "3.11")
  set(preCICE_DUMMY "${CMAKE_CURRENT_BINARY_DIR}/dummy.cpp")
  if(NOT EXISTS "${preCICE_DUMMY}")
    file(WRITE "${preCICE_DUMMY}" "")
  endif()
endif()

# Add precice as an empty target
add_library(precice ${preCICE_DUMMY})
set_target_properties(precice PROPERTIES
  # precice is a C++14 project
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED Yes
  CXX_EXTENSIONS No
  VERSION ${preCICE_VERSION}
  SOVERSION ${preCICE_SOVERSION}
  )

# Setup Boost
target_compile_definitions(precice PRIVATE BOOST_ALL_DYN_LINK BOOST_ASIO_ENABLE_OLD_SERVICES BOOST_GEOMETRY_DISABLE_DEPRECATED_03_WARNING)
target_link_libraries(precice PRIVATE
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
  target_compile_definitions(precice PRIVATE _GNU_SOURCE)
  target_link_libraries(precice PRIVATE ${CMAKE_DL_LIBS})
endif()

# Setup Eigen3
target_link_libraries(precice PRIVATE Eigen3::Eigen)
target_compile_definitions(precice PRIVATE "$<$<CONFIG:DEBUG>:EIGEN_INITIALIZE_MATRICES_BY_NAN>")

# Setup LIBXML2
target_include_directories(precice PRIVATE ${LIBXML2_INCLUDE_DIR})
target_link_libraries(precice PRIVATE ${LIBXML2_LIBRARIES})

# Setup Prettyprint
target_link_libraries(precice PRIVATE prettyprint)

# Setup JSON
target_link_libraries(precice PRIVATE JSON)

# Setup MPI
if (PRECICE_MPICommunication)
  target_link_libraries(precice PRIVATE MPI::MPI_CXX)
else()
  target_compile_definitions(precice PRIVATE PRECICE_NO_MPI)
endif()

# Setup PETSC
if (PRECICE_PETScMapping AND PRECICE_MPICommunication)
  target_link_libraries(precice PRIVATE PETSc::PETSc)
else()
  target_compile_definitions(precice PRIVATE PRECICE_NO_PETSC)
endif()

# Option Python
if (PRECICE_PythonActions)
  target_link_libraries(precice PRIVATE NumPy::NumPy)
  target_include_directories(precice PRIVATE ${PYTHON_INCLUDE_DIRS})
  target_compile_definitions(precice PRIVATE NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION)
  target_link_libraries(precice PRIVATE ${PYTHON_LIBRARIES})
else()
  target_compile_definitions(precice PRIVATE PRECICE_NO_PYTHON)
endif()


# File Configuration
include(GenerateVersionInformation)
include(${CMAKE_CURRENT_LIST_DIR}/cmake/DetectGitRevision.cmake)
configure_file("${PROJECT_SOURCE_DIR}/src/precice/impl/versions.hpp.in" "${PROJECT_BINARY_DIR}/src/precice/impl/versions.hpp" @ONLY)

# Includes Configuration
target_include_directories(precice PUBLIC
  $<BUILD_INTERFACE:${preCICE_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${preCICE_BINARY_DIR}/src>
  $<INSTALL_INTERFACE:include>
  )

# Sources Configuration
include(${CMAKE_CURRENT_LIST_DIR}/src/sources.cmake)


#
# Configuration of Target binprecice
#

add_executable(binprecice "src/drivers/main.cpp")
target_link_libraries(binprecice
  PRIVATE
  Threads::Threads
  precice
  prettyprint
  Eigen3::Eigen
  Boost::boost
  Boost::filesystem
  Boost::log
  Boost::log_setup
  Boost::program_options
  Boost::system
  Boost::thread
  Boost::unit_test_framework
  )
set_target_properties(binprecice PROPERTIES
  # precice is a C++14 project
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED Yes
  CXX_EXTENSIONS No
  )
# Copy needed properties from the lib to the executatble. This is necessary as
# this executable uses the library source, not only the interface.
copy_target_property(precice binprecice COMPILE_DEFINITIONS)
copy_target_property(precice binprecice COMPILE_OPTIONS)
if(PRECICE_MPICommunication)
  target_link_libraries(binprecice PRIVATE MPI::MPI_CXX)
endif()
if(PRECICE_MPICommunication AND PRECICE_PETScMapping)
  target_link_libraries(binprecice PRIVATE PETSc::PETSc)
endif()


#
# Configuration of Target testprecice
#
IF (BUILD_TESTING)
  add_executable(testprecice "src/testing/main.cpp")
  target_link_libraries(testprecice
    PRIVATE
    Threads::Threads
    precice
    Eigen3::Eigen
    prettyprint
    Boost::boost
    Boost::filesystem
    Boost::log
    Boost::log_setup
    Boost::program_options
    Boost::system
    Boost::thread
    Boost::unit_test_framework
    )
  set_target_properties(testprecice PROPERTIES
    # precice is a C++14 project
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED Yes
    CXX_EXTENSIONS No
    )
  # Copy needed properties from the lib to the executatble. This is necessary as
  # this executable uses the library source, not only the interface.
  copy_target_property(precice testprecice COMPILE_DEFINITIONS)
  copy_target_property(precice testprecice COMPILE_OPTIONS)

  # Testprecice fully depends on MPI and PETSc.
  if(PRECICE_MPICommunication)
    target_link_libraries(testprecice PRIVATE MPI::MPI_CXX)
  endif()
  if(PRECICE_MPICommunication AND PRECICE_PETScMapping)
    target_link_libraries(testprecice PRIVATE PETSc::PETSc)
  endif()

  message(STATUS "Including test sources")
  # Test Sources Configuration
  include(${CMAKE_CURRENT_LIST_DIR}/src/tests.cmake)
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
# binprecice - the precice binary
include(GNUInstallDirs)
install(TARGETS precice binprecice
  EXPORT preciceTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/precice
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/precice
  )

if(BUILD_TESTING AND PRECICE_InstallTest)
  # Install the testprecice target
  install(TARGETS testprecice
    EXPORT preciceTargets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )

  # Install the resources necessary for the tests
  install(DIRECTORY src
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/precice
    FILES_MATCHING
    PATTERN "*.xml"
    )
  install(DIRECTORY src/action/tests/
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/precice/src/action/tests
    FILES_MATCHING
    PATTERN "*.py"
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
  file(COPY docs/man/man1/binprecice.1 DESTINATION packaging/man1)
  if(PRECICE_InstallTest)
    file(COPY docs/man/man1/testprecice.1 DESTINATION packaging/man1)
  endif()
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

  if(BUILD_TESTING AND PRECICE_InstallTest)
    install(FILES docs/man/man1/testprecice.1
      DESTINATION share/man
      )
  endif()
endif()

# Configure a pkg-config file
configure_file(
  "${PROJECT_SOURCE_DIR}/tools/releasing/packaging/debian/precice.pc.in"
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
