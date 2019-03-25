cmake_minimum_required (VERSION 3.9.6)

project(preCICE VERSION 1.4.0 LANGUAGES CXX)

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

include(CopyTargetProperty)

# CMake Policies

# CMP0074: find_package() uses <PackageName>_ROOT variables.
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()
# CMP0075 Include file check macros honor CMAKE_REQUIRED_LIBRARIES
if(POLICY CMP0075)
  cmake_policy(SET CMP0075 NEW)
endif()

# Enforce a build type
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()
message(STATUS "Build configuration: " ${CMAKE_BUILD_TYPE})


#
# Setup Options
#

set(Platform "none" CACHE STRING "Special configurations for certain platforms")
set_property(CACHE Platform PROPERTY STRINGS none "SuperMUC")
message(STATUS "Special platform settings: " ${Platform})

option(MPI "Enables MPI-based communication and running coupling tests." ON)
option(PETSC "Enable use of the PETSc linear algebra library." ON)
option(PYTHON "Python support" ON)
option(PRECICE_Packages "Configure package generation." ON)
option(PRECICE_InstallTest "Add test binary and necessary files to install target." ON)

# Print information about this configuration
include(PrintHelper)
print_configuration(
  ADDITIONAL
  "MPI;Build with MPI"
  "PETSC;Build with PETSC"
  "PYTHON;Build with PYTHON"
  "Platform;Special platform settings"
  "PRECICE_Packages;Configure package generation"
  "PRECICE_InstallTest;Install tests/testfiles"
  )


#
# Find Mandatory Dependencies
#

find_package (Threads REQUIRED)

find_package(Boost 1.65.1 REQUIRED
  COMPONENTS filesystem log log_setup program_options system thread unit_test_framework)

find_package(Eigen3 3.2 REQUIRED)

find_package(LibXml2 REQUIRED)


#
# Find Configurable Dependencies
#

# Option: MPI
if (MPI)
  find_package(MPI REQUIRED)
endif()

# Option: PETSC
if (PETSC AND MPI)
  find_package(PETSc 3.6 REQUIRED)
  set(PETSC_VERSIONS "")
  set(PETSC_VERSION_MAJOR "")
  set(PETSC_VERSION_MINOR "")
  string(REGEX MATCHALL "[0-9]+" PETSC_VERSIONS ${PETSC_VERSION})
  list(GET PETSC_VERSIONS 0 PETSC_VERSION_MAJOR)
  list(GET PETSC_VERSIONS 1 PETSC_VERSION_MINOR)
else()
  message(STATUS "PETSc support disabled")
endif()

# Option Python
if (PYTHON)
  set(Python_ADDITIONAL_VERSIONS "2.7")
  find_package(PythonLibs 2.7 REQUIRED)
  find_package(NumPy REQUIRED)
else()
  message(STATUS "Python support disabled")
endif()


#
# Configuration of Target precice
#

# Add a dummy to silence add_library warning for cmake < 3.11.
if(CMAKE_VERSION VERSION_LESS "3.11")
  set(preCICE_DUMMY "${CMAKE_CURRENT_BINARY_DIR}/dummy.cpp")
  if(NOT EXISTS "${preCICE_DUMMY}")
    file(WRITE "${preCICE_DUMMY}" "")
  endif()
endif()

# Add precice as an empty target 
add_library(precice ${preCICE_DUMMY})
set_target_properties(precice PROPERTIES
  # precice is a C++11 project
  CXX_STANDARD 11
  SOVERSION ${preCICE_VERSION}
  )

# Setup Boost
target_compile_definitions(precice PRIVATE BOOST_ALL_DYN_LINK BOOST_ASIO_ENABLE_OLD_SERVICES)
target_link_libraries(precice PRIVATE ${Boost_LIBRARIES})
target_include_directories(precice PRIVATE ${Boost_INCLUDE_DIRS})
if(UNIX OR APPLE OR MINGW)
  target_link_libraries(precice PRIVATE ${CMAKE_DL_LIBS})
endif()

# Setup Eigen3
target_link_libraries(precice PRIVATE Eigen3::Eigen)
target_compile_definitions(precice PRIVATE "$<$<CONFIG:DEBUG>:EIGEN_INITIALIZE_MATRICES_BY_NAN>")

# Setup LIBXML2
target_include_directories(precice PRIVATE ${LIBXML2_INCLUDE_DIR})
target_link_libraries(precice PRIVATE ${LIBXML2_LIBRARIES})

if (Platform STREQUAL "SuperMUC")
  target_compile_definitions(precice PRIVATE SuperMUC_WORK)
endif()

# Setup MPI
if (MPI)
  target_link_libraries(precice PRIVATE MPI::MPI_CXX)
else()
  target_compile_definitions(precice PRIVATE PRECICE_NO_MPI)
endif()

# Setup PETSC
if (PETSC AND MPI)
  target_include_directories(precice PRIVATE ${PETSC_INCLUDES})
  target_link_libraries(precice PRIVATE ${PETSC_LIBRARIES})
else()
  target_compile_definitions(precice PRIVATE PRECICE_NO_PETSC)
endif()

# Option Python
if (PYTHON)
  target_link_libraries(precice PRIVATE NumPy::NumPy)
  target_include_directories(precice PRIVATE ${PYTHON_INCLUDE_DIRS})
  target_compile_definitions(precice PRIVATE NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION)
  target_link_libraries(precice PRIVATE ${PYTHON_LIBRARIES})
else()
  target_compile_definitions(precice PRIVATE PRECICE_NO_PYTHON)
endif()

# File Configuration
configure_file("${PROJECT_SOURCE_DIR}/src/versions.hpp.in" "${PROJECT_BINARY_DIR}/src/versions.hpp" @ONLY)

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
  Eigen3::Eigen
  ${Boost_LIBRARIES}
  )
set_target_properties(binprecice PROPERTIES
  # precice is a C++11 project
  CXX_STANDARD 11
  )
# Copy needed properties from the lib to the executatble. This is necessary as
# this executable uses the library source, not only the interface.
copy_target_property(precice binprecice COMPILE_DEFINITIONS)
copy_target_property(precice binprecice COMPILE_OPTIONS)
if(MPI)
  target_link_libraries(binprecice PRIVATE MPI::MPI_CXX)
endif()
if(MPI AND PETSC)
  target_include_directories(binprecice PRIVATE ${PETSC_INCLUDES})
  target_link_libraries(binprecice PRIVATE ${PETSC_LIBRARIES})
endif()


#
# Configuration of Target testprecice
#

add_executable(testprecice "src/testing/main.cpp")
target_link_libraries(testprecice
  PRIVATE
  Threads::Threads
  precice
  Eigen3::Eigen
  ${Boost_LIBRARIES}
  )
set_target_properties(testprecice PROPERTIES
  # precice is a C++11 project
  CXX_STANDARD 11
  )
# Copy needed properties from the lib to the executatble. This is necessary as
# this executable uses the library source, not only the interface.
copy_target_property(precice testprecice COMPILE_DEFINITIONS)
copy_target_property(precice testprecice COMPILE_OPTIONS)

# Testprecice fully depends on MPI and PETSc.
if(MPI)
  target_link_libraries(testprecice PRIVATE MPI::MPI_CXX)
endif()
if(MPI AND PETSC)
  target_include_directories(testprecice PRIVATE ${PETSC_INCLUDES})
  target_link_libraries(testprecice PRIVATE ${PETSC_LIBRARIES})
endif()

# Test Sources Configuration
include(${CMAKE_CURRENT_LIST_DIR}/src/tests.cmake)


#
# Install Targets for precice
#

# Setup General Install for:
# precice - the library
# binprecice - the precice binary
install(TARGETS precice binprecice
  EXPORT preciceTargets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include/precice
  INCLUDES DESTINATION include/precice
  )

if(PRECICE_InstallTest)
  # Install the testprecice target
  install(TARGETS testprecice
    EXPORT preciceTargets
    RUNTIME DESTINATION bin
    )

  # Install the resources necessary for the tests
  install(DIRECTORY src
    DESTINATION share/precice
    FILES_MATCHING
    PATTERN "*.xml"
    PATTERN "*.py"
    )
endif()

# Export the Targets to install
install(EXPORT preciceTargets
  FILE preciceTargets.cmake
  NAMESPACE precice::
  DESTINATION lib/cmake/precice
  )

# Generate a Package Config File for precice
include(CMakePackageConfigHelpers)
write_basic_package_version_file("preciceConfigVersion.cmake"
  VERSION ${preCICE_VERSION}
  COMPATIBILITY SameMajorVersion
  )

# Install the Config and the ConfigVersion files
install(FILES "cmake/preciceConfig.cmake" "${preCICE_BINARY_DIR}/preciceConfigVersion.cmake"
  DESTINATION lib/cmake/precice
  )

# Setup the config in the build directory
export(EXPORT preciceTargets
  NAMESPACE precice::
  FILE "preciceTargets.cmake")
file(COPY "cmake/preciceConfig.cmake"
  DESTINATION "${preCICE_BINARY_DIR}")


# Set the directory used to prepare files for packaging
set(PRECICE_PACKAGING_DIR "${CMAKE_CURRENT_BINARY_DIR}/packaging")
mark_as_advanced(PRECICE_PACKAGING_DIR)

# Compress and install changelog and manpages
find_program(GZIP_EXE gzip DOC "The gzip executable")
if(GZIP_EXE)
  # Process the changelog
  message(STATUS "Compressing changelog")
  file(COPY CHANGELOG.md DESTINATION ${PRECICE_PACKAGING_DIR})
  execute_process(COMMAND "${GZIP_EXE}" "-9nf" "${PRECICE_PACKAGING_DIR}/CHANGELOG.md")
  # Install compressed changelog
  install(FILES ${PRECICE_PACKAGING_DIR}/CHANGELOG.md.gz
    DESTINATION share/doc/libprecice${preCICE_VERSION}
    RENAME changelog.gz
    )

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
    DESTINATION share/man
    )
else()
  message(WARNING "Installing uncompressed changelog and manpages")
  # Install uncompressed changelog
  install(FILES CHANGELOG.md
    DESTINATION share/doc/libprecice${preCICE_VERSION}
    RENAME changelog
    )
  # Install uncompressed manpages
  install(DIRECTORY docs/man/man1
    DESTINATION share/man
    )
endif()

# Configure a pkg-config file
configure_file(
  "${PROJECT_SOURCE_DIR}/tools/packaging/debian/precice.pc.in"
  "lib/pkgconfig/libprecice.pc"
  @ONLY
  )
install(DIRECTORY "${preCICE_BINARY_DIR}/lib/pkgconfig" 
  DESTINATION lib
)


#
# CTest
#

include(${CMAKE_CURRENT_LIST_DIR}/cmake/CTestConfig.cmake)


#
# CPack
#

if (PRECICE_Packages)
  include(${CMAKE_CURRENT_LIST_DIR}/cmake/CPackConfig.cmake)
endif()
