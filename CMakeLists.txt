cmake_minimum_required (VERSION 3.9.6)

project(preCICE VERSION 1.3.0 LANGUAGES CXX)

#
# PREAMBLE
#

# Make our cmake modules visible to CMake
list (APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/tools/cmake-modules")

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

    # Create Target for PythonLibs
    add_library(PythonLibs INTERFACE IMPORTED)
    set_property(TARGET PythonLibs PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${PYTHON_INCLUDE_DIRS})
    set_property(TARGET PythonLibs PROPERTY INTERFACE_COMPILE_DEFINITIONS NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION)
    set_property(TARGET PythonLibs PROPERTY INTERFACE_LINK_LIBRARIES ${PYTHON_LIBRARIES})
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
target_compile_definitions(precice PUBLIC BOOST_ALL_DYN_LINK BOOST_ASIO_ENABLE_OLD_SERVICES)
target_link_libraries(precice PUBLIC ${Boost_LIBRARIES})

# Setup Eigen3
target_link_libraries(precice PRIVATE Eigen3::Eigen)
target_compile_definitions(precice PRIVATE "$<$<CONFIG:DEBUG>:EIGEN_INITIALIZE_MATRICES_BY_NAN>")

# Setup LIBXML2
target_include_directories(precice PRIVATE ${LIBXML2_INCLUDE_DIR})
target_link_libraries(precice PUBLIC ${LIBXML2_LIBRARIES})

if (Platform STREQUAL "SuperMUC")
    target_compile_definitions(precice PRIVATE SuperMUC_WORK)
endif()

# Setup MPI
if (MPI)
    target_link_libraries(precice PUBLIC MPI::MPI_CXX)
else()
    target_compile_definitions(precice PUBLIC PRECICE_NO_MPI)
endif()

# Setup PETSC
if (PETSC AND MPI)
    target_include_directories(precice PUBLIC ${PETSC_INCLUDES})
    target_link_libraries(precice PUBLIC ${PETSC_LIBRARIES})
else()
    target_compile_definitions(precice PUBLIC PRECICE_NO_PETSC)
endif()

# Option Python
if (PYTHON)
    target_link_libraries(precice PRIVATE PythonLibs)
    target_link_libraries(precice PRIVATE NumPy::NumPy)
else()
    target_compile_definitions(precice PUBLIC PRECICE_NO_PYTHON)
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
    )
set_target_properties(binprecice PROPERTIES
    # precice is a C++11 project
    CXX_STANDARD 11
    )
copy_target_property(precice binprecice COMPILE_DEFINITIONS)
copy_target_property(precice binprecice COMPILE_OPTIONS)

#
# Configuration of Target testprecice
#

add_executable(testprecice "src/testing/main.cpp")
target_link_libraries(testprecice
    PRIVATE
    Threads::Threads
    precice
    Eigen3::Eigen
    )
set_target_properties(testprecice PROPERTIES
    # precice is a C++11 project
    CXX_STANDARD 11
    )
copy_target_property(precice testprecice COMPILE_DEFINITIONS)
copy_target_property(precice testprecice COMPILE_OPTIONS)

# Add the tests
include(${CMAKE_CURRENT_LIST_DIR}/src/tests.cmake)

# Register the test using CTest
enable_testing()
set(PRECICE_TEST_DIR "${preCICE_BINARY_DIR}/TestOutput")
file(MAKE_DIRECTORY ${PRECICE_TEST_DIR})
add_test(NAME precice.Base
    COMMAND "$<TARGET_FILE:testprecice>"
    WORKING_DIRECTORY "${PRECICE_TEST_DIR}"
    )
if(MPI AND MPIEXEC_EXECUTABLE)
  file(MAKE_DIRECTORY "${PRECICE_TEST_DIR}/MPI2")
  add_test(NAME precice.MPI2
    COMMAND "${MPIEXEC_EXECUTABLE}" -np 2 "$<TARGET_FILE:testprecice>"
    WORKING_DIRECTORY "${PRECICE_TEST_DIR}/MPI2"
    )
  file(MAKE_DIRECTORY "${PRECICE_TEST_DIR}/MPI4")
  add_test(NAME precice.MPI4
    COMMAND "${MPIEXEC_EXECUTABLE}" -np 4 "$<TARGET_FILE:testprecice>"
    WORKING_DIRECTORY "${PRECICE_TEST_DIR}/MPI4"
    )
endif()

#
# Install Targets for precice
#

# Setup General Install for:
# precice - the library
# binprecice - the precice binary
# testprecice - the test for smoke testing
install(TARGETS precice binprecice testprecice
  EXPORT preciceTargets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include/precice
  INCLUDES DESTINATION include/precice
  )

# Install the resources necessary for the tests
install(DIRECTORY src
    DESTINATION share/precice
    FILES_MATCHING
    PATTERN "*.xml"
    PATTERN "*.py"
    )

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
install(FILES "preciceConfig.cmake" "${preCICE_BINARY_DIR}/preciceConfigVersion.cmake"
    DESTINATION lib/cmake/precice
    )

# Install doc files
install(FILES tools/packaging/debian/copyright
    DESTINATION share/doc/libprecice${preCICE_VERSION}
  )

# Compress and install changelog and manpages
find_program(GZIP_EXE gzip DOC "The gzip executable")
if(GZIP_EXE)
  set(PRECICE_PACKAGING_DIR "${CMAKE_CURRENT_BINARY_DIR}/packaging")
  message(STATUS "Compressing changelog")
  file(COPY CHANGELOG.md DESTINATION packaging)
  execute_process(COMMAND "${GZIP_EXE}" "-9nf" "${PRECICE_PACKAGING_DIR}/CHANGELOG.md")
  install(FILES ${PRECICE_PACKAGING_DIR}/CHANGELOG.md.gz
    DESTINATION share/doc/libprecice${preCICE_VERSION}
    RENAME changelog.gz
    )

  file(COPY tools/packaging/man1 DESTINATION packaging)
  file(GLOB PRECICE_MAN_PAGES "${PRECICE_PACKAGING_DIR}/man1/*.1")
  foreach(manpage ${PRECICE_MAN_PAGES})
    message(STATUS "Compressing manpage: ${manpage}")
    execute_process(COMMAND "${GZIP_EXE}" "-9nf" "${manpage}")
  endforeach()
  install(DIRECTORY ${PRECICE_PACKAGING_DIR}/man1
    DESTINATION share/man
  )
else()
  message(WARNING "Installing uncompressed changelog and manpages")
  install(FILES CHANGELOG.md
    DESTINATION share/doc/libprecice${preCICE_VERSION}
    RENAME changelog
    )
  install(DIRECTORY tools/packaging/man1
    DESTINATION share/man
    )
endif()


#
# Packaging Settings for CPack
#
# Variable References:
#  General https://cmake.org/cmake/help/latest/module/CPack.html#module:CPack
#  Debian Packages https://cmake.org/cmake/help/latest/cpack_gen/deb.html#cpack_gen:CPack%20DEB%20Generator
#

# Detect the system name
if(WIN32)
    set(CPACK_SYSTEM_NAME "win32")
elseif(WIN64)
    set(CPACK_SYSTEM_NAME "win64")
else()
    set(CPACK_SYSTEM_NAME "${CMAKE_SYSTEM_NAME}")
endif()

# General
set(CPACK_PACKAGE_NAME "libprecice${preCICE_VERSION}")
set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_SYSTEM_NAME}")
set(CPACK_PACKAGE_VENDOR "precice.org")
set(CPACK_PACKAGE_CONTACT "The precice developers <precice@mailman.informatik.uni-stuttgart.de>")
set(CPACK_PACKAGE_MAINTAINER "The precice developers <precice@mailman.informatik.uni-stuttgart.de>")
set(CPACK_PACKAGE_DESCRIPTION "preCICE (Precise Code Interaction Coupling Environment) is a coupling library for partitioned multi-physics simulations, including, but not restricted to fluid-structure interaction and conjugate heat transfer simulations. Partitioned means that preCICE couples existing programs (solvers) capable of simulating a subpart of the complete physics involved in a simulation. This allows for the high flexibility that is needed to keep a decent time-to-solution for complex multi-physics scenarios.")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Precise Code Interaction Coupling Environment")
set(CPACK_PACKAGE_EXECUTABLES "testprecice;binprecice")
set(CPACK_PACKAGE_HOMEPAGE_URL "www.precice.org")
#set(CPACK_PACKAGE_ICON "")
set(CPACK_PACKAGE_CHECKSUM "SHA256")
set(CPACK_RESOURCE_FILE_LICENSE "${preCICE_SOURCE_DIR}/LICENSE")
set(CPACK_MONOLITHIC_INSTALL TRUE)
set(CPACK_GENERATOR "DEB;TGZ;ZIP")
set(CPACK_STRIP_FILES TRUE)

#set(CPACK_SOURCE_PACKAGE_FILE_NAME "")
set(CPACK_SOURCE_GENERATOR ${CPACK_GENERATOR})
set(CPACK_SOURCE_IGNORE_FILES
  "/build/"
  "/.git/"
  ".gitignore"
  )

set(CPACK_DEBIAN_PACKAGE_DEPENDS "libc6, petsc-dev (>= 3.8), libboost-dev (>= 1.68), libboost-log-dev (>= 1.68), libboost-thread-dev (>= 1.68), libboost-system-dev (>= 1.68), libboost-filesystem-dev (>= 1.68), libboost-program-options-dev (>= 1.68), libboost-test-dev (>= 1.68), libeigen3-dev, libxml2-dev, python-dev, python-numpy")
set(CPACK_DEBIAN_PACKAGE_SECTION "devel")
set(CPACK_DEBIAN_PACKAGE_DESCRIPTION "\
Precise Code Interaction Coupling Environment\n\
 preCICE (Precise Code Interaction Coupling Environment) is a coupling library\n\
 for partitioned multi-physics simulations, including, but not restricted to\n\
 fluid-structure interaction and conjugate heat transfer simulations.\n\
 Partitioned means that preCICE couples existing programs (solvers) capable of\n\
 simulating a subpart of the complete physics involved in a simulation.\n\
 This allows for the high flexibility that is needed to keep a decent\n\
 time-to-solution for complex multi-physics scenarios.\
")
set(CPACK_DEBIAN_PACKAGE_CONTROL_EXTRA "${preCICE_SOURCE_DIR}/tools/packaging/debian/postinst;${preCICE_SOURCE_DIR}/tools/packaging/debian/postrm")
set(CPACK_DEBIAN_PACKAGE_CONTROL_STRUCT_PERMISSION TRUE)
set(CPACK_DEBIAN_PACKAGE_GENERATE_SHLIBS TRUE)
set(CPACK_DEBIAN_PACKAGE_GENERATE_SHLIBS_POLICY "=")

include(CPack)
