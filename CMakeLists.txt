cmake_minimum_required (VERSION 3.9.6)

project(preCICE VERSION 1.3.0 LANGUAGES CXX)


#
# PREAMBLE
#

# Make our cmake modules visible to CMake
list (APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/tools/cmake-modules")

include(CopyTargetProperty)

# Enforce a build type
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()
message(STATUS "Build configuration: " ${CMAKE_BUILD_TYPE})


#
# Prepare Source and Binary Directories
#

# Create a symlink to the last build from build/last directory.
file(MAKE_DIRECTORY "${PROJECT_SOURCE_DIR}/build")
add_custom_target(symlink ALL
  COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_BINARY_DIR} ${PROJECT_SOURCE_DIR}/build/last)


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

find_package(Boost 1.60.0 REQUIRED
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
  string(REGEX MATCHALL "[0-9]+" PETSC_VERSIONS ${PETSC_VERSION})
  list(GET PETSC_VERSIONS 0 PETSC_MAJOR)
  list(GET PETSC_VERSIONS 1 PETSC_MINOR)
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
file(WRITE dummy.cpp "")

# Add precice as an empty target 
add_library(precice dummy.cpp)
set_target_properties(precice PROPERTIES
    # precice is a C++11 project
    CXX_STANDARD 11
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
add_test(NAME precice.All COMMAND testprecice)

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
