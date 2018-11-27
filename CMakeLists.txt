cmake_minimum_required (VERSION 3.5)

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
else()
    message(STATUS "Python support disabled")
endif()


#
# Configuration of Target PRECICE
#

# Add PRECICE as an empty target 
add_library(PRECICE)
set_target_properties(PRECICE PROPERTIES
    # PRECICE is a C++11 project
    CXX_STANDARD 11
    # Set the binary name
    OUTPUT_NAME precice
    )

# Setup Boost
target_compile_definitions(PRECICE PRIVATE BOOST_ALL_DYN_LINK BOOST_ASIO_ENABLE_OLD_SERVICES)
target_link_libraries(PRECICE PUBLIC ${Boost_LIBRARIES})

# Setup Eigen3
target_link_libraries(PRECICE PRIVATE Eigen3::Eigen)
target_compile_definitions(PRECICE PRIVATE "$<$<CONFIG:DEBUG>:EIGEN_INITIALIZE_MATRICES_BY_NAN>")

# Setup LIBXML2
target_include_directories(PRECICE PRIVATE ${LIBXML2_INCLUDE_DIR})
target_link_libraries(PRECICE PUBLIC ${LIBXML2_LIBRARIES})

if (Platform STREQUAL "SuperMUC")
    target_compile_definitions(PRECICE PRIVATE SuperMUC_WORK)
endif()

# Setup MPI
if (MPI)
    target_link_libraries(PRECICE PRIVATE MPI::MPI_CXX)
    target_compile_options(PRECICE PRIVATE ${MPI_COMPILE_FLAGS})
    target_link_libraries(PRECICE PRIVATE ${MPI_LINK_FLAGS})
    target_link_libraries(PRECICE PUBLIC ${MPI_LIBRARIES})
else()
    target_compile_definitions(PRECICE PRIVATE PRECICE_NO_MPI)
endif()

# Setup PETSC
if (PETSC AND MPI)
    target_include_directories(PRECICE PUBLIC ${PETSC_INCLUDES})
    target_link_libraries(PRECICE PUBLIC ${PETSC_LIBRARIES})
else()
    target_compile_definitions(PRECICE PRIVATE PRECICE_NO_PETSC)
endif()

# Option Python
if (PYTHON)
    target_include_directories(PRECICE PRIVATE ${PYTHON_INCLUDE_DIRS})
    target_include_directories(PRECICE PRIVATE ${NumPy_INCLUDE_DIRS})
    target_compile_definitions(PRECICE PRIVATE NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION)
    target_link_libraries(PRECICE PRIVATE ${PYTHON_LIBRARIES})
else()
    target_compile_definitions(PRECICE PRIVATE PRECICE_NO_PYTHON)
endif()

# File Configuration
configure_file("${PROJECT_SOURCE_DIR}/src/versions.hpp.in" "${PROJECT_BINARY_DIR}/src/versions.hpp" @ONLY)

# Includes Configuration
target_include_directories(PRECICE PUBLIC 
    $<BUILD_INTERFACE:${preCICE_SOURCE_DIR}/src>
    $<BUILD_INTERFACE:${preCICE_BINARY_DIR}>
    $<INSTALL_INTERFACE:include>
    )

# Sources Configuration
include(${CMAKE_CURRENT_LIST_DIR}/src/sources.cmake)


#
# Configuration of Target PRECICEBIN
#

add_executable(PRECICEBIN "src/drivers/main.cpp")
target_link_libraries(PRECICEBIN 
    PRIVATE
    Threads::Threads
    PRECICE
    Eigen3::Eigen
    )
set_target_properties(PRECICEBIN PROPERTIES
    # PRECICE is a C++11 project
    CXX_STANDARD 11
    # Set the binary name
    OUTPUT_NAME binprecice
    )
copy_target_property(PRECICE PRECICEBIN COMPILE_DEFINITIONS)

#
# Configuration of Target PRECICETEST
#

add_executable(PRECICETEST "src/testing/main.cpp")
target_link_libraries(PRECICETEST
    PRIVATE
    Threads::Threads
    PRECICE
    Eigen3::Eigen
    )
set_target_properties(PRECICETEST PROPERTIES
    # PRECICE is a C++11 project
    CXX_STANDARD 11
    # Set the binary name
    OUTPUT_NAME testprecice
    )
copy_target_property(PRECICE PRECICETEST COMPILE_DEFINITIONS)
copy_target_property(PRECICE PRECICETEST COMPILE_OPTIONS)

# Add the tests
include(${CMAKE_CURRENT_LIST_DIR}/src/tests.cmake)

# Register the test using CTest
enable_testing()
add_test(NAME PRECICE.All COMMAND PRECICETEST)

#
# Install Targets for PRECICE
#

# Setup General Install for:
# PRECICE - the library
# binprecice - the precice binary
# testprecice - the test for smoke testing
install(TARGETS PRECICE PRECICEBIN PRECICETEST
  EXPORT PRECICETargets
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
install(EXPORT PRECICETargets
    FILE PRECICETargets.cmake
    NAMESPACE PRECICE::
    DESTINATION lib/cmake/precice
    )

# Generate a Package Config File for PRECICE
include(CMakePackageConfigHelpers)
write_basic_package_version_file("PRECICEConfigVersion.cmake"
    VERSION ${PRECICE_VERSION}
    COMPATIBILITY SameMajorVersion
    )

# Install the Config and the ConfigVersion files
install(FILES "PRECICEConfig.cmake" "${preCICE_BINARY_DIR}/PRECICEConfigVersion.cmake"
    DESTINATION lib/cmake/precice
    )
