cmake_minimum_required (VERSION 3.5)

project(preCICE VERSION 1.2.0 LANGUAGES CXX)

list (APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/tools/cmake-modules")

# Add preCICE as an empty target 
add_library(precice)

set_property(TARGET precice PROPERTY CXX_STANDARD 11)

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()
message(STATUS "Build configuration: " ${CMAKE_BUILD_TYPE})


if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unused-parameter")
endif()


set(Platform "none" CACHE STRING "Special configurations for certain platforms")
set_property(CACHE Platform PROPERTY STRINGS none "SuperMUC")
message(STATUS "Special platform settings: " ${Platform})
if (Platform STREQUAL "SuperMUC")
    target_compile_definitions(precice PRIVATE SuperMUC_WORK)
endif()


option(BUILD_SHARED_LIBS "Build a shared library" ON)


option(MPI "Enables MPI-based communication and running coupling tests." ON)
if (MPI)
    find_package(MPI REQUIRED)
    target_link_libraries(precice PRIVATE MPI::MPI_CXX)
    #  include_directories(${MPI_INCLUDE_PATH})
    target_compile_options(precice PRIVATE ${MPI_COMPILE_FLAGS})
    #  set(COMPILE_FLAGS  ${COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
    target_link_libraries(precice PRIVATE ${MPI_LINK_FLAGS})
    #  set(LINK_FLAGS ${LINK_FLAGS} ${MPI_LINK_FLAGS})
target_link_libraries(precice PUBLIC ${MPI_LIBRARIES})
else()
    target_compile_definitions(precice PRIVATE PRECICE_NO_MPI)
    message(STATUS "MPI support disabled")
endif()

option(PETSC "Enable use of the PETSc linear algebra library." ON)
if (PETSC AND MPI)
  find_package(PETSc 3.6 REQUIRED)
  string(REGEX MATCHALL "[0-9]+" PETSC_VERSIONS ${PETSC_VERSION})
  list(GET PETSC_VERSIONS 0 PETSC_MAJOR)
  list(GET PETSC_VERSIONS 1 PETSC_MINOR)
  target_include_directories(precice PUBLIC ${PETSC_INCLUDES})
  target_link_libraries(precice PUBLIC ${PETSC_LIBRARIES})
else()
    target_compile_definitions(precice PRIVATE PRECICE_NO_PETSC)
    message(STATUS "PETSc support disabled")
endif()

find_package (Threads REQUIRED)


find_package(Boost 1.60.0 REQUIRED
  COMPONENTS filesystem log log_setup program_options system thread unit_test_framework)
target_compile_definitions(precice PRIVATE BOOST_ALL_DYN_LINK BOOST_ASIO_ENABLE_OLD_SERVICES)
target_link_libraries(precice PUBLIC ${Boost_LIBRARIES})


find_package(Eigen3 3.2 REQUIRED)
target_include_directories(precice PRIVATE ${EIGEN3_INCLUDE_DIR})
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_definitions(precice PRIVATE EIGEN_INITIALIZE_MATRICES_BY_NAN)
endif()


find_package(LibXml2 REQUIRED)
target_include_directories(precice PRIVATE ${LIBXML2_INCLUDE_DIR})
target_link_libraries(precice PUBLIC ${LIBXML2_LIBRARIES})


option(PYTHON "Python support" ON)
if (PYTHON)
    set(Python_ADDITIONAL_VERSIONS "2.7")
    find_package(PythonLibs 2.7 REQUIRED)
    find_package(NumPy REQUIRED)
    target_include_directories(precice PRIVATE ${PYTHON_INCLUDE_DIRS})
    target_include_directories(precice PRIVATE ${NumPy_INCLUDE_DIRS})
    target_compile_definitions(precice PRIVATE NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION)
    target_link_libraries(precice PRIVATE ${PYTHON_LIBRARIES})
    # Mute register warning
    target_compile_options(precice PRIVATE -Wno-register)
else()
    target_compile_definitions(precice PRIVATE PRECICE_NO_PYTHON)
    message(STATUS "Python support disabled")
endif()


# Create versions.hpp and add binary dir to include dirs
configure_file("${PROJECT_SOURCE_DIR}/src/versions.hpp.in" "${PROJECT_BINARY_DIR}/versions.hpp" @ONLY)
target_include_directories(precice PUBLIC ${PROJECT_BINARY_DIR})

add_subdirectory("src")
target_include_directories(precice PUBLIC "src")

target_sources(precice PRIVATE ${sourcesAllNoMain})
set(headers
  src/precice/SolverInterface.hpp
  src/precice/MeshHandle.hpp
  src/precice/Constants.hpp
  src/precice/bindings/c/SolverInterfaceC.h
  src/precice/bindings/c/Constants.h)
set_target_properties(precice PROPERTIES PUBLIC_HEADER "${headers}")

add_executable(binprecice "src/drivers/main.cpp")
target_link_libraries(binprecice Threads::Threads)
target_link_libraries(binprecice precice)

add_executable(testprecice
  "src/testing/main.cpp" ${sourcesTests})
target_link_libraries(testprecice Threads::Threads)
target_link_libraries(testprecice precice)

install(TARGETS precice binprecice
  LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION include/precice
  RUNTIME DESTINATION bin)
