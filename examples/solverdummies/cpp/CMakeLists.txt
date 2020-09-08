cmake_minimum_required(VERSION 3.10.2)

list (APPEND CMAKE_MODULE_PATH "../../../cmake")
list (APPEND CMAKE_MODULE_PATH "../../../cmake/modules")

project(SolverdummyCPP LANGUAGES CXX DESCRIPTION "preCICE CPP solver dummy")

find_package(precice REQUIRED CONFIG)

find_package(Boost 1.65.1 REQUIRED
  COMPONENTS filesystem log log_setup program_options system thread unit_test_framework
  )

find_package(Eigen3 3.2 REQUIRED)
if (PRECICE_MPICommunication)
  find_package(MPI REQUIRED)
endif()

if (PRECICE_PythonActions)
  find_package(PythonLibs REQUIRED)
  find_package(NumPy REQUIRED)
endif()

if (PRECICE_PETScMapping)
  find_package(PETSc 3.12 REQUIRED)
endif()

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../../thirdparty/json
                 ${CMAKE_CURRENT_BINARY_DIR}/thirdparty/json)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../../thirdparty/prettyprint
                 ${CMAKE_CURRENT_BINARY_DIR}/thirdparty/prettyprint)

add_executable(solverdummy solverdummy.cpp)
set_target_properties(solverdummy PROPERTIES CXX_STANDARD 11)
target_link_libraries(solverdummy PRIVATE precice::precice)


enable_testing()
add_test(NAME precice.solverdummy
  COMMAND ${CMAKE_COMMAND}
  -DDUMMY_EXE=$<TARGET_FILE:solverdummy>
  -DDUMMY_CONFIG=${CMAKE_CURRENT_SOURCE_DIR}/../precice-config.xml
  -V
  -P ${CMAKE_CURRENT_SOURCE_DIR}/test.cmake
  )
